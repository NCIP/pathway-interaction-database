#!/usr/bin/perl


# Copyright SRA International
#
# Distributed under the OSI-approved BSD 3-Clause License.
# See http://ncip.github.com/pathway-interaction-database/LICENSE.txt for details.


use strict;

BEGIN {
  my @path_elems = split("/", $0);
  pop @path_elems;
  push @INC, join("/", @path_elems);
}

use CGI;
use DBI;
use PWAppConfig;
use PathwayDB;
use Pathway;
use PWLabel;
use EquationOutput;

if (-d "/app/oracle/product/dbhome/current") {
  $ENV{'ORACLE_HOME'} = "/app/oracle/product/dbhome/current";
} elsif (-d "/app/oracle/product/8.1.7") {
  $ENV{'ORACLE_HOME'} = "/app/oracle/product/8.1.7";
} elsif (-d "/app/oracle/product/8.1.6") {
  $ENV{'ORACLE_HOME'} = "/app/oracle/product/8.1.6";
}

#my $geneid = shift @ARGV;

my $query = new CGI;
my $geneid = $query->param("geneid");

my ($db_inst, $db_user, $db_pass, $pid_schema, $cgap_schema) =
    (DB_INST, DB_USER, DB_PASS, SCHEMA, "cgap");

print "Content-type: text/html\n\n";

$geneid =~ s/\s+//g;
if ($geneid !~ /^\d+$/) {
  print STDERR "illegal value for parameter geneid";
  exit;
}

my $db = DBI->connect("DBI:Oracle:" . $db_inst, $db_user, $db_pass);
if (not $db or $db->err()) {
  print STDERR "Cannot connect to " . $db_user . "@" . $db_inst . "\n";
  exit;
}

my ($lv, $pw);

my ($mol2name, $mol2ext_id);

my $mols = FindMolIdsForGene($db, $geneid);
if (keys %{ $mols } > 0) {
  ($mol2name, $mol2ext_id) = FindNamesForMols($db, $mols);
} else {
  print "No entries for Gene id = $geneid\n";
}
$db->disconnect();

print FormatOutput($mols, $mol2name, $mol2ext_id);

exit;

###################

if (keys %{ $mols } > 0) {
  $lv   = new PWLabel($db, $pid_schema);
  $pw   = new Pathway($lv);
  FindReactionsForMolId($db, $lv, $pw, [ keys %{ $mols } ]);
  if (@{ $pw->Atoms() } > 0) {
#    EquationOutput();
    my $eq = new EquationOutput($pw, $lv);
    $eq->EquationOutput();
    print @{ $eq->Lines() };
    print "\n";
  } else {
    print "No entries for Gene id = $geneid\n";
  }
} else {
  print "No entries for Gene id = $geneid\n";
}
$db->disconnect();


######################################################################
sub numerically { $a <=> $b; }

######################################################################
sub MOLPAGE_URL {
  my ($molid) = @_;
#  return "http://pid.nci.nih.gov/search/MoleculePage?molid=$molid";
  return "/search/MoleculePage?molid=$molid";
}

######################################################################
sub ATOMPAGE_URL {
  my ($atomid) = @_;
#  return "http://pid.nci.nih.gov/search/InteractionPage?atomid=$atomid";
  return "/search/InteractionPage?atomid=$atomid";
}

######################################################################
sub FormatOutput {
  my ($mols, $mol2name, $mol2ext_id) = @_;

  my @lines;
  for my $mol (keys %{ $mols }) {
    my (@lls, @ups, @names);
    for my $idtype (keys %{ $$mol2ext_id{$mol} }) {
      for my $id (keys %{ $$mol2ext_id{$mol}{$idtype} }) {
        if ($idtype eq "LL") {
          push @lls, "Entrez Gene: $id";
        } elsif ($idtype eq "UP") {
          push @ups, "UniProt: $id";
        }
      }
    }
    for my $nametype (keys %{ $$mol2name{$mol} }) {
      for my $name (keys %{ $$mol2name{$mol}{$nametype} }) {
        push @names, $name;
      }
    }
    push @lines, "<a href=\"" . MOLPAGE_URL($mol) . "\">" .
        join(", ", @names, @lls, @ups) . "</a> ";
  }
  return join("<br>\n", @lines, "");
}


######################################################################
sub FindReactionsForMolId {
  my ($db, $lv, $pw, $mols) = @_;

  use constant SOURCE_LIST          => "5";
  use constant EVIDENCE_CODE_LIST   => "";
  use constant INCLUDE_COMPLEX_USES => "1";

  my $pdb = new PathwayDB($db, $pid_schema,
      SOURCE_LIST,
      EVIDENCE_CODE_LIST, $pw, $lv);

  $pdb->AtomsOfMols(join(",", @{ $mols }), INCLUDE_COMPLEX_USES);

  $pw->MergeMolecules();
  $pw->BuildMolInstCache();
  $pw->BuildGenericLabelCache();
  $pw->PruneDuplicateAtoms();
  $pw->IdentifyMacroProcesses();
}

######################################################################
sub FindNamesForMols {
  my ($db, $mols) = @_;

  my $mollist = join(",", keys %{ $mols });
  my (%mol2name, %mol2ext_id);

  my ($sql, $stm);
  my ($molid, $extid, $id_type, $name_type, $name);

  $sql = qq!
select
  e.mol_id,
  e.id_type,
  e.ext_mol_id
from
  $pid_schema.pw_ext_mol_id e
where
     e.mol_id in ($mollist)
!;  

  $stm = $db->prepare($sql);
  if(not $stm) {
    print STDERR "prepare failed: $sql\n";
    return;
  }
  if(!$stm->execute()) {
    print STDERR "execute failed: $sql\n";
    return;
  }
  while (($molid, $id_type, $extid) = $stm->fetchrow_array()) {
    $mol2ext_id{$molid}{$id_type}{$extid} = 1;
  }
  
  $sql = qq!
select
  n.mol_id,
  n.name_type,
  n.mol_name
from
  $pid_schema.pw_mol_name n
where
     n.mol_id in ($mollist)
!;  

  $stm = $db->prepare($sql);
  if(not $stm) {
    print STDERR "prepare failed: $sql\n";
    return;
  }
  if(!$stm->execute()) {
    print STDERR "execute failed: $sql\n";
    return;
  }
  while (($molid, $name_type, $name) = $stm->fetchrow_array()) {
    $mol2name{$molid}{$name_type}{$name} = 1;
  }

  return (\%mol2name, \%mol2ext_id);

}

######################################################################
sub FindMolIdsForGene {
  my ($db, $geneid) = @_;

  my ($sql, $stm);
  my ($molid, $extid);
  my (%mols);

  ## Find direct matches of ext ids

  $sql = qq!
select
  e.mol_id,
  e.ext_mol_id
from
  $pid_schema.pw_ext_mol_id e,
  $pid_schema.pw_mol m
where
      e.ext_mol_id = '$geneid'
  and e.mol_id = m.mol_id
  and m.basic_mol_type = 'PR'
!;

  $stm = $db->prepare($sql);
  if(not $stm) {
    print STDERR "prepare failed: $sql\n";
    return;
  }
  if(!$stm->execute()) {
    print STDERR "execute failed: $sql\n";
    return;
  }
  while (($molid, $extid) = $stm->fetchrow_array()) {
    $mols{$molid} = 1;
  }
  
  ## Find matches of gene ids to proteins in PID

    $sql = qq!
select
  e.mol_id,
  g.ll_id
from
  $pid_schema.pw_ext_mol_id e,
  $pid_schema.pw_mol m,
  $cgap_schema.ll2sp s,
  $cgap_schema.sp_primary p,
  $cgap_schema.ll_gene g
where
      p.sp_id_or_secondary = e.ext_mol_id
  and p.sp_id_or_secondary = s.sp_primary
  and s.organism = 'Hs'
  and s.ll_id = g.ll_id
  and g.ll_id ='$geneid'
  and e.mol_id = m.mol_id
  and m.basic_mol_type = 'PR'
  !;

  $stm = $db->prepare($sql);
  if(not $stm) {
    print STDERR "prepare failed: $sql\n";
    return;
  }
  if(!$stm->execute()) {
    print STDERR "execute failed: $sql\n";
    return;
  }
  while (($molid, $extid) = $stm->fetchrow_array()) {
    $mols{$molid} = 1;
  }
  return \%mols;
}
