#!/usr/local/bin/perl


# Copyright SRA International
#
# Distributed under the OSI-approved BSD 3-Clause License.
# See http://ncip.github.com/pathway-interaction-database/LICENSE.txt for details.


######################################################################
# MoleculePage.pl
#
# internal molecule id
# external molecule ids
# molecule names
# used-in-complex table:
# used-in-atom table:
#   check to include
#   role (input, output, agent, inhibitor)
#   location
#   state
#   atom (hyperlink)
#   pathway (hyperlink)
#   source
#

BEGIN {
  my @path_elems = split("/", $0);
  pop @path_elems;
  push @INC, join("/", @path_elems);
}

use strict;
use CGI;
use Pathway;
use PWLabel;
use PathwayDB;
use DBI;

if (-d "/app/oracle/product/dbhome/current") {
  $ENV{'ORACLE_HOME'} = "/app/oracle/product/dbhome/current";
} elsif (-d "/app/oracle/product/8.1.7") {
  $ENV{'ORACLE_HOME'} = "/app/oracle/product/8.1.7";
} elsif (-d "/app/oracle/product/8.1.6") {
  $ENV{'ORACLE_HOME'} = "/app/oracle/product/8.1.6";
}

my $query     = new CGI;
my $base      = $query->param("BASE");
my $molid     = $query->param("molid");

my $BASE;

##!!!!! for now:
my $LIMIT_SOURCES = "1,2,3,5";

my (%atom_uses, %complex_uses, %location, %state,
    %atom2pathway, %atom2source, %pathway2name, %pathway2source);

my ($db_inst, $db_user, $db_pass,   $schema) =
   ("cgprod", "web",    "readonly", "cgap");

my $db = DBI->connect("DBI:Oracle:" . $db_inst, $db_user, $db_pass);
if (not $db or $db->err()) {
  print STDERR "Cannot connect to " . $db_user . "@" . $db_inst . "\n";
  exit;
}

my $lv   = new PWLabel($db, $schema);
my $pw   = new Pathway($lv);
my $pdb  = new PathwayDB($db, $schema, $LIMIT_SOURCES, $pw, $lv);

my $COMPLEX_TYPE  = $lv->StringToLabelValue("molecule-type", "complex");
my $PROTEIN_TYPE  = $lv->StringToLabelValue("molecule-type", "protein");
my $RNA_TYPE      = $lv->StringToLabelValue("molecule-type", "rna");
my $COMPOUND_TYPE = $lv->StringToLabelValue("molecule-type", "compound");

print "Content-type: text/html\n\n";

print "<html>\n";
print "<head><title>Molecule Page</title></head>\n";
print "<body>\n";

MoleculePage_1 ($base, $db, $schema, $molid);

$db->disconnect();

######################################################################
sub LL_URL {
  my ($locus) = @_;
  return "http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?" .
      "db=gene&cmd=Retrieve&dopt=Graphics&" .
      "list_uids=$locus";
}

######################################################################
sub MOLPAGE_URL {
  my ($molid) = @_;
  return "$BASE/cgi-bin/MoleculePage?" .
      "molid=$molid";
}

######################################################################
sub ATOM_URL {
  my ($atomid) = @_;
  return "$BASE/DrawPathway?" .
      "atom_id=$atomid&what=graphic&svg=on";
}

######################################################################
sub PATHWAY_URL {
  my ($pathwayid) = @_;
  return "$BASE/DrawPathway?" .
      "pathway_id=$pathwayid&what=graphic&svg=on";
}

######################################################################
sub MoleculePage_1 {
  my ($base, $db, $schema, $molid) = @_;

  $BASE = $base;

  $pdb->AtomsOfMols($molid);

  SimpleOrComplex($molid);
  PathwayAndSource($db, $schema);
  IdStuff($molid);
  GOStuff($molid);
  ComplexTable($molid);
  AtomTable($molid);
}

######################################################################
sub GOStuff {
  my ($molid) = @_;

  my ($db_inst, $db_user, $db_pass,   $schema) =
     ("lpgdev", "web",    "readonly", "schaefec");
  my ($locus);
  my ($sql, $stm);
  my ($go_id, $go_name, $evidence, @lines);
  
  my $db = DBI->connect("DBI:Oracle:" . $db_inst, $db_user, $db_pass);
  if (not $db or $db->err()) {
    print STDERR "Cannot connect to " . $db_user . "@" . $db_inst . "\n";
    exit;
  }

  my $ehash = $pw->MolExId($molid);
  for my $id_type (keys %{ $ehash }) {
    for my $id (keys %{ $$ehash{$id_type} }) {
      if ($id_type eq "LL") {
        $locus = $id;
      }
    }
  }

  if (! defined $locus || ! $locus) {
    return;
  }

  $sql = qq!
select unique
  n.go_id,
  n.go_name,
  s.evidence
from
  $schema.go_name n,
  $schema.sptr_goa s,
  $schema.ll2sp l
where
      l.ll_id = $locus
  and l.sp_primary = s.sp_id
  and s.go_id = n.go_id
  and n.go_class = 'MF'
order by
  s.evidence, n.go_name
!;

  $stm = $db->prepare($sql);
  if(not $stm) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "prepare call failed\n";
    exit;
  }
  if(!$stm->execute()) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "execute call failed\n";
    exit;
  }
  while (($go_id, $go_name, $evidence) = $stm->fetchrow_array()) {
    push @lines, "<tr><td>GO:$go_id</td><td>$go_name</td><td>$evidence</td></tr>";
  }

  if (@lines) {
    print "<p>GO Molecular Function Annotations:\n";
    print "<blockquote>\n";
    print "<p><table border=1 cellspacing=1 cellpadding=4>\n";
    print "<tr><td><b>GO Id</b></td><td><b>GO Name</b></td><td><b>Evidence</b></td></tr>\n";
    print join("\n", @lines) . "\n";
    print "</table>\n";
    print "</blockquote>\n";
  }

  $db->disconnect();
}

######################################################################
sub PathwayAndSource {
  my ($db, $schema) = @_;

  my ($sql, $stm);
  my (%pathways);
  my ($atom_id, $pathway_id, $pathway_name, $source);

  $sql = "select pa.atom_id, pa.pathway_id " .
    "from $schema.pw_pathway_atom pa, $schema.pw_atom a " .
    "where pa.atom_id in (" . join(",", keys %atom_uses) . ") " .
    "and a.atom_source_id in ($LIMIT_SOURCES)";
  $stm = $db->prepare($sql);
  if(not $stm) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "prepare call failed\n";
    exit;
  }
  if(!$stm->execute()) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "execute call failed\n";
    exit;
  }
  while (($atom_id, $pathway_id) = $stm->fetchrow_array()) {
    $atom2pathway{$atom_id}{$pathway_id} = 1;
    $pathways{$pathway_id} = 1;
  }

  $sql = "select a.atom_id, s.source_name from $schema.pw_atom a, " .
    "$schema.pw_source s " .
    "where a.atom_id in (" . join(",", keys %atom_uses) . ") " .
    "and a.atom_source_id = s.source_id";
  $stm = $db->prepare($sql);
  if(not $stm) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "prepare call failed\n";
    exit;
  }
  if(!$stm->execute()) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "execute call failed\n";
    exit;
  }
  while (($atom_id, $source) = $stm->fetchrow_array()) {
    $atom2source{$atom_id} = $source;
  }

  $sql = "select p.pathway_id, p.pathway_name, s.source_name " .
    "from $schema.pw_pathway p, $schema.pw_source s " .
    "where p.pathway_id in (" . join(",", keys %pathways) . ") " .
    "and p.pathway_source_id = s.source_id";
  $stm = $db->prepare($sql);
  if(not $stm) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "prepare call failed\n";
    exit;
  }
  if(!$stm->execute()) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "execute call failed\n";
    exit;
  }
  while (($pathway_id, $pathway_name, $source) = $stm->fetchrow_array()) {
    $pathway2source{$pathway_id} = $source;
    $pathway2name{$pathway_id}   = $pathway_name;
  }

}

######################################################################
sub IdStuff {
  my ($molid) = @_;

  my (@names, @ids, $locus);
  my $nhash = $pw->MolName($molid);
  for my $name_type (keys %{ $nhash }) {
    for my $name (keys %{ $$nhash{$name_type} }) {
      push @names, $name;
    }
  }
  my ($label, $moltype) = $lv->LabelValueToString($pw->MolType($molid));
  if ($moltype eq "complex") {
    unshift @names, "&ltcx_$molid&gt;";
  }
  my $ehash = $pw->MolExId($molid);
  for my $id_type (keys %{ $ehash }) {
    for my $id (keys %{ $$ehash{$id_type} }) {
      if ($id_type eq "LL") {
        $locus = $id;
      } else {
        push @ids, $id;
      }
    }
  }
  print "<p>Basic identification:\n";
  print "<blockquote>\n";
  print "<p><table border=1 cellspacing=1 cellpadding=4>\n";
  print "<tr><td>Molecule Id:</td><td>$molid</td></tr>\n";
  print "<tr><td>Molecule type:</td><td>$moltype</td></tr>\n";
  print "<tr><td>Names:</td><td>" . join(", ", sort @names) . "</td></tr>\n";
  if ($locus) {
    print "<tr><td>Entrez Gene:</td><td><a href=\"" . LL_URL($locus) . "\"" .
        ">$locus</a></td></tr>\n";   
  }
  if (@ids) {
    print "<tr><td>Other ids:</td><td>" . join(", ", sort @ids) . "</td></tr>\n";
  }
  print "</table>\n";
  print "</blockquote>\n";
}

######################################################################
sub ComplexTable {
  my ($molid) = @_;

  my ($label, $value);
  my ($cx, $m, @temp);

  if (scalar(keys %complex_uses) == 0) {
    return;
  }

  print "<p>Uses in complexes:\n";
  print "<blockquote>\n";
  print "<p><table border=1 cellspacing=1 cellpadding=4>\n";
  print "<tr><td>Complex</td><td>Other Components</td></tr>\n";
  for $cx (keys %complex_uses) {
    undef @temp;
    print "<tr><td><a href=\"" .
      MOLPAGE_URL($cx) . "\">&lt;cx_$cx&gt;</a></td><td>";
    for my $m (keys %{ $complex_uses{$cx} }) {
      if ($m != $molid) {
        push @temp, ("<a href=\"" . MOLPAGE_URL($m) . "\">" .
         $pw->PickMolName($m) . "</a>");
      }
    }
    print join(", ", @temp);
    print "</td></tr>\n";
  }
  print "</table>\n";
  print "</blockquote>\n";

}

######################################################################
sub AtomTable {
  my ($molid) = @_;

  my ($label, $value);
  my ($atom, $edge);
  my (@pids, $pathwayid, $pcell, $atomtype);
  my ($location, $state, $role, $source);
  my %tmp;

  print "<p>Interactions:\n";
  print "<blockquote>\n";
  print "<p><table border=1 cellspacing=1 cellpadding=4>\n";
  print "<tr>" .
      "<td><b>Role</b></td>" .
      "<td><b>Location</b></td>" .
      "<td><b>State</b></td>" .
      "<td><b>Interaction</b></td>" .
      "<td><b>Pathway</b></td>" .
      "<td><b>Source</b></td>" .
      "</tr>\n";
  for $atom (keys %atom_uses) {
    ($label, $atomtype) =
        $lv->LabelValueToString($pw->AtomType($atom));
    @pids = keys %{ $atom2pathway{$atom} };
    if (@pids == 0) {
      $pids[0] = "";
    }
    for $edge (keys %{ $atom_uses{$atom} }) {
      for $pathwayid (@pids) {
        if ($pathwayid) {
          $pcell = 
            "<a href=\"" . PATHWAY_URL($pathwayid) .
            "\">$pathway2name{$pathwayid}</a>";
          $source = $pathway2source{$pathwayid};
        } else {
          $pcell = "&nbsp;";
          $source = $atom2source{$atom};
        }
        $role              = $atom_uses{$atom}{$edge};
        $location          = $location{$atom}{$edge};
        $location or $location = "&nbsp;";
        $state             = $state{$atom}{$edge};
        $state or $state   = "&nbsp;";
        $source or $source = "&nbsp;";
        push @{ $tmp{"$role|$location|$state|"} }, "<tr><td>" .
          join("</td>\n<td>",
            $role,
            $location,
            $state,
            "$atomtype (<a href=\"" . ATOM_URL($atom)    . "\">$atom</a>)",
            $pcell,
            $source
          ) . "</td></tr>\n";
      }
    }
  }
  for my $i (sort keys %tmp) {
    for (@{ $tmp{$i} }) {
      print $_;
    }
  }
  print "</table>\n";
  print "<p><a href=\"" . ATOM_URL(join(",", keys %atom_uses)) .
      "\">Display all interactions</a>";
  print "</blockquote>\n";
}

######################################################################
sub SimpleOrComplex {
  my ($molid) = @_;

  my ($atom, $edge, $m, $lvid, $label, $value, $edgetype);

  for $atom (@{ $pw->Atoms() }) {
    for $edge (@{ $pw->Edges($atom) }) {
      $m = $pw->EdgeMol($atom, $edge);
      if ($m == $molid) {
        $edgetype = $pw->EdgeType($atom, $edge);
        ($label, $value) = $lv->LabelValueToString($edgetype);
        $atom_uses{$atom}{$edge} = $value;
        for $lvid (@{ $pw->MolLabel($atom, $edge) }) {
          ($label, $value) = $lv->LabelValueToString($lvid);
          if ($label eq "location") {
            $location{$atom}{$edge} = $value;
          } elsif ($label eq "activity-state") {
            $state{$atom}{$edge} = $value;
          }
        }
      }
      if ($pw->MolType($m) == $COMPLEX_TYPE) {
        my %accum;
        LeafComponents($m, \%accum);
        if (defined $accum{$molid}) {
          for my $c (keys %accum) {
            $complex_uses{$m}{$c} = 1;
          }
        }
      }
    }
  }

}

######################################################################
sub LeafComponents {
  my ($molid, $accum) = @_;

  my ($seq, $m);

  if ($pw->MolType($molid) != $COMPLEX_TYPE) {
    return;
  }
  for $seq (@{ $pw->Components($molid) }) {
    $m = $pw->ComponentMol($molid, $seq);
    if (! defined $$accum{$m}) {
      $$accum{$m} = 1;
      if ($pw->MolType($m) == $COMPLEX_TYPE) {
        LeafComponents($m, $accum);
      }
    }
  }
}

