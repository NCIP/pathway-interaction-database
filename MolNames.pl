#!/usr/local/bin/perl

use strict;
use DBI;
use CGI;

if (-d "/app/oracle/product/dbhome/current") {
  $ENV{'ORACLE_HOME'} = "/app/oracle/product/dbhome/current";
} elsif (-d "/app/oracle/product/8.1.7") {
  $ENV{'ORACLE_HOME'} = "/app/oracle/product/8.1.7";
} elsif (-d "/app/oracle/product/8.1.6") {
  $ENV{'ORACLE_HOME'} = "/app/oracle/product/8.1.6";
}
  use constant NPG_MOL_START => 200000;

  my $query = new CGI;

  my (
    $db_inst,
    $db_user,
    $db_pass,
    $schema,
  ) = ("cgprod", "web", "readonly", "pid");

  my $source_ids;
  my $mol_types = "PR,CX,RN,CM";

  my ($source_id, $mol_id, $mol_type, $name_type, $name,
      $ext_id_type, $ext_id);
  my ($component_mol_id, $component_seq_id, $complex_mol_id);
  my ($db, $sql, $stm);
  my (%id2ext, %id2src, %id2type, %id2name, %mol_types, %source_ids,
      %component);

  print "Content-type: text/plain\n\n";

  $db = DBI->connect("DBI:Oracle:" . $db_inst, $db_user, $db_pass);
  if (not $db or $db->err()) {
    print STDERR "Cannot connect to " . $db_user . "@" . $db_inst . "\n";
    exit;
  }

  for my $t (split ",", $mol_types) {
    $t =~ s/\s+//g;
    $t = uc($t);
    if ($t ne "PR" && $t ne "CM" && $t ne "CX" && $t ne "RN") {
      die "Bad mol type $t";
    }
    $mol_types{$t} = 1;
  }
  $mol_types = "'" . join("','", keys %mol_types) . "'";

#  for my $t (split ",", $source_ids) {
#    $t =~ s/\s+//g;
#    if ($t !~ /^\d+$/){
#      die "Bad source id $t";
#    }
#    $source_ids{$t} = 1;
#  }
#  $source_ids = join(",", keys %source_ids);


######################################################################
## mol ids, basic types, and external ids
######################################################################

  $sql = qq!
select
  m.mol_id,
  m.basic_mol_type,
  e.id_type,
  e.ext_mol_id
from
  $schema.pw_mol m,
  $schema.pw_ext_mol_id e
where
      m.mol_id = e.mol_id (+)
  and m.basic_mol_type in ($mol_types)
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
  while (($mol_id, $mol_type, $ext_id_type, $ext_id) =
      $stm->fetchrow_array()) {
    $id2type{$mol_id} = $mol_type;
    $id2ext{$mol_id}{$ext_id_type}{$ext_id} = 1;
  }

######################################################################
## mol ids and their sources
######################################################################

#  $sql = qq!
#select
#  s.mol_id,
#  s.mol_source_id
#from
#  $schema.pw_mol_source s
#  !;
#
#  $stm = $db->prepare($sql);
#  if(not $stm) {
#    print STDERR "$sql\n";
#    print STDERR "$DBI::errstr\n";
#    print STDERR "prepare call failed\n";
#    exit;
#  }
#  if(!$stm->execute()) {
#    print STDERR "$sql\n";
#    print STDERR "$DBI::errstr\n";
#    print STDERR "execute call failed\n";
#    exit;
#  }
#  while (($mol_id, $source_id) = $stm->fetchrow_array()) {
#    $id2src{$mol_id}{$source_id} = 1;
#  }
#
######################################################################
## mol ids and their names
######################################################################

  $sql = qq!
select
  n.mol_id,
  n.name_type,
  n.mol_name
from
  $schema.pw_mol_name n
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
  while (($mol_id, $name_type, $name) = $stm->fetchrow_array()) {
    $id2name{$mol_id}{$name_type}{$name} = 1;
  }


######################################################################
## complex_components
######################################################################

  $sql = qq!
select unique
  c.complex_mol_id,
  c.component_mol_id,
  c.component_seq_id
from
  $schema.pw_complex_component c
order by
  c.complex_mol_id,
  c.component_seq_id
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
  while (($complex_mol_id, $component_mol_id, $component_seq_id) =
      $stm->fetchrow_array()) {
    push @{ $component{$complex_mol_id} }, $component_mol_id;
  }

  $db->disconnect();

######################################################################
## put it together
######################################################################

  for $mol_id (keys %id2type) {
    my (@ext_ids, @names, $ext_ids, $names, $components);
    $mol_type = $id2type{$mol_id};
#    if (defined $id2src{$mol_id}) {
#      my ($yes);
#      for (keys %{ $id2src{$mol_id} }) {
#        if (defined $source_ids{$_}) {
#          $yes++;
#        }
#      }
#      if (! $yes) {
#        next;
#      } else {
#        $source_id = join(";", keys %{ $id2src{$mol_id} });
#      }
#    } else {     ## if there is no entry in the source table
#                 ## print anyway; probably want to know about this
#      $source_id = "-";
#    }
    for $ext_id_type (keys %{ $id2ext{$mol_id} }) {
      for $ext_id (keys %{ $id2ext{$mol_id}{$ext_id_type} }) {
        push @ext_ids, "$ext_id_type $ext_id";
      }
    }
    $ext_ids = join(";", @ext_ids);
    for $name_type (keys %{ $id2name{$mol_id} }) {
      for $name (keys %{ $id2name{$mol_id}{$name_type} }) {
        push @names, "$name_type $name";
      }
    }
    $names = join(";", @names);
    if (defined $component{$mol_id}) {
	$components = join(";", @{ $component{$mol_id} });
    } else {
      $components = "";
    }
    if ($mol_id < NPG_MOL_START) {
      next;
    }
    print join("\t",
      $mol_id,
 #     $source_id,
      $mol_type,
      $ext_ids,
      $names,
      $components
    ) . "\n";
  }
