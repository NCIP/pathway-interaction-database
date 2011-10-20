#!/usr/local/bin/perl

BEGIN {
  my @path_elems = split("/", $0);
  pop @path_elems;
  push @INC, join("/", @path_elems);
}

use strict;
use Scan;
use Pathway;
use PathwayXML;
use PathwayDB;
use PathParams;
use PWLabel;
use PathwayGene;
use FileHandle;
use DOT;
use Math::Trig;
use DBI;

use constant MAX_LONG_LEN       => 16384;

my $BASE;

my ($COMPLEX_TYPE, $PROTEIN_TYPE, $COMPOUND_TYPE, $RNA_TYPE);
#my $LIMIT_SOURCES = "1,2";
my $LIMIT_SOURCES = "";
my $LIMIT_EVIDENCE = "";

my @CCF = (    ## common co-factors
    'water',
    'atp',
    'nad',
    'nadh',
    'nadph',
    'nadp',
    'oxygen',
    'adp',
    'orthophosphate',
    'coa',
    'carbon dioxide',
    'diphosphate',
    'ammonia',
    'h+'
  );

my ($db_inst, $db_user, $db_pass, $schema) =
    (DB_INST, DB_USER, DB_PASS, SCHEMA);

use constant MAX_DOT_INTERACTIONS => 800;

my $lv;
my %mol2atom;
my %atom2atom;
my %cx2cp;

## dot2xaml constants
my $tangent = tan(deg2rad(60));
my $canvas_n = 0;
my $molCount = 0;
my $keyCount = 1;

my $PIX_PER_INCH = 72;
my $ROUND_CORNER_CLIP = 12;
my $HALF_CORNER_CLIP  = 6;

my %TEXT_LINE_INC_Y = (
  "10" => 8
);
my %TEXT_LINE_HEAD_ROOM = (
  "10" => 13
);
my %TEXT_LINE_FOOT_ROOM = (
  "10" => 7
);

my %hex_color = (
  "black"     => "000000",
  "white"     => "ffffff",
  "red"       => "ff0000",
  "green"     => "00ff00",
  "blue"      => "0000ff",
  "firebrick" => "b22222",
  "lightgray" => "736F6E"
);


######################################################################
sub r_numerically { $b <=> $a };
sub   numerically { $a <=> $b };

######################################################################
sub SetSchema {
  my ($the_schema) = @_;
# as things stand, global $schema acquires value PWAppConfig::SCHEMA
# this sub is to override that value
  $schema = $the_schema;
}

######################################################################
sub InitializeLV {

  my $db = DBI->connect("DBI:Oracle:" . $db_inst, $db_user, $db_pass);
  if (not $db or $db->err()) {
    print STDERR "Cannot connect to " . $db_user . "@" . $db_inst . "\n";
    die;
  }
  $lv   = new PWLabel($db, $schema);
  $db->disconnect();

  $COMPLEX_TYPE  = $lv->StringToLabelValue("molecule-type", "complex");
  $PROTEIN_TYPE  = $lv->StringToLabelValue("molecule-type", "protein");
  $RNA_TYPE      = $lv->StringToLabelValue("molecule-type", "rna");
  $COMPOUND_TYPE = $lv->StringToLabelValue("molecule-type", "compound");

}

######################################################################
sub InitializeComplex2Component {

  my $db = DBI->connect("DBI:Oracle:" . $db_inst, $db_user, $db_pass);
  if (not $db or $db->err()) {
    print STDERR "Cannot connect to " . $db_user . "@" . $db_inst . "\n";
    die;
  }
  my ($sql, $stm);
  my ($cx_id, $cp_id);

  $sql = qq!
select
  c.complex_mol_id,
  c.component_mol_id
from
  $schema.pw_complex_component c
  !;
  $stm = $db->prepare($sql);
  if(not $stm) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "prepare call failed\n";
    die;
  }
  if(!$stm->execute()) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "execute call failed\n";
    die;
  }
  while (($cx_id, $cp_id) = $stm->fetchrow_array()) {
    $cx2cp{$cx_id}{$cp_id} = 1;
  }
  $stm->finish(); 
  $db->disconnect();
}

######################################################################
sub InitializeAtom2Atom {

  my $db = DBI->connect("DBI:Oracle:" . $db_inst, $db_user, $db_pass);
  if (not $db or $db->err()) {
    print STDERR "Cannot connect to " . $db_user . "@" . $db_inst . "\n";
    die;
  }
  my ($sql, $stm);
  my ($a1, $a2);

  $sql = "select * from $schema.pw_connect";

  $stm = $db->prepare($sql);
  if(not $stm) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "prepare call failed\n";
    die;
  }
  if(!$stm->execute()) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "execute call failed\n";
    die;
  }
  while (($a1, $a2) = $stm->fetchrow_array()) {
    $atom2atom{$a1}{$a2} = 1;
    $atom2atom{$a2}{$a1} = 1;
    $atom2atom{$a1}{$a1} = 1;
    $atom2atom{$a2}{$a2} = 1;
  }

#  my %atom2mol;
#  my $atom_sources  = "";
#  my $atom_evidence = "";
#  my %ccf;
#  CommonCoFactors($db, \%ccf);
#  my $prune_mols = join(",", keys %ccf);
#  SetUpAtomMatrix($db, $atom_sources, $prune_mols, \%atom2mol, \%mol2atom,
#      \%atom2atom, $atom_evidence);
  $stm->finish();
  $db->disconnect();
}

######################################################################
sub MoleculeName {
  my ($pw, $molid) = @_;
  my $name = $pw->PickMolName($molid);
#  if ($pw->MolType($molid) == $COMPLEX_TYPE) {
#    $name = "&lt;" . "cx_$molid" . "&gt;";
#  }
  return $name;
}

######################################################################
sub Value2Color {
  my ($v) = @_;

  my @color_scale = (
    "0000FF",    ## blue = low
    "3399FF",
    "66CCFF",
    "99CCFF",
    "CCCCFF",    ## skip this -- too close
    "",
    "FFCCFF",    ## skip this -- too close
    "FF99FF", 
    "FF66CC", 
    "FF6666",
    "FF0000"     ## red = high
  );
  my $LO_VAL = 1.0;
  my $HI_VAL = 9.0;
  my $NORMAL_VAL = ($LO_VAL + $HI_VAL) / 2;
  my $delta1 = .25;
  my $nlo_colors = 5;
  my $nhi_colors = 5;
  my $lo_delta = (($NORMAL_VAL - $delta1) - $LO_VAL) / $nlo_colors;
  my $hi_delta = ($HI_VAL - ($NORMAL_VAL + $delta1)) / $nhi_colors;
  my ($i);
  if (($v <= $NORMAL_VAL + $delta1) && ($v >= $NORMAL_VAL - $delta1)) {
    return "";
##
## in what follows, use 4 rather than 5 to skip middle colors
##
  } elsif ($v < $NORMAL_VAL) {
    for ($i = $nlo_colors; $i > 0 ; $i--) {
      if ($v < $NORMAL_VAL - $delta1 - (($i-1) * $lo_delta)) {
        return "#" . $color_scale[$nlo_colors - $i];
      }
    }
    return "#" . $color_scale[0];
  } else {
    for ($i = $nhi_colors; $i > 0 ; $i--) {
      if ($v > $NORMAL_VAL + $delta1 + (($i-1) * $hi_delta)) {
        return "#" . $color_scale[$#color_scale - $nhi_colors + $i];
      }
    }
    return "#" . $color_scale[$#color_scale];
  }
}

######################################################################
sub BackTrack {
  my ($mazeparent, $endpoint) = @_;

  my ($x, @path, %visited);

## NOTE: $endpoint is not part of returned path
## so null path means $endpoint is not in a walk

  $x = $endpoint;
  while (1) {
    ## should not be necessary, but ...
    if (defined $visited{$x}) {
      last;
    } else {
      $visited{$x} = 1;
    }
    if (! defined $$mazeparent{$x}) {
      last;
    } else {
      $x = $$mazeparent{$x};
      push @path, $x;
    }
  }
  return \@path;
}

######################################################################
sub SetUpAtomMatrix {
  my ($db, $atom_sources, $prune_mols, $atom2mol, $mol2atom,
      $atom2atom, $atom_evidence) = @_;

  my ($sql, $stm);
  my ($source_clause, $prune_clause, $evidence_clause);
  my ($atom_id, $mol_id);

  if ($prune_mols) {
    $prune_clause = " and e.mol_id not in ($prune_mols) ";
  } else {
    $prune_clause = "";
  }

  if ($atom_sources) {
    $source_clause = " and a.atom_source_id in ($atom_sources) ";
  } else {
    $source_clause = "";
  }

  if ($atom_evidence) {
    $evidence_clause = " and v.evidence_code in ($atom_evidence) ";
  } else {
    $evidence_clause = "";
  }

  $sql = qq!
select
  e.atom_id,
  e.mol_id
from
  $schema.pw_edge e,
  $schema.pw_atom a,
  $schema.pw_evidence v
where
      e.atom_id = a.atom_id
and
      e.atom_id = v.atom_id
  $prune_clause
  $source_clause
  $evidence_clause
  !;

  $stm = $db->prepare($sql);
  if(not $stm) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "prepare call failed\n";
    die;
  }
  if(!$stm->execute()) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "execute call failed\n";
    die;
  }
  while (($atom_id, $mol_id) = $stm->fetchrow_array()) {
    $$atom2mol{$atom_id}{$mol_id} = 1;
    $$mol2atom{$mol_id}{$atom_id} = 1;
  }
  $stm->finish();
  for my $a1 (keys %{ $atom2mol}) {
    for my $m (keys %{ $$atom2mol{$a1} }) {
      for my $a2 (keys %{ $$mol2atom{$m} }) {
        $$atom2atom{$a1}{$a2} = 1;
        $$atom2atom{$a2}{$a1} = 1;
      }
    }
  }

}

######################################################################
sub WalkMaze {
  my ($atom2atom, $src, $target, $connection) = @_;
  my (%mazeparent);

#print STDERR "WalkMaze: src = " . join(",", keys %{ $src }) . "; " .
#  "target = " . join(",", keys %{ $target }) . "\n";

  my (%new, %family, $a, $b, $found);

  for $a (keys %{ $src }) {
    $new{$a} = 1;
  }
  while (1) {
    undef %family;
    for $a (keys %new) {
      $family{$a} = 1;
    }
    undef %new;
    for $a (keys %family) {
      for $b (keys %{ $$atom2atom{$a} }) {
        if (defined $$target{$b}) {
          $found = $b;
        }
        if ( !defined $mazeparent{$b} &&
             !defined $family{$b} ) {      ## no incest
          $mazeparent{$b} = $a;
          if (! defined $$target{$b}) {
            $new{$b} = 1;
          }
        }
        if ($found) {
          last;
        }
      }
      if ($found) {
        last;
      }
    }
    if ($found || scalar(keys %new) == 0) {
      last;
    }
  }
  if ($found) {
    $$connection{$found} = 1;
    for $a (@{ BackTrack(\%mazeparent, $found) }) {
      $$connection{$a} = 1;
    }
    return 1;
  } else {
    return 0;
  }
}

######################################################################
sub CommonCoFactors {
  my ($db, $ccf) = @_;

  my ($sql, $stm);
  my ($mol_id);
  my $namelist = "'" . join("','", @CCF) . "'";

  $sql = qq!
select
  m.mol_id
from
  $schema.pw_mol m,
  $schema.pw_mol_name n
where
      m.mol_id = n.mol_id
  and m.basic_mol_type = 'CM'
  and n.mol_name in ($namelist)
  !;
  $stm = $db->prepare($sql);
  if(not $stm) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "prepare call failed\n";
    die;
  }
  if(!$stm->execute()) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "execute call failed\n";
    die;
  }
  while (($mol_id) = $stm->fetchrow_array()) {
    $$ccf{$mol_id} = 1;
  }
  $stm->finish();
}


######################################################################
sub FindConnectingPath {
  my ($db, $mols_to_join, $source, $molmap) = @_;

  my (%src, %target);  ## not really directional
  my (%connection);
  my ($status);
  my ($sql, $stm, $atom_id, $mol_id, $mol_type);

####!!!! %mol2atom, %atom2atom are globals
####!!!! %cx2cp is a global

  undef %mol2atom;
  if (defined $source && keys %{ $source } > 0) {
    $sql = "select unique e.atom_id, e.mol_id " .
        "from $schema.pw_edge e, " .
        "$schema.pw_atom a " .
        "where a.atom_id = e.atom_id " .
        "and a.atom_source_id in (" .
        join(",", keys %{ $source }) .
        ")";
  } else {
    $sql = "select unique e.atom_id, e.mol_id " .
        "from $schema.pw_edge e ";
  }

  $stm = $db->prepare($sql);
  if(not $stm) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "prepare call failed\n";
    die;
  }
  if(!$stm->execute()) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "execute call failed\n";
    die;
  }
  while (($atom_id, $mol_id) = $stm->fetchrow_array()) {
    $mol2atom{$mol_id}{$atom_id} = 1;
    if (defined $cx2cp{$mol_id}) {
      for my $cp_id (keys %{ $cx2cp{$mol_id} }) {
        $mol2atom{$cp_id}{$atom_id} = 1;
      }
    }
  }
  $stm->finish();
  ## Take only one mol that matched a query:
  ## Take the one with the lowest (numerical) id that hits the
  ## most atoms

  my @temp;
  for my $q (keys %{ $molmap }) {
    my $take_me;
    for my $mol (sort numerically keys %{ $$molmap{$q} }) {
      if (defined $mol2atom{$mol}) {
        if (! defined $take_me) {
          $take_me = $mol;
        } elsif (scalar(keys %{ $mol2atom{$take_me} }) <
                 scalar(keys %{ $mol2atom{$mol}     })) {
          $take_me = $mol;
        }
      }
    }
    if (defined $take_me) {
      push @temp, $take_me;
    }

  }

# OLD: take all mols for queried name 
#  for my $m (sort numerically keys %{ $mols_to_join }) {
#    if (defined $mol2atom{$m}) {
#      push @temp, $m;
#    }
#  }

  if (@temp < 2) {
    for my $a (sort numerically keys %{ $mol2atom{$temp[0]} }) {
      $connection{$a} = 1;
      last;
    }
    return \%connection;
  }

  ## do it pairwise

  ## src is populated with $mol2atom{$temp[$i]} only if there
  ## is not yet anything in the connection; otherwise, set src = connection

  my %hits;
  for (my $i = 0; $i < @temp - 1; $i++) {
    undef %src;
    if (keys %connection == 0) {
      for my $atom (keys %{ $mol2atom{$temp[$i]} }) {
        $src{$atom} = 1;
      }
    }
#    } else {
#      for my $atom (keys %connection) {
#        $src{$atom} = 1;
#      }
#    }
    for (my $j = $i + 1; $j < @temp; $j++) {
#      if (keys %connection != 0) {
#        undef %src;
#        for my $a (keys %connection) {
#          $src{$a} = 1;
#        }
#      }
      undef %target;
      my $already_connected;
      for my $a (keys %{ $mol2atom{$temp[$j]} }) {
        if (defined $connection{$a}) {
          $already_connected++;
          last;
        }
        $target{$a} = 1;
      }
      if ($already_connected) {
        next;
      }
      for my $a (keys %{ $mol2atom{$temp[$j]} }) {
        if (keys %connection == 0 && defined $src{$a}) {
          $connection{$a} = 1;
          $already_connected++;
          last;
        }
      }

      ## if any atom for this mol is already in the connection
      ## move on to the next mol
      if ($already_connected) {
        next;
      }

      if (keys %connection == 0) {
        $status = WalkMaze(\%atom2atom, \%src, \%target, \%connection);
      } else {
        $status = WalkMaze(\%atom2atom, \%connection, \%target, \%connection);
      }
      for my $m (@temp) {
        for my $a (keys %{ $mol2atom{$m} }) {
          if (defined $connection{$a}) {
            $hits{$m} = 1;
            last;
          }
        }
      }
      if (keys %hits == @temp) {
        return \%connection;
      }
    }
    for my $m (@temp) {
      for my $a (keys %{ $mol2atom{$m} }) {
        if (defined $connection{$a}) {
          $hits{$m} = 1;
          last;
        }
      }
    }
    if (keys %hits == @temp) {
      return \%connection;
    }
  }

  ## for things that do not connect, just pick one atom

  for (my $i = 0; $i < @temp; $i++) {
    if (! defined $hits{$temp[$i]}) {
      for my $a (sort keys %{ $mol2atom{$temp[$i]} }) {
        $connection{$a} = 1;
        last;
      }
    }
  }

  return \%connection;

}


######################################################################
#
sub TimeStamp {
  my ($sec, $min, $hr, $mday, $mon, $year, $wday, $yday, $isdst) =
      localtime(time);
  my $month = ('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug',
      'Sep', 'Oct', 'Nov', 'Dec')[$mon];
  $year = $year + 1900;
  return sprintf "%d_%2.2d_%2.2d_%2.2d_%2.2d_%2.2d",
      $year, $mon+1, $mday, $hr, $min, $sec;
}


######################################################################
sub Comment_1 {
  my ($base, $name, $email, $organization,
      $atomid, $molname, $molid, $pathname, $pathid, $comment) = @_;

  my @temp;
  push @temp, TimeStamp();

  if ( ! open(LOG, ">>" . PW_COMMENT_LOG) ) {
    SetStatus(S_RESPONSE_FAIL);
    return "Cannot access comment log";
  }
  for my $x (
      $name,
      $email,
      $organization,
      $atomid,
      $molname,
      $molid,
      $pathname,
      $pathid,
      $comment
  ) {
    $x =~ s/^\s+//;
    $x =~ s/\s+$//;
    $x =~ s/\t+/ ^TAB^ /;
    $x =~ s/\n+/ ^LF^ /;
    push @temp, $x;
  }
  print LOG join("\t", @temp) . "\n";
  close LOG;
  chmod 0666, PW_COMMENT_LOG;
  return "Comment submitted\n";
}

######################################################################
sub SearchMoleculeKeywords_1 {
  my ($base, $geneid) = @_;
  my (@lines);

  $BASE = $base;
  my ($db_inst, $db_user, $db_pass, $pid_schema, $cgap_schema) =
	(DB_INST, DB_USER, DB_PASS, SCHEMA, "cgap");
 
  $geneid =~ s/\s+//g;
  if ($geneid !~ /^\d+$/) {
	return "Illegal value for parameter geneid";
  }

  my $db = DBI->connect("DBI:Oracle:" . $db_inst, $db_user, $db_pass);
  if (not $db or $db->err()) {
    return "Cannot connect to " . $db_user . "@" . $db_inst . "\n";
  }

  my ($lv, $pw);

  my ($mol2name, $mol2ext_id);

  my $mols = FindMolIdsForGene($db, $geneid, $pid_schema, $cgap_schema);
  
  if (keys %{ $mols } > 0) {
    ($mol2name, $mol2ext_id) = FindNamesForMols($db, $mols, $pid_schema);
  } else {
    return "No entries for Gene id = $geneid\n";
  }


  push @lines, FormatOutput($db, $mols, $mol2name, $mol2ext_id, $pid_schema, $geneid);

  $db->disconnect;  
  if (@lines) {
    return join("\n", @lines, "");
  } else {
    return "";
  }

}

sub FindMolIdsForGene {
  my ($db, $geneid, $pid_schema, $cgap_schema) = @_;

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

sub FormatOutput {
  my ($db, $mols, $mol2name, $mol2ext_id, $pid_schema, $gene_id) = @_;

  my ($sql, $stm, $source );
  
  my @lines;
  push @lines, "Here are molecules corresponding to Entrez Gene ID $gene_id <P>";
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

   $sql = qq!
select
  s.source_name
from
  $pid_schema.pw_source s,$pid_schema.pw_pathway p,
  $pid_schema.pw_pathway_atom pa, $pid_schema.pw_edge e
where
     (e.mol_id = $mol 
     or e.mol_id in (select family_mol_id from $pid_schema.pw_family_member  
     where member_mol_id = $mol))
     and pa.atom_id = e.atom_id
     and p.pathway_id = pa.pathway_id
     and p.pathway_source_id = s.source_id
!;

  $stm = $db->prepare($sql);
  if(not $stm) {
    return "prepare failed: $sql\n";
  }
  if(!$stm->execute()) {
    print STDERR "execute failed: $sql\n";
    return;
  }
  while (( my $s ) = $stm->fetchrow_array()) {

    $s =~ s/NATURE/NCI-Nature curated/;
    $s =~ s/BioCarta/BioCarta Imported/;

    $source = $s;
  }

    for my $nametype (keys %{ $$mol2name{$mol} }) {
      for my $name (keys %{ $$mol2name{$mol}{$nametype} }) {
        push @names, $name;
      }
    }

    if (($source eq 'NCI-Nature curated') || ($source eq 'BioCarta Imported')) {
      push @lines, "<a href=\"" . MOLPAGE_URL($mol) . "\">" .
        join(", ", @names, $source) . "</a><BR> ";
    }
  }
  return join("<br>\n", @lines, "");

}
#############################################################
sub FindNamesForMols {
  my ($db, $mols, $pid_schema) = @_;

  my $mollist = join(",", keys %{ $mols });
  my (%mol2name, %mol2ext_id);

  my ($sql, $stm);
  my ($molid, $extid, $id_type, $name_type, $name, $source);

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
sub SearchPathwayKeywords_1 {
  my ($base, $word, $format) = @_;

  $BASE = $base;

  my ($sql, $stm);
  my ($source, $pname, $ext_id, $pid);
  my (@lines);

  $word =~ s/^\s+//;
  $word =~ s/\s+$//;
  $word =~ s/\s+/ /g;
  $word = lc($word);

  if ($format eq "html") {
    return ListAllPathways_1($base, $word);
  }
  
  my $db = DBI->connect("DBI:Oracle:" . $db_inst,
        $db_user, $db_pass);
  if (not $db or $db->err()) {
    print STDERR "#! Cannot connect to " . $db_user . "@" .
        $db_inst . "\n";
    die;
  }

  $sql = "select " .
      "s.source_name, p.pathway_name, p.ext_pathway_id, p.pathway_id " .
      "from $schema.pw_source s, $schema.pw_pathway p " .
      "where s.source_id = p.pathway_source_id " .
      "order by s.source_name, p.pathway_name";

  $stm = $db->prepare($sql);
  if(not $stm) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "prepare call failed\n";
    die;
  }
  if(!$stm->execute()) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "execute call failed\n";
    die;
  }
  while (($source, $pname, $ext_id, $pid) =
      $stm->fetchrow_array()) {
    if ((index(lc($pname), $word) > -1) ||
       (index(lc($ext_id), $word) > -1)) {
      $source =~ s/NATURE/NCI-Nature curated/;
      $source =~ s/BioCarta/BioCarta Imported/;
      if ($format eq "text") {
        push @lines, join("\t", $source, $pname, $ext_id, $pid);
      } elsif ($format eq "html") {
      }
    }
  }
  $stm->finish();
  $db->disconnect();
  if (@lines) {
    return join("\n", @lines, "");
  } else {
    return "";
  }
}

###############################################################
sub NetworkHeader_1 {
  my ($base, $moleculeName, $source_id, $format) = @_;
  
  my @lines;
  my $source;
  my $macroprocess;
  my $pid;
 
   my $scancheck = Scan($source_id) + Scan($format) + Scan($moleculeName);
  if ($scancheck > 0) {
    push @lines, "Bad data in parameters\n";
    return join("\n", @lines) . "\n";

  }
 
  push @lines, "<ul class=\"xoxo viewing-option\">";
 
  push @lines, "<li><ul class=\"links\">";

 
  push @lines, "<li><a class=\"button-style\" href=\"/search/network_mol_list_landing.shtml" .
         "?source=$source_id" . "&molecule=$moleculeName\">Molecule list</a></li>";

  push @lines, "</ul></li>"; 
  push @lines, "<li><span class=\"norm\">View graphic as:</span>";
  push @lines, "<ul class=\"links\">";
  push @lines, "<li><a class=\"button-style\" href=\"/search/network_landing.shtml" .
          "?molecule=$moleculeName" . "&macro_process=$macroprocess" . "&source_id=$source_id" . "&what=graphic&jpg=on&ppage=1\">JPG</a></li>" ;

  push @lines, "<li><a class=\"button-style\" href=\"/search/network_landing.shtml" .
       "?molecule=$moleculeName" . "&macro_process=$macroprocess" . "&source_id=$source_id" . "&what=graphic&svg=on&ppage=1\">SVG</a></li>";
  
  push @lines, "<li><a class=\"button-style\" href=\"/xaml_net_landing.shtml" .
       "?molecule=$moleculeName" . "&macro_process=$macroprocess" . "&source_id=$source_id" . "&source=NATURE" . "&what=graphic&xaml=on\">Silverlight-beta</a></li>";

  push @lines, "<li><div class=\"helpbox\">";
  push @lines, "<a href=\"#\" class=\"help\">Graphics help</a>";
  push @lines, "<div class=\"description\">";
  push @lines, "<h2 class=\"title\">Information on graphics types</h2>";
  push @lines, "<p class=\"norm\">The pathway can be viewed as a JPG, SVG or Silverlight image. All molecules and interactions within the graphics formats are interactive; clicking on them will open information windows.</p>";
  push @lines, "<p class=\"norm\"><B>JPG:</b>  All common web browsers are able to display JPG images.  However, please note that some browsers may not support very large JPG files.</p>";
  push @lines, "<p class=\"norm\"><B>SVG:</b>  The advantage of SVG images is that they can be more easily zoomed and panned on-screen.  If your browser is not completely SVG enabled, please visit the <a href=\"/PID/userguide/output_formats.shtml\">User guide</a> and follow the instructions for installing the SVG plug-in specific to your browser</p>";
  push @lines, "<p class=\"norm\"><b>SVG image navigation:</b>  Use the left mouse button to activate the image.  To <b>zoom</b> in or out, use the right mouse button and select from the options on the pop-up menu.  To <b>pan</b> the image, hold down the alt key and click-and-drag the image with the mouse.</p>";
  push @lines, "<p class=\"norm\"><b>Silverlight:  </b>Please visit http://www.microsoft.com/getsilverlight/Get-Started/Install/Default.aspx to install the plug-in for your browser.</p>";
  push @lines, "<p class=\"norm\"><b>Silverlight image navigation:  Zoom</b> is available at the top-left of the image. To <b>pan</b> merely drag the image with the mouse.  In the image you may search for specific molecules or all molecules at a certain subcellular location.  To search for molecules you may enter their names, Entrez Gene identifiers or UniProt accession numbers.</p>";

  push @lines, "</div>";
  push @lines, "</div>";
  push @lines, "</li>";
  push @lines, "</ul>";
  push @lines, "</li>";

  push @lines, "<li><span class=\"norm\">Save code as:</span>";
  push @lines, "<ul class=\"links\">";

   push @lines, "<li><a class=\"button-style\" href=\"/search/network_landing.shtml" .
          "?molecule=$moleculeName" . "&macro_process=$macroprocess" . "&source_id=$source_id" . "&what=text&xml=on&ppage=1\">XML</a></li>" ;
  push @lines, "<li><a class=\"button-style\" href=\"/search/network_landing.shtml" .
"?molecule=$moleculeName" . "&macro_process=$macroprocess" . "&source_id=$source_id"."&what=text&biopax=on&ppage=1\">BioPAX</a></li>";
  push @lines, "<li><div class=\"helpbox\">";
  push @lines, "<a href=\"#\" class=\"help\">Code help</a>";
  push @lines, "<div class=\"description\">";
  push @lines, "<h2 class=\"title\">Information on types of pathway code</h2>";
  push @lines, "<p class=\"norm\">You can view the pathway code in either XML or BioPax formats.  Please see the <a href=\"/PID/userguide/output_formats.shtml\">User guide</a> for more information.</p>";
  push @lines, "</div>";
  push @lines, "</div>";
  push @lines, "</li>";

  push @lines, "</ul>";
  push @lines, "</li>";
 
  push @lines, "</ul>";
  if ($format eq 'svg') { 
    push @lines, "<p class=\"norm\">A <a href=\"/search/network_landing.shtml" .
          "?molecule=$moleculeName" . "&macro_process=$macroprocess" . "&source_id=$source_id" . "&what=graphic&svg=on/#key\">key to the image</a> icons appears below the pathway.</p>";
  } else {
    push @lines, "<p class=\"norm\">A <a href=\"/search/network_landing.shtml" . 
      "?molecule=$moleculeName" . "&macro_process=$macroprocess" . "&source_id=$source_id" . "&what=graphic&jpg=on/#key\">key to the image</a> icons appears below the pathway.</p>"; 
  }
  push @lines, "<p class=\"norm\">Click on a biomolecule or interaction for additional information.</p>";

  return join("\n", @lines, "");
}

###############################################################
sub BatchHeader_1 {
  my ($base, $pid, $genes_a, $genes_b, $format) = @_;

  $BASE = $base;

  my @ALPHA_MONTHS =
      ("Jan", "Feb", "Mar", "Apr",
       "May", "Jun", "Jul", "Aug",
       "Sep", "Oct", "Nov", "Dec");

  my (@lines);

  my $db = DBI->connect("DBI:Oracle:" . $db_inst,
        $db_user, $db_pass);
  if (not $db or $db->err()) {
    print STDERR "#! Cannot connect to " . $db_user . "@" .
        $db_inst . "\n";
    die;
  }

  my ($sql, $stm, $sqlc, $stmc);
  my ($source, $pname, $ext_id, $curator, $role, $last_updated);
  my (@curators, @reviewers);
  my ($biocarta_url,$parentName,$parentPid);

  $parentName = 'None';
  $sql = "select " .
      "s.source_name, p.pathway_name, p.ext_pathway_id, p.pathway_id, " .
      "p.last_updated " .
      "from $schema.pw_source s, $schema.pw_pathway p " .
      "where p.pathway_id = $pid and s.source_id = p.pathway_source_id";

  $stm = $db->prepare($sql);
  if(not $stm) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "prepare call failed\n";
    die;
  }
  if(!$stm->execute()) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "execute call failed\n";
    die;
  }
  ($source, $pname, $ext_id, $pid, $last_updated) = $stm->fetchrow_array();
  $stm->finish();

  my ($sec, $min, $hr, $mday, $mon, $year, $wday, $yday, $isdst) =
      localtime($last_updated);
  my $mon_alpha = $ALPHA_MONTHS[$mon];
  $year = $year + 1900;

  my $nice_source = $source;

  if ($source eq "BioCarta") {
    $sql = "select " .
      "b.bc_ext_id from pid.pw_biocarta b " .
      "where b.pid_ext_id = '$ext_id'";

    $stm = $db->prepare($sql);
    if(not $stm) {
      print STDERR "$sql\n";
      print STDERR "$DBI::errstr\n";
      print STDERR "prepare call failed\n";
      die;
    }
    if(!$stm->execute()) {
      print STDERR "$sql\n";
      print STDERR "$DBI::errstr\n";
      print STDERR "execute call failed\n";
      die;
    }
    ($biocarta_url) = $stm->fetchrow_array();
    $biocarta_url = $biocarta_url . ".asp";
    $nice_source = "<a href=\"http://www.biocarta.com/pathfiles/$biocarta_url\">BioCarta Imported</a>";
    $stm->finish();

  }
  $nice_source =~ s/NATURE/NCI-Nature curated/;
  # $nice_source =~ s/BioCarta/<a href="http:\/\/www.biocarta.com/pathfiles/$biocarta_url">BioCarta Imported<\/a>/;
  $nice_source =~ s/Reactome/Reactome Imported/;

  undef @curators;
  undef @reviewers;

  my $sqlsub = qq!

select distinct
  p2.pathway_id,
  p2.ext_pathway_id,
  p2.pathway_name
from
  $schema.pw_pathway p1,
  $schema.pw_pathway p2,
  $schema.pw_abstraction a,
  $schema.pw_pathway_atom pa,
  $schema.pw_subnet s
where
  pa.pathway_id = p2.pathway_id
  and pa.atom_id = a.atom_id
  and p1.pathway_source_id = p2.pathway_source_id
  and (a.pathway_id = p1.pathway_id or a.ext_pathway_id = p1.ext_pathway_id)
  and p1.pathway_id = $pid
  and p1.subnet = 'Y'
  and s.subnet_pathway_id = p1.pathway_id
  and p2.pathway_id = s.parent_pathway_id
!;
  my $stmsub = $db->prepare($sqlsub);
  if (not $stmsub) {
    print STDERR "$sqlsub\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "prepare call failed\n";
    die;
  }
  if (!$stmsub->execute()) {
    print STDERR "$sqlsub\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "execute call failed\n";
    die;
  }
  my $pid1;
  my $pid2;
  my $pidname;
  while (($pid1, $pid2, $pidname) = $stmsub->fetchrow_array()) {
    $parentName = $pidname;
    $parentPid = $pid1;
  }
  $stmsub->finish();

  if (($source eq "NATURE") || ($source eq "Reactome")) {
    $sqlc = "select curator, role " .
            "from $schema.pw_curators " .
            "where pathway_id = $pid " ;
    $stmc = $db->prepare($sqlc);
    if (not $stmc) {
      print STDERR "$sqlc\n";
      print STDERR "$DBI::errstr\n";
      print STDERR "prepare call failed\n";
      die;
    }
    if (!$stmc->execute()) {
      print STDERR "$sqlc\n";
      print STDERR "$DBI::errstr\n";
      print STDERR "execute call failed\n";
      die;
    }
    while (($curator, $role) = $stmc->fetchrow_array()) {
      if ($role eq "C") {
        push @curators, $curator;
      } else {
        push @reviewers, $curator;
      }
    }
    $stmc->finish();

  }

  my $other_format_line;
  ## for some reason, can't seem to get the "format"  variable
  ## through the sequence of invocations, so for now, just list both
  ## graphic formats

  push @lines, "<h1 class=\"pagetitle\">$pname</h1>";
  push @lines, "<hr />";

  push @lines, "<p class=\"citation\"><span class=\"title\">Revision date:</span> $mday" . "-" . $mon_alpha .
      "-" . $year ."</p>";
  if (($source eq "NATURE") || ($source eq "Reactome")) {
    push @lines, "<p class=\"citation\"><span class=\"title\">Curated by:</span> " . join(", ", @curators) . "</p>";
    push @lines, "<p class=\"citation\"><span class=\"title\">Reviewed by:</span> " . join(", ", @reviewers) . "</p>";
}
    push @lines, "<p class=\"citation\"><span class=\"title\">Pathway ID:</span> $ext_id </p>";
  if ($parentName ne 'None') {
    if ($format eq 'jpg') {
      push @lines, "<p class=\"citation\"><span class=\"title\">Parent pathway:  </span>  <a href=\"/search/batch_landing.shtml?pathway_id=$parentPid&source=$source&what=graphic&jpg=on&ppage=1\">$parentName</a></p><br>";
    }
    if ($format eq 'svg') {
      push @lines, "<p class=\"citation\"><span class=\"title\">Parent pathway:  </span><a href=\"/search/batch_landing.shtml?pathway_id=$parentPid&source=$source&what=graphic&svg=on&ppage=1\">$parentName</a></p>";
      push @lines, "<p><br>";
    }
  }
  if ($source eq "NATURE") {
    push @lines, "<p class=\"feedback\"><a href=\"mailto:nci-pid\@nature.com\">Submit feedback for this pathway</a></p>";
  }
  push @lines, "<ul class=\"xoxo viewing-option\">";
  push @lines, "<li><span class=\"hidden\">Additional citation information</span>";
  push @lines, "<ul class=\"links\">";
  push @lines, "<li><a class=\"button-style\" href=\"/search/mol_list_landing.shtml" .
          "?batch_id=$pid" .  "&what=text\">Molecule list</a></li>";
  push @lines, "<li><a class=\"button-style\" href=\"/search/citation_landing.shtml" .
          "?batch_id=$pid" .  "&what=text\">References</a></li>";
  push @lines, "</ul></li>";

  push @lines, "<ul class=\"xoxo viewing-option\">";
  push @lines, "<li><span class=\"norm\">View graphic as:</span>";
  push @lines, "<ul class=\"links\">";
  push @lines, "<li><a class=\"button-style\" href=\"/search/batch_landing.shtml" .
          "?pathway_id=$pid" . "&source=$source" . "&what=graphic&jpg=on&ppage=1&genes_a=$genes_a&genes_b=$genes_b\">JPG</a></li>" ;

  push @lines, "<li><a class=\"button-style\" href=\"/search/batch_landing.shtml" .
       "?pathway_id=$pid" . "&source=$source" . "&what=graphic&svg=on&ppage=1&genes_a=$genes_a&genes_b=$genes_b\">SVG</a></li>";

   push @lines, "<li><a class=\"button-style\" href=\"/search/xaml_batch_landing.shtml" .
          "?pathway_id=$pid" . "&source=$source" . "&what=graphic&xaml=on&genes_a=$genes_a&genes_b=$genes_b\">Silverlight-beta</a></li>" ;

  push @lines, "<li><div class=\"helpbox\">";
  push @lines, "<a href=\"#\" class=\"help\">Help</a>";
  push @lines, "<div class=\"description\">";
  push @lines, "<h2 class=\"title\">Information on graphics types</h2>";
  push @lines, "<p class=\"norm\">The pathway can be viewed as a JPG, SVG or Silverlight image. All molecules and interactions within the graphics formats are interactive; clicking on them will open information windows.</p>";
  push @lines, "<p class=\"norm\"><B>JPG:</b>  All common web browsers are able to display JPG images.  However, please note that some browsers may not support very large JPG files.</p>";
  push @lines, "<p class=\"norm\"><B>SVG:</b>  The advantage of SVG images is that they can be more easily zoomed and panned on-screen.  If your browser is not completely SVG enabled, please visit the <a href=\"/PID/userguide/output_formats.shtml\">User guide</a> and follow the instructions for installing the SVG plug-in specific to your browser</p>";
  push @lines, "<p class=\"norm\"><b>SVG image navigation:</b>  Use the left mouse button to activate the image.  To <b>zoom</b> in or out, use the right mouse button and select from the options on the pop-up menu.  To <b>pan</b> the image, hold down the alt key and click-and-drag the image with the mouse.</p>";
  push @lines, "<p class=\"norm\"><b>Silverlight:  </b>Please visit http://www.microsoft.com/getsilverlight/Get-Started/Install/Default.aspx to install the plug-in for your browser.</p>";
  push @lines, "<p class=\"norm\"><b>Silverlight image navigation:  Zoom</b> is available at the top-left of the image. To <b>pan</b> merely drag the image with the mouse.  In the image you may search for specific molecules or all molecules at a certain subcellular location.  To search for molecules you may enter their names, Entrez Gene identifiers or UniProt accession numbers.</p>";
  push @lines, "</div>";
  push @lines, "</div>";
  push @lines, "</li>";
  push @lines, "</ul>";
  push @lines, "</li>";

  push @lines, "<li><span class=\"norm\">Save code as:</span>";
  push @lines, "<ul class=\"links\">";
  push @lines, "<li><a class=\"button-style\" href=\"/search/batch_landing.shtml" .
          "?pathway_id=$pid" . "&source=$source" . "&what=text&xml=on&ppage=1&genes_a=$genes_a&genes_b=$genes_b\">XML</a></li>" ;
  push @lines, "<li><a class=\"button-style\" href=\"/search/batch_landing.shtml" .
       "?pathway_id=$pid" . "&source=$source" . "&what=text&biopax=on&ppage=1&genes_a=$genes_a&genes_b=$genes_b\">BioPAX</a></li>";
  push @lines, "<li><div class=\"helpbox\">";
  push @lines, "<a href=\"#\" class=\"help\">Help</a>";
  push @lines, "<div class=\"description\">";

  push @lines, "<h2 class=\"title\">Information on types of pathway code</h2>";
  push @lines, "<p class=\"norm\">You can view the pathway code in either XML or BioPax formats.  Please see the <a href=\"/PID/userguide/output_formats.shtml\">User guide</a> for more information.</p>";
  push @lines, "</div>";
  push @lines, "</div>";
  push @lines, "</li>";
  push @lines, "</ul>";
  push @lines, "</li>";
  push @lines, "</ul>";

  if ($format eq 'jpg') {
    push @lines, "<p class=\"norm\">A <a href=\"/search/batch_landing.shtml?pathway_id=$pid" . "&source=$source" . "&jpg=true&what=graphic&ppage=1&genes_a=$genes_a&genes_b=$genes_b/#key\">key to the image</a> icons appears below the pathway.</p>";
  } else {
    push @lines, "<p class=\"norm\">A <a href=\"/search/batch_landing.shtml?pathway_id=$pid" . "&source=$source" . "&svg=true&what=graphic&ppage=1&genes_a=$genes_a&genes_b=$genes_b/#key\">key to the image</a> icons appears below the pathway.</p>";
  }
  push @lines, "<p class=\"norm\">Click on a biomolecule or interaction for additional information.</p>";
  if ($format eq "svg") {
    push @lines, "<p class=\"norm\"><a href=/userguide/network_maps.shtml#>SVG Image Navigation Help</a></p>";
  }

  my %genes_a;
  my %genes_b;

  my @genes_1 = split(/,/,$genes_a);
  my @genes_2 = split(/,/,$genes_b);
  for (my $i = 0; $i < $genes_a; $i++) {
    if ($genes_1[$i] ne '') {
      $genes_a{$genes_1[$i]} = 1;
    }
  }

  for (my $i = 0; $i < $genes_b; $i++) {
    if ($genes_2[$i] ne '') {
      $genes_b{$genes_2[$i]} = 1;
    }
  }

  my $mols_a = GetMolIdsForEntrezGene($db, $pid, \%genes_a);
  my $mols_b = GetMolIdsForEntrezGene($db, $pid, \%genes_b);
  push @lines, "<div class=\"container-cite\">";
  push @lines, "<p class=\"norm\">Molecules from Group 1 are shown in <span class=\"blue\">blue</span></p>";
  push @lines, "<p class=\"norm\">Molecules from Group 2 are shown in <span class=\"red\">red</span></p>";
  push @lines, "<p class=\"norm\">Molecules from both groups are shown in <span class=\"purple\">purple</span></p>";
  push @lines, "</div>";
 
  print STDERR "PID before CreateGraphic:  $pid\n"; 
  push @lines, CreateGraphicFile($db, $schema, $pid, $source, $mols_a, $mols_b, $format);
  $db->disconnect();
  return join("\n", @lines) . "\n";

}

###############################################################
sub ConnectHeader_1 {
  my ($base, $moleculeName, $source_id, $format) = @_;

  my @lines;
  my $source;
  my $macroprocess;
  my $pid;

   my $scancheck = Scan($base) + Scan($source_id) + Scan($format) + Scan($moleculeName);
  if ($scancheck > 0) {
    push @lines, "Bad data in parameters\n";
    return join("\n", @lines) . "\n";

  }

  push @lines, "<ul class=\"xoxo viewing-option\">";

  push @lines, "<li><ul class=\"links\">";


  push @lines, "<li><a class=\"button-style\" href=\"/search/network_mol_list_landing.shtml" .
         "?source=$source_id" . "&molecule=$moleculeName\">Molecule list</a></li>";

  push @lines, "</ul></li>";
push @lines, "<li><span class=\"norm\">View graphic as:</span>";
  push @lines, "<ul class=\"links\">";
  push @lines, "<li><a class=\"button-style\" href=\"/search/connectedmolecules_landing.shtml" .
          "?connect_molecule=$moleculeName" . "&source_id=$source_id" . "&what=graphic&jpg=on&ppage=1\">JPG</a></li>" ;

  push @lines, "<li><a class=\"button-style\" href=\"/search/connectedmolecules_landing.shtml" .
       "?connect_molecule=$moleculeName" . "&source_id=$source_id" . "&what=graphic&svg=on&ppage=1\">SVG</a></li>";
 
  push @lines, "<li><a class=\"button-style\" href=\"/search/xaml_conn_landing.shtml" .
          "?connect_molecule=$moleculeName" . "&source_id=$source_id" . "&source=NATURE" . "&what=graphic&xaml=on\">Silverlight-beta</a></li>" ;

  push @lines, "<li><div class=\"helpbox\">";
  push @lines, "<a href=\"#\" class=\"help\">Graphics help</a>";
  push @lines, "<div class=\"description\">";
  push @lines, "<h2 class=\"title\">Information on graphics types</h2>";
  push @lines, "<p class=\"norm\">The pathway can be viewed as a JPG, SVG or Silverlight image. All molecules and interactions within the graphics formats are interactive; clicking on them will open information windows.</p>";
  push @lines, "<p class=\"norm\"><B>JPG:</b>  All common web browsers are able to display JPG images.  However, please note that some browsers may not support very large JPG files.</p>";
  push @lines, "<p class=\"norm\"><B>SVG:</b>  The advantage of SVG images is that they can be more easily zoomed and panned on-screen.  If your browser is not completely SVG enabled, please visit the <a href=\"/PID/userguide/output_formats.shtml\">User guide</a> and follow the instructions for installing the SVG plug-in specific to your browser</p>";
  push @lines, "<p class=\"norm\"><b>SVG image navigation:</b>  Use the left mouse button to activate the image.  To <b>zoom</b> in or out, use the right mouse button and select from the options on the pop-up menu.  To <b>pan</b> the image, hold down the alt key and click-and-drag the image with the mouse.</p>";
  push @lines, "<p class=\"norm\"><b>Silverlight:  </b>Please visit http://www.microsoft.com/getsilverlight/Get-Started/Install/Default.aspx to install the plug-in for your browser.</p>";
  push @lines, "<p class=\"norm\"><b>Silverlight image navigation:  Zoom</b> is available at the top-left of the image. To <b>pan</b> merely drag the image with the mouse.  In the image you may search for specific molecules or all molecules at a certain subcellular location.  To search for molecules you may enter their names, Entrez Gene identifiers or UniProt accession numbers.</p>";
  push @lines, "</div>";
  push @lines, "</div>";
  push @lines, "</li>";
  push @lines, "</ul>";
  push @lines, "</li>";

  push @lines, "<li><span class=\"norm\">Save code as:</span>";
  push @lines, "<ul class=\"links\">";
   push @lines, "<li><a class=\"button-style\" href=\"/search/connectedmolecules_landing.shtml" .
          "?connect_molecule=$moleculeName" . "&source_id=$source_id" . "&what=text&xml=on&ppage=1\">XML</a></li>" ;
  push @lines, "<li><a class=\"button-style\" href=\"/search/connectedmolecules_landing.shtml" .
       "?connect_molecule=$moleculeName" . "&source_id=$source_id" . "&what=text&biopax=on&ppage=1\">BioPAX</a></li>";
  push @lines, "<li><div class=\"helpbox\">";
  push @lines, "<a href=\"#\" class=\"help\">Code help</a>";
  push @lines, "<div class=\"description\">";
  push @lines, "<h2 class=\"title\">Information on types of pathway code</h2>";
  push @lines, "<p class=\"norm\">You can view the pathway code in either XML or BioPax formats.  Please see the <a href=\"/PID/userguide/output_formats.shtml\">User guide</a> for more information.</p>";
  push @lines, "</div>";
  push @lines, "</div>";
  push @lines, "</li>";
  push @lines, "</ul>";
  push @lines, "</li>";

  push @lines, "</ul>";
  if ($format eq 'svg') {
    push @lines, "<p class=\"norm\">A <a href=\"/search/connectedmolecules_landing.shtml" .
          "?connect_molecule=$moleculeName" . "&source_id=$source_id" . "&what=graphic&svg=on/#key\">key to the image</a> icons appears below the pathway.</p>";
  push @lines, "<p class=\"norm\">Click on a biomolecule or interaction for additional information.</p>";
  } else {
    push @lines, "<p class=\"norm\">A <a href=\"/search/connectedmolecules_landing.shtml" .
          "?connect_molecule=$moleculeName" . "&source_id=$source_id" . "&what=graphic&jpg=on/#key\">key to the image</a> icons appears below the pathway.</p>";
  push @lines, "<p class=\"norm\">Click on a biomolecule or interaction for additional information.</p>";

  } 
  return join("\n", @lines, "");
}

###############################################################
sub AdvancedHeader_1 {
  my ($base, $moleculeName, $macroprocess, $source_id, $evidence_code, $forward, $back, $format) = @_;

  $BASE = $base;

  my @ALPHA_MONTHS =
      ("Jan", "Feb", "Mar", "Apr",
       "May", "Jun", "Jul", "Aug",
       "Sep", "Oct", "Nov", "Dec");

  my (@lines);

  my ($sql, $stm, $sqlc, $stmc, $list, $list2, $list3, $list4, $pid);
  my ($source, $pname, $ext_id, $curator, $role, $last_updated,$mol_source_id,$found);
  my ($mname,$alias,$alias_list,$ll_id,$ll_list,$symbol,$name_type);
  $found = 0;
  $moleculeName = lc $moleculeName;

  my $db = DBI->connect("DBI:Oracle:" . "cgprod",
        "web", "readonly");
  if (not $db or $db->err()) {
    print STDERR "#! Cannot connect to " . $db_user . "@" .
        $db_inst . "\n";
    die;
  }

  # Create molecule_list from moleculeName
  my $molecule_list = '';
  my @molList = split(/[,]/,$moleculeName);
  for (my $i = 0; $i < ($#molList); $i++) {
    $molecule_list = $molecule_list . "'" . $molList[$i] . "'" . ",";
  }
  $molecule_list = $molecule_list . "'" . $molList[$#molList] . "'";

  $molecule_list = lc $molecule_list;
  # Get official symbol
  $sql = qq!
select
  distinct mol_name,name_type
from $schema.pw_mol_srch ms, $schema.pw_mol_name mn
where lower(map_name) in ($molecule_list)
and ms.mol_id = mn.mol_id
  !;

  $stm = $db->prepare($sql);
  if(not $stm) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "prepare call failed\n";
    die;
  }
  if(!$stm->execute()) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "execute call failed\n";
    die;
  }

  $list = '';
  my $first = 1;
  while (($symbol, $name_type) =
      $stm->fetchrow_array()) {
        if ($name_type eq "OF") {
          if ($first == 1) {
            $list = $symbol;
            $first = 0;
          } else {
            $list = $list . ', ' . $symbol;
          }
        }
  }

  $stm->finish();
  $db->disconnect;
  push @lines, "<ul class=\"xoxo viewing-option\">";
   push @lines, "<ul class=\"links\">";

  push @lines, "<li><a class=\"button-style\" href=\"/search/network_mol_list_landing.shtml" .
         "?source=$source_id" . "&molecule=$moleculeName\">Molecule list</a></li>";
  push @lines, "</ul></li>";
   push @lines, "<li><span class=\"norm\">View graphic as:</span>";

  push @lines, "<ul class=\"links\">";
 push @lines, "<li><a class=\"button-style\" href=\"/search/advanced_landing.shtml" .

          "?molecule=$moleculeName"  . "&pathway=$pid" . "&macro_process=$macroprocess" . "&source_id=$source_id" . "&forward=$forward" . "&back=$back" . "&what=graphic&jpg=true&complex_uses=on&family_uses=on&degree=1&ppage=1\">JPG</a></li>" ;

  push @lines, "<li><a class=\"button-style\" href=\"/search/advanced_landing.shtml" .
       "?molecule=$moleculeName"  . "&pathway=$pid" . "&macro_process=$macroprocess" . "&source_id=$source_id" . "&forward=$forward" . "&back=$back" . "&what=graphic&svg=on&complex_uses=on&family_uses=on&degree=1&ppage=1\">SVG</a></li>";

  $source="NATURE";
  if (($source_id == 2) || ($source_id == 3)) {
    $source="BioCarta";
  }

  if ($source_id == 7) {
    $source="Reactome";
  } 

  my $forward_clause = '';

  if ($forward eq 'forward') {
    $forward_clause = "&forward=forward";
  }
  
  if ($back eq 'back') {
    $forward_clause = "&back=back";
  }

  if (($forward eq 'forward') && ($back eq 'back')) {
    $forward_clause = "&forward=forward&back=back";
  }

 
  push @lines, "<li><a class=\"button-style\" href=\"/xaml_adv_landing.shtml" .
       "?what=graphic&svg=&xaml=true&xml=&biopax=&complex_uses=on&family_uses=on&degree=1" .
       "&molecule=$moleculeName" . "&pathway=$pid" . "&macro_process=$macroprocess" . "&source_id=$source_id" .
       "$forward_clause" .
       "&evidence_code=NIL&evidence_code=IAE&evidence_code=IC&evidence_code=IDA&evidence_code=IFC&evidence_code=IGI&evidence_code=IMP&evidence_code=IOS&evidence_code=IPI&evidence_code=RCA&evidence_code=RGE&evidence_code=TAS&output-format=graphic&source=$source&Submit=Go\">Silverlight-Beta</a></li>";


  push @lines, "<li><div class=\"helpbox\">";
  push @lines, "<a href=\"#\" class=\"help\">Graphics help</a>";
  push @lines, "<div class=\"description\">";
  push @lines, "<h2 class=\"title\">Information on graphics types</h2>";
  push @lines, "<p class=\"norm\">The pathway can be viewed as a JPG, SVG or Silverlight image. All molecules and interactions within the graphics formats are interactive; clicking on them will open information windows.</p>";
  push @lines, "<p class=\"norm\"><B>JPG:</b>  All common web browsers are able to display JPG images.  However, please note that some browsers may not support very large JPG files.</p>";
  push @lines, "<p class=\"norm\"><B>SVG:</b>  The advantage of SVG images is that they can be more easily zoomed and panned on-screen.  If your browser is not completely SVG enabled, please visit the <a href=\"/PID/userguide/output_formats.shtml\">User guide</a> and follow the instructions for installing the SVG plug-in specific to your browser</p>";
  push @lines, "<p class=\"norm\"><b>SVG image navigation:</b>  Use the left mouse button to activate the image.  To <b>zoom</b> in or out, use the right mouse button and select from the options on the pop-up menu.  To <b>pan</b> the image, hold down the alt key and click-and-drag the image with the mouse.</p>";
  push @lines, "<p class=\"norm\"><b>Silverlight:  </b>Please visit http://www.microsoft.com/getsilverlight/Get-Started/Install/Default.aspx to install the plug-in for your browser.</p>";
  push @lines, "<p class=\"norm\"><b>Silverlight image navigation:  Zoom</b> is available at the top-left of the image. To <b>pan</b> merely drag the image with the mouse.  In the image you may search for specific molecules or all molecules at a certain subcellular location.  To search for molecules you may enter their names, Entrez Gene identifiers or UniProt accession numbers.</p>";
  push @lines, "</div>";
  push @lines, "</div>";
  push @lines, "</li>"; 
  push @lines, "</ul>";
  push @lines, "</li>";

  push @lines, "<li><span class=\"norm\">Save code as:</span>";
  push @lines, "<ul class=\"links\">";

   push @lines, "<li><a class=\"button-style\" href=\"/search/advanced_landing.shtml" .
          "?molecule=$moleculeName"  . "&pathway=$pid" . "&macro_process=$macroprocess" . "&source_id=$source_id" . "&forward=$forward" . "&back=$back" . "&what=text&xml=on&complex_uses=on&family_uses=on&degree=1&ppage=1\">XML</a></li>" ;
  push @lines, "<li><a class=\"button-style\" href=\"/search/advanced_landing.shtml" .
       "?molecule=$moleculeName"  . "&pathway=$pid" . "&macro_process=$macroprocess" . "&source_id=$source_id" . "&forward=$forward" . "&back=$back" . "&what=text&biopax=on&complex_uses=on&family_uses=on&degree=1&ppage=1\">BioPAX</a></li>";
  push @lines, "<li><div class=\"helpbox\">";
  push @lines, "<a href=\"#\" class=\"help\">Code help</a>";
  push @lines, "<div class=\"description\">";
  push @lines, "<h2 class=\"title\">Information on types of pathway code</h2>";
  push @lines, "<p class=\"norm\">You can view the pathway code in either XML or BioPax formats.  Please see the <a href=\"/PID/userguide/output_formats.shtml\">User guide</a> for more information.</p>";
  push @lines, "</div>";
  push @lines, "</div>";
  push @lines, "</li>";
  push @lines, "</ul>";
  push @lines, "</li>";
 
  push @lines, "</ul>";
  if ($format eq 'svg') {
    push @lines, "<p class=\"norm\">A <a href=\"/search/advanced_landing.shtml" . 
          "?molecule=$moleculeName"  . "&pathway=$pid" .  "&macro_process=$macroprocess" . "&source_id=$source_id" . "&what=graphic&svg=on&ppage=1/#key\">key to the image</a> icons appears below the pathway.</p>";
  push @lines, "<p class=\"norm\">Click on a biomolecule or interaction for additional information.</p>";            
  } else {
     push @lines, "<p class=\"norm\">A <a href=\"/search/advanced_landing.shtml" .
          "?molecule=$moleculeName"  . "&pathway=$pid" . "&macro_process=$macroprocess" . "&source_id=$source_id" . "&what=graphic&jpg=on&ppage=1/#key\">key to the image</a> icons appears below the pathway.</p>";
  push @lines, "<p class=\"norm\">Click on a biomolecule or interaction for additional information.</p>";
  } 
  return join("\n", @lines) . "\n";

}
######################################################################
sub InteractionHeader_1 {
  my ($base, $atomid, $format) = @_;

     my $pname;
   my @lines;
  my $scancheck = Scan($atomid) + Scan($format);
  
  if ($scancheck > 0) {
    push @lines, "Bad data in parameters\n";
    return join("\n", @lines) . "\n";
 
  }
  
   my ($label, $label_value, $lvid, $atomtype, $pathwayid, $edge);

  my $db = DBI->connect("DBI:Oracle:" . $db_inst, $db_user, $db_pass);
  if (not $db or $db->err()) {
    print STDERR "Cannot connect to " . $db_user . "@" . $db_inst . "\n";
    die;
  }

  my $pw   = new Pathway($lv);
  my $pdb  = new PathwayDB($db, $schema,
                           $LIMIT_SOURCES, $LIMIT_EVIDENCE, $pw, $lv);
  my ($org, $srcid, $srcname) = $pdb->GetAtomBasics($atomid);
  $srcname =~ s/NATURE/NCI-Nature curated/;
  $srcname =~ s/BioCarta/BioCarta Imported/;
  $pdb->FillAtoms($atomid);

  my $other_format_line =

      "<li><a type=button href=\"/search/interaction_landing.shtml" .
      "?atom_id=$atomid" . "&what=graphic&svg=on\">SVG</a></li>" .
      "&nbsp;&nbsp;<li><a href=\"/search/interaction_landing.shtml" .
      "?atom_id=$atomid" . "&what=graphic&jpg=on\">JPG</a></li>" .
      "&nbsp;&nbsp;<li><a href=\"/search/interaction_landing.shtml" .
      "?atom_id=$atomid" . "&what=text&xml=on\">XML</a></li>" .
      "&nbsp;&nbsp;<li><a href=\"/search/interaction_landing.shtml" .
      "?atom_id=$atomid" . "&what=text&biopax=on\">BioPAX</a></li>";

  push @lines , "<form name=\"dpform\" action=\"/search/interaction_landing.shtml\" method=\"get\" id=\"advancedsearch\">";
  push @lines , "<input id=what type=\"hidden\" name=\"what\" value=\"graphic\"/> ";
  push @lines , "<input id=svg type=hidden name=svg value=true> ";
  push @lines , "<input id=jpg type=hidden name=jpg value=\"\"> ";
  push @lines , "<input id=xml type=hidden name=xml value=\"\"> ";
  push @lines , "<input id=biopax type=hidden name=biopax value=\"\"> ";
  push @lines , "<input type=hidden name=complex_uses value=on /> ";
  push @lines , "<input type=hidden name=family_uses value=on /> ";
  push @lines , "<input type=hidden name=degree value=1 /> ";
  push @lines , "<input type=\"hidden\" name=\"atom_id\" value=$atomid /> ";
  # push @lines, "<ul class=\"xoxo links\">";

  push @lines , "<B>Include Upstream/Downstream Interactions (optional)</B><BR><BR>";
  push @lines , "<input id=\"upstream\" name=\"back\" value=\"back\" type=\"checkbox\" /> ";
  push @lines , "<label for=\"upstream\">Upstream</label> ";
  push @lines , "<input id=\"downstream\" name=\"forward\" value=\"forward\" type=\"checkbox\" /> ";
  push @lines , "<label for=\"downstream\">Downstream</label> ";
push @lines , "<span class=\"label1\">";
  push @lines , "<br><br><input id=svg2 name=output-format value=SVG type=submit class=submit checked=checked set what=graphic ";
  push @lines , "onclick=\"dpform.what.value='graphic';dpform.svg.value=true;dpform.jpg.value='';dpform.xml.value='';dpform.biopax.value=''dpform.xaml.value=''\"/> ";
  push @lines , "</span>";

  push @lines , "<span class=\"label1\">";
  push @lines , "<input id=jpg2 name=output-format value=JPG type=submit class=submit set what=graphic ";
  push @lines , "onclick=\"dpform.what.value='graphic';dpform.svg.value='';dpform.jpg.value=true;dpform.xml.value='';dpform.biopax.value=''dpform.xaml.value=''\"/> ";
  push @lines , "</span>";

  push @lines , "<span class=\"label1\">";
  push @lines , "<input id=xml2 name=output-format value=XML type=submit class=submit set what=text ";
  push @lines , "onclick=\"dpform.what.value='text';dpform.svg.value='';dpform.jpg.value='';dpform.xml.value=true;dpform.biopax.value=''dpform.xaml.value=''\"/> ";
  push @lines , "</span>";

  push @lines , "<span class=\"label1\">";
  push @lines , "<input id=BioPax2 name=output-format value=BioPAX type=submit class=submit set what=text ";
  push @lines , "onclick=\"dpform.what.value='text';dpform.svg.value='';dpform.jpg.value='';dpform.xml.value='';dpform.biopax.value=true dpform.xaml.value=''\"/> ";
  push @lines , "</span>";

  push @lines , "<span class=\"label1\">";
  push @lines , "<input id=xaml2 name=output-format value=Silverlight-beta type=submit class=submit set what=graphic ";
  push @lines , "onclick=\"dpform.what.value='graphic';dpform.svg.value='';dpform.jpg.value='';dpform.xml.value='';dpform.biopax.value='';dpform.xaml.value=true\"/> ";
  push @lines , "</span>";

  push @lines , "</form>";
  # BasicInfoForAtom($pw, $atomid, $org, $srcid, $srcname, \@lines);
  # EdgeInfoForAtom($pw, $atomid, \@lines);
  # NotesForAtom($db, $atomid, \@lines);
  # ReferencesForAtom($db, $atomid, \@lines);
  # EvidenceForAtom($db, $atomid, \@lines);
  # PathwayInfoForAtom($db, $atomid, \@lines);
  # SearchDupAtoms($db, $pdb, $pw, $atomid, \@lines);

  $db->disconnect();

  return join("\n", @lines) . "\n";

}
######################################################################
sub SearchHeader_1 {
  my ($base, $atomid, $format) = @_;

    
   my $pname;
   my @lines; 
   my $scancheck = Scan($atomid) + Scan($format);
   if ($scancheck > 0) {
    push @lines, "Bad data in parameters\n";
    return join("\n", @lines) . "\n";

  }
 
  push @lines , "<input id=what type=\"hidden\" name=\"what\" value=\"graphic\" /> ";
  push @lines , "<input id=svg type=hidden name=svg value=true> ";
  push @lines , "<input id=jpg type=hidden name=jpg value=\"\"> ";
  push @lines , "<input id=xml type=hidden name=xml value=\"\"> ";
  push @lines , "<input id=biopax type=hidden name=biopax value=\"\"> ";
  push @lines , "<input type=hidden name=complex_uses value=on /> ";
  push @lines , "<input type=hidden name=family_uses value=on /> ";
  push @lines , "<input type=hidden name=degree value=1 /> ";
  push @lines , "<input type=\"hidden\" name=\"atom_id\" value=$atomid /> ";

    if ($format eq 'svg') {
    # push @lines, "<p class=\"norm\">A <a href=\"/search/network_landing.shtml?pathway=$pid" . "&source=$source" . "&jpg=true&what=graphic&ppage=1/#key\">key to the image</a> icons appears below the pathway.</p>";
  } else {
    #push @lines, "<p class=\"norm\">A <a href=\"/search/network_landing.shtml?pathway=$pid" . "&source=$source" . "&svg=true&what=graphic&ppage=1/#key\">key to the image</a> icons appears below the pathway.</p>";
  }
  push @lines, "<table border=0>";
  push @lines, "<tr><td><b>Click on a biomolecule or interaction for additional information.</b></td></tr>";
    push @lines, "</table>";

  return join("\n", @lines) . "\n";

 
}

##################################################################
sub BatchIntermediatePage_1 {
  my ($base);
 
  my %genes_a;
  my %genes_b;
  my (@rows);

  # Query Database
  # List results and p-values
  
  my $db = DBI->connect("DBI:Oracle:" . "cgprod",
        "web", "readonly");
 
  if (not $db or $db->err()) {
    print STDERR "Cannot connect to " . $db_user . "@" . $db_inst . "\n";
    die;
  }

  print STDERR "Parameters:  $base\n";
  my $g;
  my $error;
  my $file_a = 'c:\\search1.txt';
  my $file_b = 'c:\\search2.txt';
  if ($file_a) {
    while (<$file_a>) {
      print STDERR "In while $file_a\n"; 
      s/\n//g;
      s/\r//g;
      s/^ +//;
      ($g) = split /\t/;
      $g = lc($g);
      if ($g eq "") {
        next;
      }
      $genes_a{$g} = 1;
    }
  } else {
    print STDERR "Did not find $file_a\n";
  }

  if ($file_b) {
    while (<$file_b>) {
      s/\n//g;
      s/\r//g;
      s/^ +//;
      ($g) = split /\t/;
      $g = lc($g);
      if ($g eq "") {
        next;
      }
      $genes_b{$g} = 1;
    }
  }
  my $format = 'html'; 
  @rows = RankPathwaysByGeneHits($db, $schema, $format, \%genes_a, \%genes_b); 
  $db->disconnect();
  print STDERR "RankPathwaysByGeneHits\n"; 
  return join("\n", @rows) . "\n";

}
##################################################################
sub IntermediatePage_1 {
 my ($base, $moleculeName) = @_;

 my ($orig_mol_name, $keyvalue, $stm, $map_name, $pathway_name, $pathway_source_id, $pid, $atom_id, $mol_id);
 my (@lines);
 my %natureMap;
 my %naturePathway;
 my %biocartaMap;
 my %biocartaMac;
 my %biocartaPathway;
 my %reactomeMap;
 my %reactomePathway;
 my %naturePid;
 my %biocartaPid;
 my %reactomePid;
 my $natureData;
 my $biocartaData;
 my $reactomeData;
 my $macroName;
 my %natureMac;
 my %biocartaMac;
 my %reactomeMac;
 my %natureEdge;
 my %biocartaEdge;
 my %reactomeEdge;
 my %natureMol;
 my %natureMacro;
 my $genes_a = "";

  my $scancheck = Scan($moleculeName) + Scan($base);
  if ($scancheck > 0) {
    push @lines, "Bad data in parameters\n";
    return join("\n", @lines) . "\n";

  }

 my $db = DBI->connect("DBI:Oracle:" . "cgprod",
        "web", "readonly");
  if (not $db or $db->err()) {
    print STDERR "#! Cannot connect to " . $db_user . "@" .
        $db_inst . "\n";
    die;
  }


$moleculeName = lc $moleculeName;
$moleculeName =~s/, /,/g;
$moleculeName =~s/ ,/,/g;
$moleculeName =~s/\*/%/g;

print STDERR "MoleculeName:  $moleculeName\n";
$orig_mol_name = $moleculeName;

my @moleculeArray = split(/,/,$moleculeName);
$moleculeName = '';
$macroName = '';
for (my $j = 0; $j <= $#moleculeArray; $j++) {
  # if (substr($moleculeArray[$j],0,1) eq ' ') {
  #  $moleculeArray[$j] = substr($moleculeArray[$j],1);
  # }
  my $lengthTest = $moleculeArray[$j]; 
  $lengthTest =~s/%//g;
  my $lengthValue = length($lengthTest);
  print STDERR "Length:  $lengthValue\n";
  if ($lengthValue < 3) {
    $moleculeArray[$j] =~s/%//g;
  } 
  $moleculeName = $moleculeName . "'" . $moleculeArray[$j] . "'";
  print STDERR "MoleculeName:  $moleculeName\n";
  $macroName = $macroName . "'%" . $moleculeArray[$j] . "%'";
  if ($j < $#moleculeArray) {
    $moleculeName = $moleculeName . ',';
    $macroName = $macroName . ',';
  } 
}


my $sql = qq!
 select unique
      mn.mol_name,pathway_name,pathway_source_id,p.pathway_id,s1.mol_id
    from
      $schema.pw_mol_srch s1,
      $schema.pw_mol_mol m1,
      $schema.pw_edge e1,
      $schema.pw_pathway_atom pa,
      $schema.pw_pathway p,
      $schema.pw_mol_name mn
    where
      m1.mol_id_1 in 
        (select mm.mol_id_2 from $schema.pw_mol_srch s1, $schema.pw_mol_mol mm
          where s1.map_name like ? and s1.mol_id = mm.mol_id_1)
      and s1.mol_id = m1.mol_id_1
      and m1.mol_id_2 = e1.mol_id
      and pa.pathway_id = p.pathway_id
      and m1.mol_id_2 = e1.mol_id
      and e1.atom_id = pa.atom_id
      and mn.mol_id = s1.mol_id
      and s1.map_name like ? 
!;

print STDERR "SQL:  $sql\n";

  $stm = $db->prepare($sql);
  if(not $stm) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "prepare call failed\n";
    die;
  }
  for my $id (@moleculeArray) {
    if(!$stm->execute($id,$id)) {
      print STDERR "$sql\n";
      print STDERR "$DBI::errstr\n";
      print STDERR "execute call failed\n";
      die;
    }
  
    while (($map_name, $pathway_name, $pathway_source_id, $pid, $mol_id) =
      $stm->fetchrow_array()) {
      if ($pathway_source_id == 5) {
      $naturePathway{$map_name . $pathway_name} = $pathway_name;
      $natureMap{$map_name . $pathway_name} = $map_name;
      $naturePid{$map_name . $pathway_name} = $pid;
      $natureEdge{$map_name . $pathway_name} = $mol_id;
      $natureData = 1;
      }
      if ($pathway_source_id == 7) {
      $reactomePathway{$map_name . $pathway_name} = $pathway_name;
      $reactomeMap{$map_name . $pathway_name} = $map_name;
      $reactomePid{$map_name . $pathway_name} = $pid;
      $reactomeEdge{$map_name . $pathway_name} = $mol_id;
      $reactomeData = 1; 
    }
    if (($pathway_source_id == 2) || ($pathway_source_id == 3)) {
      $biocartaPathway{$map_name . $pathway_name} = $pathway_name;
      $biocartaMap{$map_name . $pathway_name} = $map_name;
      $biocartaPid{$map_name . $pathway_name} = $pid;
      $biocartaEdge{$map_name . $pathway_name} = $mol_id;
      $biocartaData = 1;
    }

  }

  $stm->finish();
  }
  print STDERR "macroName:  $macroName\n";
  my @macro=split(/,/,$macroName);
  for (my $j = 0; $j <= $#macro; $j++) 
  {
    my $macroTemp = $macro[$j];
  $sql = qq!
    select label_value_name,pathway_name,pathway_source_id,p.pathway_id
    from $schema.pw_label_value lv, $schema.pw_atom_label al,
    $schema.pw_pathway_atom pa,
    $schema.pw_pathway p
    where lower(label_value_name) like lower($macroTemp)
    and lv.label_value_id = al.label_value_id
    and pa.pathway_id = p.pathway_id
    and pa.atom_id = al.atom_id 
  !;
  print STDERR "Macro: $sql\n";
  $stm = $db->prepare($sql);
  if(not $stm) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "prepare call failed\n";
    die;
  }
  if(!$stm->execute()) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "execute call failed\n";
    die;
  }
  while (($map_name, $pathway_name, $pathway_source_id, $pid) =
      $stm->fetchrow_array()) {
    if ($pathway_source_id == 5) {
      $naturePathway{$map_name . $pathway_name} = $pathway_name;
      $natureMac{$map_name . $pathway_name} = $map_name;
      $naturePid{$map_name . $pathway_name} = $pid;
      $natureMol{$map_name} = 1;
      $natureMacro{$map_name} = 1; 
      $natureData = 1;
    }
    if ($pathway_source_id == 7) {
      $reactomePathway{$map_name . $pathway_name} = $pathway_name;
      $reactomeMac{$map_name . $pathway_name} = $map_name;
      $reactomePid{$map_name . $pathway_name} = $pid;
      $reactomeData = 1;
    }
    if (($pathway_source_id == 2) || ($pathway_source_id == 3)) {
      $biocartaPathway{$map_name . $pathway_name} = $pathway_name;
      $biocartaMac{$map_name . $pathway_name} = $map_name;
      $biocartaPid{$map_name . $pathway_name} = $pid;
      $biocartaData = 1;
    }

 
  }
  } 
  $stm->finish(); 
  my $oldmapname = '';
  if ($natureData + $biocartaData + $reactomeData == 0) {
  push @lines, "<h3>No Results Found in PID.  Please refine your search terms</h3>";
  } else { 
  push @lines, "<h2 id=\"NCI-Nature\" class=\"pagetitle2\">NCI-Nature curated</h2>";

  if ($natureData == 1) {
    $map_name='';
    my $first = 1;
    for $keyvalue (sort keys %natureMap) {
      my $testMapName = $natureMap{$keyvalue};
      if ($testMapName ne $oldmapname) { 
        if ($first == 1) {
          $first = 0;
        } else {
          $map_name = ',' . $map_name;
        }
        $map_name = $natureMap{$keyvalue} . $map_name;
        $oldmapname = $natureMap{$keyvalue}; 
        print STDERR "Map_name:  $map_name\n";
      }
    }
    print STDERR "Map_name:  $map_name\n";
    my $mac_name = '';
    my $first = 1;
    my $testMacName = '';
    my $oldmacname = '';
    for $keyvalue (sort keys %natureMac) {
      my $testMacName = $natureMac{$keyvalue};
      if ($testMacName ne $oldmacname) {
        if ($first == 1) {
          $first = 0;
        } else {
          $mac_name = ',' . $mac_name;
        }
        $mac_name = $natureMac{$keyvalue} . $mac_name;
        $oldmacname = $natureMac{$keyvalue};
      }

    }
     push @lines, "<p class=\"norm\">View the <a href=\"/search/network_landing.shtml?molecule=$orig_mol_name" . 
     "&macro_process=$orig_mol_name&family_uses=on&complex_uses=on&source_id=5" .
    "&what=graphic&jpg=on\">" .
    "complete network maps for your search terms</a> or view the predefined pathways below:</p>";
    push @lines, "<table class=data border=0 cellpadding=0 cellspacing=0>";
    push @lines, "<thead>";
    push @lines, "<tr>";
    push @lines, "<th scope=\"colgroup\" class=\"colspan column1\">Search results for</th>";
    push @lines, "<th scope=\"colgroup\" class=\"colspan column2\">Pathway name</th>";
    push @lines, "<th scope=\"colgroup\" class=\"colspan column3\">Source</th>";
    push @lines, "</tr>";
    push @lines, "</thead><tbody>";
    my $n = 0;
    my $source = 5;
    my %gd;
    my %terms_a;
    my %terms_b;
    my $mn = $moleculeName;
    $mn=~s/'//g;
    for my $j (split(",", $mn)) {
      $terms_a{$j} = 1;
      print STDERR "Term:  $j\n";
    }
    my ($genes_a, $genes_b) = GetEntrezGeneIds($db, $schema, \%terms_a, \%terms_b, $source, \%gd);
    my $genelist='';
    my $first = 1;
    for my $gtest (keys % {$genes_a}) {
      if ($first == 1) {
        $genelist = $gtest;
        $first = 0;
      } else {
        $genelist = $genelist . ',' .$gtest;
      } 
      print STDERR "$gtest\n";
    }
 
    for $keyvalue (sort keys %natureMap) {
      $map_name = $natureMap{$keyvalue};

      $pathway_name = $naturePathway{$keyvalue};
      $pid = $naturePid{$keyvalue};
      my $instance = $natureEdge{$keyvalue};
      $n++;
      push @lines, ($n %2 == 0) ? "<tr>" : "<tr class=\"odd\">";
      push @lines, "<td><a href=\"MoleculePage?molid=$instance" .
      "\">$map_name</a></td><td><a href=\"/search/pathway_landing.shtml" .
      "?what=graphic&jpg=on" .
      "&pathway_id=$pid&source=NATURE&output-format=graphic&ppage=1" .
      "&genes_a=$genelist " .
      "\">$pathway_name</a>" .
      "</td><td>NCI-Nature curated</td></tr>";
    }
  for $keyvalue (sort keys %natureMac) {
    $map_name = $natureMac{$keyvalue};
    $pathway_name = $naturePathway{$keyvalue};
    $pid = $naturePid{$keyvalue};
    $n++; 
    push @lines, ($n %2 == 0) ? "<tr>" : "<tr class=\"odd\">";
    push @lines, "<td><a href=\"/search/advanced_landing.shtml" .
      "?pathway=&molecule=&macro_process=$map_name" . "&family_uses=on&complex_uses=on&source_id=5" .
      "&what=graphic&jpg=true\">$map_name</a></td><td><a href=\"/search/pathway_landing.shtml" .
      "?what=graphic&jpg=on" .
      "&pathway_id=$pid&source=NATURE&output-format=graphic&ppage=1 " .
      "\">$pathway_name</a>" .
      "</td><td>NCI-Nature curated</td></tr>";
  }
  push @lines,"</tbody>"; 
  push @lines,"</table>";
  } else {
    $moleculeName =~ s/%//g; 
    push @lines, "<p><br>Sorry, $moleculeName did not return any results from NCI-Nature";
}
  push @lines, "<p><p><br><h2 id=\"BioCarta\" class=\"pagetitle2\">BioCarta</h2>";

  if ($biocartaData == 1) {
    $map_name='';
    my $first = 1;
    for $keyvalue (sort keys %biocartaMap) {
      my $testMapName = $biocartaMap{$keyvalue};
      if ($testMapName ne $oldmapname) {
        if ($first == 1) {
          $first = 0;
        } else {
          $map_name = ',' . $map_name;
        }
        $map_name = $biocartaMap{$keyvalue} . $map_name;
        $oldmapname = $biocartaMap{$keyvalue};
        print STDERR "Map_name:  $map_name\n";
      }
    }
    print STDERR "Map_name:  $map_name\n";
    my $mac_name = '';
    my $first = 1;
    my $testMacName = '';
    my $oldmacname = '';
    for $keyvalue (sort keys %biocartaMac) {
      my $testMacName = $biocartaMac{$keyvalue};
      if ($testMacName ne $oldmacname) {
        if ($first == 1) {
          $first = 0;
        } else {
          $mac_name = ',' . $mac_name;
        }
        $mac_name = $biocartaMac{$keyvalue} . $mac_name;
        $oldmacname = $biocartaMac{$keyvalue};
      }

    }
 
  push @lines, "<p class=\"norm\">View the <a href=\"/search/network_landing.shtml?molecule=$orig_mol_name" .
  "&macro_process=$orig_mol_name&family_uses=on&complex_uses=on&source_id=2,3" .
  "&what=graphic&jpg=on\">" .
  "complete network maps for your search terms</a> or view the predefined pathways below:</p>";

  push @lines, "<table class=data border=0 cellpadding=0 cellspacing=0>";
  push @lines, "<thead>";
  push @lines, "<tr>";
  push @lines, "<th scope=\"colgroup\" class=\"colspan column1\">Search results for</th>";
  push @lines, "<th scope=\"colgroup\" class=\"colspan column2\">Pathway name</th>";
  push @lines, "<th scope=\"colgroup\" class=\"colspan column3\">Source</th>";
  push @lines, "</tr>";
  push @lines, "</thead><tbody>";
  my $n = 0;
  
  my $source = 2;
    my %gd;
    my %terms_a;
    my %terms_b;
    my $mn = $moleculeName;
    $mn=~s/'//g;
    for my $j (split(",", $mn)) {
      $terms_a{$j} = 1;
      print STDERR "Term:  $j\n";
    }
    my ($genes_a, $genes_b) = GetEntrezGeneIds($db, $schema, \%terms_a, \%terms_b, $source, \%gd);
    my $genelist='';
    my $first = 1;
    for my $gtest (keys % {$genes_a}) {
      if ($first == 1) {
        $genelist = $gtest;
        $first = 0;
      } else {
        $genelist = $genelist . ',' .$gtest;
      }
      print STDERR "$gtest\n";
    }
 
  for $keyvalue (sort keys %biocartaMap) {
    $map_name = $biocartaMap{$keyvalue};
    $pathway_name = $biocartaPathway{$keyvalue};
    $pid = $biocartaPid{$keyvalue};
    my $instance = $biocartaEdge{$keyvalue};
    $n++;
    push @lines, ($n %2 == 0) ? "<tr>" : "<tr class=\"odd\">";

    push @lines, "<td><a href=\"MoleculePage?molid=$instance" .
      "\">$map_name</a></td><td><a href=\"/search/pathway_landing.shtml" .
      "?what=graphic&jpg=on" .
      "&pathway_id=$pid&source=2&output-format=graphic&ppage=1" .
      "&genes_a=$genelist " .
 
      "\">$pathway_name</a>" .
      "</td><td>BioCarta</td></tr>";

     # push @lines, "<tr><td>$map_name</td><td>$pathway_name</td><td>BioCarta</td></tr>";
  }
    for $keyvalue (sort keys %biocartaMac) {
    $map_name = $biocartaMac{$keyvalue};
    $pathway_name = $biocartaPathway{$keyvalue};
    $pid = $biocartaPid{$keyvalue};
    $n++;
    push @lines, ($n %2 == 0) ? "<tr>" : "<tr class=\"odd\">";
    push @lines, "<td><a href=\"/search/advanced_landing.shtml" .
      "?pathway=&molecule=&macro_process=$map_name" . "&family_uses=on&complex_uses=on&source_id=2,3" .
      "&what=graphic&jpg=true\">$map_name</a></td><td><a href=\"/search/pathway_landing.shtml" .
      "?what=graphic&jpg=on" .
      "&pathway_id=$pid&source=2,3&output-format=graphic&ppage=1 " .
      "\">$pathway_name</a>" .
      "</td><td>BioCarta</td></tr>";
  }
    push @lines,"</tbody>"; 
    push @lines,"</table>";
  } else {
   $moleculeName =~ s/%//g;
   push @lines, "<p><br>Sorry, $moleculeName did not return any results from BioCarta";
}
  push @lines, "<p><p><br><h2 id=\"Reactome\" class=\"pagetitle2\">Reactome</h2>";

  if ($reactomeData == 1) {
    $map_name='';
    my $first = 1;
    for $keyvalue (sort keys %reactomeMap) {
      my $testMapName = $reactomeMap{$keyvalue};
      if ($testMapName ne $oldmapname) {
        if ($first == 1) {
          $first = 0;
        } else {
          $map_name = ',' . $map_name;
        }
        $map_name = $reactomeMap{$keyvalue} . $map_name;
        $oldmapname = $reactomeMap{$keyvalue};
        print STDERR "Map_name:  $map_name\n";
      }
    }
     my $mac_name = '';
    my $first = 1;
    my $testMacName = '';
    my $oldmacname = '';
    for $keyvalue (sort keys %reactomeMac) {
      my $testMacName = $reactomeMac{$keyvalue};
      if ($testMacName ne $oldmacname) {
        if ($first == 1) {
          $first = 0;
        } else {
          $mac_name = ',' . $mac_name;
        }
        $mac_name = $reactomeMac{$keyvalue} . $mac_name;
        $oldmacname = $reactomeMac{$keyvalue};
      }

    } 

   push @lines, "<p class=\"norm\">View the <a href=\"/search/network_landing.shtml?molecule=$orig_mol_name" .
   "&macro_process=$orig_mol_name&family_uses=on&complex_uses=on&source_id=7" .
  "&what=graphic&jpg=on\">" .
  "complete network maps for your search terms</a> or view the predefined pathways below:</p>";
 
  push @lines, "<table class=data border=0 cellpadding=0 cellspacing=0>";
  push @lines, "<thead>";
  push @lines, "<tr>";
  push @lines, "<th scope=\"colgroup\" class=\"colspan column1\">Search results for</th>";
  push @lines, "<th scope=\"colgroup\" class=\"colspan column2\">Pathway name</th>";
  push @lines, "<th scope=\"colgroup\" class=\"colspan column3\">Source</th>";
  push @lines, "</tr>";
  push @lines, "</thead><tbody>";
  my $n = 0;
  my $source = 7;
    my %gd;
    my %terms_a;
    my %terms_b;
    my $mn = $moleculeName;
    $mn=~s/'//g;
    for my $j (split(",", $mn)) {
      $terms_a{$j} = 1;
      print STDERR "Term:  $j\n";
    }
    my ($genes_a, $genes_b) = GetEntrezGeneIds($db, $schema, \%terms_a, \%terms_b, $source, \%gd);
    my $genelist='';
    my $first = 1;
    for my $gtest (keys % {$genes_a}) {
      if ($first == 1) {
        $genelist = $gtest;
        $first = 0;
      } else {
        $genelist = $genelist . ',' .$gtest;
      }
      print STDERR "$gtest\n";
    }
 
  for $keyvalue (sort keys %reactomeMap) {
    $map_name = $reactomeMap{$keyvalue};
    $pathway_name = $reactomePathway{$keyvalue};
    $pid = $reactomePid{$keyvalue};
    my $instance = $reactomeEdge{$keyvalue};
    $n++;
    push @lines, ($n %2 == 0) ? "<tr>" : "<tr class=\"odd\">";
      push @lines, "<td><a href=\"MoleculePage?molid=$instance" .
      "\">$map_name</a></td><td><a href=\"/search/pathway_landing.shtml" .
     "?what=graphic&jpg=on" .
      "&pathway_id=$pid&source=7&output-format=graphic&ppage=1" .
      "&genes_a=$genelist " .
 
      "\">$pathway_name</a>" .
      "</td><td>Reactome</td></tr>";

    # push @lines, "<tr><td>$map_name</td><td>$pathway_name</td><td>Reactome</td></tr>";
  }
   for $keyvalue (sort keys %reactomeMac) {
    $map_name = $reactomeMac{$keyvalue};
    $pathway_name = $reactomePathway{$keyvalue};
    $pid = $reactomePid{$keyvalue};
    $n++;
    push @lines, ($n %2 == 0) ? "<tr>" : "<tr class=\"odd\">";
    push @lines, "<td><a href=\"/search/advanced_landing.shtml" .
      "?pathway=&molecule=&macro_process=$map_name" . "&family_uses=on&complex_uses=on&source_id=7" .
      "&what=graphic&jpg=true\">$map_name</a></td><td><a href=\"/search/pathway_landing.shtml" .
      "?what=graphic&jpg=on" .
      "&pathway_id=$pid&source=7&output-format=graphic&ppage=1 " .
      "\">$pathway_name</a>" .
      "</td><td>Reactome</td></tr>";
  }
  push @lines,"</tbody>"; 
  push @lines,"</table>";
  } else {
    $moleculeName =~ s/%//g;
    push @lines, "<p><br>Sorry, $moleculeName did not return any results from Reactome";

  }
  }
  $db->disconnect();
  return join("\n", @lines) . "\n";

}
 
######################################################################
sub MoleculeHeader_1 {
  my ($base, $moleculeName, $source_id, $format) = @_;

  $BASE = $base;

  my @ALPHA_MONTHS =
      ("Jan", "Feb", "Mar", "Apr",
       "May", "Jun", "Jul", "Aug",
       "Sep", "Oct", "Nov", "Dec");

  my (@lines);
   my $scancheck = Scan($moleculeName) + Scan($format) + Scan($source_id);
  if ($scancheck > 0) {
    push @lines, "Bad data in parameters\n";
    return join("\n", @lines) . "\n";

  }

  my ($sql, $stm, $sqlc, $stmc, $list);
  my ($source, $pname, $ext_id, $curator, $role, $last_updated,$mol_source_id,$found);
  my ($mname,$alias,$alias_list,$ll_id,$ll_list,$symbol,$name_type);
  $found = 0; 
  $moleculeName = lc $moleculeName;

  my $db = DBI->connect("DBI:Oracle:" . "cgprod",
        "web", "readonly");
  if (not $db or $db->err()) {
    print STDERR "#! Cannot connect to " . $db_user . "@" .
        $db_inst . "\n";
    die;
  }

# Create molecule_list from moleculeName
my $molecule_list = '';
my @molList = split(/[,]/,$moleculeName);
for (my $i = 0; $i < ($#molList); $i++) {
        $molecule_list = $molecule_list . "'" . $molList[$i] . "'" . ",";
}
$molecule_list = $molecule_list . "'" . $molList[$#molList] . "'";
  
# Get official symbol

     # Get official symbol
  $sql = qq!
select
  distinct mol_name, name_type
from $schema.pw_mol_srch ms, $schema.pw_mol_name mn
where lower(map_name) in ($molecule_list)
and ms.mol_id = mn.mol_id
  !;

  $stm = $db->prepare($sql);
  if(not $stm) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "prepare call failed\n";
    die;
  }

  if(!$stm->execute()) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "prepare call failed\n";
    die;
  }

  $list = '';
  my $first = 1;
  while (($symbol, $name_type) =
      $stm->fetchrow_array()) {
        if ($name_type eq "OF") {
          if ($first == 1) {
            $list = $symbol;
            $first = 0;
          } else { 
            $list = $list . ', ' . $symbol;
          }
        }
  }
$stm->finish();
$db->disconnect;
  push @lines, "<ul class=\"xoxo viewing-option\">";
  push @lines, "<li><span class=\"norm\">View graphic as:</span>";
  push @lines, "<ul class=\"links\">";
  push @lines, "<li><a class=\"button-style\" href=\"/search/molecule_landing.shtml" .
          "?molecule=$moleculeName" . "&source_id=$source_id" . "&what=graphic&jpg=on&family_uses=on&complex_uses=on\">JPG</a></li>" ;

  push @lines, "<li><a class=\"button-style\" href=\"/search/molecule_landing.shtml" .
       "?molecule=$moleculeName" . "&source_id=$source_id" . "&what=graphic&svg=on&family_uses=on&complex_uses=on\">SVG</a></li>";
 
  push @lines, "</ul>";
  push @lines, "</li>";

  push @lines, "<li><span class=\"norm\">Save code as:</span>";
  push @lines, "<ul class=\"links\">";
  push @lines, "<li><a class=\"button-style\" href=\"/search/molecule_landing.shtml" .
          "?molecule=$moleculeName" . "&source_id=$source_id" . "&what=text&xml=on&family_uses=on&complex_uses=on\">XML</a></li>" ;

  push @lines, "<li><a class=\"button-style\" href=\"/search/molecule_landing.shtml" .
       "?molecule=$moleculeName" . "&source_id=$source_id" . "&what=graphic&biopax=on&family_uses=on&complex_uses=on\">BioPAX</a></li>";

  push @lines, "</ul>";
  push @lines, "</li>";

  push @lines, "</ul>";
  
   
  return join("\n", @lines) . "\n";
 
}

sub NetworkMolListHeader_1 {
  my ($base, $molecule, $source) = @_;
  my (@lines);
  my $sql;
  my $stm;
  my %molhash = ();
  my $n = 0;

 my $scancheck = Scan($molecule) + Scan($source);
  if ($scancheck > 0) {
    push @lines, "Bad data in parameters\n";
    return join("\n", @lines) . "\n";

  }

  my $db = DBI->connect("DBI:Oracle:" . $db_inst,
        $db_user, $db_pass);
  if (not $db or $db->err()) {
    print STDERR "#! Cannot connect to " . $db_user . "@" .
        $db_inst . "\n";
    die;
  }

  push @lines, "<table class=data border=0 cellpadding=0 cellspacing=0<thead>";
  push @lines, "<th>Molecules</th><th>Uses in complexes (duplicate names indicate differences in component properties)</th>";
  push @lines, "</thead><tbody>";
  my @temp = split(/,/,$molecule); 
  my $list;
  for(my $i = 0; $i < ($#temp + 1); $i += ORACLE_LIST_LIMIT) {
   $list = "'" . join("','", @temp[$i..@temp-1]) . "'";
  }

  
  # $list =~s/[\s]+//g;
  $list =~s/, /,/g;
  $list =~s/' /'/g;
 
  $sql = qq!
    select distinct mol_id
    from
      $schema.pw_mol_srch s1
      where
      s1.map_name in ($list)
  !;
  $sql = qq!
  select
  distinct mm_inner_family.mol_id_2
from
  $schema.pw_mol_srch ms, 
  $schema.pw_edge e1,
  $schema.pw_mol_mol mm_outer_family,
  $schema.pw_mol_mol mm_inner_family,
  $schema.pw_mol_mol mm_complex
where
      mm_inner_family.mol_id_1 = mm_complex.mol_id_1
  and mm_complex.mol_id_2 = mm_outer_family.mol_id_2
  and mm_complex.relation in ('c','i')
  and mm_inner_family.mol_id_2 = e1.mol_id
  and mm_outer_family.relation in ('s','m','i')
  and mm_inner_family.relation in ('s','m','i')
  and ms.mol_id = mm_outer_family.mol_id_1
  and ms.map_name in ($list)
!;
  
  print STDERR "SQL--MolList: $sql\n";
  print STDERR "List:  $list\n";
  $stm = $db->prepare($sql);
  if(not $stm) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "prepare call failed\n";
    die;
  }

  if(!$stm->execute()) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "execute call failed\n";
    die;
  }
  my $molIdList = '';
  my $first = 1;
  my $mol_id;
  my %molIdhash;

  while ($mol_id = $stm->fetchrow_array()) {
    if ($first == 1) {
      $molIdList = $mol_id;
      $first = 0;
    } else {
      if ($molIdhash{$mol_id} == 1) {
      } else { 
        $molIdList = $molIdList . ',' . $mol_id;
      }
    }
  }
  $stm->finish();
 
  $sql = qq!
    select
  distinct e1.atom_id
from
  $schema.pw_edge e1,
  $schema.pw_mol_mol mm_outer_family,
  $schema.pw_mol_mol mm_inner_family,
  $schema.pw_mol_mol mm_complex,
  $schema.pw_atom a
where
      mm_inner_family.mol_id_2 = mm_complex.mol_id_2
  and mm_complex.mol_id_1 = mm_outer_family.mol_id_1
  and e1.mol_id = mm_outer_family.mol_id_1
  and mm_complex.relation in ('c','i')
  and mm_outer_family.relation in ('s','m','i')
  and mm_inner_family.relation in ('s','m','i')
  and mm_outer_family.mol_id_1 in ($molIdList)
  and e1.atom_id = a.atom_id
  and a.atom_source_id in ($source)
  !;

  print STDERR "SQL2:  $sql\n";
  $stm = $db->prepare($sql);
  if(not $stm) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "prepare call failed\n";
    die;
  }

  if(!$stm->execute()) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "execute call failed\n";
    die;
  }
  
  my $atom_id;
  
  while (($atom_id) =
      $stm->fetchrow_array()) {
    my $sql2 = qq!
    select distinct mol_name, e.mol_id
    from 
    $schema.pw_edge e,
    $schema.pw_mol_name mn,
    $schema.pw_atom a,
    $schema.pw_mol_mol mm
    where 
    e.mol_id = mn.mol_id
    and e.atom_id = $atom_id
    and mn.mol_id = mm.mol_id_1
    and mm.mol_id_1 not in (
    select mol_id_2
    from $schema.pw_mol_mol
    where
    relation='c')
    order by mol_name

    !;

    print STDERR "SQL:  $sql2\n";

    my $stm2 = $db->prepare($sql2);

    if(not $stm2) {
      print STDERR "$sql2\n";
      print STDERR "$DBI::errstr\n";
      print STDERR "prepare call failed\n";
      die;
    }

    if(!$stm2->execute()) {
      print STDERR "$sql2\n";
      print STDERR "$DBI::errstr\n";
      print STDERR "execute call failed\n";
      die;
    }
 
    my $mol_name; 
    my $mol_id; 
    while (($mol_name, $mol_id) =
      $stm2->fetchrow_array()) {

      # Query for complexes 
      my $sql3 = qq!
      select distinct mol_id_2, mol_name
      from $schema.pw_mol_mol mm,
      $schema.pw_mol_name mn,
      $schema.pw_pathway_atom pa,
      $schema.pw_edge e
      where mol_id_1 = $mol_id
      and relation = 'c'
      and mm.mol_id_2 = mn.mol_id
      and e.mol_id = mm.mol_id_2
      and e.atom_id = $atom_id
      !;

    my $stm3 = $db->prepare($sql3);
    print STDERR "sql3:  $sql3\n";
    if(not $stm3) {
      print STDERR "$sql3\n";
      print STDERR "$DBI::errstr\n";
      print STDERR "prepare call failed\n";
      die;
    }

    if(!$stm3->execute()) {
      print STDERR "$sql3\n";
      print STDERR "$DBI::errstr\n";
      print STDERR "execute call failed\n";
      die;
    }
    my $mol_id2;
    my $molName;
    my $first = 1;
    my $complexList = '';

    while (($mol_id2, $molName) =
      $stm3->fetchrow_array()) {
      if ($first == 1) {
        $first = 0;
        $complexList = "<a href=\"MoleculePage?molid=$mol_id2\">$molName</a>";
      } else {
        $complexList = $complexList . ',' . "<a href=\"MoleculePage?molid=$mol_id2\">$molName</a>";
      }
    }
    $stm3->finish();
    $mol_name = uc $mol_name;
    if ($molhash{$mol_name} == 1) {
   
    } else {
      $n++;
      push @lines, ($n %2 == 0) ? "<tr>" : "<tr class=\"odd\">";
      push @lines, "<td><a href=\"MoleculePage?molid=$mol_id\">$mol_name</a></td><td>$complexList</a></td></tr>";
      $molhash{$mol_name} = 1; 
    }
  }
  $stm2->finish();
  }

  $stm->finish();

  push @lines, "</tbody></table>";
  $db->disconnect(); 
  return join("\n", @lines) . "\n";


}

sub MolListHeader_1 {
  my ($base, $pid, $format) = @_;

  $BASE = $base;

  my @ALPHA_MONTHS =
      ("Jan", "Feb", "Mar", "Apr",
       "May", "Jun", "Jul", "Aug",
       "Sep", "Oct", "Nov", "Dec");

  my (@lines);

  my $db = DBI->connect("DBI:Oracle:" . $db_inst,
        $db_user, $db_pass);
  if (not $db or $db->err()) {
    print STDERR "#! Cannot connect to " . $db_user . "@" .
        $db_inst . "\n";
    die;
  }
  my ($complex_id, $complex, $relation, $symbol, $mol_id); 
  my ($sql, $stm, $sqlc, $stmc);
  my ($source, $pname, $ext_id, $curator, $role, $last_updated);
  my (@curators, @reviewers);
  my ($biocarta_url);

  $sql = "select " .
      "s.source_name, p.pathway_name, p.ext_pathway_id, p.pathway_id, " .
      "p.last_updated " .
      "from $schema.pw_source s, $schema.pw_pathway p " .
      "where p.pathway_id = $pid and s.source_id = p.pathway_source_id";

  $stm = $db->prepare($sql);
  if(not $stm) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "prepare call failed\n";
    die;
  }
  if(!$stm->execute()) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "execute call failed\n";
    die;
  }
  ($source, $pname, $ext_id, $pid, $last_updated) = $stm->fetchrow_array();
  $stm->finish();

  my ($sec, $min, $hr, $mday, $mon, $year, $wday, $yday, $isdst) =
      localtime($last_updated);
  my $mon_alpha = $ALPHA_MONTHS[$mon];
  $year = $year + 1900;

  my $nice_source = $source;

  if ($source eq "BioCarta") {
    $sql = "select " .
      "b.bc_ext_id from pid.pw_biocarta b " .
      "where b.pid_ext_id = '$ext_id'";

  $stm = $db->prepare($sql);
  if(not $stm) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "prepare call failed\n";
    die;
  }
  if(!$stm->execute()) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "execute call failed\n";
    die;
  }
  ($biocarta_url) = $stm->fetchrow_array();
  $biocarta_url = $biocarta_url . ".asp";
  $nice_source = "<a href=\"http://www.biocarta.com/pathfiles/$biocarta_url\">BioCarta Imported</a>";
  $stm->finish();

  }
  $nice_source =~ s/NATURE/NCI-Nature curated/;
  # $nice_source =~ s/BioCarta/<a href="http:\/\/www.biocarta.com/pathfiles/$biocarta_url">BioCarta Imported<\/a>/;
  $nice_source =~ s/Reactome/Reactome Imported/;
  undef @curators;
  undef @reviewers;
 if (($source eq "NATURE") || ($source eq "Reactome")) {
    $sqlc = "select curator, role " .
            "from $schema.pw_curators " .
            "where pathway_id = $pid " ;
    $stmc = $db->prepare($sqlc);
    if (not $stmc) {
      print STDERR "$sqlc\n";
      print STDERR "$DBI::errstr\n";
      print STDERR "prepare call failed\n";
      die;
    }
    if (!$stmc->execute()) {
      print STDERR "$sqlc\n";
      print STDERR "$DBI::errstr\n";
      print STDERR "execute call failed\n";
      die;
    }
    while (($curator, $role) = $stmc->fetchrow_array()) {
      if ($role eq "C") {
        push @curators, $curator;
      } else {
        push @reviewers, $curator;
      }
    }
    $stmc->finish();
  }
  
  my $other_format_line;
  ## for some reason, can't seem to get the "format"  variable
  ## through the sequence of invocations, so for now, just list both
  ## graphic formats
  $other_format_line =

      "<ul class=\"xoxo links\"><li><a type=button href=\"/search/pathway_landing.shtml" .
      "?pathway_id=$pid" . "&source=$source" . "&what=graphic&svg=on&ppage=1\">SVG</a></li>" .
      "&nbsp;&nbsp;<li><a href=\"/search/pathway_landing.shtml" .
      "?pathway_id=$pid" . "&source=$source" . "&what=graphic&jpg=on&ppage=1\">JPG</a></li>" .
      "&nbsp;&nbsp;<li><a href=\"/search/pathway_landing.shtml" .
 push @lines, "<h1 class=\"pagetitle\">$pname</h1>";
  push @lines, "<hr />";

  push @lines, "<p class=\"citation\"><span class=\"title\">Revision date:</span> $mday" . "-" . $mon_alpha .
      "-" . $year . "</p>";
  push @lines, "<p class=\"citation\"><span class=\"title\">Source:</span> $nice_source</p>";
  if (($source eq "NATURE") || ($source eq "Reactome")) {
    push @lines, "<p class=\"citation\"><span class=\"title\">Curated by:</span> " . join(", ", @curators). "</p>";
    push @lines, "<p class=\"citation\"><span class=\"title\">Reviewed by:</span> " . join(", ", @reviewers). "</p>";
  }
  push @lines, "<p class=\"citation\"><span class=\"title\">Pathway ID:</span> $ext_id</p>";

  $sql = qq!
    select distinct mol_name, mn.mol_id
    from
      $schema.pw_mol_srch s1,
      $schema.pw_pathway_atom pa,
      $schema.pw_edge e,
      $schema.pw_mol_mol mm,
      $schema.pw_mol_name mn
      where
      s1.mol_id = e.mol_id
      and mn.mol_id = mm.mol_id_1
      and s1.mol_id = mm.mol_id_1
      and e.atom_id = pa.atom_id
      and pa.pathway_id = $pid
      and mm.mol_id_1 not in (
      select mol_id_2
      from $schema.pw_mol_mol
      where
      relation='c')
!;
 
  $stm = $db->prepare($sql);
  if(not $stm) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "prepare call failed\n";
    die;
  }

  if(!$stm->execute()) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "execute call failed\n";
    die;
  }
  my $symbol;
  my $mol_id;
  push @lines, "<p><a href=\"ftp://ftp1.nci.nih.gov/pub/PID/molList/$pid.mol.tab.gz\">Export molecule list</a><p>";

  push @lines, "<table class=data border=0 cellpadding=0 cellspacing=0<thead>";
  push @lines, "<th>Molecules</th><th>Uses in complexes (duplicate names indicate differences in component properties)</th>";
  push @lines, "</thead><tbody>";
  my $n = 0;
  while (($symbol, $mol_id) =
      $stm->fetchrow_array()) {
    # Get complexes 
    my $sql2 = qq!
      select distinct mol_id_2, mol_name
      from $schema.pw_mol_mol mm, 
      $schema.pw_mol_name mn,
      $schema.pw_pathway_atom pa,
      $schema.pw_edge e
      where mol_id_1 = $mol_id
      and relation = 'c'   
      and mm.mol_id_2 = mn.mol_id 
      and e.atom_id = pa.atom_id
      and e.mol_id = mm.mol_id_2
      and pa.pathway_id = $pid

    !;
    my $stm2 = $db->prepare($sql2);
    print STDERR "sql2:  $sql2\n";
    if(not $stm2) {
      print STDERR "$sql2\n";
      print STDERR "$DBI::errstr\n";
      print STDERR "prepare call failed\n";
      die;
    }

    if(!$stm2->execute()) {
      print STDERR "$sql2\n";
      print STDERR "$DBI::errstr\n";
      print STDERR "execute call failed\n";
      die;
    }
    my $mol_id2;
    my $molName;
    my $complexList = '';
    my $first = 1;
    while (($mol_id2, $molName) = $stm2->fetchrow_array()) {
      if ($first == 1) {
        $first = 0;
        $complexList = "<a href=\"MoleculePage?molid=$mol_id2\">$molName</a>";
      } else {
        $complexList = $complexList . ',' . "<a href=\"MoleculePage?molid=$mol_id2\">$molName</a>"; 
      } 
    }
    $stm2->finish();
    $n++;
    push @lines, ($n %2 == 0) ? "<tr>" : "<tr class=\"odd\">";
    push @lines, "<td><a href=\"MoleculePage?molid=$mol_id\">$symbol</a></td><td>$complexList</a></td></tr>";

  }
  $stm->finish();
  $db->disconnect();
  push @lines, "</tbody></table>";
  return join("\n", @lines) . "\n";

}

sub NetworkCitationHeader_1 {
  my ($base, $moleculeName, $format, $source) = @_;

  $BASE = $base;
 
  my (@lines);
  my $scancheck = Scan($moleculeName) + Scan($format) + Scan($source);
  if ($scancheck > 0) {
    push @lines, "Bad data in parameters\n";
    return join("\n", @lines) . "\n";

  }
 
   my $db = DBI->connect("DBI:Oracle:" . $db_inst,
        $db_user, $db_pass);
  if (not $db or $db->err()) {
    print STDERR "#! Cannot connect to " . $db_user . "@" .
        $db_inst . "\n";
    die;
  }

  my ($sql, $stm, $sqlc, $stmc);
  my ($source, $pname, $ext_id, $curator, $role, $last_updated);
  my (@curators, @reviewers);
  my ($biocarta_url);

  $sql = qq! 

  !;

  $stm = $db->prepare($sql);
  if(not $stm) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "prepare call failed\n";
    die;
  }
  if(!$stm->execute()) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "execute call failed\n";
    die;
  }
  # ($source, $pname, $ext_id, $pid, $last_updated) = $stm->fetchrow_array();
  $stm->finish();

}

sub CitationHeader_1 {
  my ($base, $pid, $format) = @_;

  $BASE = $base;

  my @ALPHA_MONTHS =
      ("Jan", "Feb", "Mar", "Apr",
       "May", "Jun", "Jul", "Aug",
       "Sep", "Oct", "Nov", "Dec");

  my (@lines);

  my $db = DBI->connect("DBI:Oracle:" . $db_inst,
        $db_user, $db_pass);
  if (not $db or $db->err()) {
    print STDERR "#! Cannot connect to " . $db_user . "@" .
        $db_inst . "\n";
    die;
  }

  my ($sql, $stm, $sqlc, $stmc);
  my ($source, $pname, $ext_id, $curator, $role, $last_updated);
  my (@curators, @reviewers);
  my ($biocarta_url);

  $sql = "select " .
      "s.source_name, p.pathway_name, p.ext_pathway_id, p.pathway_id, " .
      "p.last_updated " .
      "from $schema.pw_source s, $schema.pw_pathway p " .
      "where p.pathway_id = $pid and s.source_id = p.pathway_source_id";

  $stm = $db->prepare($sql);
  if(not $stm) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "prepare call failed\n";
    die;
  }
  if(!$stm->execute()) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "execute call failed\n";
    die;
  }
  ($source, $pname, $ext_id, $pid, $last_updated) = $stm->fetchrow_array();
  $stm->finish();

  my ($sec, $min, $hr, $mday, $mon, $year, $wday, $yday, $isdst) =
      localtime($last_updated);
  my $mon_alpha = $ALPHA_MONTHS[$mon];
  $year = $year + 1900;

  my $nice_source = $source;

  if ($source eq "BioCarta") {
    $sql = "select " .
      "b.bc_ext_id from pid.pw_biocarta b " .
      "where b.pid_ext_id = '$ext_id'";

  $stm = $db->prepare($sql);
  if(not $stm) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "prepare call failed\n";
    die;
  }
  if(!$stm->execute()) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "execute call failed\n";
    die;
  }
  ($biocarta_url) = $stm->fetchrow_array();
  $biocarta_url = $biocarta_url . ".asp";
  $nice_source = "<a href=\"http://www.biocarta.com/pathfiles/$biocarta_url\">BioCarta Imported</a>";
  $stm->finish();

  }
  $nice_source =~ s/NATURE/NCI-Nature curated/;
  $nice_source =~ s/Reactome/Reactome Imported/;
  undef @curators;
  undef @reviewers;
  if (($source eq "NATURE") || ($source eq "Reactome")) {
    $sqlc = "select curator, role " .
            "from $schema.pw_curators " .
            "where pathway_id = $pid " ;
    $stmc = $db->prepare($sqlc);
    if (not $stmc) {
      print STDERR "$sqlc\n";
      print STDERR "$DBI::errstr\n";
      print STDERR "prepare call failed\n";
      die;
    }
    if (!$stmc->execute()) {
      print STDERR "$sqlc\n";
      print STDERR "$DBI::errstr\n";
      print STDERR "execute call failed\n";
      die;
    }
    while (($curator, $role) = $stmc->fetchrow_array()) {
      if ($role eq "C") {
        push @curators, $curator;
      } else {
        push @reviewers, $curator;
      }
    }
  }

  my $other_format_line;
  ## for some reason, can't seem to get the "format"  variable
  ## through the sequence of invocations, so for now, just list both
  ## graphic formats
  $other_format_line =

      "<ul class=\"xoxo links\"><li><a type=button href=\"/search/pathway_landing.shtml" .
      "?pathway_id=$pid" . "&source=$source" . "&what=graphic&svg=on&ppage=1\">SVG</a></li>" .
      "&nbsp;&nbsp;<li><a href=\"/search/pathway_landing.shtml" .
      "?pathway_id=$pid" . "&source=$source" . "&what=graphic&jpg=on&ppage=1\">JPG</a></li>" .
      "&nbsp;&nbsp;<li><a href=\"/search/pathway_landing.shtml" .
 push @lines, "<h1 class=\"pagetitle\">$pname</h1>";
  push @lines, "<hr />";

  push @lines, "<p class=\"citation\"><span class=\"title\">Revision date:</span> $mday" . "-" . $mon_alpha .
      "-" . $year . "</p>";
  push @lines, "<p class=\"citation\"><span class=\"title\">Source:</span> $nice_source</p>";
  if (($source eq "NATURE") || ($source eq "Reactome")) {
    push @lines, "<p class=\"citation\"><span class=\"title\">Curated by:</span> " . join(", ", @curators). "</p>";
    push @lines, "<p class=\"citation\"><span class=\"title\">Reviewed by:</span> " . join(", ", @reviewers). "</p>";
  }
  push @lines, "<p class=\"citation\"><span class=\"title\">Pathway ID:</span> $ext_id</p>";

  $sql = "select distinct a.pmid, pm_text from $schema.pw_pubmed a, $schema.pw_references b, $schema.pw_pathway_atom c";
  $sql = $sql . " where c.pathway_id = $pid";
  $sql = $sql . " and a.pmid = b.pmid ";
  $sql = $sql . " and b.atom_id = c.atom_id "; 
 
  $stm = $db->prepare($sql);
  if(not $stm) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "prepare call failed\n";
    die;
  }
 
  if(!$stm->execute()) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "execute call failed\n";
    die;
  }
  my $pmid;
  my $reference; 
  my $pm_text;
  push @lines, "<p><a href=\"ftp://ftp1.nci.nih.gov/pub/PID/references/$pid.csv.gz\">Export references</a><p>";
  push @lines, "<table class=data border=0 cellpadding=0 cellspacing=0><thead>";
  push @lines, "<th>References</th><th>PubMed id</th>";
  push @lines, "</thead><tbody>";
  my $n = 0;
  while (($pmid, $pm_text) =
      $stm->fetchrow_array()) {
    if ($pmid) {
      $reference = ("<td>$pm_text</td><td><a href=\"" . PUBMED_URL($pmid) . "\" > PMID:$pmid</a></td>");
      $n++;
      push @lines, ($n %2 == 0) ? "<tr>" : "<tr class=\"odd\">";
      push @lines, "$reference</tr>";
    }
  }

  $db->disconnect();
  push @lines, "</tbody></table>";
  return join("\n", @lines) . "\n";

}
#################################################
sub XamlHeader_1 {
  my ($base, $pid, $genes_a, $genes_b, $format) = @_;

  $BASE = $base;

  my @ALPHA_MONTHS =
      ("Jan", "Feb", "Mar", "Apr",
       "May", "Jun", "Jul", "Aug",
       "Sep", "Oct", "Nov", "Dec");

  my (@lines);
  my $scancheck = Scan($base) + Scan($pid) + Scan($genes_a) + Scan($genes_b) + Scan($format);
  if ($scancheck > 0) {
    push @lines, "Bad data in parameters\n";
    return join("\n", @lines) . "\n";

  }

  my $db = DBI->connect("DBI:Oracle:" . $db_inst,
        $db_user, $db_pass);
  if (not $db or $db->err()) {
    print STDERR "#! Cannot connect to " . $db_user . "@" .
        $db_inst . "\n";
    die;
  }

  my ($sql, $stm, $sqlc, $stmc);
  my ($source, $pname, $ext_id, $curator, $role, $last_updated);
  my (@curators, @reviewers);
  my ($biocarta_url,$parentName,$parentPid);

  $parentName = 'None';
  if (($pid > 99999) && ($pid < 600000)) {
    $sql = "select " .
      "s.source_name, p.pathway_name, p.ext_pathway_id, p.pathway_id, " .
      "p.last_updated " .
      "from $schema.pw_source s, $schema.pw_pathway p " .
      "where p.pathway_id = $pid and s.source_id = p.pathway_source_id";
  } else {
    $sql = "select " .
      "s.source_name, p.pathway_name, p.ext_pathway_id, p.pathway_id, " .
      "p.last_updated " .
      "from $schema.pw_source s, $schema.pw_pathway p " .
      "where p.ext_pathway_id = '$pid' and s.source_id = p.pathway_source_id";
  }


  $stm = $db->prepare($sql);
  if(not $stm) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "prepare call failed\n";
    die;
  }
  if(!$stm->execute()) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "execute call failed\n";
    die;
  }
  ($source, $pname, $ext_id, $pid, $last_updated) = $stm->fetchrow_array();
  $stm->finish();

  my ($sec, $min, $hr, $mday, $mon, $year, $wday, $yday, $isdst) =
      localtime($last_updated);
  my $mon_alpha = $ALPHA_MONTHS[$mon];
  $year = $year + 1900;

  my $nice_source = $source;

  if ($source eq "BioCarta") {
    $sql = "select " .
      "b.bc_ext_id from pid.pw_biocarta b " .
      "where b.pid_ext_id = '$ext_id'";
  $stm = $db->prepare($sql);
  if(not $stm) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "prepare call failed\n";
    die;
  }
  if(!$stm->execute()) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "execute call failed\n";
    die;
  }
  ($biocarta_url) = $stm->fetchrow_array();
  $biocarta_url = $biocarta_url . ".asp";
  $nice_source = "<a href=\"http://www.biocarta.com/pathfiles/$biocarta_url\">BioCarta Imported</a>";
  $stm->finish();

  }
  $nice_source =~ s/NATURE/NCI-Nature curated/;
  $nice_source =~ s/Reactome/Reactome Imported/;

  undef @curators;
  undef @reviewers;

  my $sqlsub = qq!

select distinct
  p2.pathway_id,
  p2.ext_pathway_id,
  p2.pathway_name
from
  $schema.pw_pathway p1,
  $schema.pw_pathway p2,
  $schema.pw_abstraction a,
  $schema.pw_pathway_atom pa,
  $schema.pw_subnet s
where
  pa.pathway_id = p2.pathway_id
  and pa.atom_id = a.atom_id
  and p1.pathway_source_id = p2.pathway_source_id
  and (a.pathway_id = p1.pathway_id or a.ext_pathway_id = p1.ext_pathway_id)
  and p1.pathway_id = $pid
  and p1.subnet = 'Y'
  and s.subnet_pathway_id = p1.pathway_id
  and p2.pathway_id = s.parent_pathway_id

!;
   my $stmsub = $db->prepare($sqlsub);
    if (not $stmsub) {
      print STDERR "$sqlsub\n";
      print STDERR "$DBI::errstr\n";
      print STDERR "prepare call failed\n";
      die;
    }
    if (!$stmsub->execute()) {
      print STDERR "$sqlsub\n";
      print STDERR "$DBI::errstr\n";
      print STDERR "execute call failed\n";
      die;
    }
    my $pid1;
    my $pid2;
    my $pidname;
    while (($pid1, $pid2, $pidname) = $stmsub->fetchrow_array()) {
      $parentName = $pidname;
      $parentPid = $pid1;
    }
    $stmsub->finish();

  if (($source eq "NATURE") || ($source eq "Reactome")) {
    $sqlc = "select curator, role " .
            "from $schema.pw_curators " .
            "where pathway_id = $pid " ;
    $stmc = $db->prepare($sqlc);
    if (not $stmc) {
      print STDERR "$sqlc\n";
      print STDERR "$DBI::errstr\n";
      print STDERR "prepare call failed\n";
      die;
    }
    if (!$stmc->execute()) {
      print STDERR "$sqlc\n";
      print STDERR "$DBI::errstr\n";
      print STDERR "execute call failed\n";
      die;
    }
    while (($curator, $role) = $stmc->fetchrow_array()) {
      if ($role eq "C") {
        push @curators, $curator;
      } else {
        push @reviewers, $curator;
      }
    }
    $stmc->finish();
  }


  push @lines, "<h1 class=\"pagetitle\">$pname</h1>";
  push @lines, "<hr />";

  push @lines, "<p class=\"citation\"><span class=\"title\">Revision date:</span> $mday" . "-" . $mon_alpha .
      "-" . $year ."</p>";

  if (($source eq "NATURE") || ($source eq "Reactome")) {
    push @lines, "<p class=\"citation\"><span class=\"title\">Curated by:</span> " . join(", ", @curators) . "</p>";
    push @lines, "<p class=\"citation\"><span class=\"title\">Reviewed by:</span> " . join(", ", @reviewers) . "</p>";
  }


  push @lines, "<p class=\"citation\"><span class=\"title\">Pathway ID:</span> $ext_id </p>";
  if ($parentName ne 'None') {
    if ($format eq 'jpg') {
      push @lines, "<p class=\"citation\"><span class=\"title\">Parent pathway:  </span>  <a href=\"/search/pathway_landing.shtml?pathway_id=$parentPid&source=$source&what=graphic&jpg=on&ppage=1\">$parentName</a></p><br>";
    }
    if ($format eq 'svg') {
      push @lines, "<p class=\"citation\"><span class=\"title\">Parent pathway:  </span><a href=\"/search/pathway_landing.shtml?pathway_id=$parentPid&source=$source&what=graphic&svg=on&ppage=1\">$parentName</a></p>";
      push @lines, "<p><br></p>";
    }
  }

  if ($source eq "NATURE") {
     my @category_name;
     my $category_count = 0;
     my $first = 1;
     my $category_url = $BASE . "/browse_categories.shtml?CMD=open&NODE=$pid&GOIDS=";
     my $parent_id = $pid;

       $sqlc = "select pathway_parent, parent_ext_id " .
            "from $schema.pw_pathway_parent " .
            "where pathway_id = $parent_id " ;
       $stmc = $db->prepare($sqlc);
       if (not $stmc) {
         print STDERR "$sqlc\n";
         print STDERR "$DBI::errstr\n";
         print STDERR "prepare call failed\n";
         die;
       }
       if (!$stmc->execute()) {
         print STDERR "$sqlc\n";
         print STDERR "$DBI::errstr\n";
         print STDERR "execute call failed\n";
         die;
       }
       my $find_parent = 0;
       print STDERR "Before while\n";
       while ((my $parent, my $parent_name) = $stmc->fetchrow_array()) {
         $category_url = $category_url . "$parent" . ",";
         $parent_id = $parent;
         if ($parent_id < 200000) {
           $category_name[$category_count] = $parent_name;
           $category_count++;
         }
         if ($category_count == 0) {
           $find_parent = $parent_id;
         }
         if ($parent_id > 4) {
           my $new_parents = GetParents($parent_id);
           $category_url = $category_url . "$new_parents" . ",";
         }
       }
       $stmc->finish();
     $category_url = $category_url . 1;
     if ($category_count == 0 ) {
       print STDERR "Before GetParentReturn\n";
       $category_name[0] = GetParentName($find_parent);
       $category_count++;
     }
     my $category_string = "<p class=\"citation\"><b><span class\"title\">Pathway category:  </span></b>";
     for (my $k = 0; $k < $category_count; $k++) {
       my $cat_name = $category_name[$k];
       $category_string = $category_string . "<a href=\"$category_url\">$cat_name</a>";
       if ($k + 1 == $category_count) {
         $category_string = $category_string . "</p>";
       } else {
         $category_string = $category_string . ", ";
       }
     }
     push @lines, "$category_string";

 }

  if ($parentName ne 'None') {
    if ($format eq 'jpg') {
      push @lines, "<p class=\"citation\"><span class=\"title\">Parent pathway:  </span>  <a href=\"/search/pathway_landing.shtml?pathway_id=$parentPid&source=$source&what=graphic&jpg=on&ppage=1\">$parentName</a></p><br>";
    }
    if ($format eq 'svg') {
      push @lines, "<p class=\"citation\"><span class=\"title\">Parent pathway:  </span><a href=\"/search/pathway_landing.shtml?pathway_id=$parentPid&source=$source&what=graphic&svg=on&ppage=1\">$parentName</a></p>";
      push @lines, "<p><br>";
    }
  }
	 
  if ($source eq "NATURE") {
    push @lines, "<p class=\"citation\"><a href=\"mailto:nci-pid\@nature.com\">Submit feedback for this pathway</a></p>";
  }

  push @lines, "<ul class=\"xoxo viewing-option\">";
  push @lines, "<li><span class=\"hidden\">Additional citation information</span>";
  push @lines, "<ul class=\"links\">";
  push @lines, "<li><a class=\"button-style\" href=\"/search/mol_list_landing.shtml" .
          "?pathway_id=$pid" .  "&what=text\">Molecule list</a></li>";
  push @lines, "<li><a class=\"button-style\" href=\"/search/citation_landing.shtml" .
          "?pathway_id=$pid" .  "&what=text\">References</a></li>";

  push @lines, "</ul>"; 
  push @lines, "</li>";
  push @lines, "<li><span class=\"norm\">View graphic as:</span>";
  push @lines, "<ul class=\"links\">";
  push @lines, "<li><a class=\"button-style\" href=\"/search/pathway_landing.shtml" .
          "?pathway_id=$pid" . "&source=$source" . "&genes_a=$genes_a" . "&genes_b=$genes_b" . "&what=graphic&jpg=on&ppage=1\">JPG</a></li>" ;

  push @lines, "<li><a class=\"button-style\" href=\"/search/pathway_landing.shtml" .
       "?pathway_id=$pid" . "&source=$source" . "&genes_a=$genes_a" . "&genes_b=$genes_b" . "&what=graphic&svg=on&ppage=1\">SVG</a></li>";

    push @lines, "<li><a class=\"button-style\" href=\"/xaml_landing.shtml" .
       "?pathway_id=$pid" . "&source=$source" . "&what=text&xaml=true\">Silverlight-beta</a></li>";

  push @lines, "<li><div class=\"helpbox\">";
  push @lines, "<a href=\"#\" class=\"help\">Help</a>";
  push @lines, "<div class=\"description\">";
  push @lines, "<h2 class=\"title\">Information on graphics types</h2>";
  push @lines, "<p class=\"norm\">The pathway can be viewed as a JPG, SVG or Silverlight image. All molecules and interactions within the graphics formats are interactive; clicking on them will open information windows.</p>";
  push @lines, "<p class=\"norm\"><B>JPG:</b>  All common web browsers are able to display JPG images.  However, please note that some browsers may not support very large JPG files.</p>";
  push @lines, "<p class=\"norm\"><B>SVG:</b>  The advantage of SVG images is that they can be more easily zoomed and panned on-screen.  If your browser is not completely SVG enabled, please visit the <a href=\"/PID/userguide/output_formats.shtml\">User guide</a> and follow the instructions for installing the SVG plug-in specific to your browser</p>";
  push @lines, "<p class=\"norm\"><b>SVG image navigation:</b>  Use the left mouse button to activate the image.  To <b>zoom</b> in or out, use the right mouse button and select from the options on the pop-up menu.  To <b>pan</b> the image, hold down the alt key and click-and-drag the image with the mouse.</p>";
  push @lines, "<p class=\"norm\"><b>Silverlight:  </b>Please visit http://www.microsoft.com/getsilverlight/Get-Started/Install/Default.aspx to install the plug-in for your browser.</p>";
  push @lines, "<p class=\"norm\"><b>Silverlight image navigation:  Zoom</b> is available at the top-left of the image. To <b>pan</b> merely drag the image with the mouse.  In the image you may search for specific molecules or all molecules at a certain subcellular location.  To search for molecules you may enter their names, Entrez Gene identifiers or UniProt accession numbers.</p>";
  push @lines, "</div>";
  push @lines, "</div>";
  push @lines, "</li>";
  push @lines, "</ul>";
  push @lines, "</li>";

   push @lines, "<li><span class=\"norm\">Save code as:</span>";
  push @lines, "<ul class=\"links\">";
  push @lines, "<li><a class=\"button-style\" href=\"/search/pathway_landing.shtml" .
          "?pathway_id=$pid" . "&pathway_short=$ext_id" . "&source=$source" . "&what=text&xml=on&ppage=1\">XML</a></li>" ;
  push @lines, "<li><a class=\"button-style\" href=\"/search/pathway_landing.shtml" .
       "?pathway_id=$pid" . "&pathway_short=$ext_id" . "&source=$source" . "&what=text&biopax=on&ppage=1\">BioPAX</a></li>";
  push @lines, "<li><div class=\"helpbox\">";
  push @lines, "<a href=\"#\" class=\"help\">Help</a>";
  push @lines, "<div class=\"description\">";
  push @lines, "<h2 class=\"title\">Information on types of pathway code</h2>";
  push @lines, "<p class=\"norm\">You can view the pathway code in either XML or BioPax formats.  Please see the <a href=\"/PID/userguide/output_formats.shtml\">User guide</a> for more information.</p>";
  push @lines, "</div>";
  push @lines, "</div>";
  push @lines, "</li>";

  push @lines, "</ul>";
  push @lines, "</li>";

  push @lines, "</ul>";


  push @lines, "<p class=\"norm\">A <a href=\"/xaml_landing.shtml?pathway_id=$pid" . "&source=$source" . "&xaml=true&what=text/#key\">key to the image</a> icons appears below the pathway.  Click on a biomolecule or interaction for additional information.</p>";
 
  $db->disconnect();
 
  return join("\n", @lines) . "\n";

}

sub GetParents {

  my ($pid) = @_;
  my $parent_return = '';
  my ($sql, $stm, $sqlc, $stmc);

  my $db = DBI->connect("DBI:Oracle:" . $db_inst,
        $db_user, $db_pass);
  if (not $db or $db->err()) {
    print STDERR "#! Cannot connect to " . $db_user . "@" .
        $db_inst . "\n";
    die;
  }


  $sqlc = "select pathway_parent, parent_ext_id " .
            "from $schema.pw_pathway_parent " .
            "where pathway_id = $pid " ;
       $stmc = $db->prepare($sqlc);
       if (not $stmc) {
         print STDERR "$sqlc\n";
         print STDERR "$DBI::errstr\n";
         print STDERR "prepare call failed\n";
         die;
       }
       if (!$stmc->execute()) {
         print STDERR "$sqlc\n";
         print STDERR "$DBI::errstr\n";
         print STDERR "execute call failed\n";
         die;
       }
       my $comma = '';
 
       while ((my $parent, my $parent_name) = $stmc->fetchrow_array()) {
         $parent_return = $parent_return . $comma . $parent; 
         $comma = ',';
         if ($parent > 4) {
           $parent_return = $parent_return . ',' . GetParents($parent); 
         }
       }
 
       $stmc->finish();

       return($parent_return);
}

sub GetParentName {

 my ($pid) = @_;
  my $parent_return = '';
  my ($sqlp, $stmp);

  my $db = DBI->connect("DBI:Oracle:" . $db_inst,
        $db_user, $db_pass);
  if (not $db or $db->err()) {
    print STDERR "#! Cannot connect to " . $db_user . "@" .
        $db_inst . "\n";
    die;
  }

my $sqlp = "select a.ext_pathway_id " .
            "from $schema.pw_category a, $schema.pw_pathway_parent b " .
            "where b.pathway_id = $pid " .
            "and a.pathway_id = b.pathway_parent " ;
           my $stmp = $db->prepare($sqlp);
           if (not $stmp) {
             print STDERR "$sqlp\n";
             print STDERR "$DBI::errstr\n";
             print STDERR "prepare call failed\n";
             die;
           }
           if (!$stmp->execute()) {
             print STDERR "$sqlp\n";
             print STDERR "$DBI::errstr\n";
             print STDERR "execute call failed\n";
             die;
           }
           print STDERR "SQL:  $sqlp\n"; 
           while ((my $t1) = $stmp->fetchrow_array()) {
             print STDERR "In while:  $t1\n";
             $parent_return = $t1 
           }
           $stmp->finish();

return $parent_return;

}
######################################################################
sub PathwayHeader_1 {
  my ($base, $pid, $genes_a, $genes_b, $format) = @_;

  $BASE = $base;

  my @ALPHA_MONTHS =
      ("Jan", "Feb", "Mar", "Apr",
       "May", "Jun", "Jul", "Aug",
       "Sep", "Oct", "Nov", "Dec");

  my (@lines);
  my $scancheck = Scan($base) + Scan($pid) + Scan($genes_a) + Scan($genes_b) + Scan($format);
  if ($scancheck > 0) {
    push @lines, "Bad data in parameters\n";
    return join("\n", @lines) . "\n";

  }

  my $db = DBI->connect("DBI:Oracle:" . $db_inst,
        $db_user, $db_pass);
  if (not $db or $db->err()) {
    print STDERR "#! Cannot connect to " . $db_user . "@" .
        $db_inst . "\n";
    die;
  }

  my ($sql, $stm, $sqlc, $stmc);
  my ($source, $pname, $ext_id, $curator, $role, $last_updated);
  my (@curators, @reviewers);
  my ($biocarta_url,$parentName,$parentPid);

  $parentName = 'None';
  if (($pid > 99999) && ($pid < 600000)) { 
    $sql = "select " .
      "s.source_name, p.pathway_name, p.ext_pathway_id, p.pathway_id, " .
      "p.last_updated " .
      "from $schema.pw_source s, $schema.pw_pathway p " .
      "where p.pathway_id = $pid and s.source_id = p.pathway_source_id";
  } else {
    $sql = "select " .
      "s.source_name, p.pathway_name, p.ext_pathway_id, p.pathway_id, " .
      "p.last_updated " .
      "from $schema.pw_source s, $schema.pw_pathway p " .
      "where p.ext_pathway_id = '$pid' and s.source_id = p.pathway_source_id";
  }


  $stm = $db->prepare($sql);
  if(not $stm) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "prepare call failed\n";
    die;
  }
  if(!$stm->execute()) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "execute call failed\n";
    die;
  }
  ($source, $pname, $ext_id, $pid, $last_updated) = $stm->fetchrow_array();
   $stm->finish();

  my ($sec, $min, $hr, $mday, $mon, $year, $wday, $yday, $isdst) =
      localtime($last_updated);
  my $mon_alpha = $ALPHA_MONTHS[$mon];
  $year = $year + 1900;

  my $nice_source = $source;

  if ($source eq "BioCarta") {
    $sql = "select " .
      "b.bc_ext_id from pid.pw_biocarta b " .
      "where b.pid_ext_id = '$ext_id'";

  $stm = $db->prepare($sql);
  if(not $stm) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "prepare call failed\n";
    die;
  }
  if(!$stm->execute()) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "execute call failed\n";
    die;
  }
  ($biocarta_url) = $stm->fetchrow_array();
  $biocarta_url = $biocarta_url . ".asp";
  $nice_source = "<a href=\"http://www.biocarta.com/pathfiles/$biocarta_url\">BioCarta Imported</a>";
  $stm->finish();
 
  }
  $nice_source =~ s/NATURE/NCI-Nature curated/;
  # $nice_source =~ s/BioCarta/<a href="http:\/\/www.biocarta.com/pathfiles/$biocarta_url">BioCarta Imported<\/a>/;
  $nice_source =~ s/Reactome/Reactome Imported/;

  undef @curators;
  undef @reviewers;

  my $sqlsub = qq!
 
select distinct
  p2.pathway_id,
  p2.ext_pathway_id,
  p2.pathway_name
from
  $schema.pw_pathway p1,
  $schema.pw_pathway p2,
  $schema.pw_subnet s
where
  p1.pathway_source_id = p2.pathway_source_id
  and p1.pathway_id = $pid
  and p1.subnet = 'Y' 
  and s.subnet_pathway_id = p1.pathway_id
  and p2.pathway_id = s.parent_pathway_id

!;
    my $stmsub = $db->prepare($sqlsub);
    if (not $stmsub) {
      print STDERR "$sqlsub\n";
      print STDERR "$DBI::errstr\n";
      print STDERR "prepare call failed\n";
      die;
    }
    if (!$stmsub->execute()) {
      print STDERR "$sqlsub\n";
      print STDERR "$DBI::errstr\n";
      print STDERR "execute call failed\n";
      die;
    }
    my $pid1;
    my $pid2;
    my $pidname;
    while (($pid1, $pid2, $pidname) = $stmsub->fetchrow_array()) {
      $parentName = $pidname; 
      $parentPid = $pid1;
    }
    $stmsub->finish();
  
  if (($source eq "NATURE") || ($source eq "Reactome")) {
    $sqlc = "select curator, role " .
            "from $schema.pw_curators " .
            "where pathway_id = $pid " ;
    $stmc = $db->prepare($sqlc);
    if (not $stmc) {
      print STDERR "$sqlc\n";
      print STDERR "$DBI::errstr\n";
      print STDERR "prepare call failed\n";
      die;
    }
    if (!$stmc->execute()) {
      print STDERR "$sqlc\n";
      print STDERR "$DBI::errstr\n";
      print STDERR "execute call failed\n";
      die;
    }
    while (($curator, $role) = $stmc->fetchrow_array()) {
      if ($role eq "C") {
        push @curators, $curator; 
      } else {
        push @reviewers, $curator; 
      }
    }
    $stmc->finish();
  }

  my $other_format_line;
  ## for some reason, can't seem to get the "format"  variable
  ## through the sequence of invocations, so for now, just list both
  ## graphic formats
  $other_format_line =
      
      "<ul class=\"xoxo links\"><li><a type=button href=\"/search/pathway_landing.shtml" .
      "?pathway_id=$pid" . "&source=$source" . "&what=graphic&jpg=on&ppage=1\">JPG</a></li>" .
      "&nbsp;&nbsp;<li><a href=\"/search/pathway_landing.shtml" .
      "?pathway_id=$pid" . "&source=$source" . "&what=graphic&svg=on&ppage=1\">SVG</a></li>" .
      "&nbsp;&nbsp;<li><a href=\"/search/pathway_landing.shtml" .
          "?pathway_id=$pid" . "&pathway_short=$ext_id" . "&source=$source" . "&what=text&xml=on\">XML</a></li>" .
      "&nbsp;&nbsp;<li><a href=\"/search/pathway_landing.shtml" .
          "?pathway_id=$pid" . "&pathway_short=$ext_id" . "&source=$source" . "&what=text&biopax=on\">BioPAX</a></li>" .
      "&nbsp;&nbsp;<li><a href=\"/search/mol_list_landing.shtml" .
          "?pathway_id=$pid" .  "&what=text\">Molecule List</a></li>" .
      "&nbsp;&nbsp;<li><a href=\"/search/citation_landing.shtml" .
          "?pathway_id=$pid" .  "&what=text\">Citation List</a></li>";

  push @lines, "<h1 class=\"pagetitle\">$pname</h1>"; 
  push @lines, "<hr />";

  push @lines, "<p class=\"citation\"><span class=\"title\">Revision date:</span> $mday" . "-" . $mon_alpha .
      "-" . $year ."</p>";
  if (($source eq "NATURE") || ($source eq "Reactome")) {
    push @lines, "<p class=\"citation\"><span class=\"title\">Curated by:</span> " . join(", ", @curators) . "</p>"; 
    push @lines, "<p class=\"citation\"><span class=\"title\">Reviewed by:</span> " . join(", ", @reviewers) . "</p>";
  }
  push @lines, "<p class=\"citation\"><span class=\"title\">Pathway ID:</span> $ext_id </p>";
  ## Add Pathway Category
  if ($source eq "NATURE") {
     my @category_name;
     my $category_count = 0;
     my $first = 1;
     my $category_url = $BASE . "/browse_categories.shtml?CMD=open&NODE=$pid&GOIDS=";
     my $parent_id = $pid;

       $sqlc = "select pathway_parent, parent_ext_id " .
            "from $schema.pw_pathway_parent " .
            "where pathway_id = $parent_id " ;
       $stmc = $db->prepare($sqlc);
       if (not $stmc) {
         print STDERR "$sqlc\n";
         print STDERR "$DBI::errstr\n";
         print STDERR "prepare call failed\n";
         die;
       }
       if (!$stmc->execute()) {
         print STDERR "$sqlc\n";
         print STDERR "$DBI::errstr\n";
         print STDERR "execute call failed\n";
         die;
       }
       my $find_parent = 0;
       print STDERR "Before while\n";
       while ((my $parent, my $parent_name) = $stmc->fetchrow_array()) {
         $category_url = $category_url . "$parent" . ",";
         $parent_id = $parent;
         if ($parent_id < 200000) {
           $category_name[$category_count] = $parent_name;
           $category_count++;
         }
         if ($category_count == 0) {
           $find_parent = $parent_id; 
         }
         if ($parent_id > 4) {
           my $new_parents = GetParents($parent_id);
           $category_url = $category_url . "$new_parents" . ",";
         }
       }
       $stmc->finish();
     $category_url = $category_url . 1;
     if ($category_count == 0 ) {
       print STDERR "Before GetParentReturn\n";
       $category_name[0] = GetParentName($find_parent); 
       $category_count++;
     } 
     my $category_string = "<b><p class=\"citation\"><span class\"title\">Pathway category:  </span></b>";
     for (my $k = 0; $k < $category_count; $k++) {
       my $cat_name = $category_name[$k];
       $category_string = $category_string . "<a href=\"$category_url\">$cat_name</a>";
       if ($k + 1 == $category_count) {
         $category_string = $category_string . "</p>";
       } else {
         $category_string = $category_string . ", ";
       }
     }
     push @lines, "$category_string";

  }
  if ($parentName ne 'None') {
    if ($format eq 'jpg') {
      push @lines, "<p class=\"citation\"><span class=\"title\">Parent pathway:  </span>  <a href=\"/search/pathway_landing.shtml?pathway_id=$parentPid&source=$source&what=graphic&jpg=on&ppage=1\">$parentName</a></p><br>";
    }
    if ($format eq 'svg') {
      push @lines, "<p class=\"citation\"><span class=\"title\">Parent pathway:  </span><a href=\"/search/pathway_landing.shtml?pathway_id=$parentPid&source=$source&what=graphic&svg=on&ppage=1\">$parentName</a></p>";
      push @lines, "<p><br>";
    }
  }
  if ($source eq "NATURE") {
    # push @lines, "<p class=\"feedback\"><a href=\"mailto:nci-pid\@nature.com\">Submit feedback for this pathway</a></p>";
  } 
  push @lines, "<ul class=\"xoxo viewing-option\">";
  push @lines, "<li><span class=\"hidden\">Additional citation information</span>"; 
  push @lines, "<ul class=\"links\">";
  push @lines, "<li><a class=\"button-style\" href=\"/search/mol_list_landing.shtml" .
          "?pathway_id=$pid" .  "&what=text\">Molecule list</a></li>";
  push @lines, "<li><a class=\"button-style\" href=\"/search/citation_landing.shtml" .
          "?pathway_id=$pid" .  "&what=text\">References</a></li>";
  push @lines, "</ul></li>";
  
  push @lines, "<ul class=\"xoxo viewing-option\">";
  push @lines, "<li><span class=\"norm\">View graphic as:</span>"; 
  push @lines, "<ul class=\"links\">";
 
  push @lines, "<li><a class=\"button-style\" href=\"/search/pathway_landing.shtml" .
          "?pathway_id=$pid" . "&source=$source" . "&genes_a=$genes_a" . "&genes_b=$genes_b" . "&what=graphic&jpg=on&ppage=1\">JPG</a></li>" ;

  push @lines, "<li><a class=\"button-style\" href=\"/search/pathway_landing.shtml" .
       "?pathway_id=$pid" . "&source=$source" . "&genes_a=$genes_a" . "&genes_b=$genes_b" . "&what=graphic&svg=on&ppage=1\">SVG</a></li>";

  push @lines, "<li><a class=\"button-style\" href=\"/xaml_landing.shtml" .
       "?pathway_id=$pid" . "&source=$source" . "&genes_a=$genes_a" . "&genes_b=$genes_b" . "&what=text&xaml=on&ppage=1\">Silverlight-beta</a></li>";
 

  push @lines, "<li><div class=\"helpbox\">";
  push @lines, "<a href=\"#\" class=\"help\">Help</a>";
  push @lines, "<div class=\"description\">";
  push @lines, "<h2 class=\"title\">Information on graphics types</h2>";
  push @lines, "<p class=\"norm\">The pathway can be viewed as a JPG, SVG or Silverlight image. All molecules and interactions within the graphics formats are interactive; clicking on them will open information windows.</p>";
  push @lines, "<p class=\"norm\"><B>JPG:</b>  All common web browsers are able to display JPG images.  However, please note that some browsers may not support very large JPG files.</p>";
  push @lines, "<p class=\"norm\"><B>SVG:</b>  The advantage of SVG images is that they can be more easily zoomed and panned on-screen.  If your browser is not completely SVG enabled, please visit the <a href=\"/PID/userguide/output_formats.shtml\">User guide</a> and follow the instructions for installing the SVG plug-in specific to your browser</p>"; 
  push @lines, "<p class=\"norm\"><b>SVG image navigation:</b>  Use the left mouse button to activate the image.  To <b>zoom</b> in or out, use the right mouse button and select from the options on the pop-up menu.  To <b>pan</b> the image, hold down the alt key and click-and-drag the image with the mouse.</p>";
  push @lines, "<p class=\"norm\"><b>Silverlight:  </b>Please visit http://www.microsoft.com/getsilverlight/Get-Started/Install/Default.aspx to install the plug-in for your browser.</p>";
  push @lines, "<p class=\"norm\"><b>Silverlight image navigation:  Zoom</b> is available at the top-left of the image. To <b>pan</b> merely drag the image with the mouse.  In the image you may search for specific molecules or all molecules at a certain subcellular location.  To search for molecules you may enter their names, Entrez Gene identifiers or UniProt accession numbers.</p>";
  push @lines, "</div>";
  push @lines, "</div>"; 
  push @lines, "</li>";
  push @lines, "</ul>"; 
  push @lines, "</li>";
  
  push @lines, "<li><span class=\"norm\">Save code as:</span>";
  push @lines, "<ul class=\"links\">";
  push @lines, "<li><a class=\"button-style\" href=\"/search/pathway_landing.shtml" .
          "?pathway_id=$pid" . "&pathway_short=$ext_id" . "&source=$source" . "&what=text&xml=on&ppage=1\">XML</a></li>" ;
  push @lines, "<li><a class=\"button-style\" href=\"/search/pathway_landing.shtml" .
       "?pathway_id=$pid" . "&pathway_short=$ext_id" . "&source=$source" . "&what=text&biopax=on&ppage=1\">BioPAX</a></li>";
  push @lines, "<li><div class=\"helpbox\">";
  push @lines, "<a href=\"#\" class=\"help\">Help</a>";
  push @lines, "<div class=\"description\">";
  push @lines, "<h2 class=\"title\">Information on types of pathway code</h2>";
  push @lines, "<p class=\"norm\">You can view the pathway code in either XML or BioPax formats.  Please see the <a href=\"/PID/userguide/output_formats.shtml\">User guide</a> for more information.</p>";
  push @lines, "</div>";
  push @lines, "</div>";
  push @lines, "</li>";
  push @lines, "</ul>";
  push @lines, "</li>";
  push @lines, "</ul>"; 

  if ($format eq 'jpg') {
    push @lines, "<p class=\"norm\">A <a href=\"/search/pathway_landing.shtml?pathway_id=$pid" . "&source=$source" . "&jpg=true&what=graphic&ppage=1/#key\">key to the image</a> icons appears below the pathway.</p>";
  } else {
    push @lines, "<p class=\"norm\">A <a href=\"/search/pathway_landing.shtml?pathway_id=$pid" . "&source=$source" . "&svg=true&what=graphic&ppage=1/#key\">key to the image</a> icons appears below the pathway.</p>";
  } 
  push @lines, "<p class=\"norm\">Click on a biomolecule or interaction for additional information.</p>";
  if ($format eq "svg") {
    push @lines, "<p class=\"norm\"><a href=/PID/userguide/network_maps.shtml#>SVG Image Navigation Help</a></p>";
  }
  
  my %genes_a;
  my %genes_b;

  my @genes_1 = split(/,/,$genes_a);
  my @genes_2 = split(/,/,$genes_b);
  for (my $i = 0; $i < $genes_a; $i++) {
    if ($genes_1[$i] ne '') {
      $genes_a{$genes_1[$i]} = 1;
    }
  } 
  for (my $i = 0; $i < $genes_b; $i++) {
    if ($genes_2[$i] ne '') {
      $genes_b{$genes_2[$i]} = 1;
    }
  }
 
  my $mols_a = GetMolIdsForEntrezGene($db, $pid, \%genes_a);
  my $mols_b = GetMolIdsForEntrezGene($db, $pid, \%genes_b);
 
  if ($format ne 'xaml') { 
    push @lines, CreateGraphicFile($db, $schema, $pid, $source, $mols_a, $mols_b, $format);
  } 
  $db->disconnect();
  return join("\n", @lines) . "\n";
  
}

######################################################################
sub ListAllPathways_1 {
  my ($base, $word) = @_;

  $BASE = $base;


  $word =~ s/^\s+//;
  $word =~ s/\s+$//;
  $word =~ s/\s+/ /g;
  $word = lc($word);

  my (@lines, @curated_lines, @biocarta_lines, @reactome_lines);
  my ($parent_id, $subnet_id);
  my $db = DBI->connect("DBI:Oracle:" . $db_inst,
        $db_user, $db_pass);
  if (not $db or $db->err()) {
    print STDERR "#! Cannot connect to " . $db_user . "@" .
        $db_inst . "\n";
    die;
  }

  push @lines, "<h2 id=\"NCI-Nature\" class=\"pagetitle2\">NCI-Nature curated</h2>";
  push @lines, "<table class=\"data\" border=\"0\" cellpadding=\"0\" cellspacing=\"0\ id=\"NCI-Nature\" >";
  push @lines, "<colgroup span=\"1\" />";
  push @lines, "<colgroup span=\"1\" />";
			
  push @lines,"<thead>";
  push @lines,"<tr>";
  push @lines,"<th scope=\"colgroup\">Pathway</th>";
  push @lines,"<th scope=\"colgroup\" class=\"colspan column3\">Source</th>";
  push @lines, "</thead>";
  push @lines, "<tbody>";

  my ($sql, $stm, $sqlc, $stmc);
  my ($source, $pname, $ext_id, $pid, $curator, $role);
  my (@curators, @reviewers);

  $sql = "select " .
      "s.source_name, p.pathway_name, p.ext_pathway_id, p.pathway_id " .
      "from $schema.pw_source s, $schema.pw_pathway p " .
      "where s.source_id = p.pathway_source_id " .
      "and subnet='N' " .
      "and p.pathway_id > 99999 " .
      "order by s.source_name, p.pathway_name "; 


  $stm = $db->prepare($sql);

  if(not $stm) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "prepare call failed\n";
    die;
  }
  if(!$stm->execute()) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "execute call failed\n";
    die;
  }
  while (($source, $pname, $ext_id, $pid) =
      $stm->fetchrow_array()) {

    if ($word) {
      if ((index(lc($pname), $word) < 0) &&
         (index(lc($ext_id), $word) < 0)) {
        next;
      } 
    }

    undef @curators;
    undef @reviewers;
    if ($source eq "NATURE")  {
      $sqlc = "select curator, role " .
              "from $schema.pw_curators " .
              "where pathway_id = $pid " ;
      $stmc = $db->prepare($sqlc);
      if (not $stmc) {
        print STDERR "$sqlc\n";
        print STDERR "$DBI::errstr\n";
        print STDERR "prepare call failed\n";
        die;
      }
      if (!$stmc->execute()) {
        print STDERR "$sqlc\n";
        print STDERR "$DBI::errstr\n";
        print STDERR "execute call failed\n";
        die;
      }
      while (($curator, $role) = $stmc->fetchrow_array()) {
        if ($role eq "C") {
          push @curators, $curator; 
        } else {
          push @reviewers, $curator; 
        }
      }
      $stmc->finish(); 
      $source = "Curated";
    }
    my $nice_source = $source;
    $nice_source =~ s/Curated/NCI-Nature curated/;
    $nice_source =~ s/BioCarta/BioCarta Imported/;
    $nice_source =~ s/Reactome/Reactome Imported/;
    if ($source eq "Curated") {
      push @curated_lines, join("\t", $nice_source, $pname, $ext_id);
    }
    if ($source eq "BioCarta") {
      push @biocarta_lines, join("\t", $nice_source, $pname, $pid);
    }
    if ($source eq "Reactome") {
      push @reactome_lines, join("\t", $nice_source, $pname, $pid);
    }
 
    # Now get subnets

    $sqlc = "select subnet_pathway_id, b.pathway_name, b.ext_pathway_id from $schema.pw_subnet a, $schema.pw_pathway b" .
      " where parent_pathway_id = $pid and a.subnet_pathway_id = b.pathway_id order by b.pathway_name\n";
    $stmc = $db->prepare($sqlc);

    if (not $stmc) {
        print STDERR "$sqlc\n";
        print STDERR "$DBI::errstr\n";
        print STDERR "prepare call failed\n";
        die;
      }
      if (!$stmc->execute()) {
        print STDERR "$sqlc\n";
        print STDERR "$DBI::errstr\n";
        print STDERR "execute call failed\n";
        die;
      }
      my $pname2;
      my $ext_path_id;
      while (($subnet_id, $pname2, $ext_path_id) = $stmc->fetchrow_array()) {
         $pname2 = "<ul style=\"background-image: url(/PID/images/clear.gif);\">Subnetwork:  " . $pname2 . "</ul>";
         if ($source eq "Curated") {
   	   push @curated_lines, join("\t", $nice_source, $pname2, $ext_path_id);
         }
         if ($source eq "BioCarta") {
           push @biocarta_lines, join("\t", $nice_source, $pname2, $subnet_id);
         }
         if ($source eq "Reactome") {
           push @reactome_lines, join("\t", $nice_source, $pname2, $subnet_id);
	 } 
      }
      $stmc->finish();
  }

  $db->disconnect();

  my $n = 0;
  my $first = 0;
  for (@curated_lines, @biocarta_lines, @reactome_lines) {
    my ($nice_source, $pname, $pid) = split /\t/;
    $n++;
    if (($first == 1) && ($nice_source eq "Reactome Imported")) {
        $first = 2;
      push @lines, "</tbody></table>";
      push @lines, "<h2 id=\"Reactome\" class=\"pagetitle2\">Reactome</h2>";
      push @lines, "<table class=\"data\" border=\"0\" cellpadding=\"0\" cellspacing=\"0\ id=\"Reactome\" >";
      push @lines, "<colgroup span=\"1\" />";
      push @lines, "<colgroup span=\"1\" />";

      push @lines,"<thead>";
      push @lines,"<tr>";
      push @lines,"<th scope=\"colgroup\">Pathway</th>";
      push @lines,"<th scope=\"colgroup\" class=\"colspan column3\">Source</th>";
      push @lines, "</thead>";
      push @lines, "<tbody>";

      push @lines, ($n %2 == 0) ? "<tr>" : "<tr class=\"odd\">";
      push @lines, "<td id=\"biocarta\"><a href=\"/search/pathway_landing.shtml" .
        "?pathway_id=$pid" . "&source=$nice_source" . "&what=graphic&jpg=on&ppage=1\">" .
        "$pname</a></td>";

      push @lines, "<td>$nice_source</td>";
    }

    if (($first == 0) && ($nice_source eq "BioCarta Imported")) {
        $first=1;
      push @lines, "</tbody></table>";
      push @lines, "<h2 id=\"BioCarta\" class=\"pagetitle2\">BioCarta</h2>";
      push @lines, "<table class=\"data\" border=\"0\" cellpadding=\"0\" cellspacing=\"0\ id=\"BioCarta\" >";
      push @lines, "<colgroup span=\"1\" />";
      push @lines, "<colgroup span=\"1\" />";

      push @lines,"<thead>";
      push @lines,"<tr>";
      push @lines,"<th scope=\"colgroup\">Pathway</th>";
      push @lines,"<th scope=\"colgroup\" class=\"colspan column3\">Source</th>";
      push @lines, "</thead>";
      push @lines, "<tbody>";

      push @lines, ($n %2 == 0) ? "<tr>" : "<tr class=\"odd\">";
      push @lines, "<td id=\"biocarta\"><a href=\"/search/pathway_landing.shtml" .
        "?pathway_id=$pid" . "&source=$nice_source" . "&what=graphic&jpg=on&ppage=1\">" .
        "$pname</a></td>";

      push @lines, "<td>$nice_source</td>";

 
     } else {
       push @lines, ($n %2 == 0) ? "<tr>" : "<tr class=\"odd\">";
       my $pname2;  
       if (substr($pname,0,1) eq '<') {
         my @temp = split(":  ", $pname);
         my $lenp = length($pname2);
         $pname2 = substr($temp[1],0,$lenp - 5); 
       } else {
         $pname2 = $pname;
       }
       push @lines, "<td><a href=\"/search/pathway_landing.shtml" .
       "?pathway_id=$pid" . "&pathway_name=$pname2" . "&source=$nice_source" . "&what=graphic&jpg=on&ppage=1\">" .
       "$pname</a></td>";
 
       push @lines, "<td>$nice_source</td>";
     }
     push @lines, "</tr>";
  }

  push @lines, "</tbody></table>";
  return join("\n", @lines) . "\n";

}

######################################################################
sub PathwayPage_1 {
  my ($base, $pathid) = @_;

  $BASE = $base;

  my $db = DBI->connect("DBI:Oracle:" . $db_inst, $db_user, $db_pass);
  if (not $db or $db->err()) {
    print STDERR "Cannot connect to " . $db_user . "@" . $db_inst . "\n";
    die;
  }

  my @lines;

  my ($sql, $stm, $sqlc, $stmc, $m, $d, $y);
  my ($source, $pname, $ext_id, $pid, $updated, $curator, $role);
  my (@curators, @reviewers);

  $sql = 
      "select s.source_name, p.pathway_name, p.ext_pathway_id, " .
      "p.pathway_id, p.last_updated " .
      "from $schema.pw_source s, $schema.pw_pathway p " .
      "where s.source_id = p.pathway_source_id " .
      "and p.pathway_id = $pathid " ;

  $stm = $db->prepare($sql);
  if(not $stm) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "prepare call failed\n";
    die;
  }
  if(!$stm->execute()) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "execute call failed\n";
    die;
  }
  while (($source, $pname, $ext_id, $pid, $updated) =
      $stm->fetchrow_array()) {
    $ext_id =~ s/pathway/Pathway/;
    my ($sec, $min, $hr, $mday, $mon, $year, $wday, $yday, $isdst) =
        localtime($updated);
    $mon++;
    $year = $year + 1900;
    my $nice_source = $source;
    $nice_source =~ s/NATURE/NCI-Nature curated/;
    $nice_source =~ s/BioCarta/BioCarta Imported/;
    $nice_source =~ s/Reactome/Reactome Imported/;
    undef @curators;
    undef @reviewers;

    if (($source eq "NATURE") || ($source eq "Reactome")) {
      $sqlc = "select curator, role " .
              "from $schema.pw_curators " .
              "where pathway_id = $pid " ;
      $stmc = $db->prepare($sqlc);
      if (not $stmc) {
        print STDERR "$sqlc\n";
        print STDERR "$DBI::errstr\n";
        print STDERR "prepare call failed\n";
        die;
      }
      if (!$stmc->execute()) {
        print STDERR "$sqlc\n";
        print STDERR "$DBI::errstr\n";
        print STDERR "execute call failed\n";
        die;
      }
      while (($curator, $role) = $stmc->fetchrow_array()) {
        if ($role eq "C") {
          push @curators, $curator;
        } else {
          push @reviewers, $curator; 
        }
      }
      $stmc->finish();
    }
    $stm->finish();
    push @lines,
        "<center><span style=\"font-family:verdana; font-size:10pt; font-weight:bold; color:#36648B; background-color:#c1cdcd;\">$pname</span></center>" .
        "<P>" .
        "<table style=\"font-family:verdana; font-size:9pt; color:#36648B\">" .
        "<tr><td><b>Pathway_id:</b></td><td>$ext_id</td></tr>" .
        "<tr><td><b>Source:</b></td><td>$nice_source</td></tr>" .
        "<tr><td><b>Last Revised:</b></td><td>$mon/$mday/$year</td></tr>" .
        ((@curators) ?
        "<tr><td valign=top><b>Curators:</b></td><td>" . join("<br>", @curators) . "</td></tr>" : "") .
        ((@reviewers) ?
        "<tr><td valign=top><b>Reviewers:</b></td><td>" . join("<br>", @reviewers) . "</td></tr>" : "") .
        "<tr><td><b>Graphical<br>Representation:</b></td><td valign=bottom>" .
        "<a href=\"/search/pathway_landing.shtml" .
            "?pathway_id=$pid" . "&source=$nice_source" . "&what=graphic&jpg=on\" style=\"color:#551A8B\">JPG</a> " .
        "<a href=\"/search/pathway_landing.shtml" .
            "?pathway_id=$pid" . "&source=$nice_source" . "&what=graphic&svg=on\" style=\"color:#551A8B\">SVG</a> " .
        "<a href=\"/search/pathway_landing.shtml" .
            "?pathway_id=$pid" . "&source=$nice_source" . "&what=text&biopax=on\" style=\"color:#551A8B\">BioPAX</a> " .
        "<a href=\"/search/pathway_landing.shtml" .
            "?pathway_id=$pid" . "&source=$nice_source" . "&what=text&xml=on\" style=\"color:#551A8B\">XML</a> " .
#        (($ext_id =~ /hsa/) ? "" :
#        "<a href=\"$BASE/Pathways/BioCarta/h_".$ext_id."\" style=\"color:#551A8B\">BioCarta</a> ") .
        "</td></tr>" .
        "</table>" ;
  }

  
  $db->disconnect();

  return join("\n", @lines) . "\n";
}

######################################################################
sub GetPwImage_1 {
  my ($base, $id) = @_;

  $BASE = $base;

  if (Scan($id) > 0) {
    return "Call failed...illegal characters in URL";
  }
  my (@buf);
  my $i = 0;
  my $fn = CACHE_ROOT . "/" . PW_CACHE_PREFIX . ".$id" ;
  open(INF, $fn) or return "Failed to open cache file $fn";
  while (read INF, $buf[$i], 16384) {
    $i++;
  }
  close INF;
  return join("", @buf);
  
}

######################################################################
sub GraphicsByClan {
  my ($pw, $lv, $parm, $graphic_type, $mol_vg, $atom_vg, $coloring) = @_;

  use Clan;
  use Cache;

  my ($size, $clan, @clans, %sizeof, @cache_ids, @map_files);

  for $clan (@{ $pw->Clans }) {
    $size = $clan->ClanSize();
    push @{ $sizeof{$size} }, $clan;
  }
  for $size (sort r_numerically keys %sizeof) {
    for $clan (@{ $sizeof{$size} }) {
      push @clans, $clan;
    }
  }
  my @dot2;
  my $counter = 0;
  for $clan (@clans) {
    $size = $clan->ClanSize();
    if ($size > MAX_DOT_INTERACTIONS  && ! $parm->{collapse}{processes}) {
      # SetStatus(S_BAD_REQUEST);
      return "Number of interactions $size exceeds maximum " .
          MAX_DOT_INTERACTIONS;
    }
    my $cache = new Cache(CACHE_ROOT, PW_CACHE_PREFIX);
    my ($cache_id, $filename) = $cache->MakeCacheFile();
    print STDERR "Filename:  $filename\n";
    
    if ($cache_id != $CACHE_FAIL) {
      print STDERR "Cache didn't fail\n";
      open(OUTF, ">$filename.dot");
      my $dot = new DOToutput($pw, $lv, *OUTF);
      $dot->SetGraphicFormat($graphic_type);
      if ($parm->{show}{atom_id}) {
        $dot->ShowAtomIds(1);
      }
      if ($parm->{collapse}{molecules}) {
        $dot->CollapseMolecules(1);
      }
      if ($parm->{show}{subtype_line}) {
        $dot->ShowSubTypeLines(1);
      }
      if ($parm->{collapse}{processes}) {
        $dot->CollapseAtoms(1);
      }
      if (defined $mol_vg) {
        $dot->SetMolValueGenerator($mol_vg);
        if (defined $coloring) {
          $dot->SetSimColoring($coloring);
        }
      } else {
        $dot->ColorMols($parm->{mol_id});
        if (defined $parm->{connect}) {
          $dot->ColorMols($parm->{connect});
          my %temp;
          for my $mol (@{ $pw->Mols() }) {
            if ($pw->MolType($mol) eq $COMPLEX_TYPE) {
              for my $comp (@{ $pw->Components($mol) }) {
                my $m = $pw->ComponentMol($mol, $comp);
                if (defined $parm->{connect}{$m}) {
                  $temp{$mol} = 1;
                }
              }
            }
          }
          $dot->ColorMols(\%temp);
        }
      }
      if (defined $coloring) {
        $dot->SetSimColoring($coloring);
      }
      if (defined $atom_vg) {
        $dot->SetAtomValueGenerator($atom_vg);
      }
      $dot->DOTGraph($clan);
      close(OUTF);
      chmod 0666, "$filename.dot";
      push @cache_ids, $cache_id;
      my $graphic_hold = $graphic_type;
      if ($graphic_type eq 'xaml') {
      #   $graphic_type = 'svg';
      }
      if ($graphic_type ne 'xaml') {
        my $cmd = DOT_PGM .  " $filename.dot -T$graphic_type > $filename";
        system($cmd);
      }
      my $cmd = DOT_PGM .  " $filename.dot -Tdot > $filename.coords";
      system($cmd);
      print STDERR "Coords created\n";
      chmod 0666, $filename;
      $cmd = "rm $filename.meta";
      system($cmd);
      $cmd = "/share/content/PID/run/dot2xaml.pl $filename.coords $filename.xamltemp $filename.meta";
      system($cmd);

      my $coordfile = $filename . ".coords";
      my $xamlfile = $filename . ".xamltemp"; 
      my $metafile = $filename . ".meta";
      
      open (META, "< $filename.meta");
      open (XAMLTEMP, "< $filename.xamltemp");
      open (XAML, "> $filename.xaml");
   
      my $count = 1;
      while (<XAMLTEMP>) {
        chomp;
        if ($count == 4) {
          while (<META>) {
            chomp;
            my $temp2 = $_;
            print XAML "$temp2\n";
          }
        }
        my $temp = $_;
        print XAML "$temp\n";
        $count++;
      }
     
      close(META);
      close(XAMLTEMP);
      close(XAML);
      my $xamlname = $filename . ".xaml";
      chmod 0666, "$xamlname";
      if ($graphic_type eq "xaml") {
        my $cmd = "cp $xamlname $filename";
        system($cmd);
        chmod 0666, "$filename";
      } 
      if ($graphic_type eq "jpg") {
        my $cmd = DOT_PGM .  " $filename.dot -Tcmap > $filename.map";
        system($cmd);
        chmod 0666, "$filename.map";
        push @map_files, "$filename.map"
      }
      unlink("$filename.dot");
    } else {
      return 0;
    }
  }
  print STDERR "Before return\n";
  return (join(",", @cache_ids), join(",", @map_files) );
}

######################################################################
sub MakeCommandFile_1 {
  my ($base, $graphic_type, $parm_string) = @_;

  $BASE = $base;

  $parm_string = "$parm_string\n";    ## add just in case

  my $parm = new PathParams($lv);
  if (! $parm->ReadParamString($parm_string)) {
    return join("\n", @{ $parm->{errors} }) . "\n";
  };

  my $db = DBI->connect("DBI:Oracle:" . $db_inst, $db_user, $db_pass);
  if (not $db or $db->err()) {
    print STDERR "Cannot connect to " . $db_user . "@" . $db_inst . "\n";
    die;
  }
  if (! $parm->CheckParams($db) ) {
    return join("\n", @{ $parm->{errors} }) . "\n";
  }

  $db->disconnect();

  my @cmds;
  for (split("\n", $parm_string)) {
    if (/^\s+$/) {
      next;
    }
    if (/^(db_user|db_pass|db_inst|db_schema)/) {
      next;
    }
    if (/^print\tcmd/) {
      next;
    }
    if (/^molecule/) {
	
    }
    push @cmds, $_;
  }
  for my $mm_locus (keys %{ $parm->{mousetrap} }) {
    for my $hs_locus (@{ $parm->{mousetrap}{$mm_locus} }) {
      push @cmds, "## using Hs LL $hs_locus for Mm LL $mm_locus";
    }
  }
  return join("\n", @cmds) . "\n";
}

######################################################################
sub MakeColorScale {
  my ($colorspec) = @_;

  my ($number, $color);
  my (@number_breaks, @color_scale);

  for $number (sort numerically keys %{ $colorspec }) {
    push @number_breaks, $number;
    push @color_scale, $$colorspec{$number};
  }
  return (\@number_breaks, \@color_scale);
}

######################################################################
sub PrPathReturn {
  my ($status, $data_type, $data, $map_files, $atom_ids,
      $parm, $pw, $interactionpage) = @_;
  print STDERR "In PrPathReturn\n";
  if ($data_type eq "xaml") {
    print STDERR "In xaml parm:  $parm\n";
    my @lines;
    my @clans     = split(",", $data);
    my $n = 0;
    my $linecount = 0;
    my $mingid = 1000;
    my $maxgid = 0;
    for my $gid (@clans) {
      if ($gid < $mingid) {
        $mingid = $gid;
      }
      if ($gid > $maxgid) {
        $maxgid = $gid;
      } 
    }
    my @nolines;
    MatchResults($parm, $pw, "html");
    my $locationString = '';
    for (my $gid = $mingid; $gid <= $maxgid; $gid++) {
      my $xamlfile = "/share/content/PID/data/cache/PW." . $gid . ".xaml";
      open (XAMLFILE, $xamlfile);
      while (<XAMLFILE>) {
        my @line = split(/=/,$_);
        if (($line[0] eq 'Location') && ($line[1] ne '""') && ($line[1] ne '')) {
          $line[1] =~ s/"//g;
          $line[1] =~ s/\n//g;
          if ($line[1] ne '') {
            $locationString = $locationString . ':' . $line[1];
          } 
        }
      }
      close(XAMLFILE);
    }
    push @lines, "<body onload=\"WaitToCallXaml($mingid, $maxgid);\" id=\"article-related\" class=\"www-nature-com-nci fullstretch\"  > ";
     push @lines, "<form id=\"formtest\" style=\"height: 10%\" >";
push @lines, "<input class=\"submit\" id=\"ZoomInButton\" type=\"button\" value=\"+\" alt=\"Zoom in\" onclick=\"CallSL_Zoom(1.5, $mingid, $maxgid)\" /> ";
    push @lines, "<input class=\"submit\" id=\"ZoomOutButton\" type=\"button\" value=\"-\" alt = \"Zoom out\" onclick=\"CallSL_Zoom(.75, $mingid, $maxgid)\" /> ";

    push @lines, "Molecule name/ID search:<input id=\"SearchXAML\" type=\"text\" value=\"\" /> ";
    push @lines, "<input id=\"searchTypeList\" value=\"entitysearch\" type=\"hidden\" /> ";
			
    push @lines, "<input class=\"submit\" id=\"SearchXAMLButton\" type=\"button\" value=\"Go\" onclick=\"callSL_SearchPathwayXAML($mingid, $maxgid)\" /> ";
    push @lines, "Location:  ";
    push @lines, "<select id=\"searchLocationList\" > ";
    push @lines, "</select> ";
    push @lines, "<script language=\"javascript\"> ";
    push @lines, "addLocation(\"$locationString\")";
    push @lines, "</script>";
    push @lines, "<input class=\"submit\" id=\"searchByLocation\" type=\"button\" value=\"Go\" onclick=\"callSL_SearchPathwayByLocation($mingid, $maxgid)\" /> ";
    push @lines, "<BR><BR><input class=\"submit\" id=\"LoadXAMLButton\" type=\"button\" value=\"Reload\" onclick=\"LoadSL_GetPathwayXAML($mingid, $maxgid)\" /> ";

    push @lines, "</form>";
			 
    for my $gid (@clans) {
      if (@clans > 1) {
        $n++; 
         push @lines, "<p><table><tr bgcolor=firebrick>" .
            "<td><b><font color=white>Network Map $n</font></b></td>" .
            "</tr></table>"; 
        
      }
       
      push @lines, "<form id=\"form$gid\" style=\"height: 100%\" >"; 
      push @lines, "<div id=\"silverlightControlHost\">";
      push @lines, "<div class=\"container-image\">";
      push @lines, "</div>";
      push @lines, "<object id=\"SLP$gid\" name=\"SLP$gid\" data=\"data:application/x-silverlight-2,\" type=\"application/x-silverlight-2\" ";
      push @lines, "width=\"100%\" height=\"90%\"> ";
      push @lines, "<param name=\"source\" value=\"/ClientBin/PathwayViewer.xap\" /> ";
      push @lines, "<param name=\"onError\" value=\"onSilverlightError\" /> ";
      push @lines, "<param name=\"background\" value=\"white\" /> ";
      push @lines, "<param name=\"minRuntimeVersion\" value=\"3.0.40624.0\" /> ";
      push @lines, "<param name=\"autoUpgrade\" value=\"true\" /> ";
      push @lines, "<param name=\"windowless\" value=\"true\" /> ";
      push @lines, "<a href=\"http://go.microsoft.com/fwlink/?LinkID=149156&v=3.0.40624.0\" style=\"text-decoration: none\"> ";
      push @lines, "<img src=\"http://go.microsoft.com/fwlink/?LinkId=108181\" alt=\"Get Microsoft Silverlight\" ";
      push @lines, "style=\"border-style: none\" />";
      push @lines, "</a>";
      push @lines, "</object>";
      push @lines, "</div>";
      push @lines, "</form>"; 
      # my $xamlfile = "/share/content/PID/data/cache/PW." . $gid . ".xaml";
      # print STDERR "Xamlfile:  $xamlfile\n"; 
      # open (XAMLFILE, $xamlfile);
      #while (<XAMLFILE>) {
      #  push @lines, $_;
      #  $linecount++;
      #}
      #print STDERR "Linecount:  $linecount\n";
      #close(XAMLFILE);
      # Add stuff for Silverlight
    }
    my $atom_ids;
    return join("\001", $status, $data_type,
        join("\n", @lines), $atom_ids);
 
  }

  if ($data_type eq "text") {
    if (($status eq S_NO_DATA) && (defined $pw)) {
      my $match = MatchResults($parm, $pw, "html");
      return join("\001", $status, $data_type,
          ($data . "<br>\n$match\n"),
          $atom_ids);
    } else {
      return join("\001", $status, $data_type, $data, $atom_ids);
    }
  } else {
    my @lines;
    my @clans     = split(",", $data);
    my @map_files = split(",", $map_files);
    my $n = 0;

    push @lines, "<blockquote>";
    push @lines, MatchResults($parm, $pw, "html");
    push @lines, "</blockquote>";

    for my $gid (@clans) {
      $n++;
      if (@clans > 1) {
        push @lines, "<p><table><tr bgcolor=firebrick>" .
            "<td><b><font color=white>Network Map $n</font></b></td>" .
            "</tr></table>";
      }
      if ($data_type eq "svg") {
        push @lines, "<embed type=\"image/svg-xml\" height=800 width=1000 ";
        push @lines, "src=\"/pidcgi/GetPwImage.pl?GID=$gid&FORMAT=SVG\" alt = \"Image of Network Search Map\" >";
      } elsif ($data_type eq "jpg") {
        push @lines, "<img src=\"/pidcgi/GetPwImage.pl?GID=$gid&FORMAT=JPG\" " .
            "usemap=#map_$gid alt = \"Image of Network Search Map\" >";
        push @lines, "<map name=map_$gid>";
        my $map_file = $map_files[$n-1];
        open(INF, "$map_file");
        while (<INF>) {
          s/[\n\r]+//;
          push @lines, $_;
        }
        close INF;
        push @lines, "</map>";
        unlink($map_file);
      }
    }
 
    $interactionpage = 1;
 
    if ((@clans > 0) && ($interactionpage < 1)) {
      push @lines, "</div>";
      push @lines, "<div id=\"key-constrain\"> ";
      push @lines, "<div class=\"key-constrain\"> ";
push @lines, "<h2 class=\"pagetitle\" id=\"key\">Key to icons</h2> ";
push @lines, "<ul class=\"key\" ";
push @lines, "><li class=\"process\"><h3 class=\"title\">Process types</h3> ";
push @lines, "<ul";
push @lines, "><li><img src=\"/images/img_transcription.gif\" class=\"process\" alt=\"transcription image\" /> Transcription</li";
push @lines, "><li><img src=\"/images/img_modification.gif\" class=\"process\" alt=\"Modification image\" /> Modification <span class=\"description\">(interaction)</span></li";
push @lines, "><li><img src=\"/images/img_reaction.gif\" class=\"process\" alt=\"reaction image\" /> Reaction</li";
push @lines, "><li><img src=\"/images/img_translocation.gif\" class=\"process\" alt=\"translocation image\" /> Translocation</li";
push @lines, "><li><img src=\"/images/img_macroprocess.gif\" class=\"process\" alt=\"Biological process image\" /> Biological process <span class=\"description\">(multi-step event)</span></li";
push @lines, "></ul";
push @lines, "></li";
push @lines, "><li class=\"edge\"><h3 class=\"title\">Edge types</h3>";
push @lines, "<ul";
push @lines, "><li class=\"input-biomolecule\">Input biomolecule <span class=\"description\">(protein/compound/complex/RNA)</span></li";
push @lines, "><li class=\"inhibitor\">Negative regulator</li";
push @lines, "><li class=\"agent\">Positive regulator</li";
push @lines, "><li class=\"output-biomolecule\">Output biomolecule <span class=\"description\">(protein/compound/complex/RNA)</span></li";
push @lines, "></ul";
push @lines, "></li";
push @lines, "><li class=\"common-locations clearfix\"><h3 class=\"title\">Common subcellular locations</h3>";
push @lines, "<ul";
push @lines, "><li>m - transmembrane</li";
push @lines, "><li>cy - cytoplasm</li";
push @lines, "><li>mi - mitochondria</li";
push @lines, "><li>n - nucleus</li";
push @lines, "><li>ex - extracellular region</li";
push @lines, "><li>v - vesicle</li";
push @lines, "><li>cs - calcium store</li";
push @lines, "><li>en - endosome</li";
push @lines, "><li>er - endoplasmic reticulum</li";
push @lines, "><li>g - Golgi apparatus</li";
push @lines, "><li>l -  lysosome</li";
push @lines, "></ul";
push @lines, "></li";
push @lines, "><li class=\"other-notations\">";
push @lines, "<table class=\"data\">";
push @lines, "<caption>Other common notation</caption>";
push @lines, "<colgroup span=\"1\" />";
push @lines, "<colgroup span=\"1\" />";
push @lines, "<thead>";
push @lines, "<tr>";
push @lines, "<th scope=\"colgroup\" class=\"colspan\">Notation</th>";
push @lines, "<th scope=\"colgroup\" class=\"colspan\">Description</th>";
push @lines, "</tr>";
push @lines, "</thead>";
push @lines, "<tbody>";
push @lines, "<tr class=\"odd\">";
push @lines, "<td><img src=\"/images/img_complex.gif\" class=\"process\" alt=\"Complex name\" /></td>";
push @lines, "<td>Complex with more than five components; click on the complex to obtain information on components</td>";
push @lines, "</tr>";
push @lines, "<tr>";
push @lines, "<td>Biomolecule name +</td>";
push @lines, "<td rowspan=\"5\">Active states of a biomolecule</td>";
push @lines, "</tr>";
push @lines, "<tr>";
push @lines, "<td>Biomolecule&nbsp;name&nbsp;+1</td>";
push @lines, "</tr>";
push @lines, "<tr>";
push @lines, "<td>Biomolecule&nbsp;name&nbsp;+2</td>";
push @lines, "</tr>";
push @lines, "<tr>";
push @lines, "<td>Biomolecule&nbsp;name&nbsp;+3</td>";
push @lines, "</tr>";
push @lines, "<tr>";
push @lines, "<td>Biomolecule&nbsp;name&nbsp;+4</td>";
push @lines, "</tr>";
push @lines, "<tr class=\"odd\">";
push @lines, "<td>Biomolecule&nbsp;name&nbsp;-</td>";
push @lines, "<td>Inactive&nbsp;biomolecule</td>";
push @lines, "</tr>";
push @lines, "<tr>";
push @lines, "<td><img src=\"/images/img_pathwaySubnetwork.gif\" class=\"process\" alt=\"Pathway subnetwork image\" /></td>";
push @lines, "<td>Hyperlink to pathway or subnetwork</td>";
push @lines, "</tr>";
push @lines, "<tr class=\"odd\">";
push @lines, "<td><img src=\"/images/img_familyname.gif\" class=\"process\" alt=\"family members image\" /></td>";
push @lines, "<td>Instances of family members in pathways</td>";
push @lines, "</tr>";
push @lines, "</tbody>";
push @lines, "</table>";
push @lines, "</li";
push @lines, "></ul>";
push @lines, "</div>";
push @lines, "</div>";
 
      push @lines, "<p class=\"back-to-top hidden\"><a href=\"#top\">Top<span class=\"hidden\"> of page</span></a></p>";
    }
    return join("\001", $status, $data_type,
        join("\n", @lines), $atom_ids);
  }
}

######################################################################
sub MatchResults {
  my ($parm, $pw, $format) = @_;

  my (%used_mols, @not_matched, @matched_not_used, @matched_used, @lines);
  my (@pmatched, @not_pmatched);
  my (@macromatched, @not_macromatched);

  my $source_id = $parm->{source_id} ;
  my $mol_id;
  my $db = DBI->connect("DBI:Oracle:" . $db_inst, $db_user, $db_pass);
  if (not $db or $db->err()) {
    print STDERR "Cannot connect to " . $db_user . "@" . $db_inst . "\n";
    die;
  }

  for my $m (@{ $pw->UsedMols() }) {
    $used_mols{$m} = 1;
  }
  
  for my $m (keys %{ $parm->{molmap} }) {
    my $used;
    for my $m1 (keys %{ $parm->{molmap}{$m} }) {
      if (defined $used_mols{$m1}) {
        push @matched_used, "$m\t$m1\t";
        $used++;
        last;
      }
    }
    if (! $used) {
      push @matched_not_used, $m;

    }
  }
  for my $m (keys %{ $parm->{no_molmap} }) {
    push @not_matched, $m;
    push @matched_used,  "$m\t$m";

  }

  for my $p (keys %{ $parm->{pathwaymap} }) {
    for my $p1 (keys %{ $parm->{pathwaymap}{$p} }) {
      push @pmatched, "$p\t$p1";
    }
  }
  for my $p (keys %{ $parm->{no_pathwaymap} }) {
    push @not_pmatched, $p;
  }

  for my $mp (keys %{ $parm->{macroprocessmap} }) {
    for my $mp1 (keys %{ $parm->{macroprocessmap}{$mp} }) {
      push @macromatched, "$mp\t$mp1";
    }
  }
  for my $mp (keys %{ $parm->{no_macroprocessmap} }) {
    push @not_macromatched, $mp;
  }
 if ($format eq "html") {
    if (@matched_used ||
        @pmatched ||
        @macromatched) {
      my $first = 1; 
      for my $m (@matched_used) {
        my ($m1, $m2) = split(/\t/, $m);
        my $m3 = $pw->PickMolName($m2);
        my $mc = int($m2);
        if ($mc > 0) {
 
        if ($first == 1) {
          $first = 0;
          push @lines, "<p class=\"norm\">Molecule search terms found:  ";
        } else {
          push @lines, ",";
        }
        my $sql;
        my $stm;
        my $tmplocus;
        my $locus;
        my $symbol;
        my $name_type;
        my $officialCount = 0;
        my $official = '';
        my $aliasCount = 0;
        my $alias = '';
        
         
        $sql = "select ll_id from $schema.pw_ext_mol_id a, cgap.ll2sp b";
        $sql = $sql . " where mol_id = $m2 and a.ext_mol_id = b.sp_primary ";
        $sql = $sql . " and organism = 'Hs' ";
       
        print STDERR "SQL2:  $sql\n"; 
        $stm = $db->prepare($sql);
        if(not $stm) {
          print STDERR "$sql\n";
          print STDERR "$DBI::errstr\n";
          print STDERR "prepare call failed\n";
          die;
        }
        if(!$stm->execute()) {
          print STDERR "$sql\n";
          print STDERR "$DBI::errstr\n";
          print STDERR "execute call failed\n";
          die;
        }
        
        while (($tmplocus) = $stm->fetchrow_array()) {
          $locus = $tmplocus;
        }
        $stm->finish();
        if ($locus ne '') {
 # Get Official Symbol

  $sql = qq!
  select distinct b.symbol, b.name_type
  from $schema.pw_molecule_search b
  where ll_id = $locus
  and name_type in ('OF','AS')
  !;

  print STDERR "SQL:  $sql\n";
  $stm = $db->prepare($sql);
  if(not $stm) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "prepare call failed\n";
    die;
  }
  if(!$stm->execute()) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "prepare call failed\n";
    die;
  }
 while (($symbol, $name_type) =
    $stm->fetchrow_array()) {
      if ($name_type eq "OF") {
        $officialCount++;
        $official = $official . ' ' . $symbol;
      } else {
        $aliasCount++;
        $alias = ", " . $symbol . $alias ;
      }
  }

  $alias = substr($alias,2);
  }

        $stm->finish();
        my $hugo = "";
        if ($officialCount > 0) {
          $hugo = $hugo .  "(HUGO symbol:" . $official . ")";
        }
        my $hyplink="<a href=\"MoleculePage?molid=$m2\">$m1</a>";
        push @lines, "$hyplink $hugo";
      }
      }
      push @lines, "</p>";
      # for my $m (@matched_not_used) {
      #  my ($m1, $m2) = split(/\t/, $m);
      #  push @lines, "<tr>" .
      #      "<td>$m1</td>" .
      #      "<td>[no match]</td>" .
      #      "</tr>";
      #}
      # for my $m (@not_matched) {
      #  push @lines, "<tr>" .
      #      "<td>$m</td>" .
      #      "<td>[no match]</td>" .
      #      "</tr>";
      #}
      $first = 1; 
      for my $p (@pmatched) {
        my ($p1, $p2) = split(/\t/, $p);
        if ($first == 1) {
          $first = 0;
          push @lines, "<p class=\"norm\">Pathways terms found:  ";
        } else {
          push @lines, ",";
        }
        push @lines, "$p2";
      }
      push @lines, "</p>";
      #for my $p (@not_pmatched) {
      #  push @lines, "<tr>" .
      #      "<td>$p</td>" .
      #      "<td>[no match]</td>" .
      # "</tr>";
      #}
      # if ($#macromatched > 0) {
        $first = 1;
        for my $mp (@macromatched) {
           
          my ($mp1, $mp2) = split(/\t/, $mp);
          if ($first == 1) {
            push @lines, "<p class=\"norm\">Biological processes terms found:  ";
            $first = 0;
          } else {
            push @lines, ",";
          }

          push @lines, "$mp2 ";
        }
        push @lines, "</p>";
      # }
      #for my $mp (@not_macromatched) {
      #  push @lines, "<tr>" .
      #      "<td>$mp</td>" .
      #      "<td>[no match]</td>" .
      #      "</tr>";
      # }
    }
}
  $db->disconnect(); 
  return join("\n", @lines, "");
}

######################################################################
sub PrBatch_1 {
  my ($base, $graphic_type, $parm_string) = @_;
  my @lines;
  my (@rows);
  my %genes_a = ();
  my %genes_b = ();

  $BASE = $base;
  $parm_string = "$parm_string\n";    ## add just in case
  
  my $scancheck = Scan($parm_string);
  print STDERR "After Scan: parmstring $scancheck\n";
  if ($scancheck > 0) {
    return PrPathReturn(S_BAD_REQUEST, "text", "Bad characters in parameter string");
  }
  print STDERR "After check in PrBatch\n";

  my $parm = new PathParams($lv);
  if (! $parm->ReadParamString($parm_string)) {
    return PrPathReturn(S_BAD_REQUEST, "text",
        join("\n", @{ $parm->{errors} }) . "\n", "", "", $parm, undef);
  };

  print STDERR "After ReadparamString\n";
  my $db = DBI->connect("DBI:Oracle:" . $db_inst, $db_user, $db_pass);
  if (not $db or $db->err()) {
    print STDERR "Cannot connect to " . $db_user . "@" . $db_inst . "\n";
    die;
  }
  if (! $parm->CheckParams($db) ) {
    return PrPathReturn(S_BAD_REQUEST, "text",
        join("\n", @{ $parm->{errors} }) . "\n", "", "", $parm, undef);
  }
  print STDERR "After CheckParams\n";
 
  my $format = 'html';
  for my $mol_names(keys %{ $parm->{mol_name}}) {
    print STDERR "Mol_name:  $mol_names\n"; 
    for my $channel (keys %{ $parm->{molecule_channel}{$mol_names} }) {
      if ($channel eq 'A') {
        $genes_a{$mol_names} = 1;
      }
      if ($channel eq 'B') {
        $genes_b{$mol_names} = 1;
      }
 
    }
  }
  my $source = ''; 
  for my $test (keys %{ $parm->{source_id} }) {
    $source = $source . $test; 
  }
  print STDERR "Before RankPathwaysByGeneHits\n"; 
  @lines = RankPathwaysByGeneHits($db, $schema, $format, $source, \%genes_a, \%genes_b);
  # push @lines, "<p>";
  $db->disconnect();
 
  return join("\n", @lines) . "\n";

}

######################################################################
sub PrPath_1 {
  my ($base, $graphic_type, $parm_string) = @_;

  $BASE = $base;
  my $pathway;
  my $source;
  my $atomidlist;
  $parm_string = "$parm_string\n";    ## add just in case

  my $scancheck = Scan($parm_string);
  print STDERR "After Scan: parmstring $parm_string\n";
  if ($scancheck > 0) {
    return PrPathReturn(S_BAD_REQUEST, "text", "Bad characters in parameter string");
  } 
  my $parm = new PathParams($lv);
  print STDERR "Before ReadParam\n";
  if (! $parm->ReadParamString($parm_string)) {
    return PrPathReturn(S_BAD_REQUEST, "text",
        join("\n", @{ $parm->{errors} }) . "\n", "", "", $parm, undef);
  };
  my $atomcheck = $parm->{atom}{100758};

  my $db = DBI->connect("DBI:Oracle:" . $db_inst, $db_user, $db_pass);
  if (not $db or $db->err()) {
    print STDERR "Cannot connect to " . $db_user . "@" . $db_inst . "\n";
    die;
  }
  if (! $parm->CheckParams($db) ) {
    return PrPathReturn(S_BAD_REQUEST, "text",
        join("\n", @{ $parm->{errors} }) . "\n", "", "", $parm, undef);
  }
  #  return PrPathReturn(S_BAD_REQUEST, "text", "xamlreturn $parm_string");
  my $transcriptionflag = $parm->{transcription};
  if ($transcriptionflag eq 'transcription') {
    $transcriptionflag = 1;
  }

  if (! defined $parm->{print} &&
      ($graphic_type eq "svg" || $graphic_type eq "jpg" || $graphic_type eq "xaml" )) {
    $parm->{print}{"dot"} = 1;
  }

  my $pw   = new Pathway($lv);
  my $pdb;
  if ($parm->{buildxml}) {
      $pdb = new PathwayDB($db, $parm->{db}{schema},
      join(",", keys %{ $parm->{source_id} }),
      join(",", keys %{ $parm->{evidence_code} }), $pw, $lv);
  } else {
      $pdb = new PathwayXML($db, $parm->{db}{schema},
      $pw, $lv);

  }
  my $blanks = "                                        ";
  my $indent_level = 0;
  my $indent_step = 2;

  if (defined $parm->{pathway}) {
    for my $pid (keys %{ $parm->{pathway} }) {
      if ($pid > 0) {
        $pw->AddPathwayId($pid);
        $pw->SetPathwayName($pid, $parm->{pathway}{$pid}{name});
        $pathway = $pathway . " " . $pid;

        $pw->SetPathwayExId($pid, $parm->{pathway}{$pid}{ext_id});
        $pw->SetPathwayOrg($pid,  $parm->{pathway}{$pid}{org});
        $pw->SetPathwaySrcId($pid, $parm->{pathway}{$pid}{src_id});
        $source = $source . " " . $parm->{pathway}{$pid}{src_id};
      } 
      if ($pid > 0) {
        my $retval = $pdb->FillPathway($pid);
      }
    }
  }

#  if (defined $parm->{macro_process}) {
#    $pdb->FillMacroProcesses([ keys %{ $parm->{macro_process} } ]);
#  }

  my $source_id_list = [keys %{ $parm->{source_id} }];
  my $evidence_code_list = [keys %{ $parm->{evidence_code} }];
  if (defined $parm->{macroprocess_id}) {
    my $retval = $pdb->FillMacroProcesses(join(",", keys %{ $parm->{macroprocess_id} }), $evidence_code_list, $source_id_list);
    # return PrPathReturn(S_BAD_REQUEST, "text", "After Macro $retval");
 
  }

  if (defined $parm->{connect}) {
    for my $a (keys %{ FindConnectingPath($db, $parm->{connect},
        $parm->{source_id}, $parm->{molmap}) }) {
      $parm->{atom}{$a} = 1;
    }
  }

  if (defined $parm->{atom}) {
    my $atom_list = join(",", keys %{ $parm->{atom} });
    my $atom_id;

    for $atom_id (keys %{ $parm->{atom} }) {
      print STDERR "atom_id:  $atom_id\n";
    } 
    print STDERR "atom_list:  $parm->{atom}\n"; 
    my $retval = $pdb->FillAtoms($atom_list, $evidence_code_list);
    # return PrPathReturn(S_BAD_REQUEST, "text", "After FillAtom $retval");
 
    print STDERR "return from FillAtoms:  $atom_list\n";
  }

  ##
  ## First fetch basic
  ##
  my $list;
  if (defined $parm->{mol_id}) {
    if ($parm->{buildxml}) {
 
      my $include_complex_uses = $parm->{include_complex_uses};
      $pdb->AtomsOfMols(join(",", keys %{ $parm->{mol_id} }), $include_complex_uses);
    } else {
      $list = join(",", keys %{ $parm->{mol_id} });
      if (!$list) {
        $list = "'0'" . $list;
      }

my $sql = qq!
      select
  mm_inner_family.mol_id_2, mn.mol_name
from
  $schema.pw_edge e1,
  $schema.pw_mol_mol mm_outer_family,
  $schema.pw_mol_mol mm_inner_family,
  $schema.pw_mol_mol mm_complex,
  $schema.pw_mol_name mn
where
      mn.mol_id = mm_inner_family.mol_id_2
  and mm_inner_family.mol_id_2 = mm_complex.mol_id_2
  and mm_complex.mol_id_1 = mm_outer_family.mol_id_1
  and mm_complex.relation in ('c','i')
  and mm_inner_family.mol_id_2 = e1.mol_id
  and mm_outer_family.relation in ('s','m','i')
  and mm_inner_family.relation in ('s','m','i')
  and mm_outer_family.mol_id_1 in ($list)
!;
  # return $sql;
  my $stm = $db->prepare($sql);
  if (not $stm) {
    print STDERR "prepare call failed\n";
    $db->disconnect();
    die "$stm query Failed!";
  }
  if (! $stm->execute()) {
    print STDERR "execute call failed\n";
    $db->disconnect();
    die "$stm query Failed!";
  }
  my $molid;
  my $molname;
    while (($molid, $molname)
      = $stm->fetchrow_array()) {
      $parm->{mol_id}{$molid} = 1;
 
      if (defined $parm->{mol_name}{$molname}) {
        for my $what (keys %{ $parm->{mol_name}{$molname} }) {
          if ($what eq "prune") {
            $parm->{prune_mol_id}{$molid} = 1;
          } elsif ($what eq "get") {
            $parm->{mol_id}{$molid} = 1;
          } elsif ($what eq "connect") {
            $parm->{connect}{$molid} = 1;
          } elsif ($what eq "value") {
            # handled immediately below
          } elsif ($what eq "channel") {
            # handled immediately below
          } 
           if (defined $parm->{mol_name_value}{$molname}) {
            for my $subtype (keys %{ $parm->{mol_name_value}{$molname} }) {
              $parm->{mol_id_value}{$molid}{$subtype} =
                  $parm->{mol_name_value}{$molname}{$subtype};
            }
          }
          if (defined $parm->{molecule_channel}{$molname}) {
            for my $channel (keys %{ $parm->{molecule_channel}{$molname} }) {
              $parm->{mol_id_channel}{$molid}{$channel} = 1;
            }
          }

        }
      }     
      # $parm->{mol_id}{$molid} = 1;
  }
  $stm->finish();
          my $retval = $pdb->MolSearch($db, $schema, $transcriptionflag, [keys %{ $parm->{mol_id} }], [keys %{ $parm->{source_id} }], [keys %{ $parm->{evidence_code}}]);
#    # return PrPathReturn(S_BAD_REQUEST, "text", "After MolSearch $retval");
 
    }
}

  ##
  ## Next prune molecules
  ##

  if (defined $parm->{prune_mol_id}) {
    $pw->PruneMol($parm->{prune_mol_id});
    print STDERR "return from PruneMol\n";
  }

  if (defined $parm->{prune_atom_id}) {
    $pw->PruneAtoms([ keys %{ $parm->{prune_atom_id} }]);
  }

  ##
  ## Next fetch by degree (each iteration within
  ## AtomsOfDegree  will recall prune mols)
  ##

  if (defined $parm->{degree}) {

    my $degree;
    my $direction;
    if (defined $parm->{degree}{both}) {
      $direction = "both";
      $degree = $parm->{degree}{both}
    } else {
      if (defined $parm->{degree}{back}) {
        $direction = "back";
        $degree = $parm->{degree}{back}
      } elsif (defined $parm->{degree}{forward}) {
        $direction = "forward";
        $degree = $parm->{degree}{forward}
      }
    }
    $atomidlist = $pdb->AtomsOfDegree($direction, $degree, $parm->{prune_mol_id},$parm->{prune_atom_id});
    # return PrPathReturn(S_BAD_REQUEST, "text", "After Degree:  $direction Mols_in:  $atomidlist");



  }
###!!!!!

  # if (defined $parm->{privateedges}) {
  #  use PrivateData;
  #  my $pd = new PrivateData($pw, $lv, $db, $parm->{db}{schema});
  #  for my $edge_data (@{ $parm->{privateedges} }) {
  #    $pd->AddRawEdge(split("\t", $edge_data));
  #  }
  #  $pd->MapRawEdges($parm);
  #}

  if (!($parm->{buildxml})) {
    $pw->SplitTranscription();
  }
  
  $db->disconnect();

  if ($parm->{collapse}{molecules}) {
    $pw->CollapseMols();
  }

  # Testing Merging of Reactome Mols with Nature
  # $pw->MergeMolecules();
  # print STDERR "return from MergeMolecules\n";

  $pw->BuildMolInstCache();
  print STDERR "return from BuildMolInstCache\n";

  $pw->BuildGenericLabelCache();
  print STDERR "return from BuildGenericLabelCache\n";

  if (defined $parm->{print}{"dot"} ||
     $graphic_type eq "svg" || $graphic_type eq "jpg" || $graphic_type eq "xaml") {
    $pw->PruneDuplicateAtoms();
    print STDERR "return from PruneDuplicateAtoms\n";
  }

  $pw->IdentifyMacroProcesses();
  print STDERR "return from IdentifyMacroProcesses\n";

  $pw->BuildClanList();

  if (@{ $pw->Atoms() } == 0) {
    SetStatus(S_NO_DATA);
    return PrPathReturn(S_NO_DATA, "text", "No interactions found ", "",
        "", $parm, $pw);
  }

  
  # Stopped here return PrPathReturn(S_BAD_REQUEST, "text", "After no atoms");

  if (defined $parm->{sim} && $parm->{sim}{method} eq "boolean") {
 
    use PWBoolean;
    my @lines;
    my $trace_on; 
    my $output_type          = $parm->{sim}{output};
    if ($output_type eq "text") {
      $trace_on = 1;
    } else {
      $trace_on = 0;
    }
    if ($output_type eq "text") {
      for my $line (@{ $parm->ListParams() }) {
        push @lines, "##$line\n";
      }
    }
# we are going to dictate the color scale 
    my %boolean_coloring = (
      1  => "FF0000",
      0  => "0000FF",
      -1 => "000000"
    );
    use Coloring;
    my $coloring = new Coloring(MakeColorScale(\%boolean_coloring));
    my $default_value = "-1";
    if (defined $parm->{sim}{booleandefault}) {
      $default_value = $parm->{sim}{booleandefault};
    }
    my $sim = new PWBoolean($pw, $lv, $default_value ,$trace_on, "");
    $sim->SetUpSimTopology();
    $sim->InitializeMolValues($parm->{mol_id_value});
    $sim->Execute($parm->{sim}{ncycle});
    if ($output_type eq "svg" || $output_type eq "jpg" || $output_type eq "xaml") {
      my $ni = @{ $pw->Atoms };
      # FFFFF 
      my $mdi = MAX_DOT_INTERACTIONS; 
      return PrPathReturn(S_BAD_REQUEST, "text", "Before interactions:  $ni $mdi ");

      if ($ni > MAX_DOT_INTERACTIONS && ! $parm->{collapse}{processes}) {
        # SetStatus(S_BAD_REQUEST);
        return PrPathReturn(S_BAD_REQUEST, "text",
            "Number of interactions $ni exceeds maximum (" .
            MAX_DOT_INTERACTIONS . ")", "", "", $parm, $pw);
      }
      print STDERR "Going to GraphicsByClan\n";
      return PrPathReturn(S_OK, $output_type,
          GraphicsByClan($pw, $lv, $parm, $output_type, $sim, $sim,
          $coloring), join(",", @{ $pw->Atoms() }), $parm, $pw, 1);
    } else {
      return PrPathReturn(S_OK, "text",
          join("\n", @lines, @{ $sim->Lines() }) . "\n",
          "", join(",", @{ $pw->Atoms() }), $parm, $pw);
    } 

  } elsif (defined $parm->{sim}) {
 
    use PWSim;
    my @lines;
    my $method               = $parm->{sim}{method};
    my $output_type          = $parm->{sim}{output};
    my $all_or_nothing       = $parm->{sim}{simple};
    my $adjust_input_weights = $parm->{sim}{compete};
    my $trace_on;
    if ($output_type eq "text") {
      $trace_on = 1;
    } else {
      $trace_on = 0;
    }
    if ($output_type eq "text") {
      for my $line (@{ $parm->ListParams() }) {
        push @lines, "##$line\n";
      }
    }
    my $coloring;
    if ($parm->{color}) {
      use Coloring;
      $coloring = new Coloring(MakeColorScale($parm->{color}));
    }
    my $sim = new PWSim($pw, $lv, $method, $all_or_nothing,
        $adjust_input_weights, $trace_on, "");
    $sim->SetUpSimTopology();
##    $sim->InitializeDeviations($parm->{sim}{mol});
    $sim->InitializeDeviations($parm->{mol_id_value});
    $sim->Execute($parm->{sim}{ncycle});
    if ($output_type eq "svg" || $output_type eq "jpg" || $output_type eq "xaml") {
      my $ni = @{ $pw->Atoms };
      if ($ni > MAX_DOT_INTERACTIONS && ! $parm->{collapse}{processes}) {
        # SetStatus(S_BAD_REQUEST);
        return PrPathReturn(S_BAD_REQUEST, "text",
            "Number of interactions $ni exceeds maximum (" .
            MAX_DOT_INTERACTIONS . ")", "", "", $parm, $pw);
      }
      return PrPathReturn(S_OK, $output_type,
          GraphicsByClan($pw, $lv, $parm, $output_type, $sim, $sim,
          $coloring), join(",", @{ $pw->Atoms() }), $parm, $pw);
    } else {
      return PrPathReturn(S_OK, "text",
          join("\n", @lines, @{ $sim->Lines() }) . "\n",
          "", join(",", @{ $pw->Atoms() }), $parm, $pw);
    } 
  }
  

  if (defined $parm->{print}{"xml"}) {
    use XMLOutput;
    my @lines;
    my $xml = new XMLOutput($pw, $lv);
    $xml->PrXML();
    return PrPathReturn(S_OK, "text",
        join("\n", @lines, @{ $xml->Lines() }) . "\n",
        "", join(",", @{ $pw->Atoms() }), $parm, $pw);
  }

 
  if (defined $parm->{print}{"xaml"}) {
 
    # Create dot output with coords
    # run dot2xaml on created file
    # load xaml file
    # Return as @lines

    # use XMLOutput;
    # my @lines;
    # my $xml = new XMLOutput($pw, $lv);
    #$xml->PrXML();
    my @lines;
    # dot->DOTGraph();

      # Write dot file
      # Convert to xaml
      # Load xaml file
      # return as lines

    return PrPathReturn(S_OK, "text",
        join("\n", @lines ) . "\n",
        "", join(",", @{ $pw->Atoms() }), $parm, $pw);
  }
 
  if (defined $parm->{print}{"biopax"}) {
    use BioPAXOutput;
    my @lines;
    my $biopax = new BioPAXOutput($pw, $lv);
    $biopax->PrOWL();
    return PrPathReturn(S_OK, "text",
        join("\n", @lines, @{ $biopax->Lines() }) . "\n",
        "", join(",", @{ $pw->Atoms() }), $parm, $pw);
  }
  
  if (defined $parm->{print}{"sql"}) {
    use SQLOutput;
    my @lines;
    my $sql = new SQLOutput($pw, $lv);
    $sql->PrSQL();
    return PrPathReturn(S_OK, "text",
        join("\n", @lines, @{ $sql->Lines() }) . "\n",
        "", join(",", @{ $pw->Atoms() }), $parm, $pw);
  }
  
  # if (defined $parm->{print}{"sbml"}) {
  #  use SBMLOutput;
  #  my @lines;
  #  my $sbml = new SBMLOutput($pw, $lv);
  #  $sbml->PrSBML();
  #  return PrPathReturn(S_OK, "text",
  #      join("\n", @lines, @{ $sbml->Lines() }) . "\n",
  #      "", join(",", @{ $pw->Atoms() }), $parm, $pw);
  # }
  
  if (defined $parm->{print}{"dot"}) {
    use DOToutput;
    my @lines;
 
    for my $line (@{ $parm->ListParams() }) {
      push @lines, "//$line";
    }
#    print join("\n", @lines) . "\n";

    my $coloring;
    my $transform;
    my ($mol_vg, $atom_vg);
    if ($parm->{color}) {
      use Coloring;
      $coloring = new Coloring(MakeColorScale($parm->{color}));
    }

    if (defined $parm->{mol_id_value}) {
      use ValueGenerator;
      if (defined $coloring) {
        $transform = $coloring->NumericBreaks;
      }
      $mol_vg = new ValueGenerator(5, 5, undef, $parm->{mol_id_value},
          undef);
      for my $mol_id (keys %{ $parm->{mol_id_value} }) {
        for my $subtype (keys %{ $parm->{mol_id_value}{$mol_id} }) {
          $mol_vg->SetMolVal(undef, $mol_id, $parm->{mol_id_value}{$mol_id}{$subtype});
        }
      }
#          $transform);
    }
    
    ## propagate to containing complexes
    if (defined $parm->{mol_id_channel}) {
     for my $mol (@{ $pw->Mols() }) {
        if ($pw->MolType($mol) eq $COMPLEX_TYPE) {
          for my $comp (@{ $pw->Components($mol) }) {
            my $m = $pw->ComponentMol($mol, $comp);
            if (defined $parm->{mol_id_channel}{$m}) {
              for my $c (keys %{ $parm->{mol_id_channel}{$m} }) {
               $parm->{mol_id_channel}{$mol}{$c} = 1;
              }
            }
          }
        }
      }

      ## propagate to family ancestors
      for my $m (keys %{ $parm->{mol_id_channel} }) {
        for my $mol (keys %{ $pw->FamilyAncestors($m) }) {
          for my $c (keys %{ $parm->{mol_id_channel}{$m} }) {
            $parm->{mol_id_channel}{$mol}{$c} = 1;
          }
        }
      }

      my $CHANNEL_3 = [
          "0000FF",
          "CC00FF",
          "FF0000" 
      ];
      $coloring = new Coloring
          ([1, 2, 3], $CHANNEL_3,);
      $mol_vg = new ValueGenerator(5, 5, undef, $parm->{mol_id_value},
          undef);
      for my $mol (keys %{ $parm->{mol_id_channel} }) {
        if (defined $parm->{mol_id_channel}{$mol}{A} &&
            defined $parm->{mol_id_channel}{$mol}{B}) {
          $mol_vg->SetMolVal(undef, $mol, 2);
        } elsif (defined $parm->{mol_id_channel}{$mol}{A}) {
          $mol_vg->SetMolVal(undef, $mol, 1);
        } elsif (defined $parm->{mol_id_channel}{$mol}{B}) {
          $mol_vg->SetMolVal(undef, $mol, 3);
        }
      }
    }

    if ($graphic_type eq "" ||
        ($graphic_type ne "xaml" && $graphic_type ne "svg" && $graphic_type ne "jpg")) {
#      my $dot = new DOToutput($pw, $lv, *STDOUT);
      my $dot = new DOToutput($pw, $lv, "");
      if ($parm->{collapse}{molecules}) {
        $dot->CollapseMolecules(1);
      }
      if ($parm->{show}{atom_id}) {
        $dot->ShowAtomIds(1);
      }
      if ($parm->{show}{subtypelines}) {
        $dot->ShowSubTypeLines(1);
      }
      if ($parm->{collapse}{process}) {
        $dot->CollapseAtoms(1);
      }

      if (defined $mol_vg) {
        $dot->SetMolValueGenerator($mol_vg);
        if (defined $coloring) {
          $dot->SetSimColoring($coloring);
        }
      } else {
        $dot->ColorMols($parm->{mol_id});
        if (defined $parm->{connect}) {
          $dot->ColorMols($parm->{connect});
          my %temp;
          for my $mol (@{ $pw->Mols() }) {
            if ($pw->MolType($mol) eq $COMPLEX_TYPE) {
              for my $comp (@{ $pw->Components($mol) }) {
                my $m = $pw->ComponentMol($mol, $comp);
                if (defined $parm->{connect}{$m}) {
                  $temp{$mol} = 1;
                }
              }
            }
          }
          $dot->ColorMols(\%temp);
        }
      }
      if (defined $atom_vg) {
        $dot->SetAtomValueGenerator($atom_vg);
      }
      $dot->DOTGraph();

      # Write dot file
      # Convert to xaml
      # Load xaml file
      # return as lines
      print STDERR "Returning DOT output\n";
      return PrPathReturn(S_OK, "text",
          join("\n", @lines, @{ $dot->Lines() }) . "\n",
          "", join(",", @{ $pw->Atoms() }), $parm, $pw);
     
    } else {
      my $ni = @{ $pw->Atoms };
      if ($ni > MAX_DOT_INTERACTIONS && ! $parm->{collapse}{processes}) {
        SetStatus(S_BAD_REQUEST);
        return PrPathReturn(S_BAD_REQUEST, "text",
            "Number of interactions $ni exceeds maximum (" .
            MAX_DOT_INTERACTIONS . ")", "", "", $parm, $pw);
      }

      print STDERR "Hopefully returning jpg and svg output $graphic_type\n"; 
      return PrPathReturn(S_OK, $graphic_type,
          GraphicsByClan($pw, $lv, $parm, $graphic_type,
              $mol_vg, $atom_vg, $coloring),
          join(",", @{ $pw->Atoms() }), $parm, $pw);
    }
  }
  
}

######################################################################
# BEGIN AtomPage
######################################################################

######################################################################
sub PUBMED_URL {
  my ($pubmed_id) = @_;
  return "http://www.ncbi.nlm.nih.gov:80/entrez/query.fcgi?".
      "cmd=Retrieve\&db=PubMed\&list_uids=$pubmed_id\&dopt=Abstract";
}

######################################################################
sub ATOMPAGE_URL {
  my ($atomid) = @_;
  return "/InteractionPage?atomid=$atomid&what=graphic&jpg=on";
}

######################################################################
sub PathwayInfoForAtom {
  my ($db, $atomid, $lines) = @_;

  my ($sql, $stm);
  my ($pid, $p_name, $p_src, $p_ex_id);

  push @{ $lines }, "<p><b>Uses in pathways:</b>";
  push @{ $lines }, "<blockquote>";
  push @{ $lines }, "<table border=1 cellspacing=1 cellpadding=4>";
  push @{ $lines }, "<tr><td><b>Source</b></td>" .
        "<td><b>Pathway name</b></td></tr>";

  $sql = qq!
select
  p.pathway_id,
  p.pathway_name,
  s.source_name,
  p.ext_pathway_id
from
  $schema.pw_pathway_atom pa,
  $schema.pw_pathway p,
  $schema.pw_source s
where
      p.pathway_id = pa.pathway_id
  and pa.atom_id = $atomid
  and s.source_id = p.pathway_source_id
  !;

  $stm = $db->prepare($sql);
  if(not $stm) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "prepare call failed\n";
    die;
  }
  if(!$stm->execute()) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "execute call failed\n";
    die;
  }
  while (($pid, $p_name, $p_src, $p_ex_id)
      = $stm->fetchrow_array()) {
    $p_src =~ s/NATURE/NCI-Nature curated/;
    $p_src =~ s/BioCarta/BioCarta Imported/;
    push @{ $lines }, "<tr><td>$p_src</td>" .
        "<td>$p_name&nbsp;&nbsp<a href=\"" . PATHWAY_URL($pid, "jpg",$p_src) .
        "\" >JPG</a>&nbsp;&nbsp;<a href=\"" .
        PATHWAY_URL($pid, "svg",$p_src) . "\" >SVG</a></td>" .
        "</tr>";
  }

  push @{ $lines }, "</table>";
  push @{ $lines }, "</blockquote>";

}

######################################################################
sub SearchDupAtoms {
  my ($db, $pdb, $pw, $atomid, $lines) = @_;

  my ($sql, $stm);
  my ($a_id, $m_id);
  my ($one_mol_id);
  my (@candidates, @temp, %mols, %summary, $standard);

  $one_mol_id = 'none';

  for (@{ $pw->Edges($atomid) }) {
    $one_mol_id = $pw->EdgeMol($atomid, $_);
    last;
  }
  if ($one_mol_id eq 'none') {}else{
  $sql = qq!
select unique
  a.atom_id,
  e.mol_id
from
  $schema.pw_atom a,
  $schema.pw_edge e,
  $schema.pw_atom_label al
where
      e.mol_id = $one_mol_id
  and a.atom_id = al.atom_id
  and a.atom_id = e.atom_id
  !;

  $stm = $db->prepare($sql);
  if(not $stm) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "prepare call failed\n";
    die;
  }
  if(!$stm->execute()) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "execute call failed\n";
    die;
  }
  while (($a_id, $m_id) = $stm->fetchrow_array()) {
    $mols{$a_id}{$m_id}++;
  }
  $stm->finish();
  }
  for $a_id (keys %mols) {
    undef @temp;
    for $m_id (sort keys %{ $mols{$a_id} }) {
      push @temp, "$m_id:$mols{$a_id}{$m_id}";
    }
    $summary{$a_id} = join("::", @temp);
  }
  $standard = $summary{$atomid};
  for $a_id (keys %summary) {
    if ($summary{$a_id} eq $standard &&
        $a_id != $atomid) {
      push @candidates, $a_id;
    }
  }

  undef @temp;
  if (@candidates) {
    $pdb->FillAtoms(join(",", @candidates));
    $pw->BuildMolInstCache();
    $pw->BuildGenericLabelCache();
    $standard = $pw->NormalAtomString($atomid);
    for $a_id (@candidates) {
      if ($pw->NormalAtomString($a_id) eq $standard) {
        push @temp, "<tr><td><a href= " .
            ATOMPAGE_URL($a_id) . " >$a_id</a></td></tr>";
      }
    }
  }

  if (@temp) {
    push @{ $lines }, "<p>Duplicate interactions:";
    push @{ $lines }, "<blockquote>";
    push @{ $lines }, "<table border=1 cellspacing=1 cellpadding=4>";
    push @{ $lines }, @temp;
    push @{ $lines }, "<table>";
    push @{ $lines }, "</blockquote>";
  }
}

######################################################################
sub BasicInfoForAtom {
  my ($pw, $atomid, $org, $srcid, $srcname, $lines) = @_;

  my ($label, $atomtype, $lvid, $label_value);

  push @{ $lines }, "<b>Identification:</b>";
  push @{ $lines }, "<blockquote>";
  push @{ $lines }, "<table border=1 cellspacing=1 cellpadding=4>";
  push @{ $lines }, "<tr><td><b>Interaction id:</b></td>" .
      "<td>$atomid&nbsp;&nbsp; ";

  ($label, $atomtype) = $lv->LabelValueToString($pw->AtomType($atomid));
  if ($atomtype eq "binding") {
    $atomtype = "modification";
  }
  if ($label eq 'process-type') {
    $label='Interaction type';
  }
  # push @{ $lines }, "<tr><td><b>$label:</b></td><td>$atomtype" .
  #    " ($org)</td><tr>";
  
  push @{ $lines }, "<tr><td><b>$label:</b></td><td>$atomtype</td></tr>";
  push @{ $lines }, "<tr><td><b>Data source:</b></td><td>$srcname</td><tr>";
  for $lvid (@{ $pw->AtomLabel($atomid) }) {
    if ($lvid == $pw->AtomType($atomid)) {
      next;
    }
    ($label, $label_value) = $lv->LabelValueToString($lvid);
    push @{ $lines }, "<tr><td><b>$label:</b></td><td>$label_value</td><tr>";
  }
  for $lvid (@{ $pw->AtomCondition($atomid) }) {
    ($label, $label_value) = $lv->LabelValueToString($lvid);
    push @{ $lines }, "<tr><td><b>Positive Condition:</b></td>" .
        "<td>$label_value</td><tr>";
  }
  for $lvid (@{ $pw->AtomNegativeCondition($atomid) }) {
    ($label, $label_value) = $lv->LabelValueToString($lvid);
    push @{ $lines }, "<tr><td><b>Negative condition:</b></td>" .
        "<td>$label_value</td><tr>";
  }

  push @{ $lines }, "<p><b>Display interaction:</b> <a href=\"" .
      INTERACTION_URL($atomid, "jpg") .
      "\" >JPG</a>&nbsp;&nbsp;<a href=\"" .
      INTERACTION_URL($atomid, "svg") .
      "\" >SVG</a>&nbsp;&nbsp;<a href=\"" .
      INTERACTION_TEXT_URL($atomid, "xml") .
      "\" >XML</a>&nbsp;&nbsp;<a href=\"" .
      INTERACTION_TEXT_URL($atomid, "biopax") .
      "\" >BioPAX</a>";
 
  push @{ $lines }, "</table>";
  push @{ $lines }, "</blockquote>";

}

######################################################################
sub ReferencesForAtom {
  my ($db, $atomid, $lines) = @_;

  my ($sql, $stm);
  my ($pathway_id, $edge_seq_id, $mol_id, $pmid, $reference, $source_id, $pm_text);
  my (@temp);

  $db->{LongReadLen} = MAX_LONG_LEN;
  
  $sql = qq!
select
  pathway_id, edge_seq_id, mol_id, a.pmid, reference, source_id, pm_text
from
  $schema.pw_references a,$schema.pw_pubmed b
where
  atom_id = $atomid
  and a.pmid = b.pmid
  !;

  $stm = $db->prepare($sql);
  if(not $stm) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "prepare call failed\n";
    die;
  }
  if(!$stm->execute()) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "execute call failed\n";
    die;
  }
  while (($pathway_id, $edge_seq_id, $mol_id, $pmid, $reference, $source_id, $pm_text) =
      $stm->fetchrow_array()) {
    if ($pmid) {
      $reference .= ("$pm_text <a href=\"" . PUBMED_URL($pmid) . "\" > PMID:$pmid</a>");
    }
    push @temp, "<li>$reference";
  }
  $stm->finish();
  if (@temp) {
    push @{ $lines }, "<p><b>References:</b>";
    push @{ $lines }, "<blockquote>";
    push @{ $lines }, "<ul>";
    push @{ $lines }, @temp;
    push @{ $lines }, "</ul>";
    push @{ $lines }, "</blockquote>";
  }

}

######################################################################
sub NotesForAtom {
  my ($db, $atomid, $lines) = @_;

  my ($sql, $stm);
  my ($pathway_id, $edge_seq_id, $mol_id, $note, $source_id);
  my (@temp);

  $sql = qq!
select
  pathway_id, edge_seq_id, mol_id, note, source_id
from
  $schema.pw_notes
where
  atom_id = $atomid
  !;

  $stm = $db->prepare($sql);
  if(not $stm) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "prepare call failed\n";
    die;
  }
  if(!$stm->execute()) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "execute call failed\n";
    die;
  }
  while (($pathway_id, $edge_seq_id, $mol_id, $note, $source_id) =
      $stm->fetchrow_array()) {
    push @temp, "<p>$note\n";
  }
  $stm->finish();
  if (@temp) {
    # push @{ $lines }, "<p><b>Notes:</b>";
    # push @{ $lines }, "<blockquote>";
    #push @{ $lines }, @temp;
    #push @{ $lines }, "</blockquote>";
  }

}

######################################################################
sub EvidenceForAtom {
  my ($db, $atomid, $lines) = @_;

  my ($sql, $stm);
  my ($evidence_code, $evidence_descr);
  my (@temp);

  $sql = qq!
 select
  e.evidence_code, d.evidence_descr
 from
  $schema.pw_evidence e,
  $schema.pw_evidence_descr d
 where
  e.atom_id = $atomid
 and
  e.evidence_code = d.evidence_code
  !;

  $stm = $db->prepare($sql);
  if(not $stm) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "prepare call failed\n";
    die;
  }
  if(!$stm->execute()) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "execute call failed\n";
    die;
  }
  while (($evidence_code, $evidence_descr) =
      $stm->fetchrow_array()) {
    push @temp, "<li>$evidence_code - <i>$evidence_descr</i>\n";
  }
  $stm->finish();
  if (@temp) {
    push @{ $lines }, "<p><b>Evidence:</b>";
    push @{ $lines }, "<blockquote>";
    push @{ $lines }, "<ul>";
    push @{ $lines }, @temp;
    push @{ $lines }, "</ul>";
    push @{ $lines }, "</blockquote>";
  }

}

######################################################################
sub EdgeInfoForAtom {
  my ($pw, $atomid, $lines) = @_;

  my ($lvid, $label, $value, $edge, $molid, $molname, $role,
      $location, $state, %tmp);

  $molid = 'none';
  for $edge (@{ $pw->Edges($atomid) }) {
    $molid = $pw->EdgeMol($atomid, $edge);
  } 
  if ($molid eq 'none') {}else{ 
  push @{ $lines }, "<p><b>Participants:</b>";
  push @{ $lines }, "<blockquote>";
  push @{ $lines }, "<table border=1 cellspacing=1 cellpadding=4>";

  push @{ $lines }, "<tr>" .
      "<td><b>Role</b></td>" .
      "<td><b>Molecule</b></td>" .
      "<td><b>Location</b></td>" .
      "<td><b>State</b></td>" .
      "<td><b>Modifications<b></td>" .
      "</tr>";

  for $edge (@{ $pw->Edges($atomid) }) {
    undef $state;
    undef $location;
    $molid = $pw->EdgeMol($atomid, $edge);
    $molname = MoleculeName($pw, $molid);
    for $lvid (@{ $pw->EdgeLabel($atomid, $edge) }) {
      ($label, $value) = $lv->LabelValueToString($lvid);
      if ($label eq "edge-type") {
        $role = $value;
      }
    }
    for $lvid (@{ $pw->MolLabel($atomid, $edge) }) {
      ($label, $value) = $lv->LabelValueToString($lvid);
      if ($label eq "location") {
        $location = $value;
      } elsif ($label eq "activity-state") {
        $state = $value;
      }
    }
    my @ptms = @{ PTMStuff($pw, $atomid, $edge) };
    my $ptms = join("<br>", @ptms);
    $ptms or $ptms = "&nbsp;";
    $location or $location = "&nbsp;";
    $state or $state   = "&nbsp;";
    $role =~ s/agent/positive regulator/;
    $role =~ s/inhibitor/negative regulator/;

    push @{ $tmp{"$role|$location|$state|"} }, "<tr><td>" .
      join("</td>\n<td>",
        $role,
        "<a href=" . MOLPAGE_URL($molid) . " >$molname</a>",
        $location,
        $state,
        $ptms
      ) . "</td></tr>";
  }
  for my $i (sort keys %tmp) {
    for (@{ $tmp{$i} }) {
      push @{ $lines }, $_;
    }
  }

  push @{ $lines }, "</table>";
  push @{ $lines }, "</blockquote>";
  }
}

######################################################################
sub AtomPage_1 {
  my ($base, $atomid) = @_;

  $BASE = $base;

   my (@lines);
  my ($label, $label_value, $lvid, $atomtype, $pathwayid, $edge);

  my $db = DBI->connect("DBI:Oracle:" . $db_inst, $db_user, $db_pass);
  if (not $db or $db->err()) {
    print STDERR "Cannot connect to " . $db_user . "@" . $db_inst . "\n";
    die;
  }

  my $pw   = new Pathway($lv);
  my $pdb  = new PathwayDB($db, $schema,
                           $LIMIT_SOURCES, $LIMIT_EVIDENCE, $pw, $lv);
  my ($org, $srcid, $srcname) = $pdb->GetAtomBasics($atomid);
  $srcname =~ s/NATURE/NCI-Nature curated/;
  $srcname =~ s/BioCarta/BioCarta Imported/;
  $pdb->FillAtoms($atomid);

  BasicInfoForAtom($pw, $atomid, $org, $srcid, $srcname, \@lines);
  EdgeInfoForAtom($pw, $atomid, \@lines);
  NotesForAtom($db, $atomid, \@lines);
  ReferencesForAtom($db, $atomid, \@lines);
  EvidenceForAtom($db, $atomid, \@lines);
  PathwayInfoForAtom($db, $atomid, \@lines);
  SearchDupAtoms($db, $pdb, $pw, $atomid, \@lines);

  $db->disconnect();

  return join("\n", @lines) . "\n";
}

######################################################################
# BEGIN MoleculePage
######################################################################

#!!! ouch... global variables for MoleculePage

my (%atom_uses, %complex_uses, %location, %state,
    %atom2pathway, %atom2source, %pathway2name, %pathway2source);


######################################################################
sub LL_URL {
  my ($locus) = @_;
  return "http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?" .
      "db=gene&cmd=Retrieve&dopt=full_report&" .
      "list_uids=$locus";
}

######################################################################
sub UNIPROT_URL {
  my ($uniprot) = @_;
  return "http://www.uniprot.org/uniprot/$uniprot"
}

######################################################################
sub MOLPAGE_URL {
  my ($molid) = @_;
  return "/MoleculePage?molid=$molid";
}

######################################################################
sub MOLINST_URL {
  my ($atom, $edge) = @_;
  return "/MoleculeInstance?inst=$atom,$edge";
}

######################################################################
sub CGAP_URL {
  my ($org, $cid) = @_;
  if ($cid =~ /,/) {
    my @temp;
    for my $c (split(/,/, $cid)) {
      push @temp, "$org.$c";
    }
    my $term = join(",", @temp);
    return "http://cgap.nci.nih.gov/Genes/RunUniGeneQuery?PAGE=1" .
        "&ORG=$org&TERM=$term";
  } else {
    return "http://cgap.nci.nih.gov/Genes/GeneInfo?ORG=$org&CID=$cid";
  }
}

######################################################################
sub INTERACTION_URL {
  my ($atomid, $what) = @_;
  return "/search/interaction_landing.shtml?" .
      "atom_id=$atomid&what=graphic&$what=on";
}
######################################################################
sub INTERACTION_TEXT_URL {
  my ($atomid, $what) = @_;
  return "/search/interaction_landing.shtml?" .
      "atom_id=$atomid&what=text&$what=on";
}

######################################################################
sub ATOM_URL {
  my ($atomid, $what) = @_;
  return "/search/search_landing.shtml?" .
      "atom_id=$atomid&what=graphic&$what=on";
}
######################################################################
sub ATOM_TEXT_URL {
  my ($atomid, $what) = @_;
  return "/search/search_landing.shtml?" .
      "atom_id=$atomid&what=text&$what=on";
}

######################################################################
sub PATHWAY_URL {
  my ($pathwayid, $what, $source) = @_;
  return "/search/pathway_landing.shtml?" .
      "pathway_id=$pathwayid&source=$source&ppage=1&what=graphic&$what=on";
}

######################################################################
sub MoleculeInstance_1 {
  my ($base, $inst) = @_;

  $BASE = $base;
  ## apparently can't get a '&' through dot
  my ($atom, $edge) = split(",", $inst);
  my $db = DBI->connect("DBI:Oracle:" . $db_inst, $db_user, $db_pass);
  if (not $db or $db->err()) {
    print STDERR "Cannot connect to " . $db_user . "@" . $db_inst . "\n";
    die;
  }
  my $pw   = new Pathway($lv);
  my $pdb  = new PathwayDB($db, $schema,
                           $LIMIT_SOURCES, $LIMIT_EVIDENCE, $pw, $lv);
  $pdb->FillAtoms($atom);
  my $molid = $pw->EdgeMol($atom, $edge);

  my @lines;

  SimpleOrComplex($pw, $molid);
  if (keys %atom_uses) {
    PathwayAndSource($db, $schema);
  }

  my $edgetype = $pw->EdgeType($atom, $edge);
  my ($label, $value) = $lv->LabelValueToString($edgetype);
  $value =~ s/agent/positive regulator/;
  $value =~ s/inhibitor/negative regulator/;

  push @lines, "<p><b>Role:</b> $value";

  IdStuff($db, $pw, $molid, \@lines);
  if ($pw->MolType($molid) == $COMPLEX_TYPE) {
    ComplexInfo($db, $pw, $molid, \@lines);
  }

  my @location;
  my @activity;
  for my $lvid (@{ $pw->MolLabel($atom, $edge) }) {
    my ($label, $value) = $lv->LabelValueToString($lvid);
    if ($label eq "location") {
      push @location, $value;
    } elsif ($label eq "activity-state") {
      push @activity, $value;
    }
  }
  my @ptms = @{ PTMStuff($pw, $atom, $edge) };

  if (@location > 0) {
    push @lines, "<p><b>Location:</b> " . join(",", @location);
  }

  if (@activity > 0) {
    push @lines, "<p><b>Activity state:</b> " . join(",", @activity);
  }

  if (@ptms > 0) {
    push @lines, "<p><b>Post-Translational modifications:</b>";
    push @lines, "<blockquote>";
    push @lines, join("<br>", @ptms);
    push @lines, "</blockquote>";
  }
  my @names; 
  my $nhash = $pw->MolName($molid);
  for my $name_type (keys %{ $nhash }) {
    for my $name (keys %{ $$nhash{$name_type} }) {
      push @names, $name;
    }
  }
  my $molNames = join(", ", sort @names);
  push @lines, "<p>View the <a href=\"" . MOLPAGE_URL($molid) . "\">$molNames Molecule Page</a> for more information";
  $db->disconnect();

  return join("\n", @lines) . "\n";

}

######################################################################
sub MoleculePage_1 {
  my ($base, $molid) = @_;

  $BASE = $base;

  ## clean out those globals

  undef %atom_uses;
  undef %complex_uses;
  undef %location;
  undef %state;
  undef %atom2pathway;
  undef %atom2source;
  undef %pathway2name;
  undef %pathway2source;

  my ($sql, $stm, $symbol, $name_type, $official, $alias);
  
  my $db = DBI->connect("DBI:Oracle:" . $db_inst, $db_user, $db_pass);
  if (not $db or $db->err()) {
    print STDERR "Cannot connect to " . $db_user . "@" . $db_inst . "\n";
    die;
  }

  my $pw   = new Pathway($lv);
  my $pdb  = new PathwayDB($db, $schema,
                           $LIMIT_SOURCES, $LIMIT_EVIDENCE, $pw, $lv);
  my @lines;

  my $INCLUDE_COMPLEX_USES = 1;
  ## elsewhere in $parm->{include_complex_uses}

  $pdb->AtomsOfMols($molid, $INCLUDE_COMPLEX_USES);
                                   
  SimpleOrComplex($pw, $molid);
  if (keys %atom_uses) {
    PathwayAndSource($db, $schema);
  }
  IdStuff($db, $pw, $molid, \@lines);
  if ($pw->MolType($molid) == $COMPLEX_TYPE) {
    ComplexInfo($db, $pw, $molid, \@lines);
  }
  GOStuff($db, $pw, $molid,\@lines);
  ComplexTable($pw, $molid, \@lines);
  if (keys %atom_uses) {
    AtomTable($pw, $molid, \@lines);
  }

  $db->disconnect();

  return join("\n", @lines) . "\n";
}

######################################################################
sub PTMStuff {
  my ($pw, $atom, $edge) = @_;

  my @ptms;
  for my $x (@{ $pw->EdgePTM($atom, $edge) }) {
    my ($protein_id, $position, $amino_acid, $modification_label_value,
        $modification_label_name) = @{ $x };
    push @ptms, "$protein_id \[$amino_acid $position $modification_label_name\]";
  }
  return \@ptms;
}

######################################################################
sub GOStuff {
  my ($db, $pw, $molid, $lines) = @_;


##!!!!!!
my $schema = "cgap";

  my ($locus, $uniprot);
  my @temp;
  my ($sql, $stm);
  my ($go_id, $go_name, $evidence, @lines);
  
  my $ehash = $pw->MolExId($molid);
  for my $id_type (keys %{ $ehash }) {
    for my $id (keys %{ $$ehash{$id_type} }) {
      if ($id_type eq "LL") {
        $locus = $id;
      }
    }
  }

  for my $id_type (keys %{ $ehash }) {
    for my $id (keys %{ $$ehash{$id_type} }) {
      if ($id_type eq "UP") {
        $uniprot = $id;
      }
    }
  }

  
if (! defined $locus || ! $locus) {} else {
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
    die;
  }
  if(!$stm->execute()) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "execute call failed\n";
    die;
  }
  while (($go_id, $go_name, $evidence) = $stm->fetchrow_array()) {
    push @temp, "<tr><td><a href=http://www.godatabase.org/cgi-bin/amigo/go.cgi?view=details\&query=GO:" . "$go_id" . ">GO:$go_id</a></td><td>$go_name</td></tr>";
  }

  $stm->finish();
}

 if (! defined $uniprot || ! $uniprot) {} else {

   $sql = qq!
select unique
  n.go_id,
  n.go_name
from
  $schema.ll_go g,
  $schema.go_name n,
  $schema.sp_primary sp,
  $schema.ll2sp l
where
      sp.sp_id_or_secondary = '$uniprot'
  and l.sp_primary = sp.sp_primary
  and g.ll_id = l.ll_id
  and g.go_id = n.go_id
  and n.go_class = 'MF'
order by
  n.go_name
!;

  $stm = $db->prepare($sql);
  if(not $stm) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "prepare call failed\n";
    die;
  }
  if(!$stm->execute()) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "execute call failed\n";
    die;
  }
  while (($go_id, $go_name) = $stm->fetchrow_array()) {
    push @temp, "<tr><td><a href=http://www.godatabase.org/cgi-bin/amigo/go.cgi?view=details\&query=GO:" . "$go_id" . ">GO:$go_id</a></td><td>$go_name</td></tr>";
  }
}
  if (@temp) {
    push @{ $lines }, "<p>Gene Ontology (GO) molecular function annotations:";
    push @{ $lines }, "<blockquote>";
    push @{ $lines }, "<p><table border=1 cellspacing=1 cellpadding=4>";
    push @{ $lines }, "<tr><td><b>GO Id</b></td><td><b>GO term</b></td></tr>";
    push @{ $lines }, join("", @temp);
    push @{ $lines }, "</table>";
    push @{ $lines }, "</blockquote>";
  }

}

######################################################################
sub PathwayAndSource {
  my ($db, $schema) = @_;

  my ($sql, $stm);
  my (%pathways);
  my ($atom_id, $pathway_id, $pathway_name, $source);

  $sql =
    "select pa.atom_id, pa.pathway_id " .
    "from $schema.pw_pathway_atom pa " .
    ($LIMIT_SOURCES ?
    ", $schema.pw_atom a " : "") .
    ($LIMIT_EVIDENCE ? 
    ", $schema.pw_evidence e " : "") .
    "where pa.atom_id in (" . join(",", keys %atom_uses) . ") " .
    ($LIMIT_SOURCES ?
     "and pa.atom_id = a.atom_id " .
     "and a.atom_source_id in ($LIMIT_SOURCES) " : "") .
    ($LIMIT_EVIDENCE ? 
      "and pa.atom_id = e.atom_id " .
      "and e.evidence_code in ($LIMIT_EVIDENCE) " : "") ;

  $stm = $db->prepare($sql);
  if(not $stm) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "prepare call failed\n";
    die;
  }
  if(!$stm->execute()) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "execute call failed\n";
    die;
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
    die;
  }
  if(!$stm->execute()) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "execute call failed\n";
    die;
  }
  while (($atom_id, $source) = $stm->fetchrow_array()) {
    $atom2source{$atom_id} = $source;
  }

  if (scalar(keys %pathways) == 0) {
    return;
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
    die;
  }
  if(!$stm->execute()) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "execute call failed\n";
    die;
  }
  while (($pathway_id, $pathway_name, $source) = $stm->fetchrow_array()) {
    $pathway2source{$pathway_id} = $source;
    $pathway2name{$pathway_id}   = $pathway_name;
  }

}

######################################################################
sub IdStuff {
  my ($db, $pw, $molid, $lines) = @_;

  my (@names, @ids, $locus, $uniprot, $go, $ca, $sql, $stm, $tmplocus, $sg);
##!! hard code for now
  my ($official, $alias, $symbol, $name_type);
  my $org = "Hs";
  my $aliasCount = 0;
  my $officialCount = 0;
  my $nhash = $pw->MolName($molid);
  my $namepush = 0;
  for my $name_type (keys %{ $nhash }) {
    for my $name (keys %{ $$nhash{$name_type} }) {
      push @names, $name;
      $namepush = 1;
    }
  }
  if ($namepush == 0) {
    $sql = "select mol_name from $schema.pw_mol_name where mol_id = $molid";
    $stm = $db->prepare($sql);
    if(not $stm) {
      print STDERR "$sql\n";
      print STDERR "$DBI::errstr\n";
      print STDERR "prepare call failed\n";
      die;
    }
    if(!$stm->execute()) {
      print STDERR "$sql\n";
      print STDERR "$DBI::errstr\n";
      print STDERR "execute call failed\n";
      die;
    }

   while ((my $newmol) = $stm->fetchrow_array()) {
     push @names, $newmol;
   }
 
  }
  my ($label, $moltype) = $lv->LabelValueToString($pw->MolType($molid));
  if ($moltype eq "complex") {
    # unshift @names, MoleculeName($pw, $molid);
  }
  my $ehash = $pw->MolExId($molid);
  for my $id_type (keys %{ $ehash }) {
    for my $id (keys %{ $$ehash{$id_type} }) {
      if ($id_type eq "LL") {
        $locus = $id;
      } elsif ($id_type eq "UP") {
        $uniprot = $id;
      } elsif ($id_type eq "CA") {
        $ca = $id; 
      } elsif ($id_type eq "GO") {
        $go = $id; 
      } elsif ($id_type eq "SG") {
        $sg = $id;
      } else {
        push @ids, $id;
      }
    }
  }

  if ($namepush == 0) {
    $sql = "select ext_mol_id, id_type from $schema.pw_ext_mol_id where mol_id = $molid";
    $stm = $db->prepare($sql);
    if(not $stm) {
      print STDERR "$sql\n";
      print STDERR "$DBI::errstr\n";
      print STDERR "prepare call failed\n";
      die;
    }
    if(!$stm->execute()) {
      print STDERR "$sql\n";
      print STDERR "$DBI::errstr\n";
      print STDERR "execute call failed\n";
      die;
    }

   while ((my $newmol, my $idtype) = $stm->fetchrow_array()) {
     if ($idtype eq "LL") {
        $locus = $newmol;
      } elsif ($idtype eq "UP") {
        $uniprot = $newmol;
      } elsif ($idtype eq "CA") {
        $ca = $newmol;
      } elsif ($idtype eq "GO") {
        $go = $newmol;
      } elsif ($idtype eq "SG") {
        $sg = $newmol;
      } else {
        push @ids, $newmol;
      }
   }

  }


  $sql = "select ll_id from $schema.pw_ext_mol_id a, cgap.ll2sp b";
  $sql = $sql . " where mol_id = $molid and b.sp_primary = substr(a.ext_mol_id,0,6)";
  $sql = $sql . " and organism = 'Hs' ";

  $stm = $db->prepare($sql);
  if(not $stm) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "prepare call failed\n";
    die;
  }
  if(!$stm->execute()) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "execute call failed\n";
    die;
  }
  
   while (($tmplocus) = $stm->fetchrow_array()) {
     $locus = $tmplocus; 
   }

   if ($locus ne '') {
 # Get Official Symbol

  $sql = qq!
  select distinct b.symbol, b.name_type
  from $schema.pw_molecule_search b
  where ll_id = $locus 
  and name_type in ('OF','AS')
  !;

  $stm = $db->prepare($sql);
  if(not $stm) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "prepare call failed\n";
    die;
  }
  if(!$stm->execute()) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "prepare call failed\n";
    die;
  }

  while (($symbol, $name_type) =
    $stm->fetchrow_array()) {
      if ($name_type eq "OF") {
        $officialCount++;
        $official = $official . $symbol . " ";
      } else {
        $aliasCount++;
        $alias = ", " . $symbol . $alias ;
      }
  }

  $alias = substr($alias,2);
  } 
  my $is_family = 0;
  my $family_kids = $pw->FamilyChildren($molid);
  if (@{ $family_kids }) {
    $is_family = 1;
    $moltype = "$moltype family";
  }

  my $is_in_family = 0;
  my $family_parents = $pw->FamilyParent($molid);
  if (@{ $family_parents }) {
    $is_in_family = 1;
  }

  my $has_subunits = 0;
  my $subunits = $pw->MolPart($molid);
  if (@{ $subunits }) {
    $has_subunits = 1;
  }

  my $is_subunit_of = 0;
  my $wholes = $pw->MolWhole($molid);
  if (@{ $wholes }) {
    $is_subunit_of = 1;
  }

  push @{ $lines }, "<p><b>Identification:</b>";
  push @{ $lines }, "<blockquote>";
  push @{ $lines }, "<p><table border=1 cellspacing=1 cellpadding=4>";
  push @{ $lines }, "<tr><td><b>Names:</b></td><td>" . join(", ", sort @names) . "</td></tr>";
  if ($officialCount > 0) { 
    push @{ $lines }, "<tr><td><b>HUGO symbol:</b></td><td>$official</td></tr
>";
  }
  if ($aliasCount > 0) {
    push @{ $lines }, "<tr><td><b>Aliases:</b></td><td>$alias</td></tr>";
  }
  push @{ $lines }, "<tr><td><b>Molecule type:</b></td><td>$moltype" .
      ($moltype ne "compound" ? " ($org)" : "") . "</td></tr>";
  
  if ($locus) {
    push @{ $lines }, "<tr><td><b>Entrez Gene:</b></td><td><a href=\"" . LL_URL($locus) . "\"" .
        ">$locus</a></td></tr>";   
    my $cid = GetCIDForLocus($db, $org, $locus);
    if ($cid) {
      push @{ $lines }, "<tr><td><b>CGAP:</b></td><td><a href=\"" .
          CGAP_URL($org, $cid) . "\"" .
          ">$org.$cid</a></td></tr>";   
    }
  }
  if ($uniprot) {
    push @{ $lines }, "<tr><td><b>UniProt:</b></td><td><a href=\"" .
        UNIPROT_URL($uniprot) . "\"" .
        ">$uniprot</a></td></tr>";   
  }

  if ($ca) {
    push @{ $lines }, "<tr><td><b>CAS ID:</b></td><td>$ca <a href=http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=pccompound&term=$ca" . ">View relevant record in PubChem</a></td></tr>";
  }

  if ($go) {
    push @{ $lines }, "<tr><td><b>GO ID:</b></td><td><a href=http://www.godatabase.org/cgi-bin/amigo/go.cgi?view=details\&query=GO:" . "$go" . ">GO:$go</a></td></tr>";
  }

  if ($sg) {
    push @{ $lines }, "<tr><td><b>Signaling Gateway Id:</b></td><td><a href=http://www.signaling-gateway.org/molecule/query?afcsid=$sg&type=orthologs&adv=latest>$sg</a></td></tr>";
  } 
  if (@ids) {
    push @{ $lines }, "<tr><td><b>Other ids:</b></td><td>" . join(", ", sort @ids) . "</td></tr>";
  }

  push @{ $lines }, "<tr><td><b>Molecule Id:</b></td><td>$molid</td></tr>";
  
  if ($is_family) {
    push @{ $lines }, "<tr><td><b>Family members:</b></td><td>";
    for my $kid (@{ $family_kids }) {
      push @{ $lines }, 
        "<a href=" . MOLPAGE_URL($kid) . " >" .
        $pw->PickMolName($kid) . "</a><br>";
    }
    push @{ $lines }, "</td></tr>";
  }
  if ($is_in_family) {
    push @{ $lines }, "<tr><td><b>Is member of family:</b></td><td>";
    for my $p (@{ $family_parents }) {
      push @{ $lines }, 
        "<a href=" . MOLPAGE_URL($p) . " >" .
        $pw->PickMolName($p) . "</a><br>";
    }
    push @{ $lines }, "</td></tr>";
  }
  if ($has_subunits) {
    push @{ $lines }, "<tr><td><b>Has subunits:</b></td><td>";
    for my $s (@{ $subunits }) {
      my $bounds;
      my ($start, $end) = split(",", $pw->PartBounds($s));
      if ($start) {
        $bounds = " ($start" . "-" . "$end aa)";
      }
      push @{ $lines }, 
        "<a href=" . MOLPAGE_URL($s) . " >" .
        $pw->PickMolName($s) . "</a><br>";
      # push @{ $lines }, "</td></tr>";

    }
     push @{ $lines }, "</td></tr>";


  }
  if ($is_subunit_of) {
    push @{ $lines }, "<tr><td><b>Is subunit of:</b></td><td>";
    my $bounds;
    my ($start, $end) = split(",", $pw->PartBounds($molid));
    if ($start) {
      $bounds = " ($start" . "-" . "$end aa)";
    }
    for my $w (@{ $wholes }) {
      push @{ $lines }, 
        "<a href=" . MOLPAGE_URL($w) . " >" .
        $pw->PickMolName($w) . "</a><br>$bounds<br>";
    }
    push @{ $lines }, "</td></tr>";
    push @{ $lines }, "<tr><td><b>Cleaved coordinates:  </b></td><td>$bounds</td></tr>";

  }
  push @{ $lines }, "</table>";
  push @{ $lines }, "</blockquote>";
}

######################################################################
sub GetCIDForLocus {
  my ($db, $org, $locus) = @_;

##!!!!!!
my $schema = "cgap";

  my ($sql, $stm);
  my ($cluster_table, $cluster_number, @temp);

  if ($org eq "Hs") {
    $cluster_table = "hs_cluster";
  } elsif ($org eq "Mm") {
    $cluster_table = "mm_cluster";
  } else {
    die "Unrecognized organism $org";
  }

  $sql = qq!
select
  cluster_number
from
  $schema.$cluster_table
where
  locuslink = $locus
  !;

  $stm = $db->prepare($sql);
  if(not $stm) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "prepare call failed\n";
    die;
  }
  if(!$stm->execute()) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "execute call failed\n";
    die;
  }
  while (($cluster_number) = $stm->fetchrow_array()) {
    push @temp, $cluster_number;
  }
  return (join",", @temp);
}


######################################################################
sub ComplexInfo {
  my ($db, $pw, $molid, $lines) = @_;

  my ($molname, $compid, $c, $state, $location, $lvid, $label, $value);
  my (@ptms, $ptms);
  
  push @{ $lines }, "<p><b>Complex components:</b>";
  push @{ $lines }, "<blockquote>";
  push @{ $lines }, "<p><table border=1 cellspacing=1 cellpadding=4>";
  push @{ $lines }, "<tr><td><b>Molecule<b></td>" .
      "<td><b>Location<b></td>" .
      "<td><b>State<b></td>" .
      "<td><b>Modifications<b></td></tr>";

  my $componentCount = 0;

  for $c (@{ $pw->Components($molid) }) {
    $componentCount++;
    $compid = $pw->ComponentMol($molid, $c);
    undef $state;
    undef $location;
    undef @ptms;
    for my $x ( @{ $pw->ComponentPTM($molid, $c) }) {
      my ($protein_id, $position, $amino_acid, $modification_label_value,
          $modification_label_name) = @{ $x };
      push @ptms,
          "$protein_id \[$amino_acid $position $modification_label_name\]";
    }
    for $lvid (@{ $pw->ComponentLabel($molid, $c) }) {
      ($label, $value) = $lv->LabelValueToString($lvid);
      if ($label eq "activity-state") {
        $state = $value;
      } elsif ($label eq "location") {
        $location = $value;
      }
    }
    if (@ptms > 0) {
      $ptms = join("<br>", @ptms);
    } else {
      $ptms = "&nbsp;";
    }
    $state or $state = "&nbsp;";
    $location or $location = "&nbsp;";
    push @{ $lines }, "<tr><td><a href=" . MOLPAGE_URL($compid) .
        " >" . MoleculeName($pw, $compid) . "</a></td>" .
        "<td>$location</td><td>$state</td><td>$ptms</td></tr>";
  }
  if ($componentCount == 0) { 
   my $sql = qq!
      select mol_id_1, mol_name
      from pid.pw_mol_mol, pid.pw_mol_name
      where mol_id_2 = $molid 
      and mol_id_1 = mol_id
      and relation = 'c'
    !;

    my $stm = $db->prepare($sql);
    if(not $stm) {
      print STDERR "$sql\n";
      print STDERR "$DBI::errstr\n";
      print STDERR "prepare call failed\n";
      die;
    }

    if(!$stm->execute()) {
      print STDERR "$sql\n";
      print STDERR "$DBI::errstr\n";
      print STDERR "execute call failed\n";
      die;
    }

    while ((my $component, my $molname) = $stm->fetchrow_array()) {
    if (@ptms > 0) {
      $ptms = join("<br>", @ptms);
    } else {
      $ptms = "&nbsp;";
    }
    $state or $state = "&nbsp;";
    $location or $location = "&nbsp;";

      push @{ $lines }, "<tr><td><a href=" . MOLPAGE_URL($component) .
        " >" . $molname . "</a></td>" .
        "<td>$location</td><td>$state</td><td>$ptms</td></tr>";
 
    }
    $stm->finish();
  } 
  push @{ $lines }, "</table>";
  push @{ $lines }, "</blockquote>";

  my ($cp_id, $cx_id, @temp);
  my $standard = join(":", sort keys %{ $cx2cp{$molid} });
  for $cx_id (keys %cx2cp) {
    if ($cx_id != $molid &&
        $standard eq join(":", sort keys %{ $cx2cp{$cx_id} })) {
      push @temp, "<tr><td><a href=" . MOLPAGE_URL($cx_id) .
          " >$cx_id</td></tr>";
    }
  }
  if (@temp) {
    push @{ $lines }, "<p><b>Other complexes with same components:</b>";
    push @{ $lines }, "<blockquote>";
    push @{ $lines }, "<p><table border=1 cellspacing=1 cellpadding=4>";
    push @{ $lines }, @temp;
    push @{ $lines }, "</table>";
    push @{ $lines }, "</blockquote>";
  }
}

######################################################################
sub ComplexTable {
  my ($pw, $molid, $lines) = @_;

  my ($label, $value);
  my ($cx, $m, @temp);

  if (scalar(keys %complex_uses) == 0) {
    return;
  }

  push @{ $lines }, "<p><b>Uses in complexes:</b>";
  push @{ $lines }, "<p>Duplicate names indicate differences in component properties</p>";
  push @{ $lines }, "<blockquote>";
  push @{ $lines }, "<p><table border=1 cellspacing=1 cellpadding=4>";
  push @{ $lines }, "<tr><td><b>Complex</b></td><td><b>Other components</b></td></tr>";
  for $cx (keys %complex_uses) {
    undef @temp;
    push @{ $lines }, "<tr><td><a href=\"" .
        MOLPAGE_URL($cx) . "\" >" .
        MoleculeName($pw, $cx) . "</a></td><td>";
    for my $m (keys %{ $complex_uses{$cx} }) {
      if ($m != $molid || keys %{ $complex_uses{$cx} } == 1) {
        ## i.e., allow for case where complex is just a dimer
        ## (so there is only one type of component)
        push @temp, ("<a href=\"" . MOLPAGE_URL($m) . "\" >" .
         MoleculeName($pw, $m) . "</a>");
      }
    }
    push @{ $lines }, join(", ", @temp);
    push @{ $lines }, "</td></tr>";
  }
  push @{ $lines }, "</table>";
  push @{ $lines }, "</blockquote>";

}

######################################################################
sub AtomTable {
  my ($pw, $molid, $lines) = @_;

  my ($label, $value);
  my ($atom, $edge);
  my (@pids, $pathwayid, $pcell, $atomtype);
  my ($location, $state, $role, $source);
  my %tmp;

  push @{ $lines }, "<p><b>Uses in interactions</b>";
  push @{ $lines }, "<blockquote>";
  push @{ $lines }, "<p><b>Display all interactions:</b> <a href=\"" .
      ATOM_URL(join(",", keys %atom_uses), "jpg") .
      "\" >JPG</a>&nbsp;&nbsp;<a href=\"" .
      ATOM_URL(join(",", keys %atom_uses), "svg") .
      "\" >SVG</a>&nbsp;&nbsp;<a href=\"" .
      ATOM_TEXT_URL(join(",", keys %atom_uses), "xml") .
      "\" >XML</a>&nbsp;&nbsp;<a href=\"" .
      ATOM_TEXT_URL(join(",", keys %atom_uses), "biopax") .
      "\" >BioPAX</a>";

  push @{ $lines }, "<p><table border=1 cellspacing=1 cellpadding=4>";
  push @{ $lines }, "<tr>" .
      "<td><b>Role</b></td>" .
      "<td><b>Location</b></td>" .
      "<td><b>State</b></td>" .
      "<td><b>Post-translational modifications</b></td>" .
      "<td><b>Interaction</b></td>" .
      "<td><b>Pathway</b></td>" .
      "<td><b>Source</b></td>" .
      "</tr>";
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
          $source = $pathway2source{$pathwayid};
          $pcell = $pathway2name{$pathwayid} .
              "&nbsp;&nbsp;<a href=\"" . PATHWAY_URL($pathwayid, "svg",$source) .
              "\" >SVG</a>&nbsp;&nbsp;" .
              "<a href=\"" . PATHWAY_URL($pathwayid, "jpg",$source) .
              "\" >JPG</a>";
          # $source = $pathway2source{$pathwayid};
        } else {
          $pcell = "&nbsp;";
          $source = $atom2source{$atom};
        }
        $role              = $atom_uses{$atom}{$edge};
        $role =~ s/agent/positive regulator/;
        $role =~ s/inhibitor/negative regulator/;
        $location          = $location{$atom}{$edge};
        $location or $location = "&nbsp;";
        $state             = $state{$atom}{$edge};
        $state or $state   = "&nbsp;";
        $source =~ s/NATURE/NCI-Nature curated/;
        $source =~ s/BioCarta/BioCarta Imported/;
        $source or $source = "&nbsp;";

        my (@ptms, $ptms);
        for my $x ( @{ $pw->EdgePTM($atom, $edge) }) {
          my ($protein_id, $position, $amino_acid,
              $modification_label_value, $modification_label_name) = @{ $x };
          push @ptms,
              "$protein_id \[$amino_acid $position $modification_label_name\]";
        }
        if (@ptms > 0) {
          $ptms = join("<br>", @ptms);
        } else {
          $ptms = "&nbsp;";
        }

        push @{ $tmp{"$role|$location|$state|"} }, "<tr><td>" .
          join("</td>\n<td>",
            $role,
            $location,
            $state,
            $ptms,
            "$atomtype (<a href=\"" . ATOMPAGE_URL($atom) .
                "\" >$atom</a>)",
            $pcell,
            $source
          ) . "</td></tr>";
      }
    }
  }
  for my $i (sort keys %tmp) {
    for (@{ $tmp{$i} }) {
      push @{ $lines }, $_;
    }
  }
  push @{ $lines }, "</table>";
  push @{ $lines }, "</blockquote>";
}

######################################################################
sub SimpleOrComplex {
  my ($pw, $molid) = @_;

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
        LeafComponents($pw, $m, \%accum);
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
  my ($pw, $molid, $accum) = @_;

  my ($seq, $m);

  if ($pw->MolType($molid) != $COMPLEX_TYPE) {
    return;
  }
  for $seq (@{ $pw->Components($molid) }) {
    $m = $pw->ComponentMol($molid, $seq);
    if (! defined $$accum{$m}) {
      $$accum{$m} = 1;
      if ($pw->MolType($m) == $COMPLEX_TYPE) {
        LeafComponents($pw, $m, $accum);
      }
    }
  }
}

######################################################################
# END MoleculePage
######################################################################
######################################################################
sub round {
  my ($x) = @_;

  if ($x >= 0) {
    return int($x + 0.5);
  } else {
    return int($x - 0.5);
  }
}

######################################################################
sub PlaceArrow {
  my ($stroke, $x0, $y0, $x1, $y1) = @_;

  ##
  ## constant dimensions for arrowhead
  ## q0 is midpoint on the short size (side connected to the edge)
  ## q1 is the tip of the arrowhead
  ##

  my ($q0x, $q0y) = (0,  0);
  my ($q1x, $q1y) = (10, 0);
  my ($q2x, $q2y) = (0,  3.5);
  my ($q3x, $q3y) = (0, -3.5);
  my ($q4x, $q4y) = (1,  3.5);
  my ($q5x, $q5y) = (1, -3.5);

  my $dx = $x1 - $x0;
  my $dy = $y1 - $y0;

  if (($dx < 0.01) && ($dx >= 0)) {
    $dx = 0.01;
  }

  my $theta = atan($dy/$dx);

 
   if ($theta > 0) {
    $q1x = -10;
  }

  if (($theta == 0) && ($dx < 0)) {
    $q1x = -10;
  }

  if ($dy > 0) {
    $q1x = 10;
  } 
 
  if (($dy > 0) && ($dx < 0)) {
    $q1x = -10;
  }

  # print STDERR "DX:  $dx DY:  $dy Theta:  $theta\n"; 
 
  ##
  ## matrix for clockwise rotation
  ##

  my $B = [ [ cos($theta), -sin($theta) ], [ sin($theta), cos($theta) ] ];

  my ($p1x, $p1y) = simple_mult($q1x, $q1y, $B);
  my ($p2x, $p2y) = simple_mult($q2x, $q2y, $B);
  my ($p3x, $p3y) = simple_mult($q3x, $q3y, $B);
  my ($p4x, $p4y) = simple_mult($q4x, $q4y, $B);
  my ($p5x, $p5y) = simple_mult($q5x, $q5y, $B);

  $p1x = sprintf("%.2f", $p1x + $x0);
  $p1y = sprintf("%.2f", $p1y + $y0);
  $p2x = sprintf("%.2f", $p2x + $x0);
  $p2y = sprintf("%.2f", $p2y + $y0);
  $p3x = sprintf("%.2f", $p3x + $x0);
  $p3y = sprintf("%.2f", $p3y + $y0);
  $p4x = sprintf("%.2f", $p4x + $x0);
  $p4y = sprintf("%.2f", $p4y + $y0);
  $p5x = sprintf("%.2f", $p5x + $x0);
  $p5y = sprintf("%.2f", $p5y + $y0);



  if ($stroke eq "ff0000") {
    # print STDERR "$p4x,$p4y ", "$p2x,$p2y ", "$p3x,$p3y", "$p4x,$p3y", "$p1x,$p2y\n";
    return [ "$p4x,$p4y", "$p2x,$p2y", "$p3x,$p3y", "$p5x,$p5y", "$p4x,$p4y" ];
  }
  return [ "$p1x,$p1y", "$p2x,$p2y", "$p3x,$p3y", "$p1x,$p1y" ];
}

######################################################################
sub simple_mult {
  my ($x, $y, $B) = @_;

  ## multiply [x y] by 2x2 rotation matrix

  my $ab0 = $x * $$B[0][0] + $y * $$B[0][1];
  my $ab1 = $x * $$B[1][0] + $y * $$B[1][1];

  return ($ab0, $ab1);
}

######################################################################
sub Decls {
  my ($root) = @_;
  return $root->{decls};
}

######################################################################
sub Color2Hex {
  my ($color_word) = @_;
  if (defined $hex_color{$color_word}) {
    return $hex_color{$color_word};
  } else {
    # print STDERR "no hex for color $color_word\n";
    return $hex_color{black};
  }
}

######################################################################
sub FindBB {
  my ($root) = @_;

  for my $decl (@{ Decls($root) }) {
    my $node_edge = $decl->{nodeoredge};
    if ($node_edge->{kind} eq "node" && $node_edge->{node} eq "graph") {
      for my $attr (@{ $decl->{attrs} }) {
        if ($attr->{label} eq "bb") {
          return split(",", $attr->{value});
        }
      }
    }
  }
}

######################################################################
sub GetDeclType {
  my ($decl) = @_;
  return $decl->{nodeoredge}->{kind};
}

######################################################################
sub GetAttrs {
  my ($decl) = @_;

  my %attrs;

  my $node_edge = $decl->{nodeoredge};
  for my $attr (@{ $decl->{attrs} }) {
    $attrs{$attr->{label}} = $attr->{value};
  }
  return \%attrs;
}

######################################################################
sub GetId {
  my ($decl) = @_;

  my $decl_type = GetDeclType($decl);
  if ($decl_type eq "node") {
    return $decl->{nodeoredge}{node};
  } else {
    return $decl->{nodeoredge}{from} . " &#45;&gt; " . $decl->{nodeoredge}{to} ;
  }
}

######################################################################
sub PlainText {
  my ($fh, $meta_f, $center_x, $center_y, $id, $width, $height, $attrs) = @_;

  my $text      = $attrs->{label};
  my $stroke    = Color2Hex($attrs->{color});
  my $delta     = 5;
  if (defined $text) {
    PrintTextBlock($fh, $meta_f, 10, $text, $attrs, $id, $center_x, $center_y, $delta);
  }
}

######################################################################
sub Box {
  my ($fh, $meta_f, $center_x, $center_y, $id, $width, $height, $attrs) = @_;

  my $text      = $attrs->{label};
  my $stroke    = Color2Hex($attrs->{color});

  my $x0 = round($center_x - $width/2);
  my $x1 = round($center_x + $width/2);
  my $y0 = round($center_y - $height/2);
  my $y1 = round($center_y + $height/2);

  my $points = join(" ",
    $x1 . "," . Invert($y1),
    $x0 . "," . Invert($y1),
    $x0 . "," . Invert($y0),
    $x1 . "," . Invert($y0),
    $x1 . "," . Invert($y1)
  );

  $stroke = "#" . $stroke;
  
  print $fh qq!
  <Polygon Points="$points" FillRule="NonZero" Fill="$stroke" Stroke="$stroke"/>
!;
  print STDERR "In box\n";
  if (defined $text) {
    PrintTextBlock($fh, $meta_f, 12, $text, $attrs, $id, $center_x, $center_y, 1, "box", $x0 );
  }
  print STDERR "Leaving box\n"; 
}

###################################################################
sub Circle {
  my ($fh, $meta_f, $x, $y, $id, $attrs) = @_;

  my $text      = $attrs->{label};
  my $stroke    = Color2Hex($attrs->{color});
   my @id_type = split(/[_]+/,$id);
  my $atom_id;
  if ($id_type[0] eq 'A') {
    $atom_id = $id_type[1];
  }
 
  $x = $x - 7;
  $y = $y + 7;

  $y = Invert($y);

  print $fh qq!
  <HyperlinkButton  NavigateUri="http://pid.nci.nih.gov/search/InteractionPage?atomid=$atom_id" >
      <HyperlinkButton.Content>
      <Canvas Tag="$keyCount" Name=\"atom$atom_id\">
      !;

  print $fh qq!
  <Ellipse Canvas.Left="$x" Width="14" Canvas.Top="$y" Height="14" Fill="#000000" Stroke="#000000" />
!;
  $x = $x - 4;
  $y = $y - 4;

  print $fh qq!
  <Ellipse Canvas.Left="$x" Width="22" Canvas.Top="$y" Height="22" Stroke="#000000" />
!;

    print $fh qq!
    </Canvas>
    </HyperlinkButton.Content>
    </HyperlinkButton>
   !;

   if (defined $text) {
    PrintTextBlock($fh, $meta_f, 10, $text, $attrs, $id, $x, $y);
  } else {
    $keyCount++;
  }
 
}

######################################################################
sub Diamond {
  my ($fh, $meta_f, $center_x, $center_y, $id, $width, $height, $attrs) = @_;

  my $text      = $attrs->{label};
  my $stroke    = Color2Hex($attrs->{color});
  
  my @id_type = split(/[_]+/,$id);
  my $atom_id; 
  if ($id_type[0] eq 'A') {
    $atom_id = $id_type[1];
  }

  my $x0 = round($center_x - $width/2);
  my $x1 = round($center_x + $width/2);
  my $y0 = round($center_y - $height/2);
  my $y1 = round($center_y + $height/2);

  my $points = join(" ",
    $center_x . "," . Invert($y1),
    $x0 . "," . Invert($center_y),
    $center_x . "," . Invert($y0),
    $x1 . "," . Invert($center_y),
    $center_x . "," . Invert($y1)
  );
  print $fh qq!
  <HyperlinkButton  NavigateUri="http://pid.nci.nih.gov/search/InteractionPage?atomid=$atom_id" >
      <HyperlinkButton.Content>
      <Canvas Tag="$keyCount" Name=\"atom$atom_id\"> 
      !;

  print $fh qq!
  <Polygon Points="$points" FillRule="NonZero" Fill="#000000" Stroke="#000000"/>
!;
 
   print $fh qq!
    </Canvas> 
    </HyperlinkButton.Content>
    </HyperlinkButton>
   !;
  
   if (defined $text) {
    PrintTextBlock($fh, $meta_f, 10, $text, $attrs, $id, $center_x, $center_y);
  } else {
    $keyCount++;
  }
}

######################################################################
sub Triangle {
  my ($fh, $meta_f, $center_x, $center_y, $id, $width, $height, $attrs) = @_;

  my $text      = $attrs->{label};
  my $stroke    = Color2Hex($attrs->{color});
  my @id_type = split(/[_]+/,$id);
  my $atom_id;
  if ($id_type[0] eq 'A') {
    $atom_id = $id_type[1];
  }

  my $x0 = round($center_x - $width/2);
  my $x1 = round($center_x + $width/2);
  my $y0 = round($center_y - $height/2);
  my $y1 = round($center_y + $height/2);

  my $points = join(" ",
    $center_x . "," . Invert($y1),
    $x0 . "," . Invert($y0),
    $x1 . "," . Invert($y0),
    $center_x . "," . Invert($y1)
  );

print $fh qq!
  <HyperlinkButton  NavigateUri="http://pid.nci.nih.gov/search/InteractionPage?atomid=$atom_id" >
      <HyperlinkButton.Content>
      <Canvas Tag="$keyCount" Name=\"atom$atom_id\">
      !;

  print $fh qq!
  <Polygon Points="$points" FillRule="NonZero" Fill="#000000" Stroke="#000000"/>
!;

  print $fh qq!
    </Canvas>
    </HyperlinkButton.Content>
    </HyperlinkButton>
    !;

  if (defined $text) {
    PrintTextBlock($fh, $meta_f, 10, $text, $attrs, $id, $center_x, $center_y);
  } else {
    $keyCount++;
  }
}

######################################################################
sub Parallelogram {
  my ($fh, $meta_f, $center_x, $center_y, $id, $width, $height, $attrs) = @_;

  my $text      = $attrs->{label};
  my $stroke    = Color2Hex($attrs->{color});
  my $delta     = 5;
  my $offset = $height/$tangent;

  my $x0 = $center_x - $width/2;
  my $x1 = $center_x + $width/2;
  my $y0 = $center_y - $height/2;
  my $y1 = $center_y + $height/2;

  my $points = join(" ",
    $x1 . "," . Invert($y1),
    $x0 + $offset . "," . Invert($y1),
    $x0 . "," . Invert($y0),
    $x1 - $offset . "," . Invert($y0),
    $x1 . "," . Invert($y1)
  );

  print $fh qq!
  <Polygon Tag="$keyCount" Points="$points" Stroke="#000000"/>
!;

  if (defined $text) {
    PrintTextBlock($fh, $meta_f, 10, $text, $attrs, $id, $center_x, $center_y, $delta, 'parallelogram');
  }
}

######################################################################
sub Mrecord {
  my ($fh, $meta_f, $center_x, $center_y, $id, $width, $height, $attrs) = @_;
  
  my $delta = 7;
  my $text      = $attrs->{label};
  my $stroke    = Color2Hex($attrs->{color});

  my $x0 = round($center_x - $width/2);
  my $x1 = round($center_x + $width/2);
  my $y0 = round($center_y - $height/2);
  my $y1 = round($center_y + $height/2);

  ## South

  my $sswx = $x0 + $ROUND_CORNER_CLIP; my $sswy = Invert($y0);
  my $ssex = $x1 - $ROUND_CORNER_CLIP; my $ssey = Invert($y0);

  ## East

  my $esex = $x1; my $esey = Invert($y0 + $ROUND_CORNER_CLIP);
  my $enex = $x1; my $eney = Invert($y1 - $ROUND_CORNER_CLIP);

  ## North

  my $nnex = $x1 - $ROUND_CORNER_CLIP; my $nney = Invert($y1);
  my $nnwx = $x0 + $ROUND_CORNER_CLIP; my $nnwy = Invert($y1);

  ## West

  my $wnwx = $x0; my $wnwy = Invert($y1 - $ROUND_CORNER_CLIP);
  my $wswx = $x0; my $wswy = Invert($y0 + $ROUND_CORNER_CLIP);

  ## SE corner

  my $se_c1x = $ssex + $HALF_CORNER_CLIP; my $se_c1y = $ssey;
  my $se_c2x = $esex;                     my $se_c2y = $esey + $HALF_CORNER_CLIP;

  ## NE corner

  my $ne_c1x = $enex;                     my $ne_c1y = $eney - $HALF_CORNER_CLIP;
  my $ne_c2x = $nnex + $HALF_CORNER_CLIP; my $ne_c2y = $nney;

  ## NW corner

  my $nw_c1x = $nnwx - $HALF_CORNER_CLIP; my $nw_c1y = $nnwy;
  my $nw_c2x = $wnwx;                     my $nw_c2y = $wnwy - $HALF_CORNER_CLIP;

  ## SW corner

  my $sw_c1x = $wswx;                     my $sw_c1y = $wswy + $HALF_CORNER_CLIP;
  my $sw_c2x = $sswx - $HALF_CORNER_CLIP; my $sw_c2y = $sswy;

  print $fh qq!
  <Polyline Points="$sswx,$sswy $ssex,$ssey" FillRule="NonZero" Stroke="#$stroke"/>
  <Path Stroke="#$stroke" Data="M$ssex $ssey\C$se_c1x $se_c1y $se_c2x $se_c2y $esex $esey"/>
  <Polyline Points="$esex,$esey $enex,$eney" FillRule="NonZero" Stroke="#$stroke"/>
  <Path Stroke="#$stroke" Data="M$enex $eney\C$ne_c1x $ne_c1y $ne_c2x $ne_c2y $nnex $nney"/>
  <Polyline Points="$nnex,$nney $nnwx,$nnwy" FillRule="NonZero" Stroke="#$stroke"/>
  <Path Stroke="#$stroke" Data="M$nnwx $nnwy\C$nw_c1x $nw_c1y $nw_c2x $nw_c2y $wnwx $wnwy"/>
  <Polyline Points="$wnwx,$wnwy $wswx,$wswy" FillRule="NonZero" Stroke="#$stroke"/>
  <Path Stroke="#$stroke" Data="M$wswx $wswy\C$sw_c1x $sw_c1y $sw_c2x $sw_c2y $sswx $sswy"/> 
!;

  if (defined $text) {
    PrintTextBlock($fh, $meta_f, 10, $text, $attrs, $id, $center_x, $center_y, $delta, "mrecord", $x0, $x1, $y0, $y1);
  }
}

######################################################################
sub Invert {
  my ($y) = @_;
  return -1*$y;
}

######################################################################
sub Invert_Y {
  my ($xy) = @_;
  my ($x, $y) = split(",", $xy);
  return ("$x," . -1*$y);
}

######################################################################
sub PrintTextBlock {
  my ($fh, $meta_f, $font_size, $text, $attrs, $id, $x, $y, $delta, $type, $x1, $x2, $y1, $y2) = @_;

  my $db = DBI->connect("DBI:Oracle:" . 'cgprod', 'web', 'readonly');
  if (not $db or $db->err()) {
    print STDERR "Cannot connect to " . 'pid' . "@" . 'cgprod' . "\n";
    die;
  }

  my $font_family = "Arial";
  my $sevenadjust = 5;
  my $label = "";
  my $prevlabel = "";
  my $shape = "";
  my $url = "";
  my $atom_id = "";
  my $mol_id = "";
  my $molType = "";
  my $state = "";
  my $ptms = "";
  my $macroprocess = "";
  my $entityType = "none";
  my $location = "";
  my $uniprot = "";
  my $entrez = "";
  my $molType = "";
  my $intType = "";
  my $nolink = 1;
  my $ltext = length($text);
  my $tmove = 0;
  my $entityVal = '';

  for (my $count = 0; $count < $ltext; $count++) {
    if ((substr($text,$count,1) eq '<') || (substr($text,$count,1) eq '>') || (substr($text,$count,1) eq '{') || (substr($text,$count,1) eq '}')) {
      $x++; 
    }
  }
  
  $text =~ s/{//g;
  $text =~ s/}//g;
  $text =~ s/\\<//g;
  $text =~ s/\\>//g;
 
  $text =~ s/[\<]+//g;
  $text =~ s/[\>]+//g;

  my @entitysearch2 = split(/[\[]/,$text);
  my $entityVal2 = $entitysearch2[0];

  my @locationCheck = split(/[\[\]]/,$text);
  
  if ($#locationCheck > 1) {
    $location = $locationCheck[1];
    # print STDERR "Location: $location\n"; 
    if ($location eq 'n') {
      $location = 'nucleus';
    }
    if ($location eq 'cy') {
      $location = 'cytoplasm';
    }
    if ($location eq 'm') {
      $location = 'transmembrane';
    }
  }
  my $statelen = length($entityVal2);
  $state = '';
  
  if (substr($entityVal2,($statelen-1),1) eq '+') {
    $state = 'active';
  }
  if (substr($entityVal2,($statelen-2),2) eq '+1') {
    $state = 'active 1';
  }
  if (substr($entityVal2,($statelen-2),2) eq '+2') {
    $state = 'active 2';
  }
  if (substr($entityVal2,($statelen-2),2) eq '+3') {
    $state = 'active 3';
  }

  if (substr($entityVal2,($statelen-1),1) eq '-') {
    $state = 'inactive';
  }

  if (substr($entityVal2,($statelen-2),2) eq '-1') {
    $state = 'inactive 1';
  }
  if (substr($entityVal2,($statelen-2),2) eq '-2') {
    $state = 'inactive 2';
  }
  if (substr($entityVal2,($statelen-2),2) eq '-3') {
    $state = 'inactive 3';
  }

  if (defined $attrs->{label}) {
    $label = $attrs->{label};
  } 

  if (defined $attrs->{shape}) {
    $shape = $attrs->{shape};
  }

  if (defined $attrs->{URL}) {
    $url = $attrs->{URL};
    $nolink=0;
  } 

  # print STDERR "Id:  $id\n";
  my @id_type = split(/[_]+/,$id);
  # print STDERR "ID:  $id_type[1]\n";
  if ($id_type[0] eq 'A') {
    $atom_id = $id_type[1];
    $entityType = "interaction";
  }
  
  my @ptmtemp; 
  if ($id_type[0] eq 'M') {
    $mol_id = $id_type[1];
    $entityType = "molecule";
    @ptmtemp = split(/__/,$id);
    $ptms = $ptmtemp[2];
  
    my $sql = qq!
      select ext_mol_id
      from pid.pw_ext_mol_id
      where id_type in ('UP')
      and mol_id = $mol_id
    !;

    my $stm = $db->prepare($sql);
    if(not $stm) {
      print STDERR "$sql\n";
      print STDERR "$DBI::errstr\n";
      print STDERR "prepare call failed\n";
      die;
    }
 
    if(!$stm->execute()) {
      print STDERR "$sql\n";
      print STDERR "$DBI::errstr\n";
      print STDERR "execute call failed\n";
      die;
    }
  
    while ((my $tempuniprot) = $stm->fetchrow_array()) {
        $uniprot = $tempuniprot; 
    }
    $stm->finish();
  
    $sql = qq!
      select
        g.ll_id
      from
        cgap.ll_gene g,
        pid.pw_ext_mol_id e,
        pid.pw_mol m,
        cgap.ll2sp s,
        cgap.sp_primary p

      where
        m.mol_id = $mol_id 
        and p.sp_id_or_secondary = e.ext_mol_id
        and p.sp_id_or_secondary = s.sp_primary
        and s.organism = 'Hs'
        and s.ll_id = g.ll_id
        and e.mol_id = m.mol_id
        and m.basic_mol_type = 'PR'

    !;

    my $stm = $db->prepare($sql);
    if(not $stm) {
      print STDERR "$sql\n";
      print STDERR "$DBI::errstr\n";
      print STDERR "prepare call failed\n";
      die;
    }

    if(!$stm->execute()) {
      print STDERR "$sql\n";
      print STDERR "$DBI::errstr\n";
      print STDERR "execute call failed\n";
      die;
    }

    while ((my $tempentrez) = $stm->fetchrow_array()) {
        $entrez = $tempentrez;
    }
    $stm->finish();



  }
 
  if ($id_type[0] eq 'P') {
    $macroprocess = $id_type[1];
    $entityType = "macroprocess";
  }

  $label =~ s/\n//g;
   
  if ($font_family = "arial") {
    $font_family = "Arial";
  }
  
  if (defined $attrs->{fontsize}) {
    $font_size = $attrs->{fontsize}
  }
  my $testlen = 0;

  my @textarray; 
  if ($shape eq "Mrecord") {
    @textarray = split(/\|/, $text);
  } else {
    @textarray = split(/\\n/, $text);
  }
  
  for (my $count = 0; $count < $#textarray; $count++) {
    my $templen = length($textarray[$count]);
    if ($templen > $testlen) {
      $testlen = $templen;
    }
  }
  
  if ($testlen == 0) {
    $testlen = length($text);
  }
  $x = $x - 2.5 * $testlen;
 
  # my $y0 = $y;
  my $start = 0;
  if ($type eq 'mrecord') {
    $start = 2;
  }
  my $y0 = $y + ((scalar(@textarray) - 1) * ($start + $TEXT_LINE_INC_Y{$font_size}));
  $y0 = $y0 + 5; 
  my $count = 0;
  my $elts = 1;
  if ($#textarray > 0) {
    $elts = $#textarray;
    if ($elts > 2) {
      $y0 = $y0 - 2;
    }
  } 
  
  my $move_y = ($y2 - $y1)/($elts+1);
  if ($move_y == 0) {
    $move_y = ($delta + $TEXT_LINE_INC_Y{$font_size});
  } 
  $y1 = Invert($y1);
  $y1 = $y1 - $move_y;
  my $textCount = 0;
  for my $t (@textarray) {
    
    my $ltext = length($t); 
    for (my $count = 0; $count < $ltext; $count++) {
      if ((substr($t,$count,1) eq '<') || (substr($t,$count,1) eq '>') || (substr($t,$count,1) eq '{') || (substr($t,$count,1) eq '}')) {
        $x++;
      }
    }
 
    $t =~ s/{//g;
    $t =~ s/}//g;
    $t =~ s/\\<//g;
    $t =~ s/\\>//g;

    $t =~ s/[\<]+//g;
    $t =~ s/[\>]+//g;

    # $t =~ s/[\\]+//g;
    # $t =~ s/[\\]+//g;

    $textCount++;
    $y = Invert($y0);
    my $foreground="#000000";

    if ($type eq "box") {
      $foreground="#ffffff";
      $x = $x - $sevenadjust;
      $sevenadjust = 0;
    }
    
    if (($shape ne 'box') && ($shape ne 'parallelogram')) {
      print $fh qq!
      <HyperlinkButton  NavigateUri="http://pid.nci.nih.gov/search/MoleculePage?molid=$mol_id" >
      <HyperlinkButton.Content>
      !;
    }
  
    if ($shape eq 'parallelogram') {
      print $fh qq!
      <HyperlinkButton  NavigateUri="http://pid.nci.nih.gov/search/intermediate_landing.shtml?molecule=$text" >
      <HyperlinkButton.Content>
      !;

    }
 
    if (($shape eq 'box') && ($nolink == 0)) {
      $url =~s/svg/xaml/g;
      $url =~s/pathway_landing/xaml_landing/g;
      $url =~s/graphic/text/g;
      $url =~s/search\///g;
      $url =~s/source%3D5/source=NATURE/g;
      
      print $fh qq!
      <HyperlinkButton  NavigateUri="http://pid-stage.nci.nih.gov$url" >
      <HyperlinkButton.Content>
      !;
    }

  # } 
    print $fh "<Canvas Name=\"mol$molCount\"> \n";
    $molCount++;

    if ($shape eq 'Mrecord') {
      my @tarray = split(/\\n/, $t);
      my $tarraylength = $#tarray + 1;
      if ($tarraylength > 3) {
        $tarraylength--;
        $tmove = 1;
      }
      $y = $y - (($tarraylength) * 8);
      my $x_temp = $x1 + 3;
      for (my $tcount = 0; $tcount <= $#tarray; $tcount++) {
        $y = $y + 8;
        if ($tcount == $#tarray) {
          $y = $y - 1;
        }

        print $fh qq!
<TextBlock Tag="$keyCount" FontSize="$font_size" FontFamily="$font_family" HorizontalAlignment="Center" Canvas.Left="$x_temp" Canvas.Top="$y" Foreground="$foreground">$tarray[$tcount]</TextBlock>
        
!;

      }
    } else {
      # $x++;
      if (($entityType ne 'none') && ($entityType ne 'macroprocess')) {    
    print $fh qq!
<TextBlock Tag="$keyCount" FontSize="$font_size" FontFamily="$font_family" HorizontalAlignment="Center" Canvas.Left="$x" Canvas.Top="$y" Foreground="$foreground">$t</TextBlock>
!;
      } else {
        my $x_temp = $x;
        if ($shape eq 'box') {
          $x_temp = $x1 + 3;
        }
        print $fh qq!
<TextBlock FontSize="$font_size" FontFamily="$font_family" HorizontalAlignment="Center" Canvas.Left="$x_temp" Canvas.Top="$y" Foreground="$foreground">$t</TextBlock>
!;
      }
    }
    if (($entityType ne 'none') && ($shape ne 'parallelogram')) {
      my $locationState = 0;
      $location = '';
      for (my $count = 0; $count <= length($t); $count++) {
        if (substr($t,$count,1) eq ']') {
          $locationState = 0;
        }
 
        if ($locationState == 1) {
          $location = $location . substr($t,$count,1);
        } 
        
        if (substr($t,$count,1) eq '[') {
          $locationState = 1;
        }

      }
      if ($location eq 'n') {
        $location = 'nucleus';
      }
      if ($location eq 'cy') {
        $location = 'cytoplasm';
      }
      if ($location eq 'm') {
        $location = 'transmembrane';
      }
      
      my @entitysearch = split(/[\[]/,$t); 
      my $entityVal = lc $entitysearch[0];

      my $statelen = length($entityVal);
      $state = '';
      if (substr($entityVal,($statelen-1),1) eq '+') {
        $state = 'active';
        $entityVal = substr($entityVal,0,$statelen-1);
      } 
  
      if (substr($entityVal2,($statelen-2),2) eq '+1') {
        $state = 'active 1';
      }
      if (substr($entityVal2,($statelen-2),2) eq '+2') {
        $state = 'active 2';
      }
      if (substr($entityVal2,($statelen-2),2) eq '+3') {
        $state = 'active 3';
      }

      if (substr($entityVal,($statelen-1),1) eq '-') {
        $state = 'inactive';
        $entityVal = substr($entityVal,0,$statelen-1);
      }

      if (substr($entityVal2,($statelen-2),2) eq '-1') {
        $state = 'inactive 1';
      }
      if (substr($entityVal2,($statelen-2),2) eq '-2') {
        $state = 'inactive 2';
      }
      if (substr($entityVal2,($statelen-2),2) eq '-3') {
        $state = 'inactive 3';
      }
      my $sql = qq!
      select ext_mol_id, b.mol_id
      from pid.pw_ext_mol_id a, pid.pw_mol_srch b
      where id_type = 'UP'
      and a.mol_id = b.mol_id
      and b.map_name = '$entityVal'
      and a.mol_id between 200000 and 299999
      !;
      
      my $stm = $db->prepare($sql);
      if(not $stm) {
      print STDERR "$sql\n";
      print STDERR "$DBI::errstr\n";
      print STDERR "prepare call failed\n";
      die;
      }
      if(!$stm->execute()) {
        print STDERR "$sql\n";
        print STDERR "$DBI::errstr\n";
        print STDERR "execute call failed\n";
        die;
      }

      my $inmol_id = 0;
      while ((my $tempuniprot, my $tempmol_id) = $stm->fetchrow_array()) {
        $uniprot = $tempuniprot;
        $inmol_id = $tempmol_id
      }
      $stm->finish();
  
      $sql = qq!
      select
        g.ll_id
      from
        cgap.ll_gene g,
        pid.pw_ext_mol_id e,
        pid.pw_mol m,
        cgap.ll2sp s,
        cgap.sp_primary p

      where
        m.mol_id = $inmol_id
        and p.sp_id_or_secondary = e.ext_mol_id
        and p.sp_id_or_secondary = s.sp_primary
        and s.organism = 'Hs'
        and s.ll_id = g.ll_id
        and e.mol_id = m.mol_id
        and m.basic_mol_type = 'PR'

    !;

    my $stm = $db->prepare($sql);
    if(not $stm) {
      print STDERR "$sql\n";
      print STDERR "$DBI::errstr\n";
      print STDERR "prepare call failed\n";
      die;
    }

    if(!$stm->execute()) {
      print STDERR "$sql\n";
      print STDERR "$DBI::errstr\n";
      print STDERR "execute call failed\n";
      die;
    }

    while ((my $tempentrez) = $stm->fetchrow_array()) {
        $entrez = $tempentrez;
    }
    $stm->finish();
 
      open(META_FILE, ">> $meta_f") || die "Cannot open $meta_f file\n";
      print META_FILE "<pd:PIDEntity\n";
      print META_FILE "x:Key=\"$keyCount\"\n";
      print META_FILE "EntityName=\"$entityVal\"\n";
      print META_FILE "EntitySearch=\"$entityVal,$uniprot,$entrez\"\n";
      print META_FILE "EntityType=\"$entityType\"\n";
      print META_FILE "MoleculeType=\"$molType\"\n";
      print META_FILE "IntType=\"$intType\"\n";
      print META_FILE "Location=\"$location\"\n";
      print META_FILE "MacroProcess=\"$macroprocess\"\n";
      print META_FILE "PTMs=\"$ptms\"\n";
      print META_FILE "ActivityState=\"$state\"\n";
      print META_FILE "MoleId=\"$inmol_id\"\n";
      print META_FILE "AtomId=\"$atom_id\"\n";
      print META_FILE "UniprotID=\"$uniprot\"\n";
      print META_FILE "Entrez=\"$entrez\"\n";
      print META_FILE "/>\n";
      close(META_FILE);
      $keyCount++;
    }
    print $fh "</Canvas>\n";
    if (($shape ne 'box') && ($shape ne 'parallelogram')) {
      print $fh qq!
      </HyperlinkButton.Content>
      </HyperlinkButton>
      !;
    }

    if ($shape eq 'parallelogram') {
      print $fh qq!
      </HyperlinkButton.Content>
      </HyperlinkButton>
      !;
    }

    if (($shape eq 'box') && ($nolink == 0)) {
      print $fh qq!
      </HyperlinkButton.Content>
      </HyperlinkButton>
      !;
    }

    $y0 -= $move_y;
    $y = Invert($y0);
    
    if (($type eq 'mrecord') && (($count) < $#textarray)) {
      # print STDERR "<Polyline Points= $x1,$y1 $x2,$y1 FillRule=NonZero Stroke=#000000 />\n";
      if ($tmove == 1) {
        $y1 = $y1 + 10; 
      } else {
        $y1 = $y1 + 2;
      }
      print $fh qq!
<Polyline Points="$x1,$y1 $x2,$y1 " FillRule="NonZero" Stroke="#000000" />
!;
    }
    $count++;
    $y1 = $y1 - $move_y;
    # $y0 -= ($delta + $TEXT_LINE_INC_Y{$font_size}) / 2;
  } 
 
    if ($shape eq "parallelogram") {
      $entityType = "macroprocess";
    }
    if ($shape eq "plaintext") {
      $entityType = "molecule";
      $molType = "protein";
    }
    if ($shape eq "Mrecord") {
      $entityType = "molecule";
      $molType = "complex";
    }
    if (($shape eq "diamond") || ($shape eq "doublecircle") || ($shape eq "triangle")) {
      $prevlabel = 1;
      $label = '';
      $entityType = "interaction";
      $intType = "transcription";
      if ($shape eq "diamond") {
        $intType = "modification";
      }
      if ($shape eq "triangle") {
        $intType = "translocation";
      }
 
    }
    $label =~ s/{//g;
    $label =~ s/}//g;
    $label =~ s/[\<]+//g;
    $label =~ s/[\>]+//g;
    $label =~ s/[\\]+//g;
    $label =~ s/[\\]+//g;

    if (($prevlabel ne $label) && ($entityType ne "interaction") && ($textCount > 1) && ($entityType ne 'none')) {
       
      open(META_FILE, ">> $meta_f") || die "Cannot open $meta_f file\n";
      print META_FILE "<pd:PIDEntity\n";
      print META_FILE "x:Key=\"$keyCount\"\n";
      print META_FILE "EntityName=\"$entityVal\"\n";
      print META_FILE "EntitySearch=\"$entityVal,$uniprot,$entrez\"\n";
      print META_FILE "EntityType=\"$entityType\"\n";
      print META_FILE "MoleculeType=\"$molType\"\n";
      print META_FILE "IntType=\"$intType\"\n";
      print META_FILE "Location=\"$location\"\n";
      print META_FILE "MacroProcess=\"$macroprocess\"\n";
      print META_FILE "PTMs=\"$ptms\"\n";
      print META_FILE "ActivityState=\"$state\"\n";
      print META_FILE "MoleId=\"$mol_id\"\n";
      print META_FILE "AtomId=\"$atom_id\"\n";
      print META_FILE "UniprotID=\"$uniprot\"\n";
      print META_FILE "Entrez=\"$entrez\"\n";
      print META_FILE "/>\n";
      close(META_FILE);
      $keyCount++;
      $prevlabel = $label;
    } 
    if ($entityType eq 'interaction') {
      print META_FILE "Interaction\n"; 
      open(META_FILE, ">> $meta_f") || die "Cannot open $meta_f file\n";
      print META_FILE "<pd:PIDEntity\n";
      print META_FILE "x:Key=\"$keyCount\"\n";
      print META_FILE "EntityName=\"$entityVal\"\n";
      print META_FILE "EntitySearch=\"$entityVal,$uniprot,$entrez\"\n";
      print META_FILE "EntityType=\"$entityType\"\n";
      print META_FILE "MoleculeType=\"$molType\"\n";
      print META_FILE "IntType=\"$intType\"\n";
      print META_FILE "Location=\"$location\"\n";
      print META_FILE "MacroProcess=\"$macroprocess\"\n";
      print META_FILE "PTMs=\"$ptms\"\n";
      print META_FILE "ActivityState=\"$state\"\n";
      print META_FILE "MoleId=\"$mol_id\"\n";
      print META_FILE "AtomId=\"$atom_id\"\n";
      print META_FILE "UniprotID=\"$uniprot\"\n";
      print META_FILE "Entrez=\"$entrez\"\n";
      print META_FILE "/>\n";
      close(META_FILE);
      $keyCount++;
    }
}

######################################################################
sub PrintEdge {
  my ($fh, $decl, $meta_f) = @_;
  
  my $attrs     = GetAttrs($decl);
  my $pos       = $attrs->{pos};
  my $stroke    = Color2Hex($attrs->{color});
  my $text      = $attrs->{label};
  my $id        = GetId($decl);

  my @dot_coords = split(/ +/, $pos);
  my @xaml_coords;

  my $path_data;
  my ($arrow, $end);
  if ($dot_coords[0] =~ /^(e|s),([0-9,]+)/) {
    my $dummy = shift @dot_coords;
    $dummy =~ /^(e|s),([0-9,]+)/;
    ($end, $arrow) = ($1, $2);
  }
  for my $xy (@dot_coords) {
    push @xaml_coords, Invert_Y($xy);
  }
  my $M      = shift @xaml_coords;
  $path_data = "M$M" . "C". join(" ", @xaml_coords);
  print $fh qq!
  <Path Stroke="#$stroke" Data="$path_data"/>
!;
  if ($arrow) {
    my $tip = $arrow;
    my $base;
    if ($end eq "e") {
      $base = $dot_coords[$#dot_coords];
    } elsif ($end eq "s") {
      $base = $dot_coords[0];
    }
    my @arrow_coords;
    for my $xy (@{ PlaceArrow($stroke, split(",", $base), split(",", $tip)) }) {
      push @arrow_coords, Invert_Y($xy);
    }
    my $point_data = join(" ", @arrow_coords);
    print $fh qq!
  <Polygon Points="$point_data" FillRule="NonZero" Fill="#$stroke" Stroke="#$stroke"/>
!;
  }
  if (defined $text) {
    PrintTextBlock($fh, $meta_f, 10, $attrs, $text, $id);
  }
}

######################################################################
sub PrintNode {
  my ($fh, $decl, $meta_f) = @_;

  my $attrs     = GetAttrs($decl);
  my $shape     = $attrs->{shape};
  my $url       = $attrs->{URL};
  my $pos       = $attrs->{pos};
  my $height    = $PIX_PER_INCH*($attrs->{height});
  my $width     = $PIX_PER_INCH*($attrs->{width});
  my $id        = GetId($decl);
  my @dot_coords = split(/ +/, $pos);

  my ($x, $y) = split(",", $dot_coords[0]);

  if ($shape eq "plaintext") {
    PlainText($fh, $meta_f, $x, $y, $id, $width, $height, $attrs);
  } elsif ($shape eq "box") {
    Box($fh, $meta_f, $x, $y, $id, $width, $height, $attrs);
  } elsif ($shape eq "diamond") {
    Diamond($fh, $meta_f, $x, $y, $id, $width, $height, $attrs);
  } elsif ($shape eq "triangle") {
    Triangle($fh, $meta_f, $x, $y, $id, $width, $height, $attrs);
  } elsif ($shape eq "parallelogram") {
    Parallelogram($fh, $meta_f, $x, $y, $id, $width, $height, $attrs);
  } elsif ($shape eq "doublecircle") {
    Circle($fh, $meta_f, $x, $y, $id, $attrs);
  } elsif ($shape eq "Mrecord") {
    Mrecord($fh, $meta_f, $x, $y, $id, $width, $height, $attrs);
  }
}

######################################################################
sub PrintDecl {
  my ($fh, $decl, $meta_f) = @_;

  my $decl_type = GetDeclType($decl);
  $canvas_n++;

  my $id        = GetId($decl);



  if ($decl_type eq "node") {
    my $attrs     = GetAttrs($decl);
    my $shape     = $attrs->{shape};
    my $label     = $attrs->{label};
    my $url       = $attrs->{URL};
    $_ = $label; 
    $label =~ s/[\<]//g;
    $label =~ s/[\>]//g;
    if ($shape ne 'parallelogram') { 
    }
    print $fh qq!
    <Canvas Name="can$canvas_n"> 
    !; 
    PrintNode($fh, $decl, $meta_f);
    print $fh qq!
    </Canvas>
    !;
    if ($shape ne 'parallelogram') {
    }
  } elsif ($decl_type eq "edge") {
    print $fh qq!
    <Canvas Name="can$canvas_n">
    !; 
    PrintEdge($fh, $decl, $meta_f);
      print $fh qq!
    </Canvas>
!;

  }

}

######################################################################
sub PrintDecls {
  my ($fh, $decls, $meta_f) = @_;

  for my $decl (@{ $decls }) {
    PrintDecl($fh, $decl, $meta_f);
  }
}

#######################################################################
sub PrintHeader {
  my ($fh) = @_;
  # print $fh "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
}

######################################################################
sub PrintMainCanvas {
  my ($fh, $metafile, $root, $x0, $x1, $y0, $y1) = @_;

  my $height = $y1 - $y0;
  my $width  = $x1 - $x0;

  my $offset = 4;
  my @coords;

  
  for my $xy (
      ($x0 - $offset) . "," . ($y0 - $offset),
      ($x0 - $offset) . "," . ($y1 + $offset),
      ($x1 + $offset) . "," . ($y1 + $offset),
      ($x1 + $offset) . "," . ($y0 - $offset),
      ($x0 - $offset) . "," . ($y0 - $offset)
    ) {
    push @coords, Invert_Y($xy);
  }
  my $points = join(" ", @coords);
  print $fh qq!
<Canvas xmlns:x="http://schemas.microsoft.com/winfx/2006/xaml" xmlns="http://schemas.microsoft.com/winfx/2006/xaml/presentation" xmlns:pd="clr-namespace:PathwayViewer;assembly=PathwayViewer" Width="$width" Height="$height">
  <Canvas.Resources >
  
  </Canvas.Resources> 
  <Canvas.RenderTransform>
    <TranslateTransform X="0" Y="$height"></TranslateTransform>
  </Canvas.RenderTransform>
  <Canvas Name="graph0">
    <Canvas>
      <Canvas.RenderTransform>
        <TransformGroup></TransformGroup>
      </Canvas.RenderTransform>
      <Polygon Points="$points" FillRule="NonZero" Fill="#ffffff" Stroke="#ffffff"/>
  
!;
  PrintDecls($fh, Decls($root), $metafile);
  print $fh qq!
    </Canvas>
  </Canvas>
</Canvas>
!;

}

######################################################################
sub PrintXAML {
  my ($fh, $root, $metafile, $x0, $x1, $y0, $y1) = @_;

  PrintHeader($fh);
  PrintMainCanvas($fh, $metafile, $root, $x0, $x1, $y0, $y1);

}

######################################################################
1;
#####################################################################
