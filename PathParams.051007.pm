

# Copyright SRA International
#
# Distributed under the OSI-approved BSD 3-Clause License.
# See http://ncip.github.com/pathway-interaction-database/LICENSE.txt for details.


package PathParams;
require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(
  ReadParams
  CheckParams
);

use strict;
use PWLabel;

use constant ORACLE_LIST_LIMIT => 500;

use constant PRINT             => "print";

use constant SHOW              => "show";

use constant MACRO_PROCESS     => "macro_process";

use constant SKIP_MOL_ID       => "skip_mol_id";
use constant SKIP_COMPLEXES    => "skip_complexes";

use constant SOURCE_ID         => "source_id";

use constant COLLAPSE          => "collapse";

use constant COLOR             => "color";

use constant VALUE             => "value";

use constant CHANNEL           => "channel";

use constant CONNECT_MOLECULE  => "connect_molecule";
use constant CONNECT_MOL_NAME  => "connect_mol_name";
use constant CONNECT_MOL_ID    => "connect_mol_id";
use constant CONNECT_MOL_EXT_ID => "connect_mol_ext_id";

use constant PRUNE_MOL_NAME    => "prune_mol_name";
use constant PRUNE_MOL_ID      => "prune_mol_id";
use constant PRUNE_MOL_EXT_ID  => "prune_mol_ext_id";

use constant PATHWAY           => "pathway";
use constant PATHWAY_ID        => "pathway_id";
use constant PATHWAY_NAME      => "pathway_name";
use constant PATHWAY_EXT_ID    => "pathway_ext_id";

use constant ATOM_ID           => "atom_id";
use constant EVIDENCE_CODE     => "evidence_code";
use constant PRUNE_ATOM_ID     => "prune_atom_id";

use constant MOLECULE          => "molecule";
use constant MOL_NAME          => "mol_name";
use constant MOL_ID            => "mol_id";
use constant MOL_EXT_ID        => "mol_ext_id";

use constant INCLUDE           => "include";

use constant EDGE              => "edge";

use constant DEGREE_BACK       => "degree_back";
use constant DEGREE_FORWARD    => "degree_forward";

use constant DB_USER           => "db_user";
use constant DB_PASS           => "db_pass";
use constant DB_INST           => "db_inst";
use constant DB_SCHEMA         => "db_schema";

use constant SIM_CYCLE         => "sim_cycle";
use constant SIM_MOL_ID        => "sim_mol_id";
use constant SIM_OUTPUT        => "sim_output";
use constant SIM_METHOD        => "sim_method";
use constant SIM_SIMPLE        => "sim_simple";
use constant SIM_COMPETE       => "sim_compete";

use constant MAX_DEGREE        => 5;

my %mol_types = ("CM" => "compound", "CX" => "complex",
    "RN" => "rna", "PR" => "protein");

my %id_types = (
    "EC" => 1,
    "LL" => 1, 
    "UP" => 1, 
    "CA" => 1,
    "KG" => 1
);

my %printopts = (
    "sql"      => 1,
    "xml"      => 1,
    "biopax"   => 1,
    "sbml"     => 1,
    "lisp"     => 1,
    "petrinet" => 1,
    "sas"      => 1,
    "dot"      => 1,
    "template" => 1
);

######################################################################
sub new {
  my ($self, $lv) = @_;
  my $x = {};
  $x->{lv} = $lv;
  return bless $x;
}

######################################################################
sub Errors {
  my ($self) = @_;

  return $self->{errors};  
}

######################################################################
sub AddError {
  my ($self, $msg) = @_;

#print STDERR "$msg\n";
  push @{ $self->{errors} }, $msg;

}

######################################################################
sub Trim {
  my ($x) = @_;
  $x =~ s/^\s+//;
  $x =~ s/\s+$//;
  return $x;
}

######################################################################
sub DoDbUser {
  my ($self, $user) = @_;
  $self->{db}{user} = $user;
}

######################################################################
sub DoDbPass {
  my ($self, $pass) = @_;
  $self->{db}{pass} = $pass;
}

######################################################################
sub DoDbInst {
  my ($self, $inst) = @_;
  $self->{db}{inst} = $inst;
}

######################################################################
sub DoDbSchema {
  my ($self, $schema) = @_;
  $self->{db}{schema} = $schema;
}

######################################################################
sub DoPrint {
  my ($self, $what) = @_;
  $what = lc(Trim($what));
  if (not defined $printopts{$what}) {
    $self->AddError("unrecognized print option $what");
    return;
  }
  $self->{print}{$what} = 1;
}

######################################################################
sub DoSourceId {
  my ($self, $id) = @_;

  $id = lc(Trim($id));
  if ($id =~ /^\d+$/) {
    $self->{source_id}{$id} = 1;
  }
}

######################################################################
sub DoMacroProcess {
  my ($self, $name) = @_;

  ## Don't lower-case; PWLabel retains case-sensitivitiy

  $name = Trim($name);
  $self->{macroprocess_name}{$name} = 1;
}

######################################################################
sub DoShow {
  my ($self, $what) = @_;

  $what = lc(Trim($what));
  if (
      $what eq "atom_id" ||
      $what eq "subtype_line" ||
      0 ) {
    $self->{show}{$what} = 1;
  } else {
    $self->AddError("unrecognized option $what for show");
    return;
  }
}

######################################################################
sub DoSkip {
  my ($self, $cmd, $mol_id) = @_;

  $cmd = lc(Trim($cmd));
  if ($cmd eq "complexes") {
    $self->{skip}{complexes} = 1;
  } elsif ($cmd eq "mol_id") {
    $self->{skip}{mol_id}{$mol_id} = 1;
  }
}

######################################################################
sub DoCollapse {
  my ($self, $cmd) = @_;

  $cmd = lc(Trim($cmd));
  if ($cmd eq "processes") {
    $self->{collapse}{processes} = 1;
  } elsif ($cmd eq "molecules") {
    $self->{collapse}{molecules} = 1;
  } elsif ($cmd eq "subnets") {
    $self->{collapse}{subnets} = 1;
  } else {
    $self->AddError("unrecognized print option $cmd");
    return;
  }
}

######################################################################
sub DoColor {
  my ($self, $number, $color) = @_;

#  $boundtype = lc(Trim($boundtype));
  $number    = lc(Trim($number));
  $color     = lc(Trim($color));

#  if ($boundtype ne "upper" || $boundtype ne "lower") {
#    $self->AddError("unrecognized bound $boundtype in color command");
#    return;
#  }
  if ($number !~ /^(-|\+)?\d+(\.\d*)?$/ && $number !~ /^(-|\+)?\.\d+$/) {
    $self->AddError("number expected for $number in color command");
    return;
  }
  if ($color !~ /^#?([0-9A-Fa-f]{6})/) {
    $self->AddError("six-place hexadecimal color specification " .
        "expected for $color");
    return;
  } else {
    $color = $1;
  }
  $self->{color}{$number} = $color;
}

######################################################################
sub DoSim {
  my ($self, $cmd, @vars) = @_;

  if ($cmd eq "cycle") {
    my $ncycle = Trim(shift @vars);
    if ($ncycle !~ /^\d+$/) {
      $self->AddError("need numeric argument for 'sim'");
      return;
    }
    $self->{sim}{ncycle} = $ncycle;
  } elsif ($cmd eq "mol_id") {
    my ($mol_id, $pct) = @vars;
    $self->{sim}{mol}{$mol_id} = $pct;
  } elsif ($cmd eq "output") {
    my $output_type = lc(Trim(shift @vars));
    $self->{sim}{output} = $output_type;
  } elsif ($cmd eq "method") {
    my $meth = lc(Trim(shift @vars));
    $self->{sim}{method} = $meth;
    if ($meth eq "boolean") {
      my $default = shift @vars;
      if ($default ne "") {
        $self->{sim}{booleandefault} = $default;
      }
    }
  } elsif ($cmd eq "simple") {
    my $flag = Trim(shift @vars);
    if ($flag ne "0" && $flag ne "" && $flag ne "1") {
      $self->AddError("unrecognized sim_simple option $flag");
      return;
    }
    $self->{sim}{simple} = $flag;
  } elsif ($cmd eq "compete") {
    my $flag = Trim(shift @vars);
    if ($flag ne "0" && $flag ne "" && $flag ne "1") {
      $self->AddError("unrecognized sim_compete option $flag");
      return;
    }
    $self->{sim}{compete} = $flag;
  } else {
    $self->AddError("unrecognized sim command $cmd");
    return;
  }
 
}

######################################################################
sub DoEdge {
  my ($self, @vars) = @_;

  my ($interaction_id, $edge_type, $mol_name, $mol_ext_id,
      $activity, $location) = @vars;
  $interaction_id = Trim($interaction_id);
  $edge_type      = lc(Trim($edge_type));
  $mol_name       = Trim($mol_name);
  $mol_ext_id     = Trim($mol_ext_id);
  $activity       = lc(Trim($activity));
  $location       = lc(Trim($location));
  if ($interaction_id eq "") {
    $self->AddError("missing interaction id in edge command: " .
        join("\t", @vars) );
    return;
  }
  if ($edge_type eq "") {
    $self->AddError("missing edge type in edge command: " .
        join("\t", @vars) );
    return;
  }
  $edge_type = lc($edge_type);
  if ($edge_type ne "input" && $edge_type ne "output" &&
      $edge_type ne "agent" && $edge_type ne "inhibitor") {
    $self->AddError("illegal edge type '$edge_type' in edge command: " .
        join("\t", @vars) );
    return;
  }
  if ($mol_name eq "" && $mol_ext_id eq "") {
    $self->AddError("must specify mol name and/or mol external id in edge command: " .
        join("\t", @vars) );
    return;
  }
  if ($activity) {
    $self->{privateactivity}{$activity} = 1;
  }
  if ($location) {
    $self->{privatelocation}{$location} = 1;
  }

  if ($mol_ext_id =~ /^\d+$/) {
    $self->{privatells}{$mol_ext_id} = 1;
  }

  push @{ $self->{privateedges} }, join("\t", $interaction_id, $edge_type,
      $mol_name, $mol_ext_id, $activity, $location);
}

######################################################################
sub DoDegree {
  my ($self, $direction, @vars) = @_;

  my $n = Trim(shift @vars);

  if ($direction ne "back" && $direction ne "forward") {
    $self->AddError("unrecognized degree direction");
    return;
  }
  if ($n !~ /^\d+$/) {
    $self->AddError("illegal value $n for degree $direction");
    return;
  } elsif ($n > MAX_DEGREE) {
    $self->AddError("value for degree $n exceeds maximum " . MAX_DEGREE);
    return;
  } else {
    $self->{degree}{$direction} = $n;
  }
}

######################################################################
sub DoInclude {
  my ($self, $what) = @_;

  $what = lc(Trim($what));
  if ($what eq "complex_uses") {
    $self->{include_complex_uses} = 1;
  } elsif ($what eq "family_uses") {
    $self->{include_family_uses} = 1;
  } else {
    $self->AddError("Illegal option $what for include");
    return;
  }
}

######################################################################
sub  MouseTrap {
  my ($self) = @_;

##! Traps only LocusLink ids for mouse

  my (%temp);

  for $_ (keys %{ $self->{mol_ext_id} }) {
    if (/^\d+$/) {
      $temp{$_} = 1;
    }
  }
  for $_ (keys %{ $self->{mol_ext_id_value} }) {
    if (/^\d+$/) {
      $temp{$_} = 1;
    }
  }
  for $_ (keys %{ $self->{privatells} }) {
    $temp{$_} = 1;
  }
  my @temp = keys %temp;
  if (@temp == 0) {
    return;
  }

  my ($db_inst, $db_user, $db_pass, $schema) =
      ("cgprod", "web", "readonly", "cgap");
  my ($db, $sql, $stm, $i, $list);
  my ($hs_locus, $mm_locus);


  $db = DBI->connect("DBI:Oracle:" . $db_inst, $db_user, $db_pass);
  if (not $db or $db->err()) {
    print STDERR "Cannot connect to " . $db_user . "@" . $db_inst . "\n";
    exit;
  }

  for($i = 0; $i < @temp; $i += ORACLE_LIST_LIMIT) {
    if(($i + ORACLE_LIST_LIMIT - 1) < @temp) {
      $list = "'" . join("','", @temp[$i..$i+ORACLE_LIST_LIMIT-1]) . "'";
    }
    else {
      $list = "'" . join("','", @temp[$i..@temp-1]) . "'";
    }

    $sql = qq!
select
  hs.locus_id,
  mm.locus_id
from
  $schema.hg_entry mm,
  $schema.hg_entry hs
where
      mm.taxon_id = 10090
  and hs.taxon_id = 9606
  and hs.entry_id = mm.entry_id
  and mm.locus_id in ($list)
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
    while (($hs_locus, $mm_locus) = $stm->fetchrow_array()) {
      push @{ $self->{mousetrap}{$mm_locus} }, $hs_locus;
    }
  }

  $db->disconnect();

  for $mm_locus (keys %{ $self->{mol_ext_id} }) {
    if (defined $self->{mol_ext_id}{$mm_locus}) {
      if (defined $self->{mousetrap}{$mm_locus}) {
        for $hs_locus (@{ $self->{mousetrap}{$mm_locus} }) {
          $self->{mol_ext_id}{$hs_locus} = $self->{mol_ext_id}{$mm_locus};
        }
        delete $self->{mol_ext_id}{$mm_locus};
        if (keys %{ $self->{mol_ext_id}{$mm_locus} } == 0) {
          delete $self->{mol_ext_id}{$mm_locus};
        }
      }
    }
  }
  
  for $mm_locus (keys %{ $self->{mol_ext_id_value} }) {
    if (defined $self->{mol_ext_id_value}{$mm_locus}) {
      if (defined $self->{mousetrap}{$mm_locus}) {
        for $hs_locus (@{ $self->{mousetrap}{$mm_locus} }) {
          for my $subtype (keys %{ $self->{mol_ext_id_value}{$mm_locus} }) {
            $self->{mol_ext_id_value}{$hs_locus}{$subtype} =
                $self->{mol_ext_id_value}{$mm_locus}{$subtype};
          }
        }
        delete $self->{mol_ext_id_value}{$mm_locus};
        if (keys %{ $self->{mol_ext_id_value}{$mm_locus} } == 0) {
          delete $self->{mol_ext_id_value}{$mm_locus};
        }
      }
    }
  }

}

######################################################################
sub DoChannel {
  my ($self, $channel, $molecule) = @_;

  my $lv = $self->{lv};
  $channel = uc(Trim($channel));
  if ($channel eq "") {
    $self->AddError("No channel specified for channel command: " .
        join("\t", $channel, $molecule));
    return;
  }
  if ($channel !~ /^(A|B)$/) {
    $self->AddError("Illegal channel specified for channel command: " .
        join("\t", $channel, $molecule));
    return;
  }

  $self->{mol_name}{lc(Trim($molecule))}{channel} = 1;
  $self->{mol_ext_id}{lc(Trim($molecule))}{channel} = 1;
  $self->{molecule_channel}{lc(Trim($molecule))}{$channel} = 1;

}

######################################################################
sub DoValue {
  my ($self, $value, $what, @vars) = @_;

  my $lv = $self->{lv};
  my ($molid, $name, $extid, $activity, $location);
  $value = Trim($value);
  if ($value eq "") {
    $self->AddError("No value specified for value command: " .
        join("\t", $what, $value, @vars));
    return;
  }
  if ($value !~ /^(-|\+)?(\d+)?(\.\d+)?$/) {
    $self->AddError("Illegal value specified for value command: " .
        join("\t", $what, $value, @vars));
    return;
  }

  if ($what eq "mol_name" || $what eq "molecule") {
    ($name, $activity, $location) = @vars;
  } elsif ($what eq "mol_id") {
    ($molid, $activity, $location) = @vars;
    $molid = Trim($molid);
    if ($molid !~ /^\d+$/) {
      $self->AddError("illegal molecule id $molid for DoMol");
      return;
    }
  } elsif ($what eq "mol_ext_id" || $what eq "molecule") {
    ($extid, $activity, $location) = @vars;
  } else {
    $self->AddError("Illegal molecular specifier for value command: " .
        join("\t", $what, $value, @vars));
    return;
  }

  $activity = lc(Trim($activity));
  $location = Trim($location);
  if ($activity ne "") {
    my $lvid = $lv->StringToLabelValue("activity-state", $activity);
    if (! $lvid) {
      $self->AddError("illegal activity value: $activity");
      return;
    } else {
#      $activity = $lvid;
    }
  }

  if ($location ne "") {
    my $lvid = $lv->StringToLabelValue("location", $location);
    if (! $lvid) {
      $self->AddError("illegal location value: $location");
      return;
    } else {
#      $location = $lvid;
    }
  }

  my $subtype = "$activity\t$location";
  if ($what eq "mol_name" || $what eq "molecule") {
    $self->{mol_name_value}{lc(Trim($name))}{$subtype} = $value;
    $self->{mol_name}{lc(Trim($name))}{value} = 1;
  } elsif ($what eq "mol_id") {
    $molid = Trim($molid);
    $self->{mol_id_value}{$molid}{$subtype} = $value;
  } elsif ($what eq "mol_ext_id" || $what eq "molecule") {
    $self->{mol_ext_id}{lc(Trim($extid))}{value}{$subtype} = $value;
    $self->{mol_ext_id_value}{lc(Trim($extid))}{$subtype} = $value;
  }

}

######################################################################
sub DoMol {
  my ($self, $op, $what, @vars) = @_;

  my ($name, $molid, $extid);

  if ($what eq "name") {
    ($name) = @vars;
    $self->{mol_name}{lc(Trim($name))}{$op} = 1;
  } elsif ($what eq "id") {
    ($molid) = @vars;
    $molid = Trim($molid);
    if ($molid !~ /^\d+$/) {
      $self->AddError("illegal molecule id $molid for DoMol");
      return;
    } else {
      if ($op eq "prune") {
        $self->{prune_mol_id}{$molid} = 1;  ## may still be an unknown molid
      } elsif ($op eq "get") {
        $self->{mol_id}{$molid} = 1;  ## may still be an unknown molid
      } elsif ($op eq "connect") {
        $self->{connect}{$molid} = 1; 
      } else {
        $self->AddError("unrecognized prune/get/connect action $op");
        return;
      }
    }
  } elsif ($what eq "ext_id") {
    ($extid) = @vars;
    $self->{mol_ext_id}{lc(Trim($extid))}{$op} = 1;
  } elsif ($what eq "general") {
    ## by name
    ($name) = @vars;
    $self->{mol_name}{lc(Trim($name))}{$op} = 1;
 
    ## by ext id
    ($extid) = @vars;
    $self->{mol_ext_id}{lc(Trim($extid))}{$op} = 1;
  } else {
    $self->AddError("unrecognized subtype $what for DoMol");
    return;
  }

}

######################################################################
sub DoAtomId {
  my ($self, $what, $atom_id) = @_;

  $atom_id = Trim($atom_id);
  if ($atom_id =~ /^\*$/) {
    if ($what eq "get") {
      $self->{atom}{$atom_id} = 1;
    }
  } elsif ($atom_id =~ /^\d+$/) {
    if ($what eq "get") {
      $self->{atom}{$atom_id} = 1;
    } elsif ($what eq "prune") {
      $self->{prune_atom_id}{$atom_id} = 1;
    }
  } else {
    $self->AddError("illegitimate atom id $atom_id");
    return;
  }

}

######################################################################
sub DoEvidenceCode {
  my ($self, $what, $evidence_code) = @_;

  $evidence_code = uc(Trim($evidence_code));
  if ($evidence_code =~ /^\*$/) {
    if ($what eq "get") {
      $self->{evidence_code}{$evidence_code} = 1;
    }
  } elsif ($evidence_code =~ /^[INRT][A-S]{1,2}$/i) {
    if ($what eq "get") {
      $self->{evidence_code}{$evidence_code} = 1;
    }
  } else {
    $self->AddError("illegitimate evidence_code $evidence_code");
    return;
  }

}

######################################################################
sub DoPathway {
  my ($self, $what, @vars) = @_;

  my ($pathway_id, $pathway_ext_id, $pathway_name);

  if ($what eq "id") {
    ($pathway_id) = @vars;
    $pathway_id = Trim($pathway_id);
    if ($pathway_id !~ /^\d+$/ && $pathway_id ne "*") {
      $self->AddError("illegitimate pathway id $pathway_id");
      return;
    }
    $self->{pathway}{id}{$pathway_id} = 1;
  } elsif ($what eq "extid") {
    ($pathway_ext_id) = @vars;
    $self->{pathway}{ext_id}{$pathway_ext_id} = 1;
  } elsif ($what eq "name") {
    ($pathway_name) = @vars;
    if ($pathway_name !~ /^\s+$/) {
      $self->{pathway}{name}{lc(Trim($pathway_name))} = 1;
    }
  } elsif ($what eq "general") {
    ($pathway_name) = @vars;
    if ($pathway_name !~ /^\s+$/) {
      $self->{pathway}{name}{lc(Trim($pathway_name))} = 1;
    }
    ($pathway_ext_id) = @vars;
    $self->{pathway}{ext_id}{$pathway_ext_id} = 1;
  } else {
    $self->AddError("unrecognized subtype $what for DoPathway");
    return;
  }

}

######################################################################
sub CheckMol {
  my ($self, $db) = @_;

  if (defined $self->{mol_name}) {
    CheckMolName($self, $db);
  }
## doing privatells to catch them in mousetrap
  if (defined $self->{mol_ext_id} || defined $self->{privatells} ||
      defined $self->{mol_ext_id_value}) {
    CheckMolExtId($self, $db);
  }
  for my $x (keys %{ $self->{mol_name} }, keys %{ $self->{mol_ext_id} }) {
    if (! defined $self->{molmap}{$x}) {
      $self->{no_molmap}{$x} = 1;
    }
  }
}

######################################################################
sub CheckMolName {
  my ($self, $db) = @_;

##!!!!!!! We aren't going to require mol type any longer

  my ($sql, $stm);
  my ($molid, $molname, $nametype, $moltype, $what, $ll_id, $ll_list);
  my ($alias, $alias_list);
  my $schema = $self->{db}{schema};

  my @temp = keys %{ $self->{mol_name} };
  my $list;
  $ll_list = 999999;
  for(my $i = 0; $i < @temp; $i += ORACLE_LIST_LIMIT) {
    if(($i + ORACLE_LIST_LIMIT - 1) < @temp) {
      $list = "'" . join("','", @temp[$i..$i+ORACLE_LIST_LIMIT-1]) . "'";
    }
    else {
      $list = "'" . join("','", @temp[$i..@temp-1]) . "'";
    }
  
  $list =~ s/,/','/;
  $list = lc $list; 
# Find alias
  $sql = qq!
select
  ll_id
from $schema.pw_molecule_search
where lower(symbol) in ($list)
  !;

print STDERR "SQL:  $sql\n";
$stm = $db->prepare($sql);
    if(not $stm) {
      $self->AddError("$sql");
      $self->AddError("$DBI::errstr");
      $self->AddError("prepare call failed");
      return;
    }

   if(!$stm->execute()) {
      $self->AddError("$sql");
      $self->AddError("$DBI::errstr");
      $self->AddError("execute call failed");
      return;
    }
   while (($ll_id) =
        $stm->fetchrow_array()) {
        $ll_list = $ll_list . ',' . $ll_id; 
   } 

$sql = qq!
select symbol from cgap.ll_gene
where ll_id in ($ll_list)
  !;

$stm = $db->prepare($sql);
    if(not $stm) {
      $self->AddError("$sql");
      $self->AddError("$DBI::errstr");
      $self->AddError("prepare call failed");
      return;
    }

   if(!$stm->execute()) {
      $self->AddError("$sql");
      $self->AddError("$DBI::errstr");
      $self->AddError("execute call failed");
      return;
    }
 while (($alias) =
        $stm->fetchrow_array()) {
        $alias = lc(Trim($alias));
        $self->{mol_name}{lc(Trim("$alias"))}{"get"} = 1;
        $alias_list = $alias_list . ',' . $alias; 
        $list = $list . "," . "'" . $alias . "'";
   }

print STDERR "List:  $list\n";
print STDERR "AliasList:  $alias_list\n";   
  ## Find exact matches to the name in PID

    $sql = qq!
select
  n.mol_id,
  lower(n.mol_name),
  n.name_type,
  m.basic_mol_type
from
  $schema.pw_mol m,
  $schema.pw_mol_name n
where
      m.mol_id = n.mol_id
  and lower(n.mol_name) in ($list)
  !;

    $stm = $db->prepare($sql);
    if(not $stm) {
      $self->AddError("$sql");
      $self->AddError("$DBI::errstr");
      $self->AddError("prepare call failed");
      return;
    }
    if(!$stm->execute()) {
      $self->AddError("$sql");
      $self->AddError("$DBI::errstr");
      $self->AddError("execute call failed");
      return;
    }
    while (($molid, $molname, $nametype, $moltype) =
        $stm->fetchrow_array()) {
      $self->{molmap}{$molname}{$molid} = 1;
      if (defined $self->{mol_name}{$molname}) {
        for $what (keys %{ $self->{mol_name}{$molname} }) {
          if ($what eq "prune") {
            $self->{prune_mol_id}{$molid} = 1;
          } elsif ($what eq "get") {
            $self->{mol_id}{$molid} = 1;
          } elsif ($what eq "connect") {
            $self->{connect}{$molid} = 1;
          } elsif ($what eq "value") {
            # handled immediately below
          } elsif ($what eq "channel") {
            # handled immediately below
          } else {
            $self->AddError("unrecognized action $what in CheckMolName");
            return;
          }
          if (defined $self->{mol_name_value}{$molname}) {
            for my $subtype (keys %{ $self->{mol_name_value}{$molname} }) {
              $self->{mol_id_value}{$molid}{$subtype} =
                  $self->{mol_name_value}{$molname}{$subtype};
	    }
          }
          if (defined $self->{molecule_channel}{$molname}) {
            for my $channel (keys %{ $self->{molecule_channel}{$molname} }) {
              $self->{mol_id_channel}{$molid}{$channel} = 1;
            }
          }
        }
      }
    }

  ## Find matches of gene names to proteins in PID

    $sql = qq!
select
  e.mol_id,
  lower(g.symbol),
  e.id_type,
  m.basic_mol_type
from
  $schema.pw_mol m,
  $schema.pw_ext_mol_id e,
  cgap.ll2sp s,
  cgap.sp_primary p,
  cgap.ll_gene g
where
      m.mol_id = e.mol_id
  and p.sp_id_or_secondary = e.ext_mol_id
  and p.sp_id_or_secondary = s.sp_primary
  and s.organism = 'Hs'
  and s.ll_id = g.ll_id
  and lower(g.symbol) in ($list)
  !;

    $stm = $db->prepare($sql);
    if(not $stm) {
      $self->AddError("$sql");
      $self->AddError("$DBI::errstr");
      $self->AddError("prepare call failed");
      return;
    }
    if(!$stm->execute()) {
      $self->AddError("$sql");
      $self->AddError("$DBI::errstr");
      $self->AddError("execute call failed");
      return;
    }
    while (($molid, $molname, $nametype, $moltype) =
        $stm->fetchrow_array()) {
      $self->{molmap}{$molname}{$molid} = 1;
      if (defined $self->{mol_name}{$molname}) {
        for $what (keys %{ $self->{mol_name}{$molname} }) {
          if ($what eq "prune") {
            $self->{prune_mol_id}{$molid} = 1;
          } elsif ($what eq "get") {
            $self->{mol_id}{$molid} = 1;
          } elsif ($what eq "connect") {
            $self->{connect}{$molid} = 1;
          } elsif ($what eq "value") {
            # handled immediately below
          } elsif ($what eq "channel") {
            # handled immediately below
          } else {
            $self->AddError("unrecognized action $what in CheckMolName");
            return;
          }
          if (defined $self->{mol_name_value}{$molname}) {
            for my $subtype (keys %{ $self->{mol_name_value}{$molname} }) {
              $self->{mol_id_value}{$molid}{$subtype} =
                  $self->{mol_name_value}{$molname}{$subtype};
	    }
          }
          if (defined $self->{molecule_channel}{$molname}) {
            for my $channel (keys %{ $self->{molecule_channel}{$molname} }) {
              $self->{mol_id_channel}{$molid}{$channel} = 1;
            }
          }
        }
      }
    }

  }
}

######################################################################
sub CheckMolExtId {
  my ($self, $db) = @_;

  ##!! convert Mm locus ids to Hs locus ids
  $self->MouseTrap();

  my ($sql, $stm);
  my ($molid, $extid, $what);

  my $schema = $self->{db}{schema};

  my @temp = keys %{ $self->{mol_ext_id} };
  my $list;

  for(my $i = 0; $i < @temp; $i += ORACLE_LIST_LIMIT) {
    if(($i + ORACLE_LIST_LIMIT - 1) < @temp) {
      $list = "'" . join("','", @temp[$i..$i+ORACLE_LIST_LIMIT-1]) . "'";
    }
    else {
      $list = "'" . join("','", @temp[$i..@temp-1]) . "'";
    }

  ## Find direct matches of ext ids

  $sql = qq!
select
  e.mol_id,
  lower(e.ext_mol_id)
from
  $schema.pw_ext_mol_id e
where
  lower(e.ext_mol_id) in ($list)
!;

    $stm = $db->prepare($sql);
    if(not $stm) {
      $self->AddError("$sql");
      $self->AddError("$DBI::errstr");
      $self->AddError("prepare call failed");
      return;
    }
    if(!$stm->execute()) {
      $self->AddError("$sql");
      $self->AddError("$DBI::errstr");
      $self->AddError("execute call failed");
      return;
    }
    while (($molid, $extid) =
        $stm->fetchrow_array()) {
      $self->{molmap}{$extid}{$molid} = 1;
      if (defined $self->{mol_ext_id}{$extid}) {
        for $what (keys %{ $self->{mol_ext_id}{$extid} }) {
          if ($what eq "prune") {
            $self->{prune_mol_id}{$molid} = 1;
          } elsif ($what eq "get") {
            $self->{mol_id}{$molid} = 1;
          } elsif ($what eq "connect") {
            $self->{connect}{$molid} = 1;
          } elsif ($what eq "value") {
            # handled below
          } elsif ($what eq "channel") {
            # handled below
          } else {
            $self->AddError("unrecognized action $what in CheckMolExtId");
            return;
          }
        }
      }
      if (defined $self->{mol_ext_id_value}{$extid}) {
        for my $subtype (keys %{ $self->{mol_ext_id_value}{$extid} }) {
          $self->{mol_id_value}{$molid}{$subtype} =
              $self->{mol_ext_id_value}{$extid}{$subtype};
        }
      }
      if (defined $self->{molecule_channel}{$extid}) {
        for my $channel (keys %{ $self->{molecule_channel}{$extid} }) {
          $self->{mol_id_channel}{$molid}{$channel} = 1;
        }
      }
    }
  }

  undef @temp;
  for my $num (keys %{ $self->{mol_ext_id} }) {
    if ($num =~ /^\d+$/) {
      push @temp, $num;
    }
  }
  
  my $list;

  for(my $i = 0; $i < @temp; $i += ORACLE_LIST_LIMIT) {
    if(($i + ORACLE_LIST_LIMIT - 1) < @temp) {
      $list = "'" . join("','", @temp[$i..$i+ORACLE_LIST_LIMIT-1]) . "'";
    }
    else {
      $list = "'" . join("','", @temp[$i..@temp-1]) . "'";
    }
  ## Find matches of gene ids to proteins in PID

    $sql = qq!
select
  e.mol_id,
  g.ll_id
from
  $schema.pw_ext_mol_id e,
  cgap.ll2sp s,
  cgap.sp_primary p,
  cgap.ll_gene g
where
      p.sp_id_or_secondary = e.ext_mol_id
  and p.sp_id_or_secondary = s.sp_primary
  and s.organism = 'Hs'
  and s.ll_id = g.ll_id
  and g.ll_id in ($list)
  !;

    $stm = $db->prepare($sql);
    if(not $stm) {
      $self->AddError("$sql");
      $self->AddError("$DBI::errstr");
      $self->AddError("prepare call failed");
      return;
    }
    if(!$stm->execute()) {
      $self->AddError("$sql");
      $self->AddError("$DBI::errstr");
      $self->AddError("execute call failed");
      return;
    }
    while (($molid, $extid) =
        $stm->fetchrow_array()) {
      if (defined $self->{mol_ext_id}{$extid}) {
      $self->{molmap}{$extid}{$molid} = 1;
        for $what (keys %{ $self->{mol_ext_id}{$extid} }) {
          if ($what eq "prune") {
            $self->{prune_mol_id}{$molid} = 1;
          } elsif ($what eq "get") {
            $self->{mol_id}{$molid} = 1;
          } elsif ($what eq "connect") {
            $self->{connect}{$molid} = 1;
          } elsif ($what eq "value") {
            # handled below
          } elsif ($what eq "channel") {
            # handled below
          } else {
            $self->AddError("unrecognized action $what in CheckMolExtId");
            return;
          }
        }
      }
      if (defined $self->{mol_ext_id_value}{$extid}) {
        for my $subtype (keys %{ $self->{mol_ext_id_value}{$extid} }) {
          $self->{mol_id_value}{$molid}{$subtype} =
              $self->{mol_ext_id_value}{$extid}{$subtype};
        }
      }
      if (defined $self->{molecule_channel}{$extid}) {
        for my $channel (keys %{ $self->{molecule_channel}{$extid} }) {
          $self->{mol_id_channel}{$molid}{$channel} = 1;
        }
      }
    }

  ## Find matches of gene ids to symbols in PID

    $sql = qq!
select
  n.mol_id,
  g.ll_id
from
  $schema.pw_mol m,
  $schema.pw_mol_name n,
  cgap.ll_gene g
where
      m.basic_mol_type = 'PR'
  and m.organism = 'Hs'
  and g.organism = 'Hs'
  and lower(g.symbol) = lower(n.mol_name)
  and m.mol_id = n.mol_id
  and g.ll_id in ($list)
  !;

    $stm = $db->prepare($sql);
    if(not $stm) {
      $self->AddError("$sql");
      $self->AddError("$DBI::errstr");
      $self->AddError("prepare call failed");
      return;
    }
    if(!$stm->execute()) {
      $self->AddError("$sql");
      $self->AddError("$DBI::errstr");
      $self->AddError("execute call failed");
      return;
    }
    while (($molid, $extid) =
        $stm->fetchrow_array()) {
      if (defined $self->{mol_ext_id}{$extid}) {
      $self->{molmap}{$extid}{$molid} = 1;
        for $what (keys %{ $self->{mol_ext_id}{$extid} }) {
          if ($what eq "prune") {
            $self->{prune_mol_id}{$molid} = 1;
          } elsif ($what eq "get") {
            $self->{mol_id}{$molid} = 1;
          } elsif ($what eq "connect") {
            $self->{connect}{$molid} = 1;
          } elsif ($what eq "value") {
            # handled below
          } elsif ($what eq "channel") {
            # handled below
          } else {
            $self->AddError("unrecognized action $what in CheckMolExtId");
            return;
          }
        }
      }
      if (defined $self->{mol_ext_id_value}{$extid}) {
        for my $subtype (keys %{ $self->{mol_ext_id_value}{$extid} }) {
          $self->{mol_id_value}{$molid}{$subtype} =
              $self->{mol_ext_id_value}{$extid}{$subtype};
        }
      }
      if (defined $self->{molecule_channel}{$extid}) {
        for my $channel (keys %{ $self->{molecule_channel}{$extid} }) {
          $self->{mol_id_channel}{$molid}{$channel} = 1;
        }
      }
    }

  }
}

######################################################################
sub CheckAtom {
  my ($self, $db) = @_;

  if (! defined $self->{atom}) {
    return;
  }

  if (defined $self->{atom}{"*"}) {
    undef %{ $self->{atom} };
    $self->{atom}{"*"} = 1;
    return;
  }

  my ($sql, $stm);
  my ($atom_id);
  my $schema = $self->{db}{schema};
  my $atomlist = join(",", keys %{ $self->{atom} });

  $sql = qq!
select
  e.atom_id
from
  $schema.pw_edge e
where
  e.atom_id in ($atomlist)
!;

  $stm = $db->prepare($sql);
  if(not $stm) {
    $self->AddError("$sql");
    $self->AddError("$DBI::errstr");
    $self->AddError("prepare call failed");
    return;
  }
  if(!$stm->execute()) {
    $self->AddError("$sql");
    $self->AddError("$DBI::errstr");
    $self->AddError("execute call failed");
    return;
  }
  while (($atom_id) = $stm->fetchrow_array()) {
    $self->{atom}{$atom_id} = 2;
  }
  for $atom_id (keys %{ $self->{atom} }) {
    if ($self->{atom}{$atom_id} < 2) {
      $self->AddError("no such atom_id $atom_id");
      delete $self->{atom}{$atom_id}
    }
  }
}

######################################################################
sub CheckEvidence {
  my ($self, $db) = @_;

  if (! defined $self->{evidence_code}) {
    return;
  }

  if (defined $self->{evidence_code}{"*"}) {
    undef %{ $self->{evidence_code} };
    $self->{evidence_code}{"*"} = 1;
    return;
  }

  my ($sql, $stm);
  my ($evidence_code);
  my $schema = $self->{db}{schema};
  my $codelist = join("','", keys %{ $self->{evidence_code} });

  $sql = qq!
select
  e.evidence_code
from
  $schema.pw_evidence_descr e
where
  e.evidence_code in ('$codelist')
!;

  $stm = $db->prepare($sql);
  if(not $stm) {
    $self->AddError("$sql");
    $self->AddError("$DBI::errstr");
    $self->AddError("prepare call failed");
    return;
  }
  if(!$stm->execute()) {
    $self->AddError("$sql");
    $self->AddError("$DBI::errstr");
    $self->AddError("execute call failed");
    return;
  }
  while (($evidence_code) = $stm->fetchrow_array()) {
    $self->{evidence_code}{$evidence_code} = 2;
  }
  for $evidence_code (keys %{ $self->{evidence_code} }) {
    if ($self->{evidence_code}{$evidence_code} < 2) {
      $self->AddError("no such evidence code $evidence_code");
      delete $self->{evidence_code}{$evidence_code}
    }
  }
}

######################################################################
sub CheckFamily {
  my ($self, $db) = @_;

  if (! defined $self->{include_family_uses}) {
    return;
  }

  my $schema = $self->{db}{schema};
  my @temp = keys %{ $self->{mol_id} };
  my ($sql, $stm, $i, $list);
  my ($member_id, $family_id);

  for($i = 0; $i < @temp; $i += ORACLE_LIST_LIMIT) {
    if(($i + ORACLE_LIST_LIMIT - 1) < @temp) {
      $list = join(",", @temp[$i..$i+ORACLE_LIST_LIMIT-1]);
    } else {
      $list = join(",", @temp[$i..@temp-1]);
    }

    $sql = qq!
select
  f.family_mol_id,
  f.member_mol_id
from
  $schema.pw_family_member f
where
      f.family_mol_id in ($list) 
  or  f.member_mol_id in ($list)
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
    while (($family_id, $member_id) = $stm->fetchrow_array()) {
      $self->{mol_id}{$family_id} = 1;
      $self->{mol_id}{$member_id} = 1;
    }
  }
}

######################################################################
sub CheckMacroprocess {
  my ($self) = @_;

  my $lv = $self->{lv};

  my (%all_mp_names, @all_mp_names, %go2lvid);

  for my $lvid ( @{ $lv->AllValuesInFamily("process-type", "macroprocess"),
      $lv->AllValuesInFamily("function", "function") }) {
    my ($dummy, $name) = $lv->LabelValueToString($lvid);
    $all_mp_names{lc($name)} = $lvid;
    push @all_mp_names, lc($name);
    for my $gid (@{ $lv->GOTermsFor($lvid) }) {
      $gid =~ s/GO:?\s*//;
      $go2lvid{$gid}{$lvid} = 1;
    }
  }
  for my $term (keys %{ $self->{macroprocess_name} }) {
    my $term1 = lc($term);
    my $hit;
    for my $mpn (@all_mp_names) {
      if (index($mpn, $term1) > -1) {
        $self->{macroprocess_id}{$all_mp_names{$mpn}} = 1;
        $self->{macroprocessmap}{$term}{$mpn} = 1;
        $hit++
      }
    }
    $term1 =~ s/go:?\s*//;
    if ($term1=~ /^\d{7}$/) {
      if (defined $go2lvid{$term1}) {
        for my $lvid (keys %{ $go2lvid{$term1} }) {
          my ($dummy, $go_name) = $lv->LabelValueToString($lvid);
          $self->{macroprocess_id}{$lvid} = 1;
          $self->{macroprocessmap}{"GO:$term"}{$go_name} = 1;
          $hit++
        }
      }
    }
    if (! $hit) {
      $self->{no_macroprocessmap}{$term} = 1;
    }
  }
}

######################################################################
sub CheckParams {
  my ($self, $db) = @_;

  my $lv = $self->{lv};
  my ($sql, $stm);
  my ($pathway_id, $pathway_ext_id, $pathway_name);

  ##
  ## DB params
  ##

  if ( (! defined $self->{db}) ||
       (! defined $self->{db}{user}) ||
       (! defined $self->{db}{pass}) ||
       (! defined $self->{db}{inst}) ||
       (! defined $self->{db}{schema}) ) {
    $self->AddError("Must specify db params [user, pass, inst, schema]");
    return;
  }
  if (! $self->{db}{user}) {
    $self->AddError("DB user not specified");
    return;
  }
  if (! $self->{db}{pass}) {
    $self->AddError("DB pass not specified");
    return;
  }
  if (! $self->{db}{inst}) {
    $self->AddError("DB instance not specified");
    return;
  }
  if (! $self->{db}{schema}) {
    $self->AddError("DB schema not specified");
    return;
  }

  my $schema = $self->{db}{schema};

  ##
  ## Get and prune mol specs
  ##

  CheckMol($self, $db);

  ##
  ## Include family members and family parents (one level only, for now)
  ## (and it won't apply to pruning or connecting mols
  ##

  CheckFamily($self, $db);


  ##
  ## Atom specs
  ##


  CheckAtom($self, $db);
  CheckEvidence($self, $db);

  ##
  ## Macroprocess names
  ##

  CheckMacroprocess($self);

  ##
  ## Pathway specs
  ##

  if (defined $self->{pathway}) {

    if (keys %{ $self->{pathway}{id} } > 0) {
      if (defined $self->{pathway}{id}{"*"}) {
        $pathway_id = "*";
      } else {
        $pathway_id = [ keys %{ $self->{pathway}{id} } ];
      }
    }
    if (keys %{ $self->{pathway}{ext_id} } > 0) {
      $pathway_ext_id = [ keys %{ $self->{pathway}{ext_id} } ];
    }
    if (keys %{ $self->{pathway}{name} } > 0) {
      $pathway_name = [ keys %{ $self->{pathway}{name} } ];
    }

    delete $self->{pathway};
    if ($pathway_id) {
      LookUpPathwayInfo($self, $db, $schema, 
          $pathway_id, undef, undef);
    }
    if ($pathway_ext_id) {
      LookUpPathwayInfo($self, $db, $schema, 
          undef, $pathway_ext_id, undef);
    }
    if ($pathway_name) {
      LookUpPathwayInfo($self, $db, $schema, 
          undef, undef, $pathway_name);
    }
    for my $p (@{ $pathway_name }, @{ $pathway_ext_id }) {
      if (! defined $self->{pathwaymap}{$p}) {
         $self->{no_pathwaymap}{$p} = 1;
      }
    }

  }

  ##
  ## Degree
  ##

  if (defined $self->{degree}) {
    if (defined $self->{degree}{back} && defined $self->{degree}{forward}) {
      if ($self->{degree}{back} == $self->{degree}{forward}) {
        $self->{degree}{both} = $self->{degree}{back};
      } else {
        $self->AddError("degree back must equal degree forward; " .
            "using lower number");
        if ($self->{degree}{back} < $self->{degree}{forward}) {
          $self->{degree}{both} = $self->{degree}{back};
          $self->{degree}{forward} = $self->{degree}{back};
        } else {
          $self->{degree}{both} = $self->{degree}{forward};
          $self->{degree}{back} = $self->{degree}{forward};
        }
      }
    }
  }

##
## Private data: check activity and location values
## Called AFTER CheckMol so that any LL ids will be caught in mousetrap
##

  if (defined $self->{privateactivity}) {
    for my $a (keys %{ $self->{privateactivity} }) {
      my $lvid = $lv->StringToLabelValue("activity-state", $a);
      if (! $lvid) {
        $self->AddError("illegal activity value: $a");
      }
    }
  }
  if (defined $self->{privatelocation}) {
    for my $loc (keys %{ $self->{privatelocation} }) {
      my $lvid = $lv->StringToLabelValue("location", $loc);
      if (! $lvid) {
        $self->AddError("illegal location value: $loc");
      }
    }
  }
  if (defined $self->{privateedges} ) {
    for my $pe (@{ $self->{privateedges} }) {
      my ($private_interaction_id,
          $edge_type,
          $mol_name,
          $mol_id,
          $activity,
          $location) = split(/\t/, $pe);
      if ($mol_id =~ /^\d+$/ && defined $self->{mousetrap}{$mol_id}) {
        for (@{ $self->{mousetrap}{$mol_id} }) {
          $mol_id = $_;
          last;          ## just take the first one
        }
      }
      $pe = join("\t", $private_interaction_id,
          $edge_type,
          $mol_name,
          $mol_id,
          $activity,
          $location);
    }
  }

##
## Must do this AFTER atoms and normal pathways have been processed
##
  my $dummy_pid;
  if (defined $self->{atom}) {
    $dummy_pid--;
    $self->{pathway}{$dummy_pid}{name} = "[synthetic pathway defined by atom list: " .
      join(",", keys %{ $self->{atom} }) . "]";
    $self->{pathway}{$dummy_pid}{ext_id} = "";
  }
  if (defined $self->{mol_id}) {
    $dummy_pid--;
    $self->{pathway}{$dummy_pid}{name} = "[synthetic pathway defined by 'get' molecule list: " .
      join(",", keys %{ $self->{mol_id} }) . "]";
    $self->{pathway}{$dummy_pid}{ext_id} = "";
  }
  if (defined $self->{prune_mol_id}) {
    $dummy_pid--;
    $self->{pathway}{$dummy_pid}{name} = "[synthetic pathway defined by 'prune' molecule list: " .
      join(",", keys %{ $self->{prune_mol_id} }) . "]";
    $self->{pathway}{$dummy_pid}{ext_id} = "";
  }
  if (defined $self->{degree}) {
    $dummy_pid--;
    my @dg;
    if (defined $self->{degree}{back}) {
      push @dg, "back=" . $self->{degree}{back};
    }
    if (defined $self->{degree}{forward}) {
      push @dg, "forward=" . $self->{degree}{forward};
    }
    $self->{pathway}{$dummy_pid}{name} = "[synthetic pathway defined by degree " .
      join(", ", @dg) . "]";
    $self->{pathway}{$dummy_pid}{ext_id} = "";
  }

  return 1;
}

######################################################################
sub LookUpPathwayInfo {
  my ($self, $db, $schema, $pathway_id, $pathway_ext_id, $pathway_name) = @_;

  my ($sql, $stm);
  my ($pathway_source_id, $pathway_organism);
  my ($is_subnet);

  my ($p);
  for ($pathway_id, $pathway_ext_id, $pathway_name) {
    if ($_ ne "") {
      $p++;
    }
  }
  if ($p == 0) {
    $self->AddError("LookUpPathwayInfo: failed to specify pathway parameter");
    return;
  }

  my $where_clause;
  if ($pathway_id) {
    if ($pathway_id eq "*") {
      $where_clause = "";
    } else {
      $where_clause = "p.pathway_id in (" .
        join(",", @{ $pathway_id }) . ")";
    }
  } elsif ($pathway_ext_id) {
    for my $ext_id (@{ $pathway_ext_id }) {
      $ext_id = lc($ext_id);
    }
    $where_clause = "(lower(p.ext_pathway_id) like '%" .
        join("%' or lower(p.ext_pathway_id) like '%", @{ $pathway_ext_id }) .
        "%')";
#    $where_clause = "lower(p.ext_pathway_id) in ('" .
#        join("','", @{ $pathway_ext_id }) . "')";        

  } elsif ($pathway_name) {
    for my $name (@{ $pathway_name }) {
      $name = lc($name);
    }
    $where_clause = "(lower(p.pathway_name) like '%" .
        join("%' or lower(p.pathway_name) like '%", @{ $pathway_name }) .
        "%')";        
  }

  if (defined $self->{source_id}) {
    if ($where_clause) {
      $where_clause .= " and ";
    }
    $where_clause .= "p.pathway_source_id in (" .
        join(",", keys %{ $self->{source_id} }) . ")";
  }
  $sql = "select p.pathway_id, lower(p.ext_pathway_id), lower(p.pathway_name), " .
    "lower(p.pathway_source_id), p.organism, p.subnet " .
    "from $schema.pw_pathway p " .
    ($where_clause ? "where $where_clause" : "");

  $stm = $db->prepare($sql);
  if(not $stm) {
    $self->AddError("$sql");
    $self->AddError("$DBI::errstr");
    $self->AddError("prepare call failed");
    return;
  }
  if(!$stm->execute()) {
    $self->AddError("$sql");
    $self->AddError("$DBI::errstr");
    $self->AddError("execute call failed");
    return;
  }
  my ($p_id, $p_ext_id, $p_name);
  while (($p_id, $p_ext_id, $p_name, $pathway_source_id,
      $pathway_organism, $is_subnet)
      = $stm->fetchrow_array()) {
    $self->{pathway}{$p_id}{ext_id} = $p_ext_id;
    $self->{pathway}{$p_id}{name}   = $p_name;
    $self->{pathway}{$p_id}{src_id} = $pathway_source_id;
    $self->{pathway}{$p_id}{org}    = $pathway_organism;
    $self->{pathway}{$p_id}{is_subnet} = $is_subnet eq "Y" ? 1 : 0;
  }
  for my $p (@{ $pathway_name }, @{ $pathway_ext_id }) {
    for my $p1 (keys %{ $self->{pathway} }) {
      if (index(lc($self->{pathway}{$p1}{name}), $p)   > -1) {
        $self->{pathwaymap}{$p}{$self->{pathway}{$p1}{name}} = 1;
      } elsif (index(lc($self->{pathway}{$p1}{ext_id}), $p) > -1) {
        $self->{pathwaymap}{$p}{$self->{pathway}{$p1}{ext_id}} = 1;
      }
    }
  }

}

######################################################################
sub ListParams {
  my ($self) = @_;

  return $self->{input_lines};
}

######################################################################
sub ReadOneParam {
  my ($self, $cmd, @vars) = @_;

    if ($cmd ne "db_pass" && $cmd ne "db_user") {
      push @{ $self->{input_lines} }, join("\t", $cmd, @vars);
    }

    if ($cmd eq PRINT) {
      DoPrint($self, @vars);
      next;
    }
    if ($cmd eq SOURCE_ID) {
      DoSourceId($self, @vars);
      next;
    }
    if ($cmd eq EVIDENCE_CODE) {
      DoEvidenceCode($self, "get", @vars);
      next;
    }
    if ($cmd eq MACRO_PROCESS) {
      DoMacroProcess($self, @vars);
      next;
    }
    if ($cmd eq SHOW) {
      DoShow($self, @vars);
      next;
    }
    if ($cmd eq SKIP_MOL_ID) {
      DoSkip($self, "mol_id", @vars);
      next;
    }
    if ($cmd eq SKIP_COMPLEXES) {
      DoSkip($self, "complexes");
      next;
    }
    if ($cmd eq COLLAPSE) {
      DoCollapse($self, @vars);
      next;
    }
    if ($cmd eq COLOR) {
      DoColor($self, @vars);
      next;
    }
    if ($cmd eq VALUE) {
      DoValue($self, @vars);
      next;
    }
    if ($cmd eq CHANNEL) {
      DoChannel($self, @vars);
      next;
    }
    if ($cmd eq EDGE) {
      DoEdge($self, @vars);
      next;
    }
    if ($cmd eq DEGREE_BACK) {
      DoDegree($self, "back", @vars);
      next;
    }
    if ($cmd eq DEGREE_FORWARD) {
      DoDegree($self, "forward", @vars);
      next;
    }
    if ($cmd eq PRUNE_MOL_ID) {
      DoMol($self, "prune", "id", @vars);
      next;
    }
    if ($cmd eq PRUNE_MOL_EXT_ID) {
      DoMol($self, "prune", "ext_id", @vars);
      next;
    }
    if ($cmd eq PRUNE_MOL_NAME) {
      DoMol($self, "prune", "name", @vars);
      next;
    }
    if ($cmd eq CONNECT_MOLECULE) {
      DoMol($self, "connect", "general", @vars);
      next;
    }
    if ($cmd eq CONNECT_MOL_ID) {
      DoMol($self, "connect", "id", @vars);
      next;
    }
    if ($cmd eq CONNECT_MOL_EXT_ID) {
      DoMol($self, "connect", "ext_id", @vars);
      next;
    }
    if ($cmd eq CONNECT_MOL_NAME) {
      DoMol($self, "connect", "name", @vars);
      next;
    }
    if ($cmd eq MOLECULE) {
      DoMol($self, "get", "general", @vars);
      next;
    }
    if ($cmd eq MOL_ID) {
      DoMol($self, "get", "id", @vars);
      next;
    }
    if ($cmd eq MOL_EXT_ID) {
      DoMol($self, "get", "ext_id", @vars);
      next;
    }
    if ($cmd eq MOL_NAME) {
      DoMol($self, "get", "name", @vars);
      next;
    }
    if ($cmd eq INCLUDE) {
      DoInclude($self, @vars);
      next;
    }
    if ($cmd eq PATHWAY) {
      DoPathway($self, "general", @vars);
      next;
    }
    if ($cmd eq PATHWAY_NAME) {
      DoPathway($self, "name", @vars);
      next;
    }
    if ($cmd eq PATHWAY_ID) {
      DoPathway($self, "id", @vars);
      next;
    }
    if ($cmd eq PATHWAY_EXT_ID) {
      DoPathway($self, "extid", @vars);
      next;
    }
    if ($cmd eq ATOM_ID) {
      DoAtomId($self, "get", @vars);
      next;
    }
    if ($cmd eq PRUNE_ATOM_ID) {
      DoAtomId($self, "prune", @vars);
      next;
    }
    if ($cmd eq DB_USER) {
      DoDbUser($self, @vars);
      next;
    }
    if ($cmd eq DB_PASS) {
      DoDbPass($self, @vars);
      next;
    }
    if ($cmd eq DB_INST) {
      DoDbInst($self, @vars);
      next;
    }
    if ($cmd eq DB_SCHEMA) {
      DoDbSchema($self, @vars);
      next;
    }
    if ($cmd eq SIM_CYCLE) {
      DoSim($self, "cycle", @vars);
      next;
    }
    if ($cmd eq SIM_OUTPUT) {
      DoSim($self, "output", @vars);
      next;
    }
    if ($cmd eq SIM_MOL_ID) {
      DoSim($self, "mol_id", @vars);
      next;
    }
    if ($cmd eq SIM_METHOD) {
      DoSim($self, "method", @vars);
      next;
    }
    if ($cmd eq SIM_SIMPLE) {
      DoSim($self, "simple", @vars);
      next;
    }
    if ($cmd eq SIM_COMPETE) {
      DoSim($self, "compete", @vars);
      next;
    }
    $self->AddError("Unrecognized command $cmd");
    return;
}

######################################################################
sub ReadParamString {
  my ($self, $pstring) = @_;

  my ($cmd, @vars);

  for (split /\n/, $pstring) {
    if (/^\s*#/) {
      next;
    }
    if (/^\s*$/) {
      next;
    }
    ($cmd, @vars) = split /\t/;
    $cmd = lc(Trim($cmd));
    if ( !$self->ReadOneParam($cmd, @vars) ) {
      return 0;
    }
  }
  return 1;
}

######################################################################
sub ReadParams {
  my ($self, $f) = @_;

  my ($cmd, @vars);

  if (! open(PARM, $f) ) {
    $self->AddError("Cannot open parameter file $f");
    return 0;
  }
  while (<PARM>) {
    chop;
    if (/^\s*#/) {
      next;
    }
    if (/^\s*$/) {
      next;
    }
    ($cmd, @vars) = split /\t/;
    $cmd = lc(Trim($cmd));
    $self->ReadOneParam($cmd, @vars);
  }
  close PARM;
  if ($self->{errors}) {
    return 0;
  } else {
    return 1;
  }
}

######################################################################
1;
######################################################################
