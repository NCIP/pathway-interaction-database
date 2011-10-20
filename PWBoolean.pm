#!/usr/local/bin/perl
package PWBoolean;
require Exporter;
@ISA = qw(Exporter);
#@EXPORT = qw(
#);

use strict;
use Pathway;
use PWLabel;

######################################################################
# 3-value logic interpretation of network.
# see subs or_3val, and_3val,not_3val
# Value of atom is and_3val of its inputs.
# Value of molecule is or_3val of the atom outputs that feed the molecule.
######################################################################

use constant NO_VAL    => "-1";
use constant NO_TYPE   => "-1";

use constant AGENT_TYPE     => "1";
use constant INHIB_TYPE     => "2";
use constant INPUT_TYPE     => "3";
use constant OUTPUT_TYPE    => "4";

##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
## "output" is an outgoing edge from an atom to a mol, but the outputs
## of a mol are all outgoing edges from all atoms to this mol
##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

## inputs: atomid -> seq of molinstid
## outputs: molinstid -> seq of atomid
## inputtype: atomid, seqid -> AGENT_TYPE | INHIB_TYPE | INPUT_TYPE
## atoms: seq of atomid
## mols: seq of molinstid
## molinstdef: molinstid -> molinststring 
## atomval: atomid -> {NO_VAL,0,1}
## molval:  molinstid -> {NO_VAL,0,1}
## molfanout: molinstid -> Seq Of atomid

######################################################################
sub new {
  my ($self, $pw, $lv, $default_val, $trace, $outfh) = @_;
  my $x = {};
  $x->{lv} = $lv;
  $x->{pw} = $pw;
  if ($trace) {
    $x->{trace} = 1;
  } else {
    $x->{trace} = 0;
  }
  $x->{output_fh} = $outfh;
  $x->{default_val} = $default_val;

  return bless $x;
}

my %location_code = (
  "transmembrane"         => "m",
  "cytoplasm"             => "cy",
  "nucleus"               => "n",
  "endoplasmic reticulum" => "er",
  "extracellular region"         => "ex",
  "calcium store"         => "cs",
  "primary endosome"      => "pe",
  "endosome"              => "e",
  "early endosome"        => "ee",
  "late endosome"         => "le",
  "recycling endosome"    => "re",
  "endosome transmembrane" => "et",
  "golgi"                 => "g",
  "lysosome"              => "l",
  "mitochondria"          => "mi",
  "vesicle"               => "v",
  "intracellular"         => "i",
  
  ""                      => ""
);

######################################################################
sub Lines {
  my ($self) = @_;
  return $self->{lines};
}

sub PrLine {
  my ($self, $x) = @_;
  my $fh = $self->{output_fh};
  if ($fh) {
    print $fh "$x\n";
  }
  push @{ $self->{lines} }, $x;
}

######################################################################
sub LocationDiacritic {
  my ($self, $label_list) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};
  my ($lvid, $label_name, $label_value);
  my $diacritic = "";
  my $location;

  for my $lvid (@{ $label_list }) {
    ($label_name, $label_value) = $lv->LabelValueToString($lvid);
    if ($label_name eq "location") {
      if (defined $location_code{$label_value}) {
        $location = $location_code{$label_value};
        if ($location ne "") {
          $diacritic .= "[$location]"
        }
      } else {
        $diacritic .= "[?]";
      }
    }
  }
  return $diacritic;
}

######################################################################
sub MolActivityDiacritic {
  my ($self, $label_list) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};
  my ($lvid, $label_name, $label_value);
  my $diacritic = "";

  for my $lvid (@{ $label_list }) {
    if ($lv->IsA($lvid, $lv->StringToLabelValue("activity-state", "active"))) {
      ($label_name, $label_value) = $lv->LabelValueToString($lvid);
      $diacritic .= "+";
      if ($label_value =~ /active(\d+)/) {
        $diacritic .= $1;
      }
    } elsif ($lv->IsA($lvid, $lv->StringToLabelValue("activity-state", "inactive"))) {
      $diacritic = "-"
    }
  }
  return $diacritic;
}

######################################################################
sub MolIsComplex {
  my ($self, $mol_id) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};

  my ($label_name, $label_value);

  ($label_name, $label_value) = $lv->LabelValueToString($pw->MolType($mol_id));
  if ($label_value eq "complex") {
    return 1;
  } else {
    return 0;
  }
}

######################################################################
sub GetMolName {
  my ($self, $molinst) = @_;

  my $pw = $self->{pw};

  my $molinststring = $self->MolInstDef($molinst);
  my ($cx_molname, $cx_diacritics, $molname);
  my ($molid, $labels) =
      Pathway::DecodeMolInstString($self->MolInstDef($molinst));

  if ($self->MolIsComplex($molid)) {
    $cx_molname = "cx_$molid";
    $cx_diacritics .= $self->MolActivityDiacritic($labels);
    $cx_diacritics .= $self->LocationDiacritic($labels);
    $molname = "<$cx_molname>$cx_diacritics";
  } else {
    $molname = $pw->PickMolName($molid);
    $molname .= $self->MolActivityDiacritic($labels);
    $molname .= $self->LocationDiacritic($labels);
  }

  return $molname;
}

######################################################################
sub EdgeType {
  my ($self, $atom, $edge) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};

  my $edgetype = $pw->EdgeType($atom, $edge);
  if ($lv->IsA($edgetype,
      $lv->StringToLabelValue("edge-type", "agent"))) {
    return AGENT_TYPE;
  } elsif ($lv->IsA($edgetype,
      $lv->StringToLabelValue("edge-type", "inhibitor"))) {
    return INHIB_TYPE;
  } elsif ($lv->IsA($edgetype,
      $lv->StringToLabelValue("edge-type", "input"))) {
    return INPUT_TYPE;
  } elsif ($lv->IsA($edgetype,
      $lv->StringToLabelValue("edge-type", "output"))) {
    return OUTPUT_TYPE;
  } else {
    print STDERR "unrecognized edge type $edgetype\n";
    return NO_TYPE;
  }
}

######################################################################
sub and_3val {
  my (@vector) = @_;
  my $default = 1;
  for my $v (@vector) {
    if ($v eq "0" || $v eq "") {
      return 0;
    } elsif ($v eq NO_VAL) {
      $default = NO_VAL;
    }
  }
  return $default;
}

######################################################################
sub or_3val {
  my (@vector) = @_;
  my $default = 0;
  for my $v (@vector) {
    if ($v eq "1") {
      return 1;
    } elsif ($v eq NO_VAL) {
      $default = NO_VAL;
    }
  }
  return $default;
}

######################################################################
sub not_3val {
  my ($v) = @_;

  if ($v eq "1") {
    return 0;
  } elsif ($v eq "0") {
    return 1;
  } else {
    return NO_VAL;
  }
}

######################################################################
sub Inputs {
  my ($self, $atom) = @_;
  if (defined $self->{inputs}{$atom}) {
    return $self->{inputs}{$atom};
  } else {
    return [];
  }
}

######################################################################
sub Outputs {
  my ($self, $mol) = @_;
  if (defined $self->{outputs}{$mol}) {
    return $self->{outputs}{$mol};
  } else {
    return [];
  }
}

######################################################################
sub Ninputs {
  my ($self, $atom) = @_;
  if (defined $self->{inputs}{$atom}) {
    return scalar(@{ $self->{inputs}{$atom} });
  } else {
    return 0;
  }
}

######################################################################
sub Noutputs {
  my ($self, $mol) = @_;
  if (defined $self->{outputs}{$mol}) {
    return scalar(@{ $self->{outputs}{$mol} });
  } else {
    return 0;
  }
}

######################################################################
sub InputType {
  my ($self, $atom, $seq) = @_;
  if (defined $self->{inputtype}{$atom}) {
    if (defined $self->{inputtype}{$atom}[$seq]) {
      return $self->{inputtype}{$atom}[$seq];
    }
  }
  return NO_TYPE;
}

######################################################################
sub Atoms {
  my ($self) = @_;
  return $self->{atoms};
}

######################################################################
sub Mols {
  my ($self) = @_;
  return $self->{mols};
}

######################################################################
sub AtomVal {
  my ($self, $atom) = @_;
  if (defined $self->{atomval}{$atom}) {
    return $self->{atomval}{$atom};
  } else {
    return NO_VAL;
  }
}

######################################################################
sub SetAtomVal {
  my ($self, $atom, $val) = @_;
  $self->{atomval}{$atom} = $val;
}

######################################################################
sub MolVal {
  my ($self, $mol, $dummy_mol_id) = @_;
  if (defined $self->{molval}{$mol}) {
    return $self->{molval}{$mol};
  } else {
    return NO_VAL;
  }
}

######################################################################
sub SetMolVal {
  my ($self, $mol, $val) = @_;
  $self->{molval}{$mol} = $val;
}


######################################################################
sub MolInstDef {
  my ($self, $mol) = @_;
  if (defined $self->{molinstdef}{$mol}) {
    return $self->{molinstdef}{$mol};
  } else {
    return "";
  }
}

######################################################################
sub InputVal {
  my ($self, $atom, $i) = @_;
  if (defined $self->{inputs}{$atom}) {
    if (defined $self->{inputs}{$atom}[$i]) {
      my $m = $self->{inputs}{$atom}[$i];
      if (defined $self->{molval}{$m}) {
        return $self->{molval}{$m};
      }
    }
  }
  return NO_VAL;
}

######################################################################
sub OutputVal {
  my ($self, $mol, $i) = @_;
  if (defined $self->{outputs}{$mol}) {
    if (defined $self->{outputs}{$mol}[$i]) {
      my $a = $self->{outputs}{$mol}[$i];
      if (defined $self->{atomval}{$a}) {
        return $self->{atomval}{$a};
      }
    }
  }
  return NO_VAL;
}

######################################################################
sub EvalAtomsPartOne {
  my ($self) = @_;
  my ($v);
  my $change = 0;
  for $a (@{ $self->Atoms }) {
    $v = $self->EvalOneAtom($a);
    if ($v ne $self->AtomVal($a)) {
      $self->SetAtomVal($a, $v);
      $change++;
    }
  }
  return $change;
}

######################################################################
sub EvalAtoms {
  my ($self) = @_;

  my $change;
  my $evalatomcycle;
  while (1) {
    $evalatomcycle++;
    $change = $self->EvalAtomsPartOne;
    if ($self->{trace}) {
      $self->PrLine("evalatoms cycle = $evalatomcycle, changed = $change");
    }
    if (! $change) {
      last;
    }
  }
}

######################################################################
sub EvalOneAtom {
  my ($self, $atom) = @_;

  my ($i, $v_sum, $mv, $v, $x);

  my ($min_input_val, $max_input_val);
  my $inputs = $self->Inputs($atom);
  my $n = $self->Ninputs($atom);

  my (@input_vector);
  for ($i = 0; $i < $n; $i++) {
    $mv = $self->InputVal($atom, $i);
    if ($self->InputType($atom, $i) == INHIB_TYPE) {
      $mv = not_3val($mv);
    }
    push @input_vector, $mv;
  }
  return and_3val(@input_vector);
}

######################################################################
sub EvalOneMol {
  my ($self, $mol) = @_;

  my ($i, $n, $v_sum, $w_sum, $mv, $v);

  my $n = $self->Noutputs($mol);
  if ($n == 0) {
    return $self->MolVal($mol);
  }
  my (@output_vector);
  for ($i = 0; $i < $n; $i++) {
    $mv = $self->OutputVal($mol, $i);
    push @output_vector, $mv;
  }
  return or_3val(@output_vector);
}

######################################################################
sub Propagate {
  my ($self) = @_;

  my ($m, $v0, $v1, $change);
  $change = 0;

  for $m (@{ $self->Mols }) {
    $v1 = $self->EvalOneMol($m);
    $v0 = $self->MolVal($m);
    if ($v1 ne $v0) {
      $change = 1;
    }
    $self->SetMolVal($m, $v1);
  }
  return $change;
}

######################################################################
sub DumpState {
  my ($self, $cycle) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};

  my ($atom, $inputs, $str,
      $av, $i, $n, $w, $molinst, $molid, $labels);

  for $molinst (@{ $self->Mols }) {
      $str = $pw->MolInstIdToString($molinst);
      ($molid) = split(":", $str);
      $str =~ s/[:,]/_/g;
      $str =~ s/_$//;
      $str = "M_$str";
      my $molname = $self->GetMolName($molinst);
      if ($self->{trace}) {
        $self->PrLine(sprintf "$cycle\tmol\t%s\t%s\tval=%d",
            $str,
            $molname,
            $self->MolVal($molinst));
      }
  }
  for my $atom (@{ $self->Atoms() }) {
    $av = $self->AtomVal($atom);
    if ($self->{trace}) {
      $self->PrLine(sprintf "$cycle\tatom\tA_$atom\t%d", $av);
    }
    $inputs = $self->Inputs($atom);
    $n = scalar(@{ $inputs });
    for ($i = 0; $i < $n; $i++) {
      $molinst = $$inputs[$i];
      $str = $pw->MolInstIdToString($molinst);
      ($molid) = split(":", $str);
      $str =~ s/[:,]/_/g;
      $str =~ s/_$//;
      $str = "M_$str";
      if ($self->{trace}) {
        $self->PrLine(sprintf "$cycle\tedge\t%s_input_%s\t%s\tval=%d",
            $atom,
            $str,
            $pw->PickMolName($molid),
            $self->MolVal($molinst))
      }
    }

  }

}

######################################################################
sub Execute {
  my ($self, $ncycle) = @_;

  for (my $i = 1; $i <= $ncycle; $i++) {
    $self->EvalAtoms;
    $self->DumpState($i);
    if (! $self->Propagate) {
      if ($self->{trace}) {
        $self->PrLine("Model quiescent");
      }
      last;
    }
  }
}

######################################################################
sub Force {
  my ($self, $v) = @_;
  if ($v > 0) {
    return 1;
  } elsif ($v < 0) {
    return NO_VAL;
  } else {
    return 0
  }
}

######################################################################
sub InitializeMolValues {
  my ($self, $mol_values) = @_;

  my $lv = $self->{lv};
  my $pw = $self->{pw};

if (0) {
  my ($mol_id, $molinstid, $molinststr, %tmp);

  for $mol_id (keys %{ $mol_values }) {

###!!! for now: "\t" denotes null-spec for activity, null-spec for location
    $tmp{$mol_id} = $$mol_values{$mol_id}{"\t"};
  }
  while (($molinstid, $molinststr) = each %{ $self->{molinstdef} }) {
    ($mol_id) = split ":", $molinststr;
    if (defined $tmp{$mol_id}) {
      my $v = $self->Force($tmp{$mol_id});
      $self->SetMolVal($molinstid, $v);
    }
  }
}

  my ($mol_id, $molinstid, $molinststr);
  my ($subtype, $activity, $location, $ls_b, @ls_b);
  while (($molinstid, $molinststr) = each %{ $self->{molinstdef} }) {
    ($mol_id, $ls_b) = split ":", $molinststr;
    if (defined $$mol_values{$mol_id}) {
      @ls_b = split ",", $ls_b;
      for $subtype (keys %{ $$mol_values{$mol_id} }) {
        my @ls_a;
        ($activity, $location) = split("\t", $subtype);
        if ($activity) {
          push @ls_a,
              $lv->StringToLabelValue("activity-state", $activity);
        }
        if ($location) {
          push @ls_a,
              $lv->StringToLabelValue("location", $location);
        }
        if (Pathway::LabelSetIsSubType(\@ls_b, \@ls_a)) {
          my $v = $self->Force($$mol_values{$mol_id}{$subtype});
          $self->SetMolVal($molinstid, $$mol_values{$mol_id}{$subtype});
        }
      }
    }
  }


}

######################################################################
sub SetUpSimTopology {
  my ($self) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};
  my $default_val = $self->{default_val};

  my ($atom, $edge, $molinst, $et);

  for $atom (@{ $pw->Atoms }) {
    push @{ $self->{atoms} }, $atom;
    $self->{atomval}{$atom} = $default_val;
    for $edge (@{ $pw->Edges($atom) }) {
      $molinst = $pw->MolInstId($atom, $edge);
      $self->{molinstdef}{$molinst} =
          $pw->NormalMolInstString($atom, $edge);
      $et = $pw->EdgeType($atom, $edge);
      if ($lv->IsA($et,
          $lv->StringToLabelValue("edge-type", "incoming-edge"))) {
        push @{ $self->{inputs}{$atom} }, $molinst;
        push @{ $self->{inputtype}{$atom} },
            $self->EdgeType($atom, $edge);
      } else {
        push @{ $self->{outputs}{$molinst} }, $atom;
        push @{ $self->{inputtype}{$molinst} },
            $self->EdgeType($atom, $edge);
      }

    }
  }
  @{ $self->{mols} } = keys %{ $self->{molinstdef} };
  for $molinst (@{ $self->Mols }) {
    $self->{molval}{$molinst} = $default_val;
  }

  $self->DumpState(0);

#test
#$self->SetMolVal(1, 1.0);

}

######################################################################
1;
######################################################################
