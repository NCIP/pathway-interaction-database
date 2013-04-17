#!/usr/local/bin/perl


# Copyright SRA International
#
# Distributed under the OSI-approved BSD 3-Clause License.
# See http://ncip.github.com/pathway-interaction-database/LICENSE.txt for details.


package PWSim;
require Exporter;
@ISA = qw(Exporter);
#@EXPORT = qw(
#);

use strict;
use Pathway;
use PWLabel;

######################################################################
# The idea...
#
# Trying to simulate deviation from some sort of normal state for
# a given network topology.
#
# At any point in the simulation, the value of a molecule (molecule
# instance) is relative to the normal value of that molecule. Thus
# molecule values cannot be compared directly to each other.
#
# At any point in the simulation, the value of an atomic process
# (atom) is relative to the normal value of that process. Thus atom
# values cannot be compared directly to each other.
#
# Values are best thought of as an ordinal scale (despite the fact
# that we are using real numbers).
#
# Assumption: In the normal state, if a molecule is an input to each of
# N processes, then the normal supply of that molecule is N times as
# great as a molecule that is input to only 1 process. Undoubtedly false,
# but we gotta simplify.
#
# The fanout of a molecule M is the set of all atoms to which M is
# an input. If A is in Fanout(M), then there is an edge (M,A).
# Each such input edge has a weight (range: 0..1).
# The weight of an input edge determines the proportion of the molecule
# value that is being used (consumed) by the given atom.
# If Fanout(M) = A, then the weight of edge (M,Ai), for Ai in A, is
# the value of Ai divided by the sum of the values of all atoms in
# A. 
#
# The fanin of a molecule M is the set of all atoms which have M as an
# output. If Fanin(M) = B, then the weight (range 0..1) of each output
# edge (B,M) determines the relative amount by which each process
# contributes to the value of M. At present, the weights of all
# output edges are fixed at 1.
#
# At the start of the simulation, the values of all molecules and all
# atoms are set to NORMAL, after which the weights of all input edges
# are computed. Then the values of individual molecules or atoms may be
# adjusted to represent deviations from the NORMAL state.
# Subsequent to this initialization, the simulation cycles until quiescence
# or until a maximum number of cycles has been executed.
#
# Each cycle proceeds as follows:
#
#   1. Determine new value for each atom
#
#      Do until no atom value changes more than EPSILON:
#       (a) Determine new value for each atom A as the weighted
#          mean of the values of all molecules Mi such that
#          A is in Fanout(Mi). [Alternative method of computation:
#          new value of A is the minimum of the values of all
#          molecules Mi such that A is in Fanout(Mi).]
#            Weighted mean =
#                sum[Mi * weight(Mi,A) * cardinality(Fanout(Mi))] / N
#                  where A is in Fanout(M1),...,Fanout(Mi),...,Fanout(MN)
#
#      (b) Compute new value for each input weight.
#
#   2. Propagate new atom values to outputs
#
#      The new value of each molecule M is the weighted mean
#      of the new values of each atom in Fanin(M).
#

use constant NO_VAL    => "-999999";
use constant NO_TYPE   => "-999999";
use constant NO_WEIGHT => "-999999";

use constant LO_VAL         => "1";
use constant HI_VAL         => "9";
use constant NORMAL_VAL     => (LO_VAL + HI_VAL)/2;
use constant EPSILON_PARTS  => "1000";
use constant EPSILON        => (HI_VAL - LO_VAL) / EPSILON_PARTS;

use constant NORMAL_WEIGHT  => "1";

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
## outputtype: molinstid, seqid -> OUTPUT_TYPE (degenerate, so far)
## inputweight: atomid, seqid -> int range 0..1
## outputweight: molinstid, seqid -> int range -1:1
## atoms: seq of atomid
## mols: seq of molinstid
## molinstdef: molinstid -> molinststring 
## atomval: atomid -> int range LO_VAL..HI_VAL
## molval:  molinstid -> int range LO_VAL..HI_VAL
## molfanout: molinstid -> Seq Of atomid

######################################################################
sub new {
  my ($self, $pw, $lv, $method, $all_or_nothing, $adjust_input_weights,
      $trace, $outfh) = @_;
  my $x = {};
  $x->{lv} = $lv;
  $x->{pw} = $pw;
  if ($method) {
    if ($method eq "mean") {
      $x->{method} = "MEAN";
    } elsif ($method eq "min") {
      $x->{method} = "MIN";
    } elsif ($method eq "max") {
      $x->{method} = "MAX";
    } else {
      die "illegal simulation method $method";
    }
  } else {
    $x->{method} = "MEAN";   ## default
  }
  if ($trace) {
    $x->{trace} = 1;
  } else {
    $x->{trace} = 0;
  }
  if ($all_or_nothing) {
    $x->{allornothing} = 1;
  } else {
    $x->{allornothing} = 0;
  }
  if ($adjust_input_weights) {
    $x->{adjustinputweights} = 1;
  } else {
    $x->{adjustinputweights} = 0;
  }
  $x->{output_fh} = $outfh;

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
sub Inputs {
  my ($self, $atom) = @_;
  if (defined $self->{inputs}{$atom}) {
    return $self->{inputs}{$atom};
  } else {
    return [];
  }
}

######################################################################
sub MolinstFanout {
  my ($self, $molinst) = @_;
  if (defined $self->{molinstfanout}{$molinst}) {
    return $self->{molinstfanout}{$molinst};
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
sub OutputType {
  my ($self, $mol, $seq) = @_;
  if (defined $self->{outputtype}{$mol}) {
    if (defined $self->{outputtype}{$mol}[$seq]) {
      return $self->{outputtype}{$mol}[$seq];
    }
  }
  return NO_TYPE;
}

######################################################################
sub SetInputWeight {
  my ($self, $atom, $molinst, $w) = @_;

  my ($i);
  my $inputs = $self->Inputs($atom);
  my $n = scalar(@{ $inputs });
  for ($i = 0; $i < $n; $i++) {
    if ($$inputs[$i] == $molinst) {
      $self->{inputweight}{$atom}[$i] = $w;
    }
  }
}

######################################################################
sub InputWeight {
  my ($self, $atom, $seq) = @_;
  if (defined $self->{inputweight}{$atom}) {
    if (defined $self->{inputweight}{$atom}[$seq]) {
      return $self->{inputweight}{$atom}[$seq];
    }
  }
  return NO_WEIGHT;
}

######################################################################
sub OutputWeight {
  my ($self, $mol, $seq) = @_;
  if (defined $self->{outputweight}{$mol}) {
    if (defined $self->{outputweight}{$mol}[$seq]) {
      return $self->{outputweight}{$mol}[$seq];
    }
  }
  return NO_WEIGHT;
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
sub MolHasValue {
  my ($self, $mol, $dummy_mol_id) = @_;
  return $self->MolVal($mol);
}

######################################################################
sub AtomHasValue {
  my ($self, $atom) = @_;
  return $self->AtomVal($atom);
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
    $v = $self->EvalOneAtom($a, $self->{method}, $self->{allornothing});
    if ($v == NO_VAL) {
      next;
    }
    if (abs($v - $self->AtomVal($a)) > EPSILON) {
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
    if ($self->{adjustinputweights}) {
      $self->AdjustInputWeights;
    }
  }

}

######################################################################
sub AdjustInputWeights {
  my ($self) = @_;

  my ($molinst, $atom, $av, $av_sum);

  for $molinst (@{ $self->Mols }) {
    $av_sum = 0;
    for $atom (@{ $self->MolinstFanout($molinst) }) {
      $av = $self->AtomVal($atom);
      if ($av == NO_VAL) {
        next;
      }
      $av_sum += $self->AtomVal($atom);
    }
    if ($av_sum > 0) {
      for $atom (@{ $self->MolinstFanout($molinst) }) {
        $av = $self->AtomVal($atom);
        if ($av == NO_VAL) {
          next;
        }
        $self->SetInputWeight($atom, $molinst, ($av / $av_sum));
      }
    }
  }

}

######################################################################
sub EvalOneAtom {
  my ($self, $atom, $method, $all_or_nothing) = @_;

  my ($i, $v_sum, $mv, $v, $x);

  my ($min_input_val, $max_input_val);
  my $inputs = $self->Inputs($atom);
  my $n = $self->Ninputs($atom);

  for ($i = 0; $i < $n; $i++) {
    $mv = $self->InputVal($atom, $i);
    if ($mv == NO_VAL) {
      return NO_VAL;
    }
    if ($self->InputType($atom, $i) == INHIB_TYPE) {
      $mv = (2 * NORMAL_VAL) - $mv;
    }
    $x = $mv * $self->InputWeight($atom, $i) *
        scalar(@{ $self->MolinstFanout($$inputs[$i]) });
    if ((! defined $min_input_val) || $x < $min_input_val) {
      $min_input_val = $x;
    }
    if ((! defined $max_input_val) || $x > $max_input_val) {
      $max_input_val = $x;
    }
    $v_sum += $x;
  }
  if ($method eq "MEAN") {
    if ($n < 1) {
      return NO_VAL;
    }
    $v = $v_sum / $n;
  } elsif ($method eq "MIN") {
    $v = $min_input_val;
  } elsif ($method eq "MAX") {
    $v = $max_input_val;
  } else {
    return NO_VAL;
  }
  if ($v > HI_VAL) {
    $v = HI_VAL;
  } elsif ($v < LO_VAL) {
    $v = LO_VAL;
  }

  if ($all_or_nothing) {
    if (NORMAL_VAL - $v > EPSILON) {
      $v = LO_VAL;
    } elsif ($v - NORMAL_VAL > EPSILON) {
      $v = HI_VAL;
    }
  }

  return $v;
  
}

######################################################################
sub EvalOneMol {
  my ($self, $mol, $all_or_nothing) = @_;

  my ($i, $n, $v_sum, $w_sum, $mv, $v);

  my $n = $self->Noutputs($mol);
  if ($n == 0) {
    return $self->MolVal($mol);
  }
  for ($i = 0; $i < $n; $i++) {
    $mv = $self->OutputVal($mol, $i);
    if ($mv == NO_VAL) {
      return NO_VAL;
    }
    $v_sum += $mv;
    $w_sum += $self->OutputWeight($mol, $i);
  }
  if ($w_sum > 0) {
    $v = $v_sum / $w_sum;
    if ($v > HI_VAL) {
      $v = HI_VAL;
    } elsif ($v < LO_VAL) {
      $v = LO_VAL;
    }
  } else {
    $v = NO_VAL;
  }

  if ($all_or_nothing) {
    if (NORMAL_VAL - $v > EPSILON) {
      $v = LO_VAL;
    } elsif ($v - NORMAL_VAL > EPSILON) {
      $v = HI_VAL;
    }
  }

  return $v;
  
}

######################################################################
sub Propagate {
  my ($self) = @_;

  my ($m, $v0, $v1, $change);
  $change = 0;

  for $m (@{ $self->Mols }) {
    $v1 = $self->EvalOneMol($m, $self->{allornothing});
    if ($v1 == NO_VAL) {
      next;
    }
    $v0 = $self->MolVal($m);
    if (abs($v1 - $v0) > EPSILON) {
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
  my ($total_mol_level, $total_proc_level,
      $avg_mol_level, $avg_proc_level,
      $n_mol, $n_proc);

  for $molinst (@{ $self->Mols }) {
      $str = $pw->MolInstIdToString($molinst);
      ($molid) = split(":", $str);
      $str =~ s/[:,]/_/g;
      $str =~ s/_$//;
      $str = "M_$str";
      my $molname = $self->GetMolName($molinst);
      if ($self->{trace}) {
        $self->PrLine(sprintf "$cycle\tmol\t%s\t%s\tval=%.2f",
            $str,
            $molname,
            $self->MolVal($molinst));
        $n_mol++;
        $total_mol_level += $self->MolVal($molinst);
      }
  }
  for my $atom (@{ $self->Atoms() }) {
    $av = $self->AtomVal($atom);
    if ($self->{trace}) {
      $self->PrLine(sprintf "$cycle\tatom\tA_$atom\t%.2f", $av);
      $n_proc++;
      $total_proc_level += $av;
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
        $self->PrLine(sprintf "$cycle\tedge\t%s_input_%s\t%s\tval=%.2f\tweight=%.2f",
            $atom,
            $str,
            $pw->PickMolName($molid),
            $self->MolVal($molinst),
            $self->InputWeight($atom, $i));
      }
    }

  }
  if ($self->{trace}) {
    $self->PrLine("Average process level = " .
        sprintf("%.2f", $total_proc_level / $n_proc) .
        " for $n_proc processes");
    $self->PrLine("Average molecule level = " .
        sprintf("%.2f", $total_mol_level / $n_mol) .
        " for $n_mol molecules");
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
sub OldInitializeDeviations {
  my ($self, $mol_deviations) = @_;

  my ($mol_id, $pct, $molinstid, $molinststr, %tmp, $x);

  for $mol_id (keys %{ $mol_deviations }) {

    if (! $mol_id) {
      next;
    }

    $pct = $$mol_deviations{$mol_id};
    if ($pct =~ /^-?\d*.?\d*$/) {
      if ($pct < -1 || $pct > 1) {
        die "illegal percent deviation '$pct' for molecule id $mol_id";
      }
      if ($pct < 0) {
        $x = NORMAL_VAL + (NORMAL_VAL - LO_VAL) * $pct;
      } else {
        $x = NORMAL_VAL + (HI_VAL - NORMAL_VAL) * $pct
      }
      if ($x < LO_VAL) {
        $x = LO_VAL;
      }
      if ($x > HI_VAL) {
        $x = HI_VAL;
      }
##      $tmp{$mol_id} = $x;
##
## Use molinstid rather than mol_id
##
      if (defined $self->{molval}{$mol_id}) {
        $self->SetMolVal($mol_id, $x);
      } else {
        die "illegal molecule instance id specified (simulation)";
      }
    } else {
      die "illegal percent deviation '$pct' for molecule id $mol_id";
    }
  }

##  while (($molinstid, $molinststr) = each %{ $self->{molinstdef} }) {
##    ($mol_id) = split ":", $molinststr;
##    if (defined $tmp{$mol_id}) {
##      $self->SetMolVal($molinstid, $tmp{$mol_id});
##    }
##  }

}

######################################################################
sub InitializeDeviations {
  my ($self, $mol_deviations) = @_;

  my $lv = $self->{lv};
  my $pw = $self->{pw};

if (0) {
  my ($mol_id, $molinstid, $molinststr, %tmp);

  for $mol_id (keys %{ $mol_deviations }) {

###!!! for now: "\t" denotes null-spec for activity, null-spec for location
    $tmp{$mol_id} = $$mol_deviations{$mol_id}{"\t"};
  }
  while (($molinstid, $molinststr) = each %{ $self->{molinstdef} }) {
    ($mol_id) = split ":", $molinststr;
    if (defined $tmp{$mol_id}) {
      $self->SetMolVal($molinstid, $tmp{$mol_id});
    }
  }
}

  my ($mol_id, $molinstid, $molinststr);
  my ($subtype, $activity, $location, $ls_b, @ls_b);
  while (($molinstid, $molinststr) = each %{ $self->{molinstdef} }) {
    ($mol_id, $ls_b) = split ":", $molinststr;
    if (defined $$mol_deviations{$mol_id}) {
      @ls_b = split ",", $ls_b;
      for $subtype (keys %{ $$mol_deviations{$mol_id} }) {
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
          $self->SetMolVal($molinstid, $$mol_deviations{$mol_id}{$subtype});
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

  my ($atom, $edge, $molinst, $et);

  for $atom (@{ $pw->Atoms }) {
    push @{ $self->{atoms} }, $atom;
    $self->{atomval}{$atom} = NORMAL_VAL;
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
        push @{ $self->{molinstfanout}{$molinst} }, $atom;
      } else {
        push @{ $self->{outputs}{$molinst} }, $atom;
        push @{ $self->{inputtype}{$molinst} },
            $self->EdgeType($atom, $edge);
        push @{ $self->{outputweight}{$molinst} }, NORMAL_WEIGHT;
      }

    }
  }
  @{ $self->{mols} } = keys %{ $self->{molinstdef} };
  for $molinst (@{ $self->Mols }) {
    $self->{molval}{$molinst} = NORMAL_VAL;
  }

  $self->AdjustInputWeights;
  $self->DumpState(0);

#test
#$self->SetMolVal(1, 1.0);

}

######################################################################
1;
######################################################################
