#!/usr/local/bin/perl
package SimDOToutput;
require Exporter;
@ISA = qw(Exporter);
#@EXPORT = qw(
#  CreateIsA
#  CreateLabelSort
#);

use strict;
use Pathway;
use PWLabel;

my (%visited);

my %layer_color = (
  1 => "#cc6600",
  2 => "#ffcc33",
  3 => "#ffcc99",
  4 => "#ffffcc",
  5 => "#cccccc"
);

######################################################################
sub new {
  my ($self, $sim, $value2color, $pw, $lv, $outfh) = @_;
  my $x = {};
  $x->{sim} = $sim;
  $x->{lv} = $lv;
  $x->{pw} = $pw;
  $x->{f_Value2Color} = $value2color;
  $x->{output_fh} = $outfh;

  ##
  ## Molecules
  ##

  $x->{shape}{$lv->StringToLabelValue("molecule-type",
      "protein")}  = "ellipse";

  $x->{color}{$lv->StringToLabelValue("molecule-type",
      "molecule-type")} = "black";                  ## default

  $x->{height}{$lv->StringToLabelValue("molecule-type",
      "molecule-type")}   = "";

  $x->{width}{$lv->StringToLabelValue("molecule-type",
      "molecule-type")}    = "";

  $x->{fontsize}{$lv->StringToLabelValue("molecule-type",
      "molecule-type")} = "14";

  $x->{style}{$lv->StringToLabelValue("molecule-type",
      "molecule-type")}    = "filled";

  ##
  ## Processes
  ##

  $x->{shape}{$lv->StringToLabelValue("process-type",
     "process-type")}  = "hexagon";
  $x->{shape}{$lv->StringToLabelValue("process-type",
      "macroprocess")}  = "parallelogram";
  $x->{shape}{$lv->StringToLabelValue("process-type",
      "macroprocess")}  = "parallelogram";
  $x->{shape}{$lv->StringToLabelValue("process-type",
      "reaction")}      = "box";
  $x->{shape}{$lv->StringToLabelValue("process-type",
      "binding")}       = "diamond";
  $x->{shape}{$lv->StringToLabelValue("process-type",
      "translocation")} = "triangle";
##  $x->{shape}{$lv->StringToLabelValue("process-type",
##      "modification")}  = "parallelogram";
  $x->{shape}{$lv->StringToLabelValue("process-type",
      "modification")}  = "diamond";
  $x->{shape}{$lv->StringToLabelValue("process-type",
      "transcription")}  = "doublecircle";

  $x->{color}{$lv->StringToLabelValue("process-type",
      "process-type")} = "black";                   ## default

  $x->{height}{$lv->StringToLabelValue("process-type",
      "process-type")}   = "0.2";

  $x->{width}{$lv->StringToLabelValue("process-type",
      "process-type")}    = "0.2";

  $x->{fontsize}{$lv->StringToLabelValue("process-type",
      "process-type")} = "10";

  $x->{style}{$lv->StringToLabelValue("process-type",
		      "process-type")}    = "filled";
  $x->{style}{$lv->StringToLabelValue("process-type",
      "macroprocess")}   = "solid";

  ##
  ## Edges
  ##

  $x->{color}{$lv->StringToLabelValue("edge-type",
      "edge-type")}    = "black";                    ## default
  $x->{color}{$lv->StringToLabelValue("edge-type",
      "agent")}    = "green";
  $x->{color}{$lv->StringToLabelValue("edge-type",
      "inhibitor")}    = "red";
  $x->{arrowhead}{$lv->StringToLabelValue("edge-type",  ## default
      "edge-type")}    = "normal";
  $x->{arrowhead}{$lv->StringToLabelValue("edge-type",
      "inhibitor")}    = "tee";

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
sub ShowAtomIds {
  my ($self, $off_or_on) = @_;

  $self->{showatomids} = $off_or_on;
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
sub ComplexName {
  my ($self, $cx_mol_id, $lines) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};
  my ($comp_seq_id, $sub_mol_id);
  my ($molname, $label_name, $label_value, $lvid, $location);

  for $comp_seq_id (@{ $pw->Components($cx_mol_id) }) {
    $sub_mol_id = $pw->ComponentMol($cx_mol_id, $comp_seq_id);
    if ($self->MolIsComplex($sub_mol_id)) {
      $self->ComplexName($sub_mol_id, $lines);
    } else {
      $molname = $pw->PickMolName($sub_mol_id);
      $molname .= $self->MolActivityDiacritic(
          $pw->ComponentLabel($cx_mol_id, $comp_seq_id));
      $molname .= $self->LocationDiacritic(
          $pw->ComponentLabel($cx_mol_id, $comp_seq_id));
    }
    push @{ $lines }, $molname;
  }
}

######################################################################
sub DOTAttr {
  my ($self, $attr_name, $attr_value) = @_;

  return "$attr_name=\"" .
      $self->AttrValueOf($attr_name, $attr_value) . "\"";

}

######################################################################
sub AttrValueOf {
  my ($self, $attr, $label_value_id) = @_;
  my $h  = $self->{$attr};
  my $lv = $self->{lv};
  if (defined $h->{$label_value_id}) {
    return $h->{$label_value_id};
  } else {
    my $v = $lv->FindLowestAncestor($label_value_id, [ keys %{ $h } ]);
    if ($v == $label_value_id) {
      ##
      ## nothing defined for this or any parent !!!
      ##
      return "";
    } else {
      $h->{$label_value_id} = $h->{$v};
      return $h->{$v};
    }
  }
}

######################################################################
my $indent_level;
my $blanks =  "                                        ";
use constant INDENT_INCR => 2;

sub Lines {
  my ($self) = @_;
  return $self->{lines};
}

sub PrLine {
  my ($self, $x) = @_;
  my $fh = $self->{output_fh};
  if ($fh) {
    print $fh substr($blanks, 0, $indent_level*INDENT_INCR);
    print $fh "$x\n";
  }
  push @{ $self->{lines} },
      substr($blanks, 0, $indent_level*INDENT_INCR) . $x;
}

sub PrLineNoIndent {
  my ($self, $x) = @_;
  my $fh = $self->{output_fh};
  if ($fh) {
    print $fh "$x\n";
  }
  push @{ $self->{lines} }, $x;
}

sub Indent {
  if ($indent_level*INDENT_INCR >= length($blanks)) {
    print STDERR "attempt to indent beyond limit\n";
  } else {
    $indent_level++;
  }
}

sub Exdent {
  if ($indent_level > 0) {
    $indent_level--;
  } else {
    print STDERR "attempt to exdent < 0\n";
  }
}
######################################################################


##
## There may be a conflict between molinst ids and atom ids,
## and since these are both nodes, we need to map these to
## something else.
##
my (%nodemap, $nextnodeid);

######################################################################
sub MapNodeId {
  my ($self, $what, $id) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};

  if ($what eq "ATOM" && $pw->IsMacroProcess($id)) {
    return $self->MapNodeId("MACROPROCESS", $pw->AtomType($id));
  }

  if ($what eq "MACROPROCESS") {
    ## $id is the label_value_id for the kind of macroprocess
    return "m_$id";
  } elsif ($what eq "ATOM" || $what eq "MOLECULE") {
    if (defined $nodemap{"$what$id"}) {
      return $nodemap{"$what$id"};
    } else {
      $nextnodeid++;
      $nodemap{"$what$id"} = $nextnodeid;
      return $nextnodeid;
    }
  } else {
    print STDERR "MapNodeId: illegal what = $what\n";
  }
}

######################################################################
sub DOTMacroProcesses {
  my ($self) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};

  my ($atomid, $processtype, $cond, %processtype2cond,
      %processtype2atom);
  my ($label, $value);
  my (%missing);
  my $macro_atom_ids = $pw->MacroProcesses();

  for $atomid (keys %{ $macro_atom_ids }) {
    $processtype = $$macro_atom_ids{$atomid};
    $processtype2atom{$processtype}{$atomid} = 1;
    for $cond (@{ $pw->AtomCondition($atomid) }) {
      $processtype2cond{$processtype}{$cond} = 1;
    }
  }

  ##
  ## Do conditions on "regular" atoms; note referenced
  ## macroprocess types that are not "declared" in this model;
  ##
  for $atomid (@{ $pw->Atoms() }) {
    for $cond (@{ $pw->AtomCondition($atomid) }) {
      if (! defined $processtype2atom{$cond}) {
        $missing{$cond} = 1;
      }
      if (defined $$macro_atom_ids{$atomid}) {
        next;     ## skip macroprocesses... do later
      }
      $self->PrLine($self->MapNodeId("MACROPROCESS", $cond) . " " .
          "-> " . 
          $self->MapNodeId("ATOM", $atomid) . " " .
          "[" .
          $self->DOTAttr("color",
                $lv->StringToLabelValue("edge-type", "agent")) .
          "];") ;
    }
  }

  ##
  ## Create nodes for "undeclared" macroprocess types
  ##
  for $processtype (keys %missing) {
    ($label, $value) = $lv->LabelValueToString($processtype);
    my @attrs;
    for my $a ("shape", "height", "width", "fontsize", "style") {
      push @attrs, DOTAttr($self, $a, $processtype);
    }
    push @attrs, "label=\"$value\"";
    $self->PrLine($self->MapNodeId("MACROPROCESS", $processtype) . " " .
        "[" .
        join(", ", @attrs) .
        "];") ;
  }

  ##
  ## Now do "declared" macroprocess types
  ## These may have agents, inhibitors, inputs, outputs and conditions.
  ##
  for $processtype (keys %processtype2atom) {
    ($label, $value) = $lv->LabelValueToString($processtype);
    my @attrs;
    for my $a ("shape", "height", "width", "fontsize", "color", "style") {
      push @attrs, DOTAttr($self, $a, $processtype);
    }
    push @attrs, "label=\"$value\"";
    $self->PrLine($self->MapNodeId("MACROPROCESS", $processtype) . " " .
        "[" .
        join(", ", @attrs) .
        "];") ;

    for $atomid (keys %{ $processtype2atom{$processtype} }) {
      ##
      ## process-conditions
      ##
      for $cond (@{ $pw->AtomCondition($atomid) }) {
        my $src  = $self->MapNodeId("MACROPROCESS", $cond);
        my $targ = $self->MapNodeId("MACROPROCESS", $processtype);
        if (! defined $self->{edgecache}{"$src,$targ"} ) {
          $self->{edgecache}{"$src,$targ"} = 1;
          $self->PrLine($src . " " .
              "-> " . 
              $targ . " " .
              "[" .
              $self->DOTAttr("color",
                    $lv->StringToLabelValue("edge-type", "agent")) .
              "];") ;
        }
      }
      ##
      ## now do regular edges
      ##
      for my $edgeid (@{ $pw->Edges($atomid) }) {
        my $molinst = $pw->MolInstId($atomid, $edgeid);
        DOTEdge($self, $atomid, $edgeid, $molinst);
        if (not $visited{$molinst}) {
          $visited{$molinst} = 1;
          DOTMol($self, $atomid, $edgeid, $molinst);
        }
      }
    }

  }

}

######################################################################
sub SubtypeEdges {
  my ($self, $molinstset) = @_;
  
  my $pw = $self->{pw};
  my $lv = $self->{lv};

  my @set = keys %{ $molinstset };
  my ($i, $j, $i2j, $j2i, $source, $target);
  for ($i = 0; $i < @set-1; $i++) {
    for ($j = $i + 1; $j < @set; $j++) {
#print STDERR "comparing $set[$i], $set[$j]\n";
      $i2j = 0; $j2i = 0;
      if ($pw->CachedMolInstSubtype($set[$i], $set[$j])) {
        $i2j = 1;
      }
      if ($pw->CachedMolInstSubtype($set[$j], $set[$i])) {
        $j2i = 1;
      }
      if ($i2j) {
        $source = $set[$i];
        $target = $set[$j];
      } elsif ($j2i) {
        $source = $set[$j];
        $target = $set[$i];
      } else {
        next;
      }
      $self->PrLine($self->MapNodeId("MOLECULE", $source) . " " .
             "-> " . 
             $self->MapNodeId("MOLECULE", $target) . " " .
               "[" .
                 (($i2j && $j2i) ? "dir=\"both\" " : "") .
                 "color=\"black\" " .
                 "style=\"dotted\" " .
               "];") ;
    }
  }

}

######################################################################
sub DOTSubtypeEdges {
  my ($self) = @_;
  
  my $pw = $self->{pw};
  my $lv = $self->{lv};

  my $cache = $self->{id2inst};

  for (keys %{ $cache } ) {
    SubtypeEdges($self, $$cache{$_});
  }

}

######################################################################
sub DOTGraph {
  my ($self) = @_;
  
  my $pw = $self->{pw};
  my $lv = $self->{lv};

  my ($atom);

  $self->PrLine("digraph G {");
  Indent();
  DOTMacroProcesses($self);
  for $atom (@{ $pw->Atoms }) {
    if (! $pw->IsMacroProcess($atom)) {
      DOTAtom($self, $atom);
    }
  }

  DOTSubtypeEdges($self);

  Exdent();

  $self->PrLine("}");

}

######################################################################
sub DOTAtom {
  my ($self, $atom) = @_;

  my ($molinst, $edge);
  my $pw = $self->{pw};
  my $lv = $self->{lv};

  DOTProcess($self, $atom);

  for $edge (@{ $pw->Edges($atom) }) {
    $molinst = $pw->MolInstId($atom, $edge);
    DOTEdge($self, $atom, $edge, $molinst);
    if (not $visited{$molinst}) {
      $visited{$molinst} = 1;
      DOTMol($self, $atom, $edge, $molinst);
    }
  }

}

######################################################################
sub DOTMol {
  my ($self, $atom, $edge, $molinst) = @_;

  my $pw = $self->{pw};
  my $sim = $self->{sim};
  my $molid = $pw->EdgeMol($atom, $edge);
  my $nt = $pw->MolType($molid);
  my @attrs;
  my ($molname, $cx_molname, $cx_diacritics);
  my $fh = $self->{output_fh};


  my $shapecolor = $self->{f_Value2Color}($sim->MolVal($molinst));
  my $fontcolor = "white";
  my $style = "filled";
  if ($shapecolor eq "white" || $shapecolor eq "#FFFFFF" ||
      $shapecolor eq "") {
    $fontcolor = "black";
    $style = "";
  }


  ## save so we can check for subtypes later
  $self->{id2inst}{$molid}{$molinst} = 1;

  for my $a ("shape", "height", "width", "fontsize") {
    push @attrs, DOTAttr($self, $a, $nt);
  }

  push @attrs, "color=\"$shapecolor\"";
  push @attrs, "fontcolor=\"$fontcolor\"";
  push @attrs, "style=\"$style\"";

  if ($self->MolIsComplex($molid)) {
    $cx_molname = "cx_$molid";
    $cx_diacritics .= $self->MolActivityDiacritic(
        $pw->MolLabel($atom, $edge));
    $cx_diacritics .= $self->LocationDiacritic(
        $pw->MolLabel($atom, $edge));
    $molname = "<$cx_molname>$cx_diacritics";
  } else {
    $molname = $pw->PickMolName($molid);
    $molname .= $self->MolActivityDiacritic(
        $pw->MolLabel($atom, $edge));
    $molname .= $self->LocationDiacritic(
        $pw->MolLabel($atom, $edge));
  }

  push @attrs, "label=\"" . $molname . "\"";
  my $node_id = $self->MapNodeId("MOLECULE", $molinst);

  $self->PrLine($self->MapNodeId("MOLECULE", $molinst) . " " .
      "[" .
      join(", ", @attrs) .
      "];") ;
}

######################################################################
sub DOTProcess {
  my ($self, $atomid) = @_;

  my $pw = $self->{pw};
  my $sim = $self->{sim};
  my $lv = $self->{lv};
  my $nt = $pw->AtomType($atomid);

  my $shapecolor = $self->{f_Value2Color}($sim->AtomVal($atomid));
  my $fontcolor = "white";
  my $style = "filled";
  if ($shapecolor eq "white" || $shapecolor eq "#FFFFFF" ||
      $shapecolor eq "") {
    $fontcolor = "black";
    $style = "";
  }

  my @attrs;
  my $layer;

  for my $a ("shape", "height", "width", "fontsize") {
    push @attrs, DOTAttr($self, $a, $nt);
  }

  push @attrs, "fontcolor=\"$fontcolor\"";
  push @attrs, "style=\"$style\"";

#  if (@{ $pw->LayersOf($atomid) }) {
#    for (@{ $pw->LayersOf($atomid) }) {
#      if (!defined $layer || $layer > $_) {
#        $layer = $_;
#      }
#    }
#  } else {
#    $layer = 0;
#  }
#  if ($layer == 0) {
#    push @attrs, DOTAttr($self, "color", $nt);
#  } else {
#    push @attrs, "color=\"" . $layer_color{$layer} . "\"";
#  }

  if ($self->{showatomids}) {
    push @attrs, "label=\"$atomid\"";
    push @attrs, "fontcolor=\"white\"";
  } else {
    push @attrs, "label=\"\"";
  }

  $self->PrLine($self->MapNodeId("ATOM", $atomid) . " " .
      "[" .
      join(", ", @attrs) .
      "];") ;
}

######################################################################
sub DOTEdgeDirection {
  my ($self, $atom, $edge) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};
  my $reversible = 0;

  for my $value (@{ $pw->AtomLabel($atom) }) {
    if ($lv->IsA($value, $lv->StringToLabelValue("reversible", "yes"))) {
      $reversible = 1;
      last;
    }
  }

  my $edgetype = $pw->EdgeType($atom, $edge);
  if ($lv->IsA($edgetype,
      $lv->StringToLabelValue("edge-type", "incoming-edge"))) {
    if ($lv->IsA($edgetype,
      $lv->StringToLabelValue("edge-type", "agent"))) {
      if ($reversible) {
        return "IN:NEITHER";
      }
    } elsif ($lv->IsA($edgetype,
      $lv->StringToLabelValue("edge-type", "inhibitor"))) {
      if ($reversible) {
        return "IN:NEITHER";
      }
    }
    if ($reversible) {
      return "IN:BOTH";
    } else {
      return "IN";
    }
  } elsif ($lv->IsA($edgetype,
      $lv->StringToLabelValue("edge-type", "outgoing-edge"))) {
    if ($reversible) {
      return "OUT:BOTH";
    } else {
      return "OUT";
    }
  } else {
    print STDERR "unrecognized edge type $edgetype\n";
    return "";
  }
}

######################################################################
sub DOTEdge {
  my ($self, $atom, $edge, $molinst) = @_;
  
  my $pw = $self->{pw};
  my $lv = $self->{lv};

  my ($edgetype, $label, @attrs);
  my $edgetype = $pw->EdgeType($atom, $edge);
  my $dir = DOTEdgeDirection($self, $atom, $edge);
  my ($x, $edgelabel) = $lv->LabelValueToString($edgetype);

  for my $a ("arrowhead", "color") {
    push @attrs, DOTAttr($self, $a, $edgetype);
  }

  for (@{ $pw->EdgeLabel($atom, $edge) }) {
    my ($label, $value) = $lv->LabelValueToString($_);
    if ($label eq "function") {
      push @attrs, "label = \"$value\"";
      last;
    }
  }

  if ($dir =~ /:BOTH$/) {
    push @attrs, "dir=\"both\"";
  } elsif ($dir =~ /:NEITHER$/) {
    push @attrs, "dir=\"none\"";
  }

  Indent();
  my ($src, $targ);
  if ($dir =~ /^IN/) {
    $src  = $self->MapNodeId("MOLECULE", $molinst);
    $targ = $self->MapNodeId("ATOM", $atom);
    if ($pw->IsMacroProcess($atom) && defined $self->{edgecache}{"$src,$targ"}) {
    } else {
      $self->{edgecache}{"$src,$targ"} = 1;
      $self->PrLine($src . " " .
             "-> " . 
             $targ . " " .
               "[" .
                 join(", ", @attrs) .
               "];") ;
    }
  } elsif ($dir =~ /^OUT/) {
    $targ = $self->MapNodeId("MOLECULE", $molinst);
    $src  = $self->MapNodeId("ATOM", $atom);
    if ($pw->IsMacroProcess($atom) && defined $self->{edgecache}{"$src,$targ"}) {
    } else {
      $self->{edgecache}{"$src,$targ"} = 1;
      $self->PrLine($src . " " .
             "-> " . 
             $targ . " " .
               "[" .
                 join(", ", @attrs) .
               "];") ;
    }
  } else {
    print STDERR "unrecognized direction: $dir\n";
  }

  Exdent();
}

######################################################################
1;
######################################################################

 
