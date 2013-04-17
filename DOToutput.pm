#!/usr/local/bin/perl


# Copyright SRA International
#
# Distributed under the OSI-approved BSD 3-Clause License.
# See http://ncip.github.com/pathway-interaction-database/LICENSE.txt for details.


package DOToutput;
require Exporter;
@ISA = qw(Exporter);
#@EXPORT = qw(
#  CreateIsA
#  CreateLabelSort
#);

use strict;
use Pathway;
use PWLabel;
use Clan;
use URI::Escape;
use Coloring;

my (%mol_visited, %edge_visited, %atom_visited);

use constant MAX_COMPONENTS_TO_SHOW => 5;

my $BLUE_RED_8 = [
    "0000FF",    ## blue = low
    "3399FF",
    "66CCFF",
    "99CCFF",
    "FFFFFF",  ## insert white neutral
#    "CCCCFF",    ## skip this -- too close
#    "FFCCFF",    ## skip this -- too close
    "FF99FF", 
    "FF66CC", 
    "FF6666",
    "FF0000"     ## red = high
  ];

my $CHANNEL_3 = [
    "0000FF",
    "CC00FF",
    "FF0000" 
];

my $BROWN_GRAY_5 = [
    "000000",
    "cc6600",
    "ffcc33",
    "ffcc99",
    "ffffcc",
    "cccccc"
  ];

my $SIM_COLORING = new Coloring
      ([1.8, 2.7, 3.6, 4.5, 5.5, 6.4, 7.3, 8.2, 9], $BLUE_RED_8);
my $DEGREE_COLORING = new Coloring
      ([0, 1, 2, 3, 4, 5], $BROWN_GRAY_5,);
my $CHANNEL_COLORING = new Coloring
      ([0, 1, 2, 3], $CHANNEL_3,);
#my $SIM_COLORING = new Coloring
#      (4.5, 5.5, [1.8,2.7,3.6,4.5,6.4,7.3,8.2,9], $BLUE_RED_8,   "FFFFFF");
#my $DEGREE_COLORING = new Coloring
#      (0,   0,   [1,2,3,4,5],                     $BROWN_GRAY_5, "000000");


######################################################################
sub r_numerically { $b <=> $a };

######################################################################
sub MOLPAGE_URL {
  my ($molid) = @_;
  
  return "/MoleculePage?molid=$molid";
}

######################################################################
sub MOLINST_URL {
  my ($atom, $edge) = @_;
  
  $_ = $atom;
  $atom =~ s/_[\d]+//g;

  return "/MoleculeInstance?inst=$atom,$edge";
}

######################################################################
sub ATOMPAGE_URL {
  my ($atomid) = @_;
  return "/InteractionPage?atomid=$atomid";
}

######################################################################
sub PATHWAYPAGE_URL {
  my ($pid, $graphic_format) = @_;

  return "/search/pathway_landing.shtml%3Fpathway_id%3D$pid" .
      "%26what%3Dgraphic%26$graphic_format%3Don%26source%3D5%26ppage%3D1";
}

######################################################################
sub SetGraphicFormat {
  my ($self, $format) = @_;
  $self->{graphicformat} = lc($format);
}

######################################################################
sub SetHTMLFileName {
  my ($self, $f) = @_;
  $self->{html_fn} = $f;
}

######################################################################
sub numerically { $a <=> $b };

######################################################################
sub BreakUpName {
  my ($name) = @_;

  my $SEGMENT = 15;
  my @chars = ("/", " ", "_");
  my @lines;

  while(length($name) > 0) {
    my %hits;
    for my $c (@chars) {
      my $idx = index($name, $c, $SEGMENT);
      if ($idx > -1) {
        $hits{$idx} = $c;
      }
    }
    if (keys %hits) {
      my @hits = sort numerically keys %hits;
      my $i = $hits[0];
      push @lines, substr($name, 0, $i+1);
      $name = substr($name, $i+1);
    } else {
      last;
    }
  }
  if (length($name) > 0) {
    push @lines, $name;
  }
  return join("\\n", @lines);
}

######################################################################
sub SourceLabel {
  my ($self, $clan,$pid) = @_;

  my $pw = $self->{pw};

  my $curated;
  my $biocarta;
  my $kegg;
  my $reactome;
  my $other;
  my $pathwayName = 'x';
  my $pidtemp = 999999;
  my $pathwayNameCount = 0;
  my $first = 1;

  $pathwayName = $pw->{pathway}{$pid}{pname};

  for my $atom (@{ $pw->Atoms }) {
    if ($clan && ! $clan->HasAtom($atom)) {
      next;
    }
    my $srcid = $pw->AtomSource($atom);
    if ($srcid eq "5") {
      $curated++;
    } elsif ($srcid eq "2" || $srcid eq "3") {
      $biocarta++;
    } elsif ($srcid eq "1") {
      $kegg++;
    } elsif ($srcid eq "7") {
      $reactome++;
    }
  }
  my $source_count;
  my $source_string;
  for my $x ($curated, $biocarta, $kegg, $reactome) {
    if ($x) {
      $source_count++;
    }
  }
  if ($source_count > 1) {
    $source_string = "";
  } else {
    if ($curated) {
      $source_string = "NCI-Nature Curated";
    } elsif ($biocarta) {
      $source_string = "BioCarta Imported"
    } elsif ($kegg) {
      $source_string = "KEGG";
    } elsif ($reactome) {
      $source_string = "Reactome";
    }
  }
  
  # if ($pathwayNameCount > 0)  {
     $source_string = $pathwayName;
  # }
  
  $self->PrLine(
     "source_label [ " .
     "shape=\"box\" " .
     "label=\"PID $source_string\" " .
     "color=\"firebrick\" " .
     "style=\"filled\" " .
     "fontcolor=\"white\" " .
     "fontname=\"arial\" " .
     "fontsize=\"12\" " .
    "];"
  );

}
######################################################################
sub new {
  my ($self, $pw, $lv, $outfh) = @_;
  my $x = {};

  undef %mol_visited;
  undef %edge_visited;
  undef %atom_visited;

  $x->{lv} = $lv;
  $x->{pw} = $pw;
  $x->{output_fh} = $outfh;


  $x->{graphicformat} = "svg";    ## default

  ##
  ## Molecules
  ##

  $x->{shape}{$lv->StringToLabelValue("molecule-type",
      "molecule-type")}  = "plaintext";
  $x->{shape}{$lv->StringToLabelValue("molecule-type",
      "protein")}  = "plaintext";
  $x->{shape}{$lv->StringToLabelValue("molecule-type",
      "compound")} = "plaintext";
  $x->{shape}{$lv->StringToLabelValue("molecule-type",
      "rna")}      = "plaintext";
  $x->{shape}{$lv->StringToLabelValue("molecule-type",
      "complex")}  = "Mrecord";

  $x->{color}{$lv->StringToLabelValue("molecule-type",
      "molecule-type")} = "black";                  ## default

  $x->{height}{$lv->StringToLabelValue("molecule-type",
      "molecule-type")}   = "";

  $x->{width}{$lv->StringToLabelValue("molecule-type",
      "molecule-type")}    = "";

  $x->{fontsize}{$lv->StringToLabelValue("molecule-type",
      "molecule-type")} = "10";

  $x->{fontsize}{$lv->StringToLabelValue("molecule-type",
      "complex")} = "10";

  $x->{style}{$lv->StringToLabelValue("molecule-type",
      "molecule-type")}    = "";

  $x->{fontname}{$lv->StringToLabelValue("molecule-type",
      "molecule-type")} = "arial";

  ##
  ## Processes
  ##

  $x->{shape}{$lv->StringToLabelValue("process-type",
     "process-type")}  = "hexagon";
  $x->{shape}{$lv->StringToLabelValue("process-type",
     "subnet")}  = "box";
  $x->{shape}{$lv->StringToLabelValue("process-type",
     "pathway")}  = "box";
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

  $x->{color}{$lv->StringToLabelValue("process-type",
      "subnet")} = "lightgray";

  $x->{color}{$lv->StringToLabelValue("process-type",
      "pathway")} = "lightgray";

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

  $x->{style}{$lv->StringToLabelValue("process-type",
      "subnet")}   = "filled";
  $x->{style}{$lv->StringToLabelValue("process-type",
      "pathway")}   = "filled";

  $x->{fontname}{$lv->StringToLabelValue("process-type",
      "process-type")} = "arial";

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

  $x->{fontname}{$lv->StringToLabelValue("edge-type",
      "edge-type")} = "arial";

  return bless $x;
}

my %location_code = (
  "transmembrane"         => "m",
  "cytoplasm"             => "cy",
  "nucleus"               => "n",
#  "endoplasmic reticulum" => "er",
  "extracellular region"         => "ex",
#  "calcium store"         => "cs",
#  "primary endosome"      => "pe",
#  "endosome"              => "e",
#  "early endosome"        => "ee",
#  "late endosome"         => "le",
#  "recycling endosome"    => "re",
#  "endosome transmembrane" => "et",
  "golgi"                 => "g",
  "lysosome"              => "l",
  "mitochondria"          => "mi",
  "vesicle"               => "v",
  "intracellular"         => "i",
  
  ""                      => ""
);

######################################################################
sub SetStandalone {
  my ($self, $standalone) = @_;
  $self->{standalone} = $standalone;
}

######################################################################
sub SetSimColoring {
  my ($self, $coloring) = @_;

  $SIM_COLORING = $coloring;
}

######################################################################
sub SetMolValueGenerator {
  my ($self, $vg) = @_;

  $self->{molvg} = $vg;
}

######################################################################
sub SetAtomValueGenerator {
  my ($self, $vg) = @_;

  $self->{atomvg} = $vg;
}

######################################################################
sub ShowAtomIds {
  my ($self, $off_or_on) = @_;

  $self->{showatomids} = $off_or_on;
}

######################################################################
sub ShowSubTypeLines {
  my ($self, $off_or_on) = @_;

  $self->{showsubtypelines} = $off_or_on;
}

######################################################################
sub CollapseAtoms {
  my ($self, $off_or_on) = @_;

  $self->{collapseatoms} = $off_or_on;
}

######################################################################
sub CollapseMolecules {
  my ($self, $off_or_on) = @_;

  $self->{collapsemolecules} = $off_or_on;
}

######################################################################
sub ColorMols {
  my ($self, $mols) = @_;

  for my $m (keys %{ $mols }) {
    $self->{colormols}{$m} = 1;
  }
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
        $label_value = CleanString($label_value);
        $diacritic .= "[$label_value]";
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
  my ($label_name, $label_value);
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

  if (scalar(@{ $pw->Components($cx_mol_id) }) > MAX_COMPONENTS_TO_SHOW) {
    push @{ $lines }, "...";
    return;
  }
  for $comp_seq_id (@{ $pw->Components($cx_mol_id) }) {
    $sub_mol_id = $pw->ComponentMol($cx_mol_id, $comp_seq_id);

#    if ($self->MolIsComplex($sub_mol_id)) {
#      $self->ComplexName($sub_mol_id, $lines);
#    } else {
#      $molname = $pw->PickMolName($sub_mol_id);
#      $molname .= $self->MolActivityDiacritic(
#          $pw->ComponentLabel($cx_mol_id, $comp_seq_id));
#      $molname .= $self->LocationDiacritic(
#          $pw->ComponentLabel($cx_mol_id, $comp_seq_id));
#    }

    if ($self->MolIsComplex($sub_mol_id)) {
##    $self->ComplexName($sub_mol_id, $lines);
      $molname = BreakUpName(CleanString($pw->PickMolName($sub_mol_id)));
#      if ($molname =~ /\//) {
#        $molname = "\\<cx_$sub_mol_id\\>";
#      } else {
        $molname = "\\<$molname\\>";
#      }
    } else {
      $molname = BreakUpName(CleanString($pw->PickMolName($sub_mol_id)));
    }
    $molname .= $self->MolActivityDiacritic(
        $pw->ComponentLabel($cx_mol_id, $comp_seq_id));
    $molname .= $self->LocationDiacritic(
        $pw->ComponentLabel($cx_mol_id, $comp_seq_id));

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
sub CleanString {
  my ($s) = @_;

  while ($s =~ /([\"\'\200-\377])/) {
    my $c = $1;
    my $x = sprintf("&#x%x;", ord($c));
    $s =~ s/$c/$x/g;
  }
  return $s;
}

######################################################################
sub MapNodeId {
  my ($self, $what, $id) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};

  if ($what eq "ATOM") {
    my $atomtype = $pw->AtomType($id);
    if ($pw->IsMacroProcess($id)) {
      return $self->MapNodeId("MACROPROCESS", $atomtype);
    } elsif (
        $atomtype eq $lv->StringToLabelValue("process-type", "subnet") ||
        $atomtype eq $lv->StringToLabelValue("process-type", "pathway")) {
      my $abs_id = $pw->AbstractionId($id);
      if ($abs_id eq "") {
        $abs_id = $pw->AbstractionExtId($id);
      }
      return ($self->MapNodeId("PATHWAY", $abs_id));
    }
  }

  if ($what eq "MACROPROCESS") {
    ## $id is the label_value_id for the kind of macroprocess
    return "P_$id";
  } elsif ($what eq "ATOM") {
    return "A_$id";
  } elsif ($what eq "PATHWAY") {
    return "N_$id";
  } elsif ($what eq "MOLECULE") {
    my $i = $pw->MolInstIdToString($id);
    $i =~ s/:/__/g;
    $i =~ s/,/_/g;
    $i =~ s/_$//;
    $i = "M_$i";
    return $i;
  } else {
    print STDERR "MapNodeId: illegal what = $what\n";
  }
}

######################################################################
sub DOTAbstractions {
  my ($self, $clan) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};

  my (%pid2atom, %extpid2atom);
  my $graphic_format = $self->{graphicformat};


  my $subnet_type = $lv->StringToLabelValue("process-type",
      "subnet");
  my $pathway_type = $lv->StringToLabelValue("process-type",
      "pathway");

  if ($subnet_type eq "0" && $pathway_type ne "0") {
    $subnet_type = $pathway_type;
  } elsif ($pathway_type eq "0" && $subnet_type ne "0") {
    $pathway_type = $subnet_type;
  }

  for my $atomid (@{ $pw->Atoms() }) {
    if ($clan && (! $clan->HasAtom($atomid))) {
      next;
    }
    my $atomtype = $pw->AtomType($atomid);
    if ($atomtype ne $subnet_type && $atomtype ne $pathway_type) {
      next;
    }

    my $abs_pid = $pw->AbstractionId($atomid);
    if ($abs_pid ne "") {
      $pid2atom{$abs_pid}{$atomid} = 1;
    } else {
      my $ext_pid = $pw->AbstractionExtId($atomid);
      $extpid2atom{$ext_pid}{$atomid} = 1;
    }
  }
  ## following situation should never arise, but we'll check
  ## to make sure
  for my $abs_pid (keys %pid2atom) {
    my $ext_pid = $pw->PathwayExId($abs_pid);
    if (defined $extpid2atom{$ext_pid}) {
      for my $atom (keys %{ $extpid2atom{$ext_pid} }) {
        $pid2atom{$abs_pid}{$atom} = 1;
      }
      delete $extpid2atom{$ext_pid};
    }
  }

  for my $abs_pid (keys %pid2atom) {
    my (@attrs);
    for my $a ("color", "shape", "height", "width", "fontsize", "style", "fontname") {
      push @attrs, DOTAttr($self, $a, $pathway_type);
    }
    push @attrs, "label=\"" . BreakUpName(CleanString($pw->PathwayName($abs_pid))) .
          "\"";
    if ($self->{standalone}) {
      push @attrs, "URL=\"" .  $pw->PathwayExId($abs_pid) .
          ".$graphic_format\"";
    } else {
        push @attrs, "URL=\"" . PATHWAYPAGE_URL($abs_pid,
            $graphic_format) . "\"";
    }
    
    $self->PrLine($self->MapNodeId("PATHWAY", $abs_pid) . " " .
        "[" .
        join(", ", @attrs) .
        "];") ;

    for my $atom (keys %{ $pid2atom{$abs_pid} }) {
      my %seen;
      for my $edge (@{ $pw->Edges($atom) }) {
        my $molinst = $pw->MolInstId($atom, $edge);
        my $edge_type = $pw->EdgeType($atom, $edge);
        if (defined $seen{$molinst}{$edge_type}) {
          next;
        } else {
          $seen{$molinst}{$edge_type} = 1;
          $self->DOTEdge($atom, $edge, $molinst);
          $self->DOTMol($atom, $edge, $molinst);
        }
      }
    }
  }

  for my $ext_pid (keys %extpid2atom) {
    my (@attrs);
    for my $a ("color", "shape", "height", "width", "fontsize", "style", "fontname") {
      push @attrs, DOTAttr($self, $a, $pathway_type);
    }
    push @attrs, "label=\"" . BreakUpName($ext_pid) .
          "\"";
    if ($self->{standalone}) {
      push @attrs, "URL=\"" .  $ext_pid .
          ".$graphic_format\"";
    } else {
      push @attrs, "URL=\"" . "/search/pathway_landing.shtml%3F" .
          "pathway_ext_id%3D" . $ext_pid .
          "%26what%3Dgraphic%26$graphic_format%3Don%26source%3D5" .
          "\"";
    }
    $self->PrLine($self->MapNodeId("PATHWAY", $ext_pid) . " " .
        "[" .
        join(", ", @attrs) .
        "];") ;

    for my $atom (keys %{ $extpid2atom{$ext_pid} }) {
      my %seen;
      for my $edge (@{ $pw->Edges($atom) }) {
        my $molinst = $pw->MolInstId($atom, $edge);
        my $edge_type = $pw->EdgeType($atom, $edge);
        if (defined $seen{$molinst}{$edge_type}) {
          next;
        } else {
          $seen{$molinst}{$edge_type} = 1;
          $self->DOTEdge($atom, $edge, $molinst);
          $self->DOTMol($atom, $edge, $molinst);
        }
      }
    }
  }

}

######################################################################
sub DOTMacroProcesses {
  my ($self, $clan) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};

  my $macroprocess_type = $lv->StringToLabelValue("process-type",
      "macroprocess");
  my ($atomid, $processtype, $cond, %processtype2cond,
      %processtype2atom);
  my ($label, $value);
  my (%missing);
  my $macro_atom_ids = $pw->MacroProcesses();

  for $atomid (keys %{ $macro_atom_ids }) {
    if ($clan && (! $clan->HasAtom($atomid))) {
      next;
    }
    $processtype = $$macro_atom_ids{$atomid};
    $processtype2atom{$processtype}{$atomid} = 1;
    for $cond (@{ $pw->AtomCondition($atomid) }) {
      $processtype2cond{$processtype}{$cond} = 1;
    }
  
    for $cond (@{ $pw->AtomNegativeCondition($atomid) }) {
      $processtype2cond{$processtype}{$cond} = 1;
    }

  }

  ##
  ## Do conditions on "regular" atoms; note referenced
  ## macroprocess types that are not "declared" in this model;
  ##
  for $atomid (@{ $pw->Atoms() }) {
    if ($clan && (! $clan->HasAtom($atomid))) {
      next;
    }
    for $cond (@{ $pw->AtomCondition($atomid) }) {
      if (! defined $processtype2atom{$cond}) {
        $missing{$cond} = 1;
      }
      if (defined $$macro_atom_ids{$atomid}) {
        next;     ## skip macroprocesses... do later
      }
      if ($self->{collapseatoms}) {
        for my $edge (@{ $pw->Edges($atomid) }) {
          my $edgetype = $pw->EdgeType($atomid, $edge);
          if ($lv->IsA($edgetype,
              $lv->StringToLabelValue("edge-type", "outgoing-edge"))) {
            my $molinst = $pw->MolInstId($atomid, $edge);
            my ($src, $targ);
            $src  = $self->MapNodeId("MACROPROCESS", $cond);
            $targ = $self->MapNodeId("MOLECULE", $molinst);
            Indent();
            if (! defined $self->{edgecache}{"$src,$targ"}) {
              $self->{edgecache}{"$src,$targ"} = 1;
              $self->PrLine($src . " " .
                     "-> " . 
                     $targ . " " .
                       "[" .
                         "color=green" .
                       "];") ;
            }
            Exdent();
          }
        }
      } else {
        $self->PrLine($self->MapNodeId("MACROPROCESS", $cond) . " " .
            "-> " . 
            $self->MapNodeId("ATOM", $atomid) . " " .
            "[" .
            $self->DOTAttr("color",
                  $lv->StringToLabelValue("edge-type", "agent")) .
            "];") ;
      }
    }
    for $cond (@{ $pw->AtomNegativeCondition($atomid) }) {
      if (! defined $processtype2atom{$cond}) {
        $missing{$cond} = 1;
      }
      if (defined $$macro_atom_ids{$atomid}) {
        next;     ## skip macroprocesses... do later
      }
      if ($self->{collapseatoms}) {
        for my $edge (@{ $pw->Edges($atomid) }) {
          my $edgetype = $pw->EdgeType($atomid, $edge);
          if ($lv->IsA($edgetype,
              $lv->StringToLabelValue("edge-type", "inhibitor"))) {
            my $molinst = $pw->MolInstId($atomid, $edge);
            my ($src, $targ);
            $src  = $self->MapNodeId("MACROPROCESS", $cond);
            $targ = $self->MapNodeId("MOLECULE", $molinst);
            Indent();
            if (! defined $self->{edgecache}{"$src,$targ"}) {
              $self->{edgecache}{"$src,$targ"} = 1;
              $self->PrLine($src . " " .
                     "-> " .
                     $targ . " " .
                       "[" .
                         "color=red" .
                       "];") ;
            }
            Exdent();
          }
        }
      } else {
        $self->PrLine($self->MapNodeId("MACROPROCESS", $cond) . " " .
            "-> " .
            $self->MapNodeId("ATOM", $atomid) . " " .
            "[" .
            $self->DOTAttr("color",
                  $lv->StringToLabelValue("edge-type", "inhibitor")) .
            "];") ;
      }
    }

  }

  ##
  ## Create nodes for "undeclared" macroprocess types
  ##
  for $processtype (keys %missing) {
    ($label, $value) = $lv->LabelValueToString($processtype);
    my @attrs;
    for my $a ("shape", "height", "width", "fontsize", "style", "fontname") {
#      push @attrs, DOTAttr($self, $a, $processtype);
# some of the "macroprocesses are actually functions
      push @attrs, DOTAttr($self, $a, $macroprocess_type);
    }
    push @attrs, "label=\"" . BreakUpName(CleanString($value)) . "\"";
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
    for my $a ("shape", "height", "width", "fontsize", "color", "style", "fontname") {
#      push @attrs, DOTAttr($self, $a, $processtype);
# some of the "macroprocesses are actually functions
      push @attrs, DOTAttr($self, $a, $macroprocess_type);
    }
    push @attrs, "label=\"" . BreakUpName(CleanString($value)) . "\"";
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
      for $cond (@{ $pw->AtomNegativeCondition($atomid) }) {
        my $src  = $self->MapNodeId("MACROPROCESS", $cond);
        my $targ = $self->MapNodeId("MACROPROCESS", $processtype);
        if (! defined $self->{edgecache}{"$src,$targ"} ) {
          $self->{edgecache}{"$src,$targ"} = 1;
          $self->PrLine($src . " " .
              "-> " .
              $targ . " " .
              "[" .
              $self->DOTAttr("color",
                    $lv->StringToLabelValue("edge-type", "inhibitor")) .
              "];") ;
        }
      }

      ##
      ## now do regular edges
      ##
      for my $edgeid (@{ $pw->Edges($atomid) }) {
        my $molinst = $pw->MolInstId($atomid, $edgeid);
        DOTEdge($self, $atomid, $edgeid, $molinst);
        DOTMol($self, $atomid, $edgeid, $molinst);
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
sub DOTSubnet {
  my ($self, $subnet, $clan) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};

  my $count = 0;
  for my $atom (@{ $pw->AtomsInSubnet($subnet) }) {
    if ($clan && (! $clan->HasAtom($atom))) {
      next;
    }
    if ($pw->IsMacroProcess($atom)) {
      next;
    }
    if (defined $atom_visited{$atom}) {
      next;
    }
    $count++;
  }
  if (! $count) {
    return;
  }
  $self->PrLine("subgraph cluster_$subnet {");
  $self->PrLine("graph [ label=\"subnet_$subnet\" ]");
  for my $atom (@{ $pw->AtomsInSubnet($subnet) }) {
    if ($clan && (! $clan->HasAtom($atom))) {
      next;
    }
    if ($pw->IsMacroProcess($atom)) {
      next;
    }
    if (defined $atom_visited{$atom}) {
      next;
    }
    if ($self->{collapseatoms}) {
      DOTCollapsedAtom($self, $atom);
    } else {
      DOTAtom($self, $atom);
    }
  }
  $self->PrLine("}");
}

######################################################################
sub DoGlobalStuff {
  my ($self, $clan) = @_;
  
  my $pw = $self->{pw};
  my $lv = $self->{lv};

  my (%atom_in_subnet, %molinst_in_subnet);
  my (%atom_in_global, %molinst_in_global);

  ## Assume connectivity of subnets. I.e., if one atom in a subnet
  ## is in the  clan, then the other atoms in the subnet are also in
  ## the clan.

  for my $subnet (@{ $pw->Subnets() }) {
    for my $atom (@{ $pw->AtomsInSubnet($subnet) }) {
      if ($clan && (! $clan->HasAtom($atom))) {
        next;
      }
      $atom_in_subnet{$atom}{$subnet} = 1;
      for my $edge (@{ $pw->Edges($atom) }) {
        my $molinst = $pw->MolInstId($atom, $edge);
        $molinst_in_subnet{$molinst}{$subnet} = 1;
      }
    }
  }
  for my $atom (@{ $pw->Atoms() }) {
    if ($clan && (! $clan->HasAtom($atom))) {
      next;
    }
    if (defined $atom_in_subnet{$atom}) {
      next;
    }
    $atom_in_global{$atom} = 1;
    for my $edge (@{ $pw->Edges($atom) }) {
      my $molinst = $pw->MolInstId($atom, $edge);
      $molinst_in_global{$molinst} = 1;
    }
  }
  for my $molinst (keys %molinst_in_subnet) {
    if (scalar(keys %{ $molinst_in_subnet{$molinst} }) > 1) {
      $molinst_in_global{$molinst} = 1;
    }
  }
  
  for my $atom (@{ $pw->Atoms() }) {
    if ($clan && (! $clan->HasAtom($atom))) {
      next;
    }
    if ($pw->IsMacroProcess($atom)) {
      next;
    }
    if (defined $atom_in_global{$atom}) {
      if ($self->{collapseatoms}) {
        DOTCollapsedAtom($self, $atom);
      } else {
        DOTAtom($self, $atom);
      }
    }
    for my $edge (@{ $pw->Edges($atom) }) {
      my $molinst = $pw->MolInstId($atom, $edge);
      if (defined $molinst_in_global{$molinst}) {
        DOTMol($self, $atom, $edge, $molinst);
        if (! $self->{collapseatoms}) {
          DOTEdge($self, $atom, $edge, $molinst);
        }
      }
    }
  }
}

######################################################################
sub DOTGraph {
  my ($self, $clan, $pid) = @_;
  
  my $pw = $self->{pw};
  my $lv = $self->{lv};

  my ($atom);
  my (%in_a_subnet, %count2subnet);

  ## As for getting molecules into/out of subnet boxes, doesn't seem
  ## to matter whether subgraph clusters appear before or after global
  ## interactions.

  ## dot will put a node inside a subgraph cluster if the node is attached
  ## to an edge that is declared a subgraph, and it will put the
  ## node in the first such subgraph that it encounters.

  ## It's possible for a given interaction to be in multiple subnets.
  ## But in the display, it's only possible to have an interaction
  ## in one subnet. Furthermore, dot will draw the interaction in the
  ## first cluster subgraph in which it is mentioned. So be greedy
  ## and specify the larger subnets first.

  for my $subnet (@{ $pw->Subnets() }) {
    my $count = 0;
    for my $atom (@{ $pw->AtomsInSubnet($subnet) }) {
      if ($clan && (! $clan->HasAtom($atom))) {
        next;
      }
      $count++;
      $in_a_subnet{$atom}++;
    }
    push @{ $count2subnet{$count} }, $subnet;
  }

  $self->PrLine("digraph G {");
  $self->SourceLabel($clan,$pid);
  Indent();

  ## Good or bad, this ends up putting all macroprocesses into
  ## the global area.

  DOTMacroProcesses($self, $clan);
  DOTAbstractions($self, $clan);
  DoGlobalStuff($self, $clan);

  for my $count (sort r_numerically keys %count2subnet) {
    for my $subnet (@{ $count2subnet{$count} }) {
      $self->DOTSubnet($subnet, $clan);
    }
  }

  if ($self->{showsubtypelines}) {
    DOTSubtypeEdges($self);
  }

  #HookUpFamilies($self, $clan);

  Exdent();

  $self->PrLine("}");

}

######################################################################
sub DOTAtom {
  my ($self, $atom) = @_;

  my ($molinst, $edge);
  my $pw = $self->{pw};
  my $lv = $self->{lv};

  if ($atom_visited{$atom}) {
    return;
  } else {
    $atom_visited{$atom} = 1;
  }

  my $atom_type = $pw->AtomType($atom);
  my $subnet_type  = $lv->StringToLabelValue("process-type", "subnet");
  my $pathway_type = $lv->StringToLabelValue("process-type", "pathway");

## actually, all abstractions are done elsewhere
  if ($atom_type eq $subnet_type || $atom_type eq $pathway_type) {
    return;
  }

  ## don't include subnet placeholders that don't connect to anything
  ## (this put in to accommodate Reactome data)
  if (($atom_type eq $subnet_type || $atom_type eq $pathway_type) &&
      @{ $pw->Edges($atom) } == 0 &&
      @{ $pw->AtomCondition($atom) } == 0) {
    return;
  }

  DOTProcess($self, $atom);

  for $edge (@{ $pw->Edges($atom) }) {
    $molinst = $pw->MolInstId($atom, $edge);
    DOTEdge($self, $atom, $edge, $molinst);
    DOTMol($self, $atom, $edge, $molinst);
  }

}

######################################################################
sub DOTCollapsedAtom {
  my ($self, $atom) = @_;

  my ($molinst, $edge, $edgetype, $dir);
  my $pw = $self->{pw};
  my $lv = $self->{lv};
  my (@attrs);
  my (%ins, %outs);

  for $edge (@{ $pw->Edges($atom) }) {
    $molinst = $pw->MolInstId($atom, $edge);
    $edgetype = $pw->EdgeType($atom, $edge);
    if ($lv->IsA($edgetype,
        $lv->StringToLabelValue("edge-type", "incoming-edge"))) {
      $ins{$edge} = $molinst;
    } else {
      $outs{$edge} = $molinst;
    }
    DOTMol($self, $atom, $edge, $molinst);
  }

  for $edge (@{ $pw->Edges($atom) }) {
    if (defined $ins{$edge}) {
      $edgetype = $pw->EdgeType($atom, $edge);
      for my $i (keys %outs) {
        if ($ins{$edge} eq $outs{$i}) {
          next; #!!! drop reflexive edges
        }
        undef @attrs;
        $dir = DOTEdgeDirection($self, $atom, $edge);
        if ($dir =~ /:BOTH$/) {
          push @attrs, "dir=\"both\"";
        } elsif ($dir =~ /:NEITHER$/) {
          push @attrs, "dir=\"none\"";
        }
        for my $a ("arrowhead", "color", "fontname") {
          push @attrs, DOTAttr($self, $a, $edgetype);
        }
        for (@{ $pw->EdgeLabel($atom, $edge) }) {
          my ($label, $value) = $lv->LabelValueToString($_);
          if ($label eq "function") {
            push @attrs, "label = \"$value\"";
            last;
          }
        }
        Indent();
        my ($src, $targ);
        $src  = $self->MapNodeId("MOLECULE", $ins{$edge});
        $targ = $self->MapNodeId("MOLECULE", $outs{$i});
        if (! defined $self->{edgecache}{"$src,$targ,$edgetype"}) {
          $self->{edgecache}{"$src,$targ,$edgetype"} = 1;
          $self->PrLine($src . " " .
                 "-> " . 
                 $targ . " " .
                   "[" .
                     join(", ", @attrs) .
                   "];") ;
        }
        Exdent();
      }
    }
  }

}

######################################################################
sub DOTMol {
  my ($self, $atom, $edge, $molinst) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};
  my $molid = $pw->EdgeMol($atom, $edge);
  my $nt = $pw->MolType($molid);
  my @attrs;
  my ($molname, $cx_molname, $cx_diacritics, @lines);
  my ($fontcolor, $shapecolor);
  my $fh = $self->{output_fh};
  my $vg;
  if (defined $self->{molvg}) {
    $vg = $self->{molvg};
  }

  if (defined $mol_visited{$molinst}) {
    return;
  } else {
    $mol_visited{$molinst} = 1;
  }

  ## save so we can check for subtypes later
  $self->{id2inst}{$molid}{$molinst} = 1;

  if (defined $self->{molvg} && $self->{molvg}->MolHasValue($molinst, $molid)) {
    $shapecolor = "#" .
        $SIM_COLORING->Value2Color($vg->MolVal($molinst, $molid));

#    $fontcolor = "#FFFFFF";
#    if (uc($shapecolor) eq "#FFFFFF") {
#      push @attrs, "color=\"#000000\"";       ## black border
#      push @attrs, "fontcolor=\"#000000\"";
#    } else {
#      push @attrs, "style=\"filled\"";
#      push @attrs, "color=\"$shapecolor\"";
#     if ($self->MolIsComplex($molid)) {
#        push @attrs, "fontcolor=\"$shapecolor\"";
#      } else {
#        push @attrs, "fontcolor=\"$fontcolor\"";
#      }
#    }
#    if ($nt == $lv->StringToLabelValue("molecule-type", "complex")) {
#      push @attrs, "shape=\"Mrecord\"";
#    } else {
#      push @attrs, "shape=\"ellipse\"";
#    }
#    for my $a ("height", "width", "fontsize") {
#      push @attrs, DOTAttr($self, $a, $nt);
#    }

## Let's try just coloring the name, not a filled ellipse
## In any case, style=filled won't actually color the background infields of
## a Mrecord
    $fontcolor;
    if (uc($shapecolor) eq "#FFFFFF") {
      $fontcolor = "#000000";
    } else {
      $fontcolor = $shapecolor;
    }
    push @attrs, "fontcolor=\"$fontcolor\"";
    for my $a ("shape", "height", "width", "fontsize", "fontname") {
      push @attrs, DOTAttr($self, $a, $nt);
    }


  } else {
    for my $a ("shape", "height", "width", "fontsize", "color", "style", "fontname") {
      push @attrs, DOTAttr($self, $a, $nt);
    }
    if (defined $self->{colormols}{$molid}) {
      push @attrs, "fontcolor=blue";
    }
  }

  if ($molid =~ /^\d+$/) {
#    push @attrs, "URL=\"" . MOLPAGE_URL($molid) . "\"";
    if ($self->{standalone}) {
      push @attrs, "URL=\"" . $self->{html_fn} .
          "#edge_$atom" . "_" . "$edge\"";
    } else {
      push @attrs, "URL=\"" . MOLINST_URL($atom, $edge) . "\"";
    }
  }

  my @comment_pairs = ("molid", $molid);

  if ($self->MolIsComplex($molid)) {
    $self->ComplexName($molid, \@lines);
    $cx_molname = BreakUpName(CleanString($pw->PickMolName($molid)));
#    if ($cx_molname =~ /\//) {
#      $cx_molname = "cx_$molid";
#    }
    $cx_diacritics .= $self->MolActivityDiacritic(
        $pw->MolLabel($atom, $edge));
    $cx_diacritics .= $self->LocationDiacritic(
        $pw->MolLabel($atom, $edge));
    if (! $self->{collapsemolecules}) {
      unshift @lines, "\\<$cx_molname\\>$cx_diacritics";
    }
    $molname = "{" . join("|", @lines) . "}";
    my @component_mol = map {$pw->ComponentMol($molid, $_)}
        @{$pw->Components($molid)};
    push @comment_pairs, "complex", join ",", @component_mol;
  } else {
    $molname = $pw->PickMolName($molid);
    $molname .= $self->MolActivityDiacritic(
        $pw->MolLabel($atom, $edge));
    $molname .= $self->LocationDiacritic(
        $pw->MolLabel($atom, $edge));
    $molname = BreakUpName(CleanString($molname));
  }

  $molname =~ s/"/\\"/g;
  push @attrs, "label=\"" . $molname . "\"";
  my $node_id = $self->MapNodeId("MOLECULE", $molinst);
  push @comment_pairs, ("node", $node_id);

  $self->PrLine($self->MapNodeId("MOLECULE", $molinst) . " " .
      "[" .
      join(", ", @attrs) .
      "];") ;

  my ($hash, $type);
  foreach $hash ($pw->MolExId($molid), $pw->MolName($molid)) {
    # external IDs, LocusLink, etc.
    foreach $type (keys %{$hash}) {
      push @comment_pairs, map {$type, $_} (keys %{$hash->{$type}});
    }
  }
  $self->PrLineNoIndent("// 200 " . join("|", @comment_pairs));
  # 200 = HTTP OK
}

######################################################################
sub DOTProcess {
  my ($self, $atomid) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};
  my $vg;
  my $nt = $pw->AtomType($atomid);
  my @attrs;
  my $layer;

  my $graphic_format = $self->{graphicformat};

  for my $a ("shape", "height", "width", "fontsize", "fontname") {
    if ($pw->IsMacroProcess($atomid)) {
      ## allow for macroprocesses that are actually functions
      push @attrs, DOTAttr($self, $a, $lv->StringToLabelValue("process-type",
          "macroprocess"));
    } else {
      push @attrs, DOTAttr($self, $a, $nt);
    }
  }

  if (defined $self->{atomvg} && $self->{atomvg}->AtomHasValue($atomid)) {
    $vg = $self->{atomvg};
    my $shapecolor = "#" . $SIM_COLORING->Value2Color(
        $vg->AtomVal($atomid));
    if (uc($shapecolor) eq "#FFFFFF") {
      push @attrs, "color=\"#000000\"";    ## put on a black border
    } else {
      push @attrs, "color=\"$shapecolor\"";
      push @attrs, "style=\"filled\"";
    } 
  } else {
    if ($nt eq $lv->StringToLabelValue("process-type", "subnet") ||
        $nt eq $lv->StringToLabelValue("process-type", "pathway")) {
      push @attrs, "fontcolor=\"black\"";
    }
    push @attrs, "style=\"filled\"";
    if (@{ $pw->LayersOf($atomid) }) {
      for (@{ $pw->LayersOf($atomid) }) {
        if (!defined $layer || $layer > $_) {
          $layer = $_;
        }
      }
    } else {
      $layer = 0;
    }
    if ($layer == 0) {
      push @attrs, DOTAttr($self, "color", $nt);
    } else {
      push @attrs, "color=\"#" . $DEGREE_COLORING->Value2Color($layer) . "\"";
    }
  }

  if ($self->{showatomids}) {
    if ($nt eq $lv->StringToLabelValue("process-type", "pathway") ||
        $nt eq $lv->StringToLabelValue("process-type", "subnet")) {
##      push @attrs, "label=\"" . $pw->PathwayExId($pw->Abstraction($atomid)) .
##          "\"";
      my $abs_pid = $pw->Abstraction($atomid);
      if ($abs_pid ne "") {
        push @attrs, "label=\"" . BreakUpName(CleanString($pw->PathwayName($abs_pid))) .
            "\"";
      } else {
        push @attrs, "label=\"" . BreakUpName($pw->AbstractionExtId($atomid)) .
            "\"";
      }
    } else {
      push @attrs, "label=\"$atomid\"";
      push @attrs, "fontcolor=\"white\"";
    }
  } elsif ($nt eq $lv->StringToLabelValue("process-type", "pathway") ||
           $nt eq $lv->StringToLabelValue("process-type", "subnet")) {
##    push @attrs, "label=\"" . $pw->PathwayExId($pw->Abstraction($atomid)) .
##        "\"";
    my $abs_pid = $pw->Abstraction($atomid);
    if ($abs_pid ne "") {
      push @attrs, "label=\"" . BreakUpName(CleanString($pw->PathwayName($abs_pid))) .
          "\"";
    } else {
      push @attrs, "label=\"" . BreakUpName($pw->AbstractionExtId($atomid)) .
          "\"";
    }
  } elsif ($nt eq $lv->StringToLabelValue("process-type", "subnet")) {
    push @attrs, "label=\"$atomid\"";
  } else {
    push @attrs, "label=\"\"";
  }

  if ($nt eq $lv->StringToLabelValue("process-type", "pathway") ||
      $nt eq $lv->StringToLabelValue("process-type", "subnet")) {
    if ($self->{standalone}) {
##      push @attrs, "URL=\"" .  $pw->PathwayExId($pw->Abstraction($atomid)) .
##          ".$graphic_format\"";
      push @attrs, "URL=\"" .  $pw->AbstractionExtId($atomid) .
          ".$graphic_format\"";
    } else {
##      push @attrs, "URL=\"" . PATHWAYPAGE_URL($pw->Abstraction($atomid),
##          $graphic_format) . "\"";
      if ($pw->AbstractionId($atomid)) {
        push @attrs, "URL=\"" . PATHWAYPAGE_URL($pw->AbstractionId($atomid),
            $graphic_format) . "\"";
      } elsif ($pw->AbstractionExtId($atomid)) {
        push @attrs, "URL=\"" . "/search/pathway_landing.shtml%3F" .
            "pathway_ext_id%3D" . $pw->AbstractionExtId($atomid) .
            "%26what%3Dgraphic%26$graphic_format%3Don%26source%3D5" .
            "\"";
      }
    
    }
  } elsif ($atomid =~ /^\d+$/) {
    if ($self->{standalone}) {
      push @attrs, "URL=\"" . $self->{html_fn} .
          "#interaction_$atomid\"";
    } else {
      push @attrs, "URL=\"" . ATOMPAGE_URL($atomid) . "\"";
    }
  } elsif ($atomid =~ /^(\d+)_\d+$/) {
## temp fix: deal with atoms resulting from the split of a
## multiple-output transcription
    if ($self->{standalone}) {
      push @attrs, "URL=\"" . $self->{html_fn} .
          "#interaction_$1\"";
    } else {
      push @attrs, "URL=\"" . ATOMPAGE_URL($1) . "\"";
    }
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
##    print STDERR "unrecognized edge type $edgetype\n";
    return ":NEITHER";
  }
}

######################################################################
sub DOTEdge {
  my ($self, $atom, $edge, $molinst) = @_;
  
  my $pw = $self->{pw};
  my $lv = $self->{lv};

  if (defined $edge_visited{"$atom,$edge"}) {
    return;
  } else {
    $edge_visited{"$atom,$edge"} = 1;
  }

  my ($edgetype, $label, @attrs);
  my $edgetype = $pw->EdgeType($atom, $edge);
  my $dir = DOTEdgeDirection($self, $atom, $edge);
  my ($x, $edgelabel) = $lv->LabelValueToString($edgetype);

  for my $a ("arrowhead", "color", "fontname") {
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
  } else {
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
##    print STDERR "unrecognized direction: $dir\n";
## go ahead and do something anyway

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
  }

  Exdent();
}

######################################################################
sub HookUpFamilies {
  my ($self, $clan) = @_;

  my $pw = $self->{pw};

  my (%temp, %tempinst);
  my @atoms;

  if (defined $clan) {
    @atoms = @{ $clan->ClanAtoms() }
  } else {
    @atoms = @{ $pw->Atoms() }
  }

  for my $atom (@atoms) {
    for my $edge (@{ $pw->Edges($atom) }) {
      my $mol     = $pw->EdgeMol($atom, $edge);
      my $molinst = $pw->MolInstId($atom, $edge);
      for my $parent (@{ $pw->FamilyParent($mol) }) {
        $temp{$parent}{$mol} = 1;
        $tempinst{$mol}{$molinst} = 1;
      }
    }
  }

  my $node_id;
  for my $parent (keys %temp) {
    if (scalar(keys %{ $temp{$parent} }) > 1) {
      $node_id++;
      my $src = "F_$node_id";
      my @node_attrs;
      push @node_attrs, "shape=\"box\"";
      push @node_attrs, "color=\"brown\"";
      push @node_attrs, "fontname=\"arial\"";
      push @node_attrs, "fontsize=\"10\"";
      push @node_attrs, "URL=\"" . MOLPAGE_URL($parent) . "\"";
      push @node_attrs, "label=" .
          "\"family:\\n" . BreakUpName(CleanString($pw->PickMolName($parent))) . "\"";
      $self->PrLine("F_$node_id" . " " .
          "[" .
             join(", ", @node_attrs) .
          "];") ;
      for my $mol (keys %{ $temp{$parent} }) {
        for my $molinst (keys  %{ $tempinst{$mol} }) {
          my $targ = $self->MapNodeId("MOLECULE", $molinst);
          my @edge_attrs;
          push @edge_attrs, "color=\"brown\"";
          push @edge_attrs, "style=\"dotted\"";
          $self->PrLine($src . " " .
                 "-> " . 
                 $targ . " " .
                   "[" .
                     join(", ", @edge_attrs) .
                   "];") ;
        }
      }
    }
  }

}

######################################################################
1;
######################################################################
