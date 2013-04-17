#!/usr/local/bin/perl


# Copyright SRA International
#
# Distributed under the OSI-approved BSD 3-Clause License.
# See http://ncip.github.com/pathway-interaction-database/LICENSE.txt for details.


package SBMLOutput;
require Exporter;
@ISA = qw(Exporter);
#@EXPORT = qw(
#);

use strict;
use Pathway;
use PWLabel;

my %basic_atom_types = (
  "reaction"      => 1,
  "binding"       => 1,         ## keep for now
  "modification"  => 1,
  "transcription" => 1,
  "translocation" => 1
#  "cell-process"  => 1
);

######################################################################
sub new {
  my ($self, $pw, $lv, $outfh) = @_;
  my $x = {};
  $x->{pw} = $pw;
  $x->{lv} = $lv;
  $x->{output_fh} = $outfh;
  return bless $x;
}

######################################################################
sub  numerically { $a <=> $b };

######################################################################
my $indent_level;
my $blanks =  "                                        ";
use constant INDENT_INCR => 2;

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

######################################################################
sub MoleculeList {
  my ($self) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};

  $self->PrLine("<listOfSpecies>");
  Indent();

  my ($a, $e, $ml, $mi, $label, $value, $location, %seen);

  for $a (@{ $pw->Atoms() }) {
    for $e (@{ $pw->Edges($a) }) {
      $mi = $pw->MolInstId($a, $e);
      $mi = "s$mi";
      if (! $seen{$mi}) {
        $seen{$mi} = 1;
        undef $location;
        $location = "cell_part";
        for $ml (@{ $pw->MolLabel($a, $e) }) {
          ($label, $value) = $lv->LabelValueToString($ml);
          if ($label eq "location") {
            $value =~ s/\s+/_/; 
            $location = $value;
          }
        }
        $self->PrLine("<species id=\"$mi\"" .
            ($location ? " compartment=\"$location\"" : "") .">");
        Indent();
        # $self->PrLine("<annotation>");
        # Indent();
        #$self->Molecule($pw->EdgeMol($a,$e));
        #for $ml (@{ $pw->MolLabel($a, $e) }) {
        #  ($label, $value) = $lv->LabelValueToString($ml);
        #  $self->Label($label, $value);
        #}
        #Exdent();
        #$self->PrLine("</annotation>");
        Exdent();
        $self->PrLine("</species>");
      }
    }
  }

  Exdent();
  $self->PrLine("</listOfSpecies>");
}

######################################################################
sub Molecule {
  my ($self, $id) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};

  if (defined $self->{molecules_to_define}{$id}) {
    delete $self->{molecules_to_define}{$id};
  }
  if (defined $self->{molecules_defined}{$id}) {
    return;
  } else {
    $self->{molecules_defined}{$id} = 1;
  }

  $self->PrLine("<Molecule id=\"$id\">");
  Indent();

  my $mol_id = $id;
  my ($label, $mol_type) = $lv->LabelValueToString($pw->MolType($mol_id));
  $self->Label($label, $mol_type);

  my $mehash = $pw->MolExId($mol_id);
  for my $id_type (keys %{ $mehash }) {
    for my $i (keys %{ $$mehash{$id_type} }) {
      $self->Name($id_type, $i);
    }
  }

  my $mnhash = $pw->MolName($mol_id);
  for my $name_type (keys %{ $mnhash }) {
    for my $i (keys %{ $$mnhash{$name_type} }) {
      $self->Name($name_type, $i);
    }
  }

  if ($mol_type eq "complex") {
    Indent();
    for my $comp (sort numerically @{ $pw->Components($mol_id) }) {
      my @labels;
      my $compmol = $pw->ComponentMol($mol_id, $comp);
      if (! defined $self->{molecules_defined}{$compmol} ) {
        $self->{molecules_to_define}{$compmol} = 1;
      }      
      for my $lvid (@{ $pw->ComponentLabel($mol_id, $comp) }) {
        my ($name, $value) = $lv->LabelValueToString($lvid);
        if ($name && $name ne "molecule-type") {
          push @labels, $lvid;
        }
      }
      $self->Component("ComplexComponent", $compmol, \@labels);
    }

    Exdent();
  }

  Exdent();
  $self->PrLine("</Molecule>");
  for my $m (keys %{ $self->{molecules_to_define} }) {
    $self->Molecule($m);
  }
}

######################################################################
sub ComponentList {
  my ($self) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};

}

######################################################################
sub Component {
  my ($self, $what, $id, $labels) = @_;

  ## what == PathwayComponent | InteractionComponent | ComplexComponent

  my $pw = $self->{pw};
  my $lv = $self->{lv};

  if ($labels && @{ $labels }) {
    $self->PrLine("<$what idref=\"$id\">");
    Indent();
    for my $lvid (@{ $labels }) {
      my ($label, $value) = $lv->LabelValueToString($lvid);
      $self->Label($label, $value);
    }
    Exdent();
    $self->PrLine("</$what>");
  } else {
    $self->PrLine("<$what idref=\"$id\" />");
  }
}

######################################################################
sub Label {
  my ($self, $label, $value) = @_;

  $self->PrLine("<Label label_type=\"$label\" value=\"$value\" />");
}

######################################################################
sub Name {
  my ($self, $nametype, $name) = @_;

  $self->PrLine("<Name name_type=\"$nametype\" value=\"$name\" />");
}

######################################################################
sub InteractionList {
  my ($self) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};

  $self->PrLine("<listOfReactions>");
  Indent();
  for my $atom (@{ $pw->Atoms() }) {
    $self->Interaction($atom);
  }
  Exdent();
  $self->PrLine("</listOfReactions>");
}

######################################################################
sub Interaction {
  my ($self, $id) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};

  my ($reversible);
  my $atom = $id;

  for my $lvid (@{ $pw->AtomLabel($atom) }) {
    my ($label, $value) = $lv->LabelValueToString($lvid);
    if ($label eq "reversible") {
      if ($value eq "yes") {
        $reversible="true";
      } elsif ($value eq "no") {
        $reversible = "false";
      }
    }
  }
  $self->PrLine("<reaction metaid=\"pid_i_$id\"" . " id=\"pid_i_$id\"" . " name=\"pid_i_$id\"" .
      ($reversible ? " reversible=\"$reversible\"" : "") . ">");

  Indent();
  # $self->PrLine("<annotation>");
  #Indent();
  #for my $lvid (@{ $pw->AtomCondition($atom) }) {
  #  my ($label, $value) = $lv->LabelValueToString($lvid);
  #  $self->PrLine("<Condition>$value</Condition>");
  #}
  #for my $lvid (@{ $pw->AtomLabel($atom) }) {
  #  my ($label, $value) = $lv->LabelValueToString($lvid);
  #  if ($label eq "revsersible") {
  #    if ($value eq "yes") {
  #      $reversible="true";
  #    } elsif ($value eq "no") {
  #      $reversible = "false";
  #    }
  #  }
  #  $self->Label($label, $value);
  #}
  #Exdent();
  #$self->PrLine("</annotation>");

  my (@reactants, @products, @modifiers);
  for my $edge (sort numerically @{ $pw->Edges($atom) }) {
    my $edgetype = $pw->EdgeType($atom, $edge);
    if ($lv->IsA($edgetype,
        $lv->StringToLabelValue("edge-type", "input"))) {
      push @reactants, $edge;
    } elsif ($lv->IsA($edgetype,
        $lv->StringToLabelValue("edge-type", "output"))) {
      push @products, $edge;
    } elsif ($lv->IsA($edgetype,
        $lv->StringToLabelValue("edge-type", "agent"))) {
      push @modifiers, $edge;
    } elsif ($lv->IsA($edgetype,
        $lv->StringToLabelValue("edge-type", "inhibitor"))) {
      push @modifiers, $edge;
    }
  }

  if (@reactants) {
    $self->PrLine("<listOfReactants>");
    Indent();
    for my $edge (@reactants) {
      my $molid = $pw->EdgeMol($atom, $edge);
      my @labels;
      for my $lv (@{ $pw->EdgeLabel($atom, $edge) }) {
        push @labels, $lv;
      }
      for my $lv (@{ $pw->MolLabel($atom, $edge) }) {
        push @labels, $lv;
      }
      $self->PrLine("<speciesReference species=\"" . "s" .
          $pw->MolInstId($atom, $edge) . "\"/>");
##      $self->Component("InteractionComponent", $molid, \@labels);
    }
    Exdent();
    $self->PrLine("</listOfReactants>");
  }

  if (@products) {
    $self->PrLine("<listOfProducts>");
    Indent();
    for my $edge (@products) {
      my $molid = $pw->EdgeMol($atom, $edge);
      my @labels;
      for my $lv (@{ $pw->EdgeLabel($atom, $edge) }) {
        push @labels, $lv;
      }
      for my $lv (@{ $pw->MolLabel($atom, $edge) }) {
        push @labels, $lv;
      }
      $self->PrLine("<speciesReference species=\"" . "s" .
          $pw->MolInstId($atom, $edge) . "\"/>");
##      $self->Component("InteractionComponent", $molid, \@labels);
    }
    Exdent();
    $self->PrLine("</listOfProducts>");
  }

  if (@products) {
    $self->PrLine("<listOfModifiers>");
    Indent();
    for my $edge (@modifiers) {
      my $molid = $pw->EdgeMol($atom, $edge);
      my @labels;
      for my $lv (@{ $pw->EdgeLabel($atom, $edge) }) {
        push @labels, $lv;
      }
      for my $lv (@{ $pw->MolLabel($atom, $edge) }) {
        push @labels, $lv;
      }
      $self->PrLine("<modifierSpeciesReference species=\"" . "s" .
          $pw->MolInstId($atom, $edge) . "\"/>");
##      $self->Component("InteractionComponent", $molid, \@labels);
    }
    Exdent();
    $self->PrLine("</listOfModifiers>");
  }

  Exdent();
  $self->PrLine("</reaction>");
}

######################################################################
sub PathwayList {
  my ($self) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};

  for my $pid (@{ $pw->PathwayId }) {
    $self->Pathway($pid);
  }

}

######################################################################
sub PathwayAtomMap {
  my ($self) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};

  for my $atom (@{ $pw->Atoms }) {
    for my $pid (@{ $pw->AtomPathway($atom) }) {
      $self->{pathwayatom}{$pid}{$atom} = 1;
    }
  }
}

######################################################################
sub Pathway {
  my ($self, $id) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};

  my ($name, $org, $extid) = 
      ($pw->PathwayName($id), $pw->PathwayOrg($id),
      $pw->PathwayExId($id));
  $self->PrLine("<Pathway id=\"$id\">");
  Indent();
  $self->PrLine("<Organism>$org</Organism>");
  $self->PrLine("<Name>$name</Name>");
  $self->PrLine("<ShortName>$extid</ShortName>");
  if (defined $self->{pathwayatom}{$id}) {
    for my $atom (keys %{ $self->{pathwayatom}{$id} }) {
      $self->Component("PathwayComponent", $atom, undef);
    }
  }
  Exdent();
  $self->PrLine("</Pathway>");
}

######################################################################
sub Locations {
  my ($self) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};

  my ($a, $e, $ml, %locations, $label, $value);

  for $a (@{ $pw->Atoms() }) {
    for $e (@{ $pw->Edges($a) }) {
      for $ml (@{ $pw->MolLabel($a, $e) }) {
        ($label, $value) = $lv->LabelValueToString($ml);
        if ($label eq "location") {
          $locations{$value} = 1;
        }
      }
    }
  }

  $self->PrLine("<listOfCompartments>");
  Indent();
  for $value (keys %locations) {
    $value =~ s/\s+/_/;
    $self->PrLine("<compartment metaid=\"$value\" id=\"$value\" name=\"$value\"/>");
  }
  $self->PrLine("<compartment metaid=\"cell_part\" id=\"cell_part\" name=\"cell_part\"/>");
  Exdent();
  $self->PrLine("</listOfCompartments>");
}

######################################################################
sub PrSBML {
  my ($self) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};

  $self->PathwayAtomMap();

  $self->PrLine("<?xml version=\"1.0\" encoding=\"UTF-8\"?>");
  $self->PrLine("<sbml xmlns=\"http://www.sbml.org/sbml/level2/version4\" level=\"2\" version=\"4\">");
  Indent(); 
  $self->PrLine("<model>");
  Indent();
  # $self->PrLine("<annotation>");
  # $self->PathwayList();
  # $self->PrLine("</annotation>");
  $self->Locations();
  $self->MoleculeList();
  $self->InteractionList();
  Exdent();
  $self->PrLine("</model>");
  Exdent();
  $self->PrLine("</sbml>");
}


######################################################################
1;
######################################################################


