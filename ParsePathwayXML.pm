

# Copyright SRA International
#
# Distributed under the OSI-approved BSD 3-Clause License.
# See http://ncip.github.com/pathway-interaction-database/LICENSE.txt for details.


package ParsePathwayXML;

require Exporter;
@ISA = qw(Exporter);
#@EXPORT = qw(
#  CreateIsA
#  CreateLabelSort
#);

use strict;
use XML::Parser;
use Pathway;
use PWLabel;
use DBI;
use XMLOutput;

my %basic_mol_type = (
  "protein"  => "PR",
  "complex"  => "CX",
  "compound" => "CM",
  "rna"      => "RN",
  "molecule-type" => "MO"
);

my %mol_name_types = (
  "AS" => 1,
  "OF" => 1,
  "PF" => 1
);

my %mol_ext_id_types = (
  "LL" => 1,
  "CA" => 1,
  "EC" => 1,
  "UP" => 1,
  "KG" => 1,
  "GO" => 1
);

my %mol_label_types = (
  "activity-state" => 1,
  "location"       => 1
);

my %edge_label_types = (
  "function"       => 1
);

my %atom_label_types = (
  "reversible"     => 1
);

my ($lv, $pw);
my (%current_elems);
my (@elem_stack, @attr_stack, %char_value, %id_def);
my (%mol_labels, %atom_labels, %edge_labels, %atom_conditions);
my (@ptm_expression);
my (%temp_label_sort, %temp_child_parent);
my $p;

######################################################################
sub parse {
  my ($self, $fh) = @_;
  $p->parse($fh);
}

######################################################################
sub attr_top {
  my ($x) = @_;
  if (!$x) {
    $x = 0;
  }
  return $attr_stack[$x];
}

######################################################################
sub elem_top {
  my ($x) = @_;
  if (!$x) {
    $x = 0;
  }
  return $elem_stack[$x];
}

######################################################################
sub handle_attrs {
  my (@attrs) = @_;
  my %pairs;
  while (@attrs) {
    my $attr_name  = shift @attrs;
    my $attr_value = shift @attrs;
    $pairs{$attr_name} = $attr_value;
  }
  return \%pairs;
}

######################################################################
sub handle_start {
  my ($p, $elem, @attrs) = @_;

  unshift @elem_stack, $elem;

  ## the business about undef'ing the upcoming char value is to
  ## get around a problem with the library xml parser: if it encounters
  ## &3x27; in a char value in an entity (but not in an attribute value)
  ## it will create multiple calls to handle_char; so we need to accumulate
  ## the string, after an initial undef

  undef $char_value{elem_top()};

  unshift @attr_stack, handle_attrs(@attrs);
  my $attr_top = attr_top();
  for my $x (keys %{ $attr_top }) {
#    print "$x=$$attr_top{$x}\n";
    if ($x eq "id") {
      $id_def{$elem} = $$attr_top{id};
    }
  }
  if ($elem eq "ComplexComponent" || $elem eq "InteractionComponent") {
    $id_def{$elem}++;
  }
}

######################################################################
sub handle_end {
  my ($p, $elem) = @_;

  my ($attr_top, $elem_top) = (attr_top(), elem_top());

  if ($elem eq "Model") {

  ############################################################

  } elsif ($elem eq "Molecule") {
    $pw->AddMol($$attr_top{id}, $lv->BasicMolTypeCode(
        $basic_mol_type{$$attr_top{molecule_type}}));

  } elsif ($elem eq "MoleculeList") {

  } elsif ($elem eq "MoleculeType") {

  } elsif ($elem eq "ComplexComponent") {
    my $newComponent = $id_def{ComplexComponent};
    my $newMolecule = $id_def{Molecule};
    $pw->AddComponent($id_def{Molecule}, $id_def{ComplexComponent},
        $$attr_top{molecule_idref});
    for my $ptmexpr (@ptm_expression) {
      my ($uniprot_acc, $pos, $residue, $mod_type, $value) =
          @{ $ptmexpr };
      $pw->AddComponentPTM($id_def{Molecule}, $id_def{ComplexComponent},
          $uniprot_acc, $pos, $residue, $value, $mod_type);
    }
    undef @ptm_expression;

  } elsif ($elem eq "ComplexComponentList") {
    undef $id_def{ComplexComponent};

  } elsif ($elem eq "Family") {
    $pw->AddFamilyMember($$attr_top{member_molecule_idref},
        $$attr_top{family_molecule_idref})

  } elsif ($elem eq "FamilyMemberList") {

  } elsif ($elem eq "Member") {
    $pw->AddFamilyMember($$attr_top{member_molecule_idref}, $id_def{Molecule});
    for my $ptmexpr (@ptm_expression) {
      my ($uniprot_acc, $pos, $residue, $mod_type, $value) =
          @{ $ptmexpr };
      $pw->AddMemberPTM($id_def{Molecule}, $attr_stack[0]{member_molecule_idref},
          $uniprot_acc, $pos, $residue, $value, $mod_type);
    }
    undef @ptm_expression;

  } elsif ($elem eq "Part") {
    $pw->AddMolPart($$attr_top{part_molecule_idref},
        $$attr_top{whole_molecule_idref}, $$attr_top{start}, $$attr_top{end});

  ############################################################

  } elsif ($elem eq "Interaction") {
    my $interaction_type = $$attr_top{interaction_type};
    my $interaction_value = $lv->StringToLabelValue("process-type",
        $interaction_type);
    if ($interaction_value == 0) {
      $interaction_value = $lv->StringToLabelValue("function",
          $interaction_type);
    }
    if ($interaction_value == 0) {
      print STDERR "no label value for interaction type $interaction_type\n";
    } else {
      $pw->AddAtom($id_def{Interaction}, $interaction_value);
      $pw->AddAtomSource($id_def{Interaction}, $attr_stack[0]{_Source});
    }
    undef %atom_labels;

  } elsif ($elem eq "PositiveCondition") {
    my $cond_type = $$attr_top{condition_type};
    my $cond_value = $lv->StringToLabelValue("process-type", $cond_type);
    if ($cond_value == 0) {
      $cond_value = $lv->StringToLabelValue("function", $cond_type);
    }
    if ($cond_value == 0) {
      print STDERR "no label value for positivecondition $cond_type\n";
    } else {
      $pw->AddAtomCondition($id_def{Interaction}, $cond_value);
    }
  } elsif ($elem eq "NegativeCondition") {
    my $cond_type = $$attr_top{condition_type};
    my $cond_value = $lv->StringToLabelValue("process-type", $cond_type);
    if ($cond_value == 0) {
      $cond_value = $lv->StringToLabelValue("function", $cond_type);
    }
    if ($cond_value == 0) {
      print STDERR "no label value for negativecondition $cond_type\n";
    } else {
      $pw->AddAtomNegativeCondition($id_def{Interaction}, $cond_value);
    }
 
  }
    elsif ($elem eq "RoleType") {

  } elsif ($elem eq "InteractionType") {

  } elsif ($elem eq "InteractionComponent") {
    $pw->AddEdge($id_def{Interaction}, $id_def{InteractionComponent},
        $lv->StringToLabelValue("edge-type",
        $$attr_top{role_type}), $$attr_top{molecule_idref});
    for my $ptmexpr (@ptm_expression) {
      my ($uniprot_acc, $pos, $residue, $mod_type, $value) =
          @{ $ptmexpr };
      $pw->AddEdgePTM($id_def{Interaction}, $id_def{InteractionComponent},
          $uniprot_acc, $pos, $residue, $value, $mod_type);
    }
    undef @ptm_expression;
    undef %mol_labels;
    undef %edge_labels;

  } elsif ($elem eq "InteractionComponentList") {
    undef $id_def{InteractionComponent};

  ############################################################

  } elsif ($elem eq "PTMExpression") {

  } elsif ($elem eq "PTMTerm") {
    my $uniprot_acc = $attr_stack[0]{"protein"};
    my $pos         = $attr_stack[0]{"position"};
    my $residue     = $attr_stack[0]{"aa"};
    my $mod_type    = $attr_stack[0]{"modification"};
    my $value  = $lv->StringToLabelValue("ptm",$mod_type);
    push @ptm_expression,
        [ $uniprot_acc, $pos, $residue, $mod_type, $value ];

  ############################################################

  } elsif ($elem eq "Subnet") {

  } elsif ($elem eq "SubnetComponentList") {

  } elsif ($elem eq "SubnetComponent") {
    $pw->AddSubnetAtom($id_def{Subnet}, $$attr_top{interaction_idref})

  ############################################################

  } elsif ($elem eq "Abstraction") {
    $pw->AddAbstraction($id_def{Interaction},
        $$attr_top{pathway_idref},
        $$attr_top{pathway_name},
        $$attr_top{external_pathway_id});

  ############################################################

  } elsif ($elem eq "Source") {
    if ($char_value{Source} eq "NCI-Nature Curated") {
      $char_value{Source} = "NATURE";
    }
    $pw->AddSource($$attr_top{id}, $char_value{Source});
    $attr_stack[1]{_Source} = $attr_stack[0]{id};
#    $pw->AddSource($char_value{Source}, $char_value{Source});

  } elsif ($elem eq "Evidence") {
    $pw->AddAtomEvidence($id_def{Interaction}, $char_value{Evidence});

  } elsif ($elem eq "Reference") {
    if ($elem_stack[2] eq "Interaction") {
      $pw->AddAtomReferences($id_def{Interaction}, $char_value{Reference});
    } elsif ($elem_stack[2] eq "Pathway") {
      $pw->AddPathwayReferences($id_def{Pathway}, $char_value{Reference});
    }

  } elsif ($elem eq "Note") {
    $pw->AddAtomNotes($id_def{Interaction}, $char_value{Note});

  } elsif ($elem eq "Curator") {
    $pw->AddCurator($id_def{Pathway}, $char_value{Curator});

  } elsif ($elem eq "Reviewer") {
    $pw->AddReviewer($id_def{Pathway}, $char_value{Reviewer});

  } elsif ($elem eq "Pathway") {
    $pw->SetPathwaySrcId($id_def{Pathway}, $attr_stack[0]{_Source});
    $pw->AddPathwaySource($id_def{Pathway}, $attr_stack[0]{_Source});
    $pw->SetPathwayName($id_def{Pathway}, $char_value{LongName});
    $pw->SetPathwayExId($id_def{Pathway}, $char_value{ShortName});
    $pw->SetPathwayOrg($id_def{Pathway}, $char_value{Organism});
    $pw->SetIsSubnet($id_def{Pathway}, $attr_stack[0]{subnet} eq "true" ?
        1 : 0);
  } elsif ($elem eq "PathwayComponent") {
    $pw->AddAtomPathway($$attr_top{interaction_idref}, $id_def{Pathway})

  } elsif ($elem eq "PathwayComponentList") {


  } elsif ($elem eq "Organism") {
  } elsif ($elem eq "LongName") {
  } elsif ($elem eq "ShortName") {

  } elsif ($elem eq "Name") {
    my $name_type  = $$attr_top{name_type};
    my $name_value = $$attr_top{value};
    if (defined $mol_name_types{$name_type}) {
      $pw->AddMolName($id_def{Molecule}, $name_type, $name_value);
    } elsif (defined $mol_ext_id_types{$name_type}) {
      $pw->AddMolExId($id_def{Molecule}, $name_type, $name_value);
    }
  } elsif ($elem eq "Label") {
    my $label_type  = $$attr_top{label_type};
    my $label_value = $$attr_top{value};
    if (defined $edge_label_types{$label_type}) {
      $pw->AddEdgeLabel($id_def{Interaction}, $id_def{InteractionComponent},
          $lv->StringToLabelValue($label_type, $label_value));
    } elsif (defined $mol_label_types{$label_type}) {
      if (elem_top(1) eq "InteractionComponent") {
        $pw->AddMolLabel($id_def{Interaction}, $id_def{InteractionComponent},
            $lv->StringToLabelValue($label_type, $label_value));
      } elsif (elem_top(1) eq "ComplexComponent") {
        $pw->AddComponentLabel($id_def{Molecule},
            $id_def{ComplexComponent},
            $lv->StringToLabelValue($label_type, $label_value));
      } elsif (elem_top(1) eq "Member") {
        $pw->AddMemberLabel($id_def{Molecule},
            $attr_stack[1]{member_molecule_idref},
            $lv->StringToLabelValue($label_type,$label_value));
      }
    } elsif (defined $atom_label_types{$label_type}) {
        $pw->AddAtomLabel($id_def{Interaction},
            $lv->StringToLabelValue($label_type, $label_value));
    }

  ############################################################

  } elsif ($elem eq "Ontology") {

    $lv->CreateLabelSort();
    $lv->CreateIsA();

  } elsif ($elem eq "LabelType") {

    $lv->AddSort($$attr_top{"name"}, $$attr_top{"id"});

  } elsif ($elem eq "LabelValueList") {

  } elsif ($elem eq "LabelValue") {

    $temp_label_sort{$$attr_top{"id"}}      = $id_def{"LabelType"};
    $temp_child_parent{$$attr_top{"id"}}{$$attr_top{"parent_idref"}} = 1;
    $lv->AddParent($$attr_top{"id"}, $$attr_top{"parent_idref"});
    my $label_type_attrs = attr_top(2);
    my $label_type_name  = $$label_type_attrs{"name"};
    $lv->AddNamePair($$attr_top{"id"}, $$attr_top{"name"},
        $label_type_name);
    if ($label_type_name eq "molecule-type") {
      my $mol_type_name = $$attr_top{"name"};
      if ($mol_type_name eq "protein") {
        $lv->AddBasicMolType("PR", $$attr_top{"id"});
      } elsif ($mol_type_name eq "compound") {
        $lv->AddBasicMolType("CM", $$attr_top{"id"});
      } elsif ($mol_type_name eq "complex") {
        $lv->AddBasicMolType("CX", $$attr_top{"id"});
      } elsif ($mol_type_name eq "rna") {
        $lv->AddBasicMolType("RN", $$attr_top{"id"});
      } elsif ($mol_type_name eq "molecule-type") {
        $lv->AddBasicMolType("MO", $$attr_top{"id"});
      }
    }

    if ($$attr_top{"GO"}) {
      for my $go (split(",", $$attr_top{"GO"})) {
        $lv->AddGOTerm($$attr_top{"id"}, $go);
      }
    }

  }

  shift @elem_stack;
  shift @attr_stack;

}

######################################################################
sub handle_char {
  my ($p, $s) = @_;

  $char_value{elem_top()} .= $s;
}

######################################################################
sub new {
  my ($self, $labeldb, $pathway) = @_;
  $lv = $labeldb;
  $pw = $pathway;
  $p = new XML::Parser(Handlers => {
    Start   => \&handle_start,
    End     => \&handle_end,
    Char    => \&handle_char
  } );
  return bless {};
}

######################################################################
1;
