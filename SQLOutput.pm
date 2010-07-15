#!/usr/local/bin/perl
package SQLOutput;
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

my %mol_type_abbrev = (
  "protein"  => "PR",
  "compound" => "CM",
  "complex"  => "CX",
  "rna"      => "RN",
  "molecule-type" => "MO"
);

######################################################################
sub new {
  my ($self, $pw, $lv, $outfh, $what) = @_;
  my $x = {};
  $x->{pw} = $pw;
  $x->{lv} = $lv;
  $x->{output_fh} = $outfh;
  $x->{what} = "table";
  $x->{what} = $what;
  $x->{date} = time();
  return bless $x;
}

######################################################################
sub SetExistingDefs {
  my ($self, $list) = @_;
  $self->{existingdef} = $list;
}

######################################################################
sub SetPTMExprId {
  my ($self, $id) = @_;
  $self->{ptm_expr_id} = $id;
}

######################################################################
sub  numerically { $a <=> $b };

######################################################################
sub TimeStamp {
  my ($sec, $min, $hr, $mday, $mon, $year, $wday, $yday, $isdst) =
      localtime(time);
  my $month = ('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug',
      'Sep', 'Oct', 'Nov', 'Dec')[$mon];
  $year = $year + 1900;
#  return sprintf "%d_%2.2d_%2.2d %2.2d:%2.2d::%2.2d",
#      $year, $mon+1, $mday, $hr, $min, $sec;
  return sprintf "%2.2d%2.2d%2.2d",
      $mon+1, $mday, $year;
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
  my ($self, $table_name, @values) = @_;
  my $fh = $self->{output_fh};
  my $line;
  if ($self->{what} eq "table") {
    $line = join("\t", $table_name, @values);
  } elsif ($self->{what} eq "sql") {
    $line = "insert into $table_name values ('" .
        join("', '", @values) . "');";
  }
  if ($fh) {
    print $fh "$line\n";
  }
  push @{ $self->{lines} }, $line;
}

######################################################################
sub Ontology {
  my ($self) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};

  my $used_labels = $pw->AllUsedLabels();

  for my $lvid (keys %{ $used_labels }) {
    for my $x (@{ $lv->AncestorsOf($lvid) }) {
      $$used_labels{$x} = 1;
    }
  }

  my %temp;

  my $lvid_index = $lv->AllNamePairs();
  for my $pair (keys %{ $lvid_index }) {
    my ($sort_name, $label_value_name) = split("\t", $pair);
    if (defined $$used_labels{$$lvid_index{$pair}}) {
      push @{ $temp{$sort_name} }, join("\t", $label_value_name,
          $$lvid_index{$pair});
    }
  }

  undef $lvid_index;
  for my $label_id (@{ $lv->AllSorts }) {
    my $sort_name = $lv->LabelString($label_id);
    if (defined $temp{$sort_name}) {
      $self->PrLine("pw_label", $label_id, CleanString($sort_name));
      for my $pair (sort @{ $temp{$sort_name} }) {
        my ($label_value_name, $lvid) = split("\t", $pair);
        my $parent_id = $lv->ParentOf($lvid);

        $self->PrLine("pw_label_value", $lvid, $label_id,
            CleanString($label_value_name));
        $self->PrLine("pw_label_value_parent", $lvid, $parent_id);
        for my $g (@{ $lv->GOTermsFor($lvid) }) {
          if ($g =~ /^\d+$/) {
            $g = "GO:$g";
          }
          $self->PrLine("pw_ext_lv_id", $lvid, $g, "GO");
        }
      }
    }
  }

}

######################################################################
sub MoleculeList {
  my ($self) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};

  my %seen;

#  for my $molid (@{ $pw->UsedMols() }) {
  for my $molid (@{ $pw->Mols() }) {
    if ($molid ne "") {
      $self->Molecule($molid);
      $seen{$molid} = 1;
    }
  }

  for my $molid (keys %{ $self->{family_member_list} }) {
    if (! $seen{$molid}) {
      $self->Molecule($molid);
      $seen{$molid} = 1;
    }
  }

  for my $molid (keys %{ $self->{whole_molecule_list} }) {
    if (! $seen{$molid}) {
      $self->Molecule($molid);
      $seen{$molid} = 1;
    }
  }

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

  if (defined $self->{existingdef}{$id}) {
    return;
  }

  my $mol_id = $id;
  my ($label, $mol_type) = $lv->LabelValueToString($pw->MolType($mol_id));

  $self->PrLine("pw_mol", $mol_id, $mol_type_abbrev{$mol_type}, "Hs");

  my $mehash = $pw->MolExId($mol_id);
  for my $id_type (keys %{ $mehash }) {
    for my $i (keys %{ $$mehash{$id_type} }) {
      $self->PrLine("pw_ext_mol_id", $mol_id, $i, $id_type);
    }
  }

  my $mnhash = $pw->MolName($mol_id);
  for my $name_type (keys %{ $mnhash }) {
    for my $i (keys %{ $$mnhash{$name_type} }) {
      $self->PrLine("pw_mol_name", $mol_id, $i, $name_type);
    }
  }

  if ($mol_type eq "complex") {
    for my $comp (sort numerically @{ $pw->Components($mol_id) }) {
      my $compmol = $pw->ComponentMol($mol_id, $comp);
      if (! defined $self->{molecules_defined}{$compmol} ) {
        $self->{molecules_to_define}{$compmol} = 1;
      }      
      for my $lvid (@{ $pw->ComponentLabel($mol_id, $comp) }) {
        my ($name, $value) = $lv->LabelValueToString($lvid);
        if ($name && $name ne "molecule-type") {
          $self->PrLine("pw_component_labeling", $mol_id, $comp, $lvid);
        }
      }
      $self->PrLine("pw_complex_component", $mol_id, $comp, $compmol);
      for my $ptm (@{ $pw->ComponentPTM($mol_id, $comp) }) {
        my ($uniprot, $pos, $residue, $ptm_lvid, $ptm_name) = @{ $ptm };
        my ($operation, $left, $right) = ("", "", "");
        $self->{ptm_expr_id}++;
        $self->PrLine("pw_component_modification",
            $mol_id, $comp, $self->{ptm_expr_id});
        $self->PrLine("pw_ptm_expression",
            $self->{ptm_expr_id}, $uniprot, $ptm_lvid, $residue,
            $pos, $operation, $left, $right);
      }

    }
  }

  if (defined $pw->FamilyChildren($id)) {
    for my $c (@{ $pw->FamilyChildren($id) }) {
      $self->PrLine("pw_family_member", $id, $c);
      $self->{family_member_list}{$c} = $id;
    }
  }

  if (defined $pw->MolWhole($id)) {
    for my $w (@{ $pw->MolWhole($id) }) {
      my ($stop, $end) = split(",", $pw->PartBounds($id));
      $self->PrLine("pw_mol_part", $id, $w, $stop, $end);
      $self->{whole_molecule_list}{$w} = $id;
    }
  }

  for my $m (keys %{ $self->{molecules_to_define} }) {
    $self->Molecule($m);
  }
}

######################################################################
sub InteractionList {
  my ($self) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};

  for my $atom (@{ $pw->Atoms() }) {
    $self->Interaction($atom);
  }
}

######################################################################
sub SubnetList {
  my ($self) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};

  if (@{ $pw->Subnets() } == 0) {
    return;
  }

  $self->PrLine("<SubnetList>");
  Indent();
  for my $subnet (@{ $pw->Subnets() }) {
    $self->PrLine("<Subnet id=\"$subnet\">");
    Indent();
    $self->PrLine("<SubnetComponentList>");
    Indent();
    for my $atom (@{ $pw->AtomsInSubnet($subnet) }) {
      $self->PrLine("<SubnetComponent idref=\"$atom\" />");
    }
    Exdent();
    $self->PrLine("</SubnetComponentList>");
    Exdent();
    $self->PrLine("</Subnet>");
  }
  Exdent();
  $self->PrLine("</SubnetList>");
}

######################################################################
sub Source {
  my ($self, $name) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};

  for my $src_id (keys %{ $self->{allsources} }) {
    my $name = $pw->SourceName($src_id);
    $self->PrLine("pw_source", $src_id, $name, "");
  }
}

######################################################################
sub Interaction {
  my ($self, $id) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};

  my $src_id = $pw->AtomSource($id);
  $self->{allsources}{$src_id} = 1;
  $self->PrLine("pw_atom", $id, "Hs", $src_id);

  my $atom = $id;
  for my $lvid (@{ $pw->AtomCondition($atom) }) {
    my ($label, $value) = $lv->LabelValueToString($lvid);
    $self->PrLine("pw_atom_condition", $atom, $lvid);
  }
  for my $lvid (@{ $pw->AtomLabel($atom) }) {
    $self->PrLine("pw_atom_label", $atom, $lvid);
  }
  for my $edge (sort numerically @{ $pw->Edges($atom) }) {
    my $molid = $pw->EdgeMol($atom, $edge);
    $self->PrLine("pw_edge", $atom, $edge, $molid);
    for my $lvid (@{ $pw->EdgeLabel($atom, $edge) }) {
      $self->PrLine("pw_edge_label", $atom, $edge, $lvid);
    }
    my $edge_type = $pw->EdgeType($atom, $edge);
    $self->PrLine("pw_edge_label", $atom, $edge, $edge_type);

    for my $lvid (@{ $pw->MolLabel($atom, $edge) }) {
      $self->PrLine("pw_mol_label", $atom, $edge, $lvid);
    }
    for my $ptm (@{ $pw->EdgePTM($atom, $edge) }) {
      my ($uniprot, $pos, $residue, $ptm_lvid, $ptm_name) = @{ $ptm };
      my ($operation, $left, $right) = ("", "", "");
      $self->{ptm_expr_id}++;
      $self->PrLine("pw_edge_modification", $atom, $edge,
          $self->{ptm_expr_id});
      $self->PrLine("pw_ptm_expression",
          $self->{ptm_expr_id}, $uniprot, $ptm_lvid, $residue,
          $pos, $operation, $left, $right);
    }
  }
  if ($pw->AtomType($atom) eq
      $lv->StringToLabelValue("process-type", "pathway") ||
      $pw->AtomType($atom) eq
      $lv->StringToLabelValue("process-type", "subnet")) {
##    my $pid = $pw->Abstraction($atom);
    my $pid = $pw->AbstractionId($atom);
    my $ext_id = $pw->AbstractionExtId($atom);
    $self->PrLine("pw_abstraction", $atom, $pid, $ext_id);
  }

  for my $evid (@{ $pw->AtomEvidence($atom) }) {
    my ($pathway, $edge, $mol, $src) = ("", "", "", "");
    $self->PrLine("pw_evidence",
        $pathway, $atom, $edge, $mol, $evid, $src);
  }

  for my $ref (@{ $pw->AtomReferences($atom) }) {
    my ($pathway, $edge, $mol, $src) = ("", "", "", "");
    my ($pmid, $text) = split(/\t/, $ref);
    $self->PrLine("pw_references",
        $pathway, $atom, $edge, $mol, $pmid, $text, $src);
  }

  for my $note (@{ $pw->AtomNotes($atom) }) {
    my ($pathway, $edge, $mol, $src) = ("", "", "", "");
    $self->PrLine("pw_notes",
        $pathway, $atom, $edge, $mol, $note, $src);
  }

}

######################################################################
sub PathwayList {
  my ($self) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};

  if (@{ $pw->PathwayId }) {
    for my $pid (@{ $pw->PathwayId }) {
      $self->Pathway($pid);
    }
  }

}

######################################################################
sub PathwayAtomMap {
  my ($self) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};

  for my $atom (@{ $pw->Atoms() }) {
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

  if (! defined $self->{pathwayatom}{$id} ||
      keys %{ $self->{pathwayatom}{$id} } < 1) {
    return;
  }
  my ($name, $org, $extid, $src_id) = 
      (CleanString($pw->PathwayName($id)), $pw->PathwayOrg($id),
      CleanString($pw->PathwayExId($id)), $pw->PathwaySource($id));
  my $last_updated = $self->{date};
  my $is_subnet = $pw->PathwayIsSubnet($id) == 1 ? "Y" : "N";
  $self->{allsources}{$src_id} = 1;
  $self->PrLine("pw_pathway", $id,
      $name, $org, $src_id, $extid, $is_subnet, $last_updated);
  if (defined $self->{pathwayatom}{$id}) {
    for my $atom (keys %{ $self->{pathwayatom}{$id} }) {
      $self->PrLine("pw_pathway_atom", $id, $atom);
    }
  }
  for my $c (@{ $pw->Curators($id) }) {
    $self->PrLine("pw_curators", $id, $c, $src_id, 'C');
  }
  ## use pw_curators for both reviewer and curator
  for my $r (@{ $pw->Reviewers($id) }) {
    $self->PrLine("pw_curators", $id, $r, $src_id, 'R');
  }

  for my $ref (@{ $pw->PathwayReferences($id) }) {
    my ($atom, $edge, $mol, $src) = ("", "", "", "");
    my ($pmid, $text) = split(/\t/, $ref);
    $self->PrLine("pw_references",
        $id, $atom, $edge, $mol, $pmid, $text, $src);
  }


}

######################################################################
sub PrSQL {
  my ($self) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};

  $self->Ontology();
  $self->PathwayAtomMap();
  $self->MoleculeList();
  $self->InteractionList();
  $self->SubnetList();
  $self->PathwayList();  
  $self->Source();
}


######################################################################
1;
######################################################################


