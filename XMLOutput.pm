#!/usr/local/bin/perl
package XMLOutput;
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

my %name_types = (

  "AS"           => "alias",
  "OF"           => "official symbol",
  "PF"           => "preferred symbol",

  "CA"           => "Chemical Abstracts",
  "CH"           => "Chemical Entities of Biological Interest",
  "EC"           => "Enzyme Consortium",
  "GO"           => "Gene Ontology",
  "KG"           => "KEGG",
  "LL"           => "EntrezGene",
  "UP"           => "UniProt"

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
sub TimeStamp {
  my ($sec, $min, $hr, $mday, $mon, $year, $wday, $yday, $isdst) =
      localtime(time);
  my $month = ('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug',
      'Sep', 'Oct', 'Nov', 'Dec')[$mon];
  $year = $year + 1900;
  return sprintf "%d_%2.2d_%2.2d %2.2d:%2.2d::%2.2d",
      $year, $mon+1, $mday, $hr, $min, $sec;
}

######################################################################
sub CleanString {
  my ($s) = @_;

  $s =~ s/<//g;
  $s =~ s/>//g;
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
sub Ontology {
  my ($self) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};

# <Ontology>
#   <LabelType name= id= >
#   <LabelValueList>
#     <LabelValue name= id= parent_idref= >
#   </LabelValueList>
# </Ontology>

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

  $self->PrLine("<Ontology>");
  Indent();
  for my $label_id (@{ $lv->AllSorts }) {
    my $sort_name = $lv->LabelString($label_id);
    if (defined $temp{$sort_name}) {
      $self->PrLine("<LabelType name=\"" . CleanString($sort_name) .
          "\" id=\"$label_id\">");
      Indent();
      $self->PrLine("<LabelValueList>");

      Indent();
      for my $pair (sort @{ $temp{$sort_name} }) {
        my ($label_value_name, $lvid) = split("\t", $pair);
        my $parent_id = $lv->ParentOf($lvid);

        my $temp = $lv->GOTermsFor($lvid);
        my $go_ids = "";
        if (@{ $temp }) {
          $go_ids = join(",", map { if (/^GO:/) { $_ } else {"GO:$_" } }
              @{ $temp });
        }

        $self->PrLine("<LabelValue name=\"" .
            CleanString($label_value_name) . "\" " .
            "id=\"$lvid\" parent_idref=\"$parent_id\" " .
            ($go_ids ? "GO=\"$go_ids\" " : "") .
            "/>")
      }
      Exdent();
      $self->PrLine("</LabelValueList>");
      Exdent();
      $self->PrLine("</LabelType>");
    }
  }
  Exdent();
  $self->PrLine("</Ontology>");
}

######################################################################
sub MoleculeList {
  my ($self) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};

  my %seen;

  $self->PrLine("<MoleculeList>");
  Indent();
## if we don't process all mols, used or not, we will not
## recover whole/part info for unused parts
  for my $molid (@{ $pw->UsedMols() }) {
#  for my $molid (@{ $pw->Mols() }) {
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

  Exdent();
  $self->PrLine("</MoleculeList>");
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

  my $mol_id = $id;
  my ($label, $mol_type) = $lv->LabelValueToString($pw->MolType($mol_id));
  $self->PrLine("<Molecule molecule_type=\"$mol_type\" id=\"$id\">");
  Indent();
##  $self->PrLine("<MoleculeType>$mol_type</MoleculeType>");
#  $self->Label($label, $mol_type);

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
    $self->PrLine("<ComplexComponentList>");
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
      my $ptms = $pw->ComponentPTM($mol_id, $comp);
      $self->ComplexComponent($compmol, \@labels, $ptms);
    }
    Exdent();
    $self->PrLine("</ComplexComponentList>");

  }

  if (@{ $pw->FamilyChildren($id) }) {
    $self->PrLine("<FamilyMemberList>");
    Indent();
    for my $m (@{ $pw->FamilyChildren($id) }) {
      my @labels;
      for my $lvid (@{ $pw->MemberLabel($id, $m) }) {
        my ($name, $value) = $lv->LabelValueToString($lvid);
        if ($name && $name ne "molecule-type") {
          push @labels, $lvid;
        }
      }
      my $ptms = $pw->MemberPTM($id, $m);
      $self->FamilyMember($m, \@labels, $ptms);
#      $self->PrLine("<Family family_molecule_idref=\"$id\" " .
#          "member_molecule_idref=\"$m\" />");
      $self->{family_member_list}{$m} = $id;
      if (! defined $self->{molecules_defined}{$m} ) {
        $self->{molecules_to_define}{$m} = 1;
      }
    }
    Exdent();
    $self->PrLine("</FamilyMemberList>");
  }

  if (defined $pw->MolWhole($id)) {
    for my $w (@{ $pw->MolWhole($id) }) {
      my ($start, $end) = split(",", $pw->PartBounds($id));
      $self->PrLine("<Part whole_molecule_idref=\"$w\" " .
          "part_molecule_idref=\"$id\" start=\"$start\" end=\"$end\" />");
      $self->{whole_molecule_list}{$w} = $id;
    }
  }

  Exdent();
  $self->PrLine("</Molecule>");
  for my $m (keys %{ $self->{molecules_to_define} }) {
    $self->Molecule($m);
  }
}

######################################################################
sub FamilyMember {
  my ($self, $id, $labels, $ptms) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};

  $self->PrLine("<Member member_molecule_idref=\"$id\">");
  Indent();
  if ($labels && @{ $labels }) {
    for my $lvid (@{ $labels }) {
      my ($label, $value) = $lv->LabelValueToString($lvid);
      $self->Label($label, $value);
    }
  }
  if ($ptms && @{ $ptms }) {
    $self->PrLine("<PTMExpression>");
    Indent();
    for my $ptm (@{ $ptms } ){
      my ($protein_id, $position, $amino_acid, $modification_label_value,
        $modification_label_name) = @{ $ptm };
      $self->PrLine("<PTMTerm protein=\"$protein_id\" position=\"$position\" " .
          "aa=\"$amino_acid\" modification=\"$modification_label_name\" />");
    }
    Exdent();
    $self->PrLine("</PTMExpression>");
  }
  Exdent();
  $self->PrLine("</Member>");
}

######################################################################
sub ComplexComponent {
  my ($self, $id, $labels, $ptms) = @_;

  my $what = "ComplexComponent";
  my $pw = $self->{pw};
  my $lv = $self->{lv};

  $self->PrLine("<$what molecule_idref=\"$id\">");
  Indent();
  if ($labels && @{ $labels }) {
    for my $lvid (@{ $labels }) {
      my ($label, $value) = $lv->LabelValueToString($lvid);
      $self->Label($label, $value);
    }
  }
  if ($ptms && @{ $ptms }) {
    $self->PrLine("<PTMExpression>");
    Indent();
    for my $ptm (@{ $ptms } ){
      my ($protein_id, $position, $amino_acid, $modification_label_value,
        $modification_label_name) = @{ $ptm };
      $self->PrLine("<PTMTerm protein=\"$protein_id\" position=\"$position\" " .
          "aa=\"$amino_acid\" modification=\"$modification_label_name\" />");
    }
    Exdent();
    $self->PrLine("</PTMExpression>");
  }
  Exdent();
  $self->PrLine("</$what>");
}

######################################################################
sub InteractionComponent {
  my ($self, $id, $edgetype, $labels, $ptms) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};

  $self->PrLine("<InteractionComponent " .
      "role_type=\"$edgetype\" molecule_idref=\"$id\">");
  Indent();
##  $self->PrLine("<RoleType>$edgetype</RoleType>");
  if ($labels && @{ $labels }) {
    for my $lvid (@{ $labels }) {
      my ($label, $value) = $lv->LabelValueToString($lvid);
      $self->Label($label, $value);
    }
  }
  if (@{ $ptms }) {
    $self->PrLine("<PTMExpression>");
    Indent();
    for my $ptm (@{ $ptms } ){
      my ($protein_id, $position, $amino_acid, $modification_label_value,
        $modification_label_name) = @{ $ptm };
      $self->PrLine("<PTMTerm protein=\"$protein_id\" position=\"$position\" " .
          "aa=\"$amino_acid\" modification=\"$modification_label_name\" />");
    }
    Exdent();
    $self->PrLine("</PTMExpression>");
  }
  Exdent();
  $self->PrLine("</InteractionComponent>");
}

######################################################################
sub PathwayComponent {
  my ($self, $id, $labels) = @_;

  my $what = "PathwayComponent";
  my $pw = $self->{pw};
  my $lv = $self->{lv};

  if ($labels && @{ $labels }) {
    $self->PrLine("<$what interaction_idref=\"$id\">");
    Indent();
    for my $lvid (@{ $labels }) {
      my ($label, $value) = $lv->LabelValueToString($lvid);
      $self->Label($label, $value);
    }
    Exdent();
    $self->PrLine("</$what>");
  } else {
    $self->PrLine("<$what interaction_idref=\"$id\" />");
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

  my $long_name_type = $name_types{$nametype};
  $self->PrLine("<Name name_type=\"$nametype\" " .
      "long_name_type=\"$long_name_type\" value=\"" .
      CleanString($name) . "\" />");
}

######################################################################
sub InteractionList {
  my ($self, $pid) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};

  $self->PrLine("<InteractionList>");
  Indent();
  if ($pid) {
    for my $atom (keys %{ $self->{pathwayatom}{$pid} } ) {
      $self->Interaction($atom);
    }
  } else {
    for my $atom (@{ $pw->Atoms() }) {
      $self->Interaction($atom);
    }
  }
  Exdent();
  $self->PrLine("</InteractionList>");
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
  my ($self, $id) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};

  my $name = $pw->SourceName($id);
  if ($name eq "NATURE") {
    $name = "NCI-Nature Curated";
  }
  $self->PrLine("<Source id=\"$id\">$name</Source>")

}

######################################################################
sub Interaction {
  my ($self, $id) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};

  my ($label, $value) = $lv->LabelValueToString($pw->AtomType($id));
  $self->PrLine("<Interaction interaction_type=\"$value\" id=\"$id\">");
  Indent();
  if ($value eq "pathway" || $value eq "subnet") {
    my $pid = $pw->AbstractionId($id);
    $self->PrLine("<Abstraction " .
        "pathway_idref=\"$pid\" " .
        "pathway_name=\"" . CleanString($pw->PathwayName($pid)) . "\" " .
        "external_pathway_id=\"" . CleanString($pw->AbstractionExtId($id)) .
        "\" />");
  }
  $self->Source($pw->AtomSource($id));

  my $atom = $id;
  for my $lvid (@{ $pw->AtomCondition($atom) }) {
    my ($label, $value) = $lv->LabelValueToString($lvid);
    $self->PrLine("<PositiveCondition condition_type=\"$value\">$value</PositiveCondition>");
  }
  for my $lvid (@{ $pw->AtomNegativeCondition($atom) }) {
    my ($label, $value) = $lv->LabelValueToString($lvid);
    $self->PrLine("<NegativeCondition condition_type=\"$value\">$value</NegativeCondition>");
  }
  for my $lvid (@{ $pw->AtomLabel($atom) }) {
    my ($label, $value) = $lv->LabelValueToString($lvid);
    if ($label ne "process-type") {
      $self->Label($label, $value);
    }
  }
  if ($pw->AtomEvidence($atom)) {
    if (@{ $pw->AtomEvidence($atom) }) {
      $self->PrLine("<EvidenceList>");
      Indent();
      for my $evidence (@{ $pw->AtomEvidence($atom) }) {
        $self->PrLine("<Evidence value=\"$evidence\">$evidence</Evidence>");
      }
      Exdent();
      $self->PrLine("</EvidenceList>");
    }
  }
  if (@{ $pw->AtomReferences($atom) }) {
    $self->PrLine("<ReferenceList>");
    Indent();
    for my $reference (@{ $pw->AtomReferences($atom) }) {
      $self->PrLine("<Reference pmid=\"$reference\">$reference</Reference>");
    }
    Exdent();
    $self->PrLine("</ReferenceList>");
  }
# At present, the notes are not for public consumption
#  if (@{ $pw->AtomNotes($atom) }) {
#    $self->PrLine("<NoteList>");
#    Indent();
#    for my $note (@{ $pw->AtomNotes($atom) }) {
#      $self->PrLine("<Note value=\"" .
#        CleanString($note) . "\">" . CleanString($note) . "</Note>");
#    }
#    Exdent();
#    $self->PrLine("</NoteList>");
#  }
  $self->PrLine("<InteractionComponentList>");
  Indent();
  for my $edge (sort numerically @{ $pw->Edges($atom) }) {
    my $molid = $pw->EdgeMol($atom, $edge);
    my @labels;
    my ($label, $edgetype) =
        $lv->LabelValueToString($pw->EdgeType($atom, $edge));
    for my $lvid (@{ $pw->EdgeLabel($atom, $edge) }) {
      my ($label, $value) = $lv->LabelValueToString($lvid);
      if ($label ne "edge-type") {
        push @labels, $lvid;
      }
    }
    for my $lvid (@{ $pw->MolLabel($atom, $edge) }) {
      push @labels, $lvid;
    }
    my $ptms = $pw->EdgePTM($atom, $edge);
    $self->InteractionComponent($molid, $edgetype, \@labels, $ptms);
  }
  Exdent();
  $self->PrLine("</InteractionComponentList>");
  Exdent();
  $self->PrLine("</Interaction>");

}

######################################################################
sub PathwayList {
  my ($self) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};

  if (@{ $pw->PathwayId }) {
    $self->PrLine("<PathwayList>");
    Indent();
    for my $pid (@{ $pw->PathwayId }) {
      $self->Pathway($pid);
    }
    Exdent();
    $self->PrLine("</PathwayList>");
  }

}

######################################################################
sub PathwayAtomMap {
  my ($self, $the_pid) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};

  for my $atom (@{ $pw->Atoms() }) {
    for my $pid (@{ $pw->AtomPathway($atom) }) {
      if (! $the_pid || $the_pid eq $pid) {
        $self->{pathwayatom}{$pid}{$atom} = 1;
      }
    }
  }
}

######################################################################
sub Pathway {
  my ($self, $id) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};

  my ($name, $org, $extid, $is_subnet) = 
      ($pw->PathwayName($id), $pw->PathwayOrg($id),
      $pw->PathwayExId($id), $pw->PathwayIsSubnet($id));
  if ($is_subnet) {
    $is_subnet = "true";
  } else {
    $is_subnet = "false";
  }
  $self->PrLine("<Pathway id=\"$id\" subnet=\"$is_subnet\">");
  Indent();
  $self->PrLine("<Organism>$org</Organism>");
  $self->PrLine("<LongName>" . CleanString($name) . "</LongName>");
  $self->PrLine("<ShortName>$extid</ShortName>");
  $self->Source($pw->PathwaySource($id));
  if (@{ $pw->Curators($id) }) {
    $self->PrLine("<CuratorList>");
    Indent();
    for my $curator (@{ $pw->Curators($id) }) {
      $self->PrLine("<Curator>" . CleanString($curator) . "</Curator>");
    }
    Exdent();
    $self->PrLine("</CuratorList>");
  }
  if (@{ $pw->Reviewers($id) }) {
    $self->PrLine("<ReviewerList>");
    Indent();
    for my $reviewer (@{ $pw->Reviewers($id) }) {
      $self->PrLine("<Reviewer>" . CleanString($reviewer) . "</Reviewer>");
    }
    Exdent();
    $self->PrLine("</ReviewerList>");
  }
  if (@{ $pw->PathwayReferences($id) }) {
    $self->PrLine("<ReferenceList>");
    Indent();
    for my $reference (@{ $pw->PathwayReferences($id) }) {
      $self->PrLine("<Reference pmid=\"$reference\">$reference</Reference>");
    }
    Exdent();
    $self->PrLine("</ReferenceList>");
  }
  if (defined $self->{pathwayatom}{$id}) {
    $self->PrLine("<PathwayComponentList>");
    Indent();
    for my $atom (keys %{ $self->{pathwayatom}{$id} }) {
      $self->PathwayComponent($atom, undef);
    }
    Exdent();
    $self->PrLine("</PathwayComponentList>");
  }
  Exdent();
  $self->PrLine("</Pathway>");
}

######################################################################
sub PrXML {
  my ($self, $pid) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};

  $self->PrLine("<?xml version=\"1.0\" encoding=\"UTF-8\"?>");
  $self->PrLine("<NCI_PID_XML>");
  $self->PrLine("<Created>" . TimeStamp() . "</Created>");
  Indent();
  $self->Ontology();
  $self->PrLine("<Model>");
  $self->Indent();
  $self->PathwayAtomMap($pid);
  $self->MoleculeList();
  $self->InteractionList($pid);
  $self->SubnetList();
  $self->PathwayList();  
  $self->Exdent();
  $self->PrLine("</Model>");
  Exdent();
  $self->PrLine("</NCI_PID_XML>");
}


######################################################################
1;
######################################################################
