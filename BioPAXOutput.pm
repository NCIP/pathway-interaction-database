#!/usr/local/bin/perl


# Copyright SRA International
#
# Distributed under the OSI-approved BSD 3-Clause License.
# See http://ncip.github.com/pathway-interaction-database/LICENSE.txt for details.


package BioPAXOutput;
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

my $GO_TRANSCRIPTION_ID = "0006350";

my %id_type_lookup = (
  "UP" => "UniProt",
  "LL" => "EntrezGene",
  "EC" => "EnzymeConsortium",
  "CA" => "ChemicalAbstracts",
  "GO" => "GO",
  "KG" => "KEGG",
  "Pubmed" => "Pubmed",
  "MOD"    => "PSI-MOD",
  "MI"     => "PSI-MI",
  "CHEBI"  => "Chemical Entities of Biological Interest",
  "ECO"    => "Evidence Codes Ontology",
## added these pid_b_ things; in a few cases, there were duplicate
## definitions of ids like "GO_12345", one as a unificationXref and
## one as a relationshipXref. So we will prefix all relationshipXrefs
## with pid_b_
  "pid_b_GO" => "pid_b_GO",
  "pid_b_LL" => "pid_b_EntrezGene"
);

#####
##### PTM types
#####

# Acetylation: MOD:00394; acetyl group (CHEBI:22190)
# Biotinylation: MOD:00126; biotin (CHEBI:15956)??
# Dephosphoryaltion
# Farnesyaltion: MOD:00437; farnesyl group (CHEBI:24017)
# Glysoaminoglycan: MOD:00224????; glycosaminoglycans (CHEBI:18085)
# Geranylgeranylation: MOD:00441; geranylgeranyl group (CHEBI:24231)
# Glycosylation: MOD:00693; glycosyl glycosides (CHEBI:24407)???
# Hydroxylation: CHEBI:43176
# Methylation: MOD:00599; methyl group (CHEBI:32875)
# Myristoylation: MOD:00438; myristoyl group (CHEBI:25456)
# Oxidation: MOD:00412; 
# Palmitoylation: MOD:00440; palmitoyl group (CHEBI:25839)
# Phosphorylation: MOD:00696; phosphate group (CHEBI:32958)??
# poly (ADP-ribosyl)ation: CHEBI:16960
# Sumoylation
# Ubiquitination

  my %ptm_type_xrefs = (
    "acetylation"           => { "MOD" => "00394", "CHEBI" => "22190" },
    "biotinylation"         => { "MOD" => "00126" },
#    "compound modification" =>
    "dephosphorylation"     => { "GO" => "0016311" },
    "farnesylation"         => { "MOD" => "00437", "CHEBI" => "24017" },
    "geranylgeranylation"   => { "MOD" => "00441", "CHEBI" => "24231" },
    "glycosaminogen"        => { "CHEBI" => "18085" },
    "glycosaminoglycan"     => { "CHEBI" => "18085" },
    "glycosylation"         => { "MOD" => "00693" }, 
    "hydroxylation"         => { "CHEBI" => "43176" },
    "methylation"           => { "MOD" => "00599", "CHEBI" => "32875" },
    "myristoylation"        => { "MOD" => "00438", "CHEBI" => "25456" },
    "oxidation"             => { "MOD" => "00412" },
    "n-palmitoylation"      => { "MOD" => "00440", "CHEBI" => "25839" },
    "s-palmitoylation"      => { "MOD" => "00440", "CHEBI" => "25839" },
    "palmitoylation"        => { "MOD" => "00440", "CHEBI" => "25839" },
    "phosphorylation"       => { "MOD" => "00696" },
    "poly (ADP-ribosyl)ation" => { "CHEBI" => "16960" },
    "sumoylation"           => { "MI" => "0554", "GO" => "0016925" },
    "ubiquitination"        => { "MI" => "0189", "GO" => "0016567" },
    "_1_4_alpha_D_glucosyl_n_group" => { "CHEBI" => "15444" },
    "L_selenocysteinyl_group" => { "CHEBI" => "30003" },
    "palmitoyl_group"         => { "CHEBI" => "25839" },
    "phosphate_group"         => { "CHEBI" => "35780" },
    "carboxyl_group"          => { "CHEBI" => "23025" },
    "glycogen_group"          => { "CHEBI" => "28087" },
    "myristoyl_group"         => { "CHEBI" => "25456" },
    "acetyl_group"            => { "CHEBI" => "22190" },
    "hydroxyl_group"          => { "CHEBI" => "24706" },
    "half_cystyl_group"       => { "CHEBI" => "30770" },
    "methyl_group"            => { "CHEBI" => "25251" },
    "L_cystinyl_group"        => { "CHEBI" => "50066" },
    "decanoyl_group"	      => { "CHEBI" => "23574" },
    "octanoyl_group"          => { "CHEBI" => "25650" },
    "N_Glycan"                => { "CHEBI" => "52472" },
    "O_Phospho_L_serine"      => { "CHEBI" => "7692" },
    "pantetheine_4__phosphate_group" =>{ "CHEBI" => "47982" }
  );

#####
##### Evidence codes
#####

my %evidence_code_xrefs = (
  "IAE" => { "ECO" => "0000055" },
  "IC"  => { "ECO" => "0000001" },
  "IDA" => { "ECO" => "0000002" },
  "IEA" => { "ECO" => "0000036" },
  "IEP" => { "ECO" => "0000008" },
  "IFC" => { "ECO" => "0000012" },
  "IGI" => { "ECO" => "0000011" },
  "IMP" => { "ECO" => "0000015" },
#  "IOS" => "Inferred from Other Species",
  "IPI" => { "ECO" => "0000021" },
  "ISS" => { "ECO" => "0000041" },
  "NAS" => { "ECO" => "0000034" },
  "ND"  => { "ECO" => "0000035" },
#  "NIL" => "No Evidence",
  "NR"  => { "ECO" => "0000037" },
  "RCA" => { "ECO" => "0000053" },
  "RGE" => { "ECO" => "0000049" },
  "TAS" => { "ECO" => "0000033" }
);

my %evidence_codes = (
  "IAE" => "Inferred from Array Experiments",
  "IC"  => "Inferred by Curator",
  "IDA" => "Inferred from Direct Assay",
  "IEA" => "Inferred from Electronic Annotation",
  "IEP" => "Inferred from Expression Pattern",
  "IFC" => "Inferred from Functional Complementation",
  "IGI" => "Inferred from Genetic Interaction",
  "IMP" => "Inferred from Mutant Phenotype",
  "IOS" => "Inferred from Other Species",
  "IPI" => "Inferred from Physical Interaction",
  "ISS" => "Inferred from Sequence or Structural Similarity",
  "NAS" => "Non-traceable Author Statement",
  "ND"  => "No biological Data available",
  "NIL" => "No Evidence",
  "NR"  => "Not Recorded",
  "RCA" => "Inferred from Reviewed Computational Analysis",
  "RGE" => "Inferred from Reporter Gene Expression",
  "TAS" => "Traceable Author Statement"
);


#####
##### Other generic stuff
#####

my $TIME_STAMP = TimeStamp();

my $BIOPAX_HEADER = qq!<?xml version="1.0" encoding="UTF-8"?>
<rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"
  xmlns:bp="http://www.biopax.org/release/biopax-level2.owl#"
  xmlns:rdfs="http://www.w3.org/2000/01/rdf-schema#"
  xmlns:owl="http://www.w3.org/2002/07/owl#"
  xmlns:xsd="http://www.w3.org/2001/XMLSchema#"
  xmlns="http://pid.nci.nih.gov/biopax#"
  xml:base="http://pid.nci.nih.gov/biopax">
  <owl:Ontology rdf:about="">
    <owl:imports rdf:resource="http://www.biopax.org/release/biopax-level2.owl" />
    <rdfs:comment rdf:datatype="http://www.w3.org/2001/XMLSchema#string">BioPAX output created $TIME_STAMP, converted from the Pathway Interaction Database, National Cancer Institute, http://pid.nci.nih.gov.</rdfs:comment>
  </owl:Ontology>!;

my $BIOPAX_TRAILER = qq!</rdf:RDF>!;

my $HOMO_SAPIENS_DEF = qq!<bp:bioSource rdf:ID="Homo_sapiens">
  <bp:NAME rdf:datatype="http://www.w3.org/2001/XMLSchema#string">Homo sapiens</bp:NAME>
  <bp:TAXON-XREF rdf:resource="#NCBI_taxonomy_9606" />
</bp:bioSource>
<bp:unificationXref rdf:ID="NCBI_taxonomy_9606">
  <bp:DB rdf:datatype="http://www.w3.org/2001/XMLSchema#string">NCBI_taxonomy</bp:DB>
  <bp:ID rdf:datatype="http://www.w3.org/2001/XMLSchema#string">9606</bp:ID>
</bp:unificationXref>!;

my $PID_DATA_SOURCES = qq!<bp:dataSource rdf:ID="PID_DataSource">
  <bp:NAME rdf:datatype="http://www.w3.org/2001/XMLSchema#string">Pathway Interaction Database</bp:NAME>
  <bp:COMMENT rdf:datatype="http://www.w3.org/2001/XMLSchema#string">http://pid.nci.nih.gov</bp:COMMENT>
</bp:dataSource>
<bp:dataSource rdf:ID="PID_Curated_DataSource">
  <bp:NAME rdf:datatype="http://www.w3.org/2001/XMLSchema#string">Pathway Interaction Database NCI-Nature Curated Data</bp:NAME>
  <bp:COMMENT rdf:datatype="http://www.w3.org/2001/XMLSchema#string">http://pid.nci.nih.gov</bp:COMMENT>
</bp:dataSource>
<bp:dataSource rdf:ID="PID_BioCarta_DataSource">
  <bp:NAME rdf:datatype="http://www.w3.org/2001/XMLSchema#string">Pathway Interaction Database BioCarta Data</bp:NAME>
  <bp:COMMENT rdf:datatype="http://www.w3.org/2001/XMLSchema#string">http://pid.nci.nih.gov</bp:COMMENT>
</bp:dataSource>
<bp:dataSource rdf:ID="PID_Reactome_DataSource">
  <bp:NAME rdf:datatype="http://www.w3.org/2001/XMLSchema#string">Pathway Interaction Database Reactome Data</bp:NAME>
  <bp:COMMENT rdf:datatype="http://www.w3.org/2001/XMLSchema#string">http://pid.nci.nih.gov</bp:COMMENT>
</bp:dataSource>
<bp:dataSource rdf:ID="PID_KEGG_DataSource">
  <bp:NAME rdf:datatype="http://www.w3.org/2001/XMLSchema#string">Pathway Interaction Database KEGG Data</bp:NAME>
  <bp:COMMENT rdf:datatype="http://www.w3.org/2001/XMLSchema#string">http://pid.nci.nih.gov</bp:COMMENT>
</bp:dataSource>!;

######################################################################
sub ExtIdRef {
  my ($id_type, $id) = @_;
  if (! defined $id_type_lookup{$id_type}) {
    return "__" . "_$id";    
  } else {
    return $id_type_lookup{$id_type} . "_$id";
  }
}

######################################################################
sub RDF_ID {
  my ($tag, $id) = @_;
  return "<$tag rdf:ID=\"$id\" >";
}

######################################################################
sub RDF_resource {
  my ($tag, $id) = @_;
  return "<$tag rdf:resource=\"#$id\" />";
}

######################################################################
sub unificationXref {
  my ($id, $db, $entry) = @_;

  return qq!<bp:unificationXref rdf:ID="$id">
  <bp:DB rdf:datatype="http://www.w3.org/2001/XMLSchema#string">$db</bp:DB>
  <bp:ID rdf:datatype="http://www.w3.org/2001/XMLSchema#string">$entry</bp:ID>
</bp:unificationXref>!;
}

######################################################################
sub relationshipXref {
  my ($id, $db, $entry, $relationship_type) = @_;

  if ($relationship_type ne "") {
    return qq!<bp:relationshipXref rdf:ID="$id">
  <bp:DB rdf:datatype="http://www.w3.org/2001/XMLSchema#string">$db</bp:DB>
  <bp:ID rdf:datatype="http://www.w3.org/2001/XMLSchema#string">$entry</bp:ID>
  <bp:RELATIONSHIP-TYPE rdf:datatype="http://www.w3.org/2001/XMLSchema#string">$relationship_type</bp:RELATIONSHIP-TYPE>
</bp:relationshipXref>!;
  } else {
    return qq!<bp:relationshipXref rdf:ID="$id">
  <bp:DB rdf:datatype="http://www.w3.org/2001/XMLSchema#string">$db</bp:DB>
  <bp:ID rdf:datatype="http://www.w3.org/2001/XMLSchema#string">$entry</bp:ID>
</bp:relationshipXref>!;
  }
}

######################################################################
sub publicationXref {
  my ($id, $db, $entry) = @_;

  return qq!<bp:publicationXref rdf:ID="$id" >
  <bp:DB rdf:datatype="http://www.w3.org/2001/XMLSchema#string">$db</bp:DB>
  <bp:ID rdf:datatype="http://www.w3.org/2001/XMLSchema#string">$entry</bp:ID>
</bp:publicationXref>!;
}

######################################################################
sub CloseTag {
  my ($tag) = @_;

  return "</$tag>";
}

######################################################################
sub String {
  my ($tag, $s) = @_;
  return "<$tag " .
      "rdf:datatype=\"http://www.w3.org/2001/XMLSchema#string\"" .
      ">$s</$tag>";
}

######################################################################
sub Integer {
  my ($tag, $s) = @_;
  return "<$tag " .
      "rdf:datatype=\"http://www.w3.org/2001/XMLSchema#integer\"" .
      ">$s</$tag>";
}

######################################################################
sub ParseLabelIds {
  my ($self, $labels) = @_;
  
  my $pw = $self->{pw};
  my $lv = $self->{lv};

####### push lvids

  my (@activity, @location);
  if (defined $labels) {
    for my $lvid (sort numerically @{ $labels }) {
      my ($label, $value) = $lv->LabelValueToString($lvid);
      if ($label eq "activity-state") {
        push @activity, $lvid;
      } elsif ($label eq "location") {
        push @location, $lvid;
        $self->{locations}{$lvid} = 1;
      }
    }
  }
  return (\@activity, \@location);
}

######################################################################
sub ParseLabelValues {
  my ($self, $labels) = @_;
  
  my $pw = $self->{pw};
  my $lv = $self->{lv};

####### push values

  my (@activity, @location);
  if (defined $labels) {
    for my $lvid (sort numerically @{ $labels }) {
      my ($label, $value) = $lv->LabelValueToString($lvid);
      if ($label eq "activity-state") {
        push @activity, $value;
      } elsif ($label eq "location") {
        push @location, $value;
        $self->{locations}{$lvid} = 1;
      }
    }
  }
  return (\@activity, \@location);
}

######################################################################
sub MoleculeInstanceIdData {
  my ($self, $mol, $labels, $ptms) = @_;
  
  my $pw = $self->{pw};
  my $lv = $self->{lv};

  my $ptm_string = Pathway::NormalPTMString($ptms);
  $ptm_string =~ s/,/_/g;
  my ($activity, $location) = $self->ParseLabelIds($labels);
  my $activity_string = join("_", @{ $activity });
  my $location_string = join("_", @{ $location });
  my $part_bounds = $pw->PartBounds($mol);
  if ($part_bounds =~ /,/) {
    my ($start, $stop) = split(",", $part_bounds);
    if ($start eq "") {
      $start = "0";
    }
    if ($stop eq "") {
      $stop = "0";
    }
    $part_bounds = join("-", $start, $stop);
  }

  my $inst_id = join("_", "pid", "x", $mol, $part_bounds,
      $activity_string, $ptm_string, $location_string);
  $inst_id =~ s/_+/_/g;
  $inst_id =~ s/_$//;

  return ($part_bounds, $activity_string, $ptm_string, $location_string);

}

######################################################################
sub MoleculeInstanceNameData {
  my ($self, $mol, $labels, $ptms) = @_;
  
  my $pw = $self->{pw};
  my $lv = $self->{lv};

  my $mol_name = $pw->PickMolName($mol);
  my $ptm_string = Pathway::NormalPTMString($ptms);
  $ptm_string =~ s/,/_/g;

  my ($activity, $location) = $self->ParseLabelValues($labels);
  my $activity_string = join("_", @{ $activity });
  my $location_string = join("_", @{ $location });
  my $part_bounds = $pw->PartBounds($mol);
  if ($part_bounds =~ /,/) {
    my ($start, $stop) = split(",", $part_bounds);
    if ($start eq "") {
      $start = "0";
    }
    if ($stop eq "") {
      $stop = "0";
    }
    $part_bounds = join("-", $start, $stop);
  }

  return ($mol_name, $part_bounds, $activity_string, $ptm_string,
      $location_string);

}

######################################################################
sub PhysicalEntityId {
  my ($self, $mol, $labels, $ptms) = @_;
  
  my $pw = $self->{pw};
  my $lv = $self->{lv};


  my ($part_bounds, $activity_string, $ptm_string, $location_string) =
    $self->MoleculeInstanceIdData($mol, $labels, $ptms);

## 01/31/07: per suggestion from Guanming Wu (CSHL) we are moving
## activity-state and ptms off of the physical entity and onto the
## physical entity participant. So we no longer need physical entities
## that include ptms and activity-state
##  my $pe_id = join("_", "pid", "m", $mol, $part_bounds,
##      $activity_string, $ptm_string);

#  my $pe_id = join("_", "pid", "m", $mol, $part_bounds);
  my $pe_id = join("_", "pid", "m", $mol);

  $pe_id =~ s/_+/_/g;
  $pe_id =~ s/_$//;
  return $pe_id;
}

######################################################################
sub PhysicalEntityParticipantId {
  my ($self, $mol, $labels, $ptms) = @_;
  
  my $pw = $self->{pw};
  my $lv = $self->{lv};

  my ($part_bounds, $activity_string, $ptm_string, $location_string) =
    $self->MoleculeInstanceIdData($mol, $labels, $ptms);

  my $pep_id = join("_", "pid", "x", $mol, $part_bounds,
      $activity_string, $ptm_string, $location_string);
  $pep_id =~ s/_+/_/g;
  $pep_id =~ s/_$//;

  $self->{moleculeinstanceids}{$pep_id} = [ $mol, $labels,
      $ptms ];

  for my $lvid (@{ $labels }) {
    my ($label, $value) = $lv->LabelValueToString($lvid);
    if ($label eq "activity-state") {
      my ($protein_id, $position, $amino_acid,
          $modification_label_value, $modification_label_name) =
          ($mol, "0", "X", $lvid, $value);
      $self->{sequencefeatures}{ join("\t",
          $protein_id, $position, $amino_acid,
          $modification_label_value, $modification_label_name) } = 1;
    }
  }
  for my $ptm (@{ $ptms }) {
    my ($protein_id, $position, $amino_acid,
        $modification_label_value, $modification_label_name) = @{ $ptm };
    $self->{sequencefeatures}{
        join( "\t", $protein_id, $position, $amino_acid,
        $modification_label_value, $modification_label_name) } = 1;
  }


## 01/31/07: per suggestion from Guanming Wu (CSHL) we are moving
## activity-state and ptms off of the physical entity and onto the
## physical entity participant. So we no longer need physical entities
## that include ptms and activity-state
##  my $pe_id = join("_", "pid", "m", $mol, $part_bounds,
  ## we need a physical entity declaration for a modified molecule
#  my $pe_id = $self->PhysicalEntityId($mol, $labels, $ptms);
#  $pe_id =~ s/_+/_/g;
#  $pe_id =~ s/_$//;
#  my ($activity, $location) = $self->ParseLabelIds($labels);
#  if (@{ $activity } || @{ $ptms }) {
#      $self->{modifiedmols}{$pe_id} = [ $mol, $activity, $ptms ];
#  }

  return $pep_id;
}

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
sub MoleculeList {
  my ($self) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};

  my %seen;

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
}

######################################################################
sub Molecule {
  my ($self, $id, $labels, $ptms) = @_;

## this has to cover plain mols and mols with activity-state
## labels and/or ptms

  my $pw = $self->{pw};
  my $lv = $self->{lv};

  my $mol_name = $pw->PickMolName($id);
  my $pe_id = $self->PhysicalEntityId($id, $labels, $ptms);

  if (defined $self->{molecules_to_define}{$pe_id}) {
    delete $self->{molecules_to_define}{$pe_id};
  }
  if (defined $self->{molecules_defined}{$pe_id}) {
    return;
  } else {
    $self->{molecules_defined}{$pe_id} = [ $id, $labels, $ptms ];
  }

## types: protein, dna, rna, complex, smallMolecule
  my $tag;
  my $mol_id = $id;
  my ($label, $mol_type) = $lv->LabelValueToString($pw->MolType($mol_id));
  if ($mol_type eq "protein") {
    $tag = "bp:protein";
  } elsif ($mol_type eq "rna") {
    $tag = "bp:rna";
  } elsif ($mol_type eq "complex") {
    $tag = "bp:complex";
  } elsif ($mol_type eq "compound") {
    $tag = "bp:smallMolecule";
  } else {
    $tag = "physicalEntity";  ## this could only happen in error (or producing bpx for Reactome)
##    print STDERR "## untyped molecule id = $mol_id\n";
##    return;
  }

  $self->PrLine(RDF_ID($tag, $pe_id));
  Indent();

  if ($mol_type eq "protein" || $mol_type eq "rna" ||
      $mol_type eq "complex") {
    $self->PrLine(RDF_resource("bp:ORGANISM", "Homo_sapiens"));
  }
  $self->PrLine(RDF_resource("bp:DATA-SOURCE", "PID_DataSource"));
  $self->PrLine(String("bp:NAME", CleanString($self->PhysicalEntityName($id,
      $labels, $ptms))));
  my $mehash = $pw->MolExId($mol_id);
  for my $id_type (keys %{ $mehash }) {
    if (defined $id_type_lookup{$id_type}) {
      for my $i (keys %{ $$mehash{$id_type} }) {
        if ($mol_type eq "protein" && $id_type eq "LL") {
          $self->PrLine(RDF_resource("bp:XREF", ExtIdRef("pid_b_$id_type", $i)));
          $self->{relationship_xrefs}{ExtIdRef("pid_b_$id_type", $i)} =
              join("\t", $id_type, $i, "gene");
        } elsif (($mol_type eq "protein" || $mol_type eq "complex") &&
            $id_type eq "GO") {
          $self->PrLine(RDF_resource("bp:XREF", ExtIdRef("pid_b_$id_type", $i)));
          $self->{relationship_xrefs}{ExtIdRef("pid_b_$id_type", $i)} =
              join("\t", $id_type, $i, "Gene Ontology term for a protein");
        } else {
          $self->PrLine(RDF_resource("bp:XREF", ExtIdRef($id_type, $i)));
          $self->{unification_xrefs}{ExtIdRef($id_type, $i)} =
              join("\t", $id_type, $i);
        }
      }
    }
  }

  if ($pw->PartBounds($mol_id) ne "") {
    my $mehash = $pw->MolExId($mol_id);
    if (! defined $$mehash{UP}) {
      for my $whole (@{ $pw->MolWhole($mol_id) }) {
        my $mehash1 = $pw->MolExId($whole);
        if (defined $$mehash1{UP}) {
          for my $uniprot (keys %{ $$mehash1{UP} }) {
            $self->PrLine(RDF_resource("bp:XREF", ExtIdRef("UP", $uniprot)));
            $self->{unification_xrefs}{ExtIdRef("UP", $uniprot)} =
                join("\t", "UP", $uniprot);
          }
          last;
        }
      }
    }
  }

  my $mnhash = $pw->MolName($mol_id);
  for my $name_type (keys %{ $mnhash }) {
    for my $i (keys %{ $$mnhash{$name_type} }) {
      if ($i ne $mol_name) {
        $self->PrLine(String("bp:SYNONYMS", CleanString($i)));
      }
    }
  }

  if ($mol_type eq "complex") {
    for my $comp (sort numerically @{ $pw->Components($mol_id) }) {
      my $labels = $pw->ComponentLabel($mol_id, $comp);
      my $ptms = $pw->ComponentPTM($mol_id, $comp);
      my $compmol = $pw->ComponentMol($mol_id, $comp);
      my $pe_id = $self->PhysicalEntityId($compmol, $labels, $ptms);
      if (! defined $self->{molecules_defined}{$pe_id} ) {
        $self->{molecules_to_define}{$pe_id} = [ $compmol, $labels, $ptms ];
      }
      $self->ComplexComponent($compmol, $labels, $ptms);
    }
  }

  my $part_bounds = $pw->PartBounds($mol_id);
  if ($part_bounds ne "") {
    my ($start, $stop) = split(",", $part_bounds);
    if ($start eq "") {
      $start = "0";
    }
    if ($stop eq "") {
      $stop = "0";
    }
    $self->PrLine("<bp:COMMENT rdf:datatype=\"http://www.w3.org/2001/XMLSchema#string\">");
    Indent();
    my $uniprots;
    my $mehash = $pw->MolExId($mol_id);
    if (defined $$mehash{UP}) {
      $uniprots = join(", ", keys %{ $$mehash{UP} });
    } else {
      for my $whole (@{ $pw->MolWhole($mol_id) }) {
        my $mehash1 = $pw->MolExId($whole);
        if (defined $$mehash1{UP}) {
          $uniprots = join(", ", keys %{ $$mehash1{UP} });
          last;
        } else {
          $uniprots = "[unspecified]";
        }
      }
    }
    $self->PrLine("protein subunit; protein id = $uniprots; " .
        "start = $start, stop = $stop");
    Exdent();
    $self->PrLine("</bp:COMMENT>");
  }

  if ($mol_type eq "protein") {
    for my $member (@{ $pw->FamilyChildren($mol_id) }) {
      my $mehash = $pw->MolExId($member);
      if (defined $$mehash{UP}) {
        for my $uniprot (keys %{ $$mehash{UP} }) {
          my $xrefid = "pid_b_$mol_id" . "_$member";
          $self->PrLine(RDF_resource("bp:XREF", $xrefid));
          $self->{relationship_xrefs}{$xrefid} =
              join("\t", "UP", $uniprot, "protein family member");
        }
      } elsif (defined $$mehash{LL}) {
        for my $eg (keys %{ $$mehash{LL} }) {
          my $xrefid = "pid_b_$mol_id" . "_$member";
          $self->PrLine(RDF_resource("bp:XREF", $xrefid));
          $self->{relationship_xrefs}{$xrefid} =
              join("\t", "LL", $eg, "protein family member");
        }
      }
    }
  }

  Exdent();
  $self->PrLine(CloseTag($tag));

  for my $m (keys %{ $self->{molecules_to_define} }) {
    $self->Molecule(@{ $self->{molecules_to_define}{$m} });
  }
}

######################################################################
sub MoleculeId {
  my ($self, $molid) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};

  my $mol_name = $pw->PickMolName($molid);
  $mol_name =~ s/ /_/g;

  return "pid_m_$molid";
}

######################################################################
sub PhysicalEntityParticipantName {
  my ($self, $molid, $labels, $ptms) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};

  my $mol_name = $pw->PickMolName($molid);
  my $ptm_string = Pathway::NormalPTMString($ptms);
  $ptm_string =~ s/,/_/g;
  my (@activity, @location);
  my (@activity_name, @location_name);
  for my $lvid (sort numerically @{ $labels }) {
    my ($label, $value) = $lv->LabelValueToString($lvid);
    if ($label eq "activity-state") {
      push @activity, $lvid;
      push @activity_name, $value;
    } elsif ($label eq "location") {
      push @location, $lvid;
      push @location_name, $value;
      $self->{locations}{$lvid} = 1;
    }
  }
  my $activity_string = join("_", sort @activity_name);
  my $location_string = join("_", sort @location_name);
  my $part_bounds = $pw->PartBounds($molid);
  if ($part_bounds =~ /,/) {
    my ($start, $stop) = split(",", $part_bounds);
    if ($start eq "") {
      $start = "0";
    }
    if ($stop eq "") {
      $stop = "0";
    }
    $part_bounds = join("-", $start, $stop);
  }

  my $pep_name = join("_", "pid", "x", $molid, $mol_name,
      $part_bounds,
      $activity_string, $ptm_string, $location_string);
  $pep_name =~ s/_+/_/g;
  $pep_name =~ s/_$//;
  return $pep_name;
}

######################################################################
sub PhysicalEntityName {
  my ($self, $molid, $labels, $ptms) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};

  my $mol_name = $pw->PickMolName($molid);
  my $part_bounds = $pw->PartBounds($molid);
  if ($part_bounds =~ /,/) {
    my ($start, $stop) = split(",", $part_bounds);
    if ($start eq "") {
      $start = "0";
    }
    if ($stop eq "") {
      $stop = "0";
    }
    $part_bounds = join("-", $start, $stop);
  }

#  my $pe_name = join("_", "pid", "m", $molid, $mol_name,
#      $part_bounds);
  my $pe_name = join("_", $mol_name, $part_bounds);
  $pe_name =~ s/_+/_/g;
  $pe_name =~ s/_$//;
  return $pe_name;
}

######################################################################
sub ComplexComponent {
  my ($self, $id, $labels, $ptms) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};

  my $participantid = $self->PhysicalEntityParticipantId($id, $labels, $ptms);
  $self->PrLine(RDF_resource("bp:COMPONENTS", $participantid));

}

######################################################################
sub InteractionComponent {
  my ($self, $id, $edgetype, $labels, $ptms) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};

  if ($edgetype eq "input") {

    my $participantid =
        $self->PhysicalEntityParticipantId($id, $labels, $ptms);
#    if (! defined $self->{moleculeinstanceids}{$participantid}) {
#      $self->{moleculeinstanceids}{$participantid} = [ $id, $labels, $ptms ];
#    }
    $self->PrLine(RDF_resource("bp:LEFT", $participantid));

  } elsif ($edgetype eq "agent" || $edgetype  eq "inhibitor") {

    my $participantid =
        $self->PhysicalEntityParticipantId($id, $labels, $ptms);
#    if (! defined $self->{moleculeinstanceids}{$participantid}) {
#      $self->{moleculeinstanceids}{$participantid} = [ $id, $labels, $ptms ];
#    }
    $self->PrLine(RDF_resource("bp:CONTROLLER", $participantid));

  } elsif ($edgetype eq "output") {

    my $participantid =
        $self->PhysicalEntityParticipantId($id, $labels, $ptms);
#    if (! defined $self->{moleculeinstanceids}{$participantid}) {
#      $self->{moleculeinstanceids}{$participantid} = [ $id, $labels, $ptms ];
#    }

    $self->PrLine(RDF_resource("bp:RIGHT", $participantid));

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

  for my $subnet (@{ $pw->Subnets() }) {
    for my $atom (@{ $pw->AtomsInSubnet($subnet) }) {

      my ($label, $value) = $lv->LabelValueToString($pw->AtomType($atom));
      if ($value eq "pathway" || $value eq "subnet") {
        my $pid = $pw->AbstractionId($atom);
        my $pext_id = $pw->PathwayExId($pid);
## but then have to define this pathway somewhere...
        $self->PrLine(RDF_resource("bp:PATHWAY-COMPONENTS",
            "pid_p_$pid" . "_$pext_id"));
      } else {
        $self->PrLine(RDF_resource("bp:PATHWAY-COMPONENTS", "pid_i_$atom"));
        for my $ci (keys %{ $self->{control_interactions}{$atom} }) {
          $self->PrLine(RDF_resource("bp:PATHWAY-COMPONENTS", $ci));
        }
      }
    }
  }
}

######################################################################
sub AnalyzeInteractionType {
  my ($self, $id) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};

## rules:
##   1. this should return a BioPAX type, so don' call it for
##      PID types like pathway, subnet, macroprocess
##   2. if reaction, declare biochemicalReaction
##   3. if transcription, declare biochemicalReaction
##   4. from here on, look only at inputs and outputs
##   5. if there is at least one complex in inputs or outputs and
##      number of protein inputs != number of protein outputs
##      or number of complex inputs != number of complex outputs
##      declare complexAssembly
##   6. if same mol (input & output) has activity (activity labels or ptms)
##      change and location change, then declare
##      transportWithBiochemicalReaction
##   7. if same mol (input & output) has activity (activity labels or ptms)
##      change and NO location change, then declare BiochemicalReaction
##   8. if same mol (input & output) has location change and NO activity
##      (activity labels or ptms) change, then declare transport
##   9. if input compound is not matched by output compound, then
##      declare biochemicalReaction
##  10. if input protein is not matched by output  protein (e.g.,
##      cleaving), declare biochemicalReaction

  if ($pw->IsMacroProcess($id)) {
    return "interaction";
  }

  my $atomtype_lvid = $pw->AtomType($id);

  my $pathway_lvid = $lv->StringToLabelValue("process-type", "pathway");
  my $subnet_lvid = $lv->StringToLabelValue("process-type", "subnet");

  my ($dummy, $atomtype) = $lv->LabelValueToString($atomtype_lvid);

  if ($lv->IsA($atomtype_lvid, $pathway_lvid) ||
      $lv->IsA($atomtype_lvid, $subnet_lvid) ) {
    print STDERR "can't call AnalyzeInteractionType with type $atomtype\n";
    return "";
  }

  if ($atomtype eq "transcription") {
    return "biochemicalReaction";
  } elsif ($atomtype eq "reaction") {
    return "biochemicalReaction";
  }

  my (%mol2edge, %mol2type, %activity_labels, %ptms, %location_labels);
  my (@inputs, @outputs, %components, %edge2type, %no_match);
  my ($pr_in, $cx_in, $pr_out, $cx_out);

  for my $edge (@{ $pw->Edges($id) }) {

    my ($dummy, $edgetype) =
        $lv->LabelValueToString($pw->EdgeType($id, $edge));
    $edge2type{$edge} = $edgetype;
    if ($edgetype eq "input") {
      push @inputs, $edge;
    } elsif ($edgetype eq "output") {
      push @outputs, $edge;
    }

    my $mol = $pw->EdgeMol($id, $edge);
    $mol2edge{$mol}{$edge} = 1;
    my ($dummy, $moltype) =
        $lv->LabelValueToString($pw->MolType($mol));
      $mol2type{$mol} = $moltype;
    ## are we doing anything with these components ???
    ## Not now, but they could be used to test for state changes
    ## within complexes
    if ($moltype eq "complex") {
      for my $comp (@{ $pw->Components($mol) }) {
        push @{ $components{$edge} }, $pw->ComponentMol($mol, $comp);
      }
    }

    if ($moltype eq "protein") {
      if ($edgetype eq "input") {
        $pr_in++;
      } elsif ($edgetype eq "output") {
        $pr_out++;
      }
    } elsif ($moltype eq "complex") {
      if ($edgetype eq "input") {
        $cx_in++;
      } elsif ($edgetype eq "output") {
        $cx_out++;
      }
    }

    for my $lvid (@{ $pw->MolLabel($id, $edge) }) {
      my ($label, $value) = $lv->LabelValueToString($lvid);
      if ($label eq "activity-state") {
        $activity_labels{$edge}{$lvid} = 1;
      } elsif ($label eq "location") {
        $location_labels{$edge}{$lvid} = 1;
      }
    }    

    my $ptms = $pw->EdgePTM($id, $edge);
    if (@{ $ptms }) {
      $ptms{$edge} = Pathway::NormalPTMString($ptms);
    }
  }

  if (($cx_in || $cx_out) && ($pr_in != $pr_out || $cx_in != $cx_out) ) {
    return "complexAssembly";
  }

  ## if an input mol matches an output mol but there are 
  ## changes in ptms or state labels, then declare biochemicalConversion

  my ($location_change, $state_change);
  for my $mol (keys %mol2edge) {
    my @temp = keys %{ $mol2edge{$mol} };
    
    if (@temp) {
      my $mol2typer=$mol2type{$mol};
      # print STDERR "Moltype:  $mol	$mol2typer\n";
      $no_match{$mol2type{$mol}}++;
    }
    for (my $i = 0; $i < @temp - 1; $i++) {
      my $mol2typer=$mol2type{$mol};
      # print STDERR "MultiMoltype:  $mol      $mol2typer\n";
 
      my $edge_i = $temp[$i];
      my $type_i = $edge2type{$edge_i};
      if ($type_i ne "input" && $type_i ne "output") {
        next;
      }
      my $activity_labels =
          join(",", sort keys %{ $activity_labels{$edge_i} });
      my $location_labels =
          join(",", sort keys %{ $location_labels{$edge_i} });
      my $ptm_string = $ptms{$edge_i};
      for (my $j = $i+1; $j < @temp; $j++) {
        my $edge_j = $temp[$j];
        my $type_j = $edge2type{$edge_j};
        if ($type_j ne "input" && $type_j ne "output") {
          next;
        }
        if ($type_i eq $type_j) {
          next;
        }
        if (join(",", sort keys %{ $activity_labels{$edge_j} }) ne
            $activity_labels) {
          $state_change++;
        }
        if ($ptms{$edge_j} ne $ptm_string) {
          $state_change++;
        }
        if ($atomtype eq "translocation") {
          ## assume there is in fact a location change
          $location_change++;
        } elsif ($location_labels) {
          ## For present purposes, don't treat a mismatch between a
          ## null location spec
          ## and a non-null location spec as implying a transport.
          my $location_labels_j =
              join(",", sort keys %{ $location_labels{$edge_j} });
          if ($location_labels_j && $location_labels_j ne $location_labels) {
            $location_change++;
          }
        }
      }
    }
  }
  if ($state_change) {
    if ($location_change) {
      return "transportWithBiochemicalReaction";
    } else {
      return "biochemicalReaction";
    }
  } elsif ($location_change) {
    return "transport";
  } else {
    my $nomatchcompound = $no_match{"compound"};
        my $nomatchprotein = $no_match{"protein"};
        my $nomatchcomplex = $no_match{"complex"};
# print STDERR "Compound:  $nomatchcompound       Protein:  $nomatchprotein       Complex:  $nomatchcomplex\n";
 
    if ($no_match{"compound"}) {
      ## metabolic
      return "biochemicalReaction";
    } elsif ($no_match{"protein"}) {
      ## signaling (e.g. cleaving)
      return "biochemicalReaction";
    } elsif ($no_match{"complex"}) {
## here we are, simplistically, assuming that there is a state change
## inside a complex (resulting in a different complex, so the input complex
## is not matched by the outside complex)
      return "biochemicalReaction";
    } else {
## presumably ill-formed stuff
      if ($atomtype eq "translocation") {
        return "transport";
      } else {
        my $nomatchcompound = $no_match{"compound"};
        my $nomatchprotein = $no_match{"protein"};
        my $nomatchcomplex = $no_match{"complex"};
# print STDERR "Compound:  $nomatchcompound	Protein:  $nomatchprotein	Complex:  $nomatchcomplex\n";
print STDERR "can't determine interaction type for $id ($atomtype)\n";
        return "interaction";
      }
    }

  }
}

######################################################################
sub Interaction {
  my ($self, $id) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};

  my (@agent_edges, @inhibitor_edges);
  my ($tag);

  my $atomtype_lvid = $pw->AtomType($id);

  my $pathway_lvid = $lv->StringToLabelValue("process-type", "pathway");
  my $subnet_lvid = $lv->StringToLabelValue("process-type", "subnet");

  if ($lv->IsA($atomtype_lvid, $pathway_lvid) ||
      $lv->IsA($atomtype_lvid, $subnet_lvid) ||
      $pw->IsMacroProcess($id) ) {
    $tag = "interaction";
  } else {
    $tag = $self->AnalyzeInteractionType($id);
  }

  $self->PrLine(RDF_ID("bp:$tag", "pid_i_$id"));
  Indent();

  my $datasource = $pw->SourceName($pw->AtomSource($id));
  if ($datasource eq "NATURE") {
    $self->PrLine(RDF_resource("bp:DATA-SOURCE", "PID_Curated_DataSource"));
  } elsif ($datasource eq "BioCarta") {
    $self->PrLine(RDF_resource("bp:DATA-SOURCE", "PID_BioCarta_DataSource"));
  } elsif ($datasource eq "Reactome") {
    $self->PrLine(RDF_resource("bp:DATA-SOURCE", "PID_Reactome_DataSource"));
  } elsif ($datasource eq "KEGG") {
    $self->PrLine(RDF_resource("bp:DATA-SOURCE", "PID_KEGG_DataSource"));
  } else {
    $self->PrLine(RDF_resource("bp:DATA-SOURCE", "PID_DataSource"));
  }

  my $atom = $id;

  if ($pw->IsMacroProcess($id)) {
    my ($dummy, $name) = $lv->LabelValueToString($atomtype_lvid);
    $self->PrLine(String("bp:NAME", CleanString($name)));
    ## sort to make predictable
    for my $go_id (sort @{ $lv->GOTermsFor($atomtype_lvid) }) {
      ## just take the first one
      $go_id =~ s/^GO:?\s*//;
      $self->PrLine(RDF_resource("bp:XREF", ExtIdRef("pid_b_GO", $go_id)));
      $self->{relationship_xrefs}{ExtIdRef("pid_b_GO", $go_id)} =
          join("\t", "GO", $go_id, "Gene Ontology biological process");
      last;
    }

  } elsif ($lv->IsA($atomtype_lvid, $pathway_lvid) ||
      $lv->IsA($atomtype_lvid, $subnet_lvid)) {
      my $pext_id = $pw->AbstractionExtId($id);
      my $pid = $pw->AbstractionId($id);
      my $participantid = "pid_p_$pid" . "_$pext_id";
      $self->PrLine(RDF_resource("bp:PARTICIPANTS", $participantid));
  } else {
    my ($dummy, $atomtype) = $lv->LabelValueToString($pw->AtomType($id));
    if ($atomtype eq "transcription") {
      my $go_id = $GO_TRANSCRIPTION_ID;
      $self->PrLine(RDF_resource("bp:XREF", ExtIdRef("pid_b_GO", $go_id)));
      $self->{relationship_xrefs}{ExtIdRef("pid_b_GO", $go_id)} =
          join("\t", "GO", $go_id, "Gene Ontology biological process");
    }
  }

  if ($pw->AtomEvidence($atom)) {
    if (@{ $pw->AtomEvidence($atom) }) {
      for my $evidence (@{ $pw->AtomEvidence($atom) }) {
## needs a bunch more stuff: bp:openControlledVocabulary with a TERM property and a unification XREF
#        $self->PrLine("<bp:EVIDENCE-CODE>$evidence</bp:EVIDENCE-CODE>");
##begin new
        $self->{evidence_codes}{$evidence} = 1;
        $self->PrLine("<bp:EVIDENCE>");
        Indent();
        $self->PrLine("<bp:evidence rdf:ID=\"pid_ev_$atom" . "_$evidence\">");
        Indent();
        $self->PrLine(RDF_resource("bp:EVIDENCE-CODE", $self->EvidenceId($evidence)));
        Exdent();
        $self->PrLine("</bp:evidence>");
        Exdent();
        $self->PrLine("</bp:EVIDENCE>");
##end new
      }
    }
  }
  if (@{ $pw->AtomReferences($atom) }) {
    for my $reference (@{ $pw->AtomReferences($atom) }) {
      $self->PrLine(RDF_resource("bp:XREF", "Pubmed_$reference"));
      $self->{publication_xrefs}{"Pubmed_$reference"} =
          join("\t", "Pubmed", $reference);
    }
  }
## do not print notes: they contain unfiltered information
#  if (@{ $pw->AtomNotes($atom) }) {
#    for my $note (@{ $pw->AtomNotes($atom) }) {
#      $self->PrLine(String("bp:COMMENT", $note));
#    }
#  }

  for my $edge (sort numerically @{ $pw->Edges($atom) }) {
    my $molid = $pw->EdgeMol($atom, $edge);
    my @labels;
    my ($label, $edgetype) =
        $lv->LabelValueToString($pw->EdgeType($atom, $edge));
    if ($edgetype eq "agent") {
      push @agent_edges, $edge;
      next;
    } elsif ($edgetype eq "inhibitor") {
      push @inhibitor_edges, $edge;
      next;
    }
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
## "interaction" has PARTICIPANTS, not LEFT or RIGHT
    if ($tag eq "interaction") {
      my $participantid =
          $self->PhysicalEntityParticipantId($molid, \@labels, $ptms);
      $self->PrLine(RDF_resource("bp:PARTICIPANTS", $participantid));
    } else {
      $self->InteractionComponent($molid, $edgetype, \@labels, $ptms);
    }
  }

  Exdent();
  $self->PrLine(CloseTag("bp:$tag"));

  ##
  ## conditions
  ##

  my $condition_num = 0;
  for my $lvid (@{ $pw->AtomCondition($atom) }) {
    $condition_num++;
    $self->BuildControlInteraction($atom, undef, "condition", $condition_num);
  }

  my $negative_condition_num = 0;
  for my $lvid (@{ $pw->AtomNegativeCondition($atom) }) {
    $condition_num++;
    $self->BuildControlInteraction($atom, undef, "negative_condition", $negative_condition_num);
  }

  ##
  ## agents
  ##

  my $agent_num = 0;
  for my $edge (@agent_edges) {
    $agent_num++;
    $self->BuildControlInteraction($atom, $edge, "agent", $agent_num);
  }

  ##
  ## inhibitors
  ##

  my $inhibitor_num = 0;
  for my $edge (@inhibitor_edges) {
    $inhibitor_num++;
    $self->BuildControlInteraction($atom, $edge, "inhibitor", $inhibitor_num);
  }

}

######################################################################
sub BuildControlInteraction {
  my ($self, $original_atom_id, $original_edge, $control_type,
      $control_num) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};

  my ($CT, $ct);

  if ($control_type eq "agent") {
    $CT = "ACTIVATION";
    $ct = "activator";
  } elsif ($control_type eq "inhibitor") {
    $CT = "INHIBITION";
    $ct = "inhibitor";
  } elsif ($control_type eq "condition") {
    $CT = "ACTIVATION";
    $ct = "condition";
  } elsif ($control_type eq "negative_condition") {
    $CT = "INHIBITION";
    $ct = "condition";
  } else {
    print STDERR "bad control type $control_type\n";
    return;
  }

  my $control_interaction_id = "pid_i_$original_atom_id" .
      "_" . $ct . "_" . $control_num;
  $self->PrLine(RDF_ID("bp:control", $control_interaction_id));
  $self->{control_interactions}{$original_atom_id}
      {$control_interaction_id} = 1;
  Indent();
  $self->PrLine(String("bp:CONTROL-TYPE", $CT));

  if ($control_type eq "agent" || $control_type eq "inhibitor") {

    my $molid = $pw->EdgeMol($original_atom_id, $original_edge);
    my @labels;
    for my $lvid (@{ $pw->EdgeLabel($original_atom_id, $original_edge) }) {
      my ($label, $value) = $lv->LabelValueToString($lvid);
      if ($label ne "edge-type") {
        push @labels, $lvid;
      }
    }
    for my $lvid (@{ $pw->MolLabel($original_atom_id, $original_edge) }) {
      push @labels, $lvid;
    }
    my $ptms = $pw->EdgePTM($original_atom_id, $original_edge);
    $self->InteractionComponent($molid, "agent", \@labels, $ptms);

  }
  if ($control_type eq "condition") {      ## condition need to add negative condition with control_type

    my $lvid = ${ $pw->AtomCondition($original_atom_id) } [$control_num-1];
    my ($dummy, $name) = $lv->LabelValueToString($lvid);
    $self->PrLine(String("bp:NAME", CleanString($name)));
    ## sort to make predictable
    for my $go_id (sort @{ $lv->GOTermsFor($lvid) }) {
      ## just take the first one
      $go_id =~ s/^GO:?\s*//;
      $self->PrLine(RDF_resource("bp:XREF", ExtIdRef("pid_b_GO", $go_id)));
      $self->{relationship_xrefs}{ExtIdRef("pid_b_GO", $go_id)} =
          join("\t", "GO", $go_id, "Gene Ontology biological process");
      last;
    }

  }
  if ($control_type eq "negative_condition") {      ## condition need to add negative condition with control_type

    my $lvid = ${ $pw->AtomNegativeCondition($original_atom_id) } [$control_num-1];
    my ($dummy, $name) = $lv->LabelValueToString($lvid);
    $self->PrLine(String("bp:NAME", CleanString($name)));
    ## sort to make predictable
    for my $go_id (sort @{ $lv->GOTermsFor($lvid) }) {
      ## just take the first one
      $go_id =~ s/^GO:?\s*//;
      $self->PrLine(RDF_resource("bp:XREF", ExtIdRef("pid_b_GO", $go_id)));
      $self->{relationship_xrefs}{ExtIdRef("pid_b_GO", $go_id)} =
          join("\t", "GO", $go_id, "Gene Ontology biological process");
      last;
    }

  }


  $self->PrLine(RDF_resource("bp:CONTROLLED", "pid_i_$original_atom_id"));

  Exdent();

  $self->PrLine(CloseTag("bp:control"));
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

  if ($id < 0) {
    return;       ## this is not a predefined pathway
  }

  my ($name, $org, $extid, $is_subnet) = 
      ($pw->PathwayName($id), $pw->PathwayOrg($id),
      $pw->PathwayExId($id), $pw->PathwayIsSubnet($id));

  if ($is_subnet) {
    $is_subnet = "true";
  } else {
    $is_subnet = "false";
  }

  my $pid = "pid_p_$id" . "_$extid";
  $self->PrLine(RDF_ID("bp:pathway", $pid));
  Indent();

  $self->PrLine(RDF_resource("bp:ORGANISM", "Homo_sapiens"));

  my $datasource = $pw->SourceName($pw->PathwaySource($id));
  if ($datasource eq "NATURE") {
    $self->PrLine(RDF_resource("bp:DATA-SOURCE", "PID_Curated_DataSource"));
  } elsif ($datasource eq "BioCarta") {
    $self->PrLine(RDF_resource("bp:DATA-SOURCE", "PID_BioCarta_DataSource"));
  } elsif ($datasource eq "Reactome") {
    $self->PrLine(RDF_resource("bp:DATA-SOURCE", "PID_Reactome_DataSource"));
  } else {
    $self->PrLine(RDF_resource("bp:DATA-SOURCE", "PID_DataSource"));
  }

  $self->PrLine(String("bp:NAME", CleanString($name)));

  if (@{ $pw->Curators($id) }) {
## AUTHORS probably only for an article citation, not for curator
    my $note = "Pathway curators: " . join(", ", @{ $pw->Curators($id) });
    $self->PrLine(String("bp:COMMENT", CleanString($note)));
  }

  if (@{ $pw->Reviewers($id) }) {
## AUTHORS probably only for an article citation, not for reviewer
    my $note = "Pathway reviewers: " . join(", ", @{ $pw->Reviewers($id) });
    $self->PrLine(String("bp:COMMENT", CleanString($note)));
  }

  if (defined $self->{pathwayatom}{$id}) {
    for my $atom (keys %{ $self->{pathwayatom}{$id} }) {
      my ($label, $value) = $lv->LabelValueToString($pw->AtomType($atom));
      if ($value eq "pathway" || $value eq "subnet") {
        my $pid = $pw->AbstractionId($atom);
        my $pext_id = $pw->PathwayExId($pid);
## but then have to define this pathway somewhere...
        $self->PrLine(RDF_resource("bp:PATHWAY-COMPONENTS",
            "pid_p_$pid" . "_$pext_id"));
      } else {
        $self->PrLine(RDF_resource("bp:PATHWAY-COMPONENTS", "pid_i_$atom"));
        for my $ci (keys %{ $self->{control_interactions}{$atom} }) {
          $self->PrLine(RDF_resource("bp:PATHWAY-COMPONENTS", $ci));
        }
      }
    }
  }

  Exdent();
  $self->PrLine(CloseTag("bp:pathway"));

}

######################################################################
sub LocationId {
  my ($self, $location) = @_;

  my $loc = $location;
  $loc =~ s/ /_/g;
  return $loc;
}

######################################################################
sub PTMtypes {
  my ($self) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};

  for my $ptm_type (keys %{ $self->{ptm_types} }) {
    if ($ptm_type eq "") {
      next;
    }
    $self->PrLine(RDF_ID("bp:openControlledVocabulary",
        $ptm_type));
    Indent();
    $self->PrLine(String("bp:TERM", $ptm_type));

    if (defined $ptm_type_xrefs{$ptm_type}) {
      for my $db (keys %{ $ptm_type_xrefs{$ptm_type} }) {
        my $id = $ptm_type_xrefs{$ptm_type}{$db};
        my $compound_id = $db . "_$id";
        $self->PrLine(RDF_resource("bp:XREF", $compound_id));
        if ($db eq "MOD") {
          $self->{unification_xrefs}{$compound_id} =
              join("\t", $db, $id);
        } else {
## openControlledVocabulary cannot have a relationshipXref
##          $self->{relationship_xrefs}{$compound_id} =
##              join("\t", $db, $id);
          $self->{unification_xrefs}{$compound_id} =
              join("\t", $db, $id);
        }
      }
    }

    Exdent();
    $self->PrLine(CloseTag("bp:openControlledVocabulary"));
  }
}

######################################################################
sub SequenceSites {
  my ($self) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};

  for my $ss (keys %{ $self->{sequence_sites} }) {
    my ($protein, $position) = split(/\t/, $ss);
    $self->SequenceSite($protein, $position);
  }
}

######################################################################
sub SequenceFeatures {
  my ($self) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};

  my $complex_type = $lv->StringToLabelValue("molecule-type", "complex");
  for my $sf (keys %{ $self->{sequencefeatures} }) {
    my ($protein, $position, $amino_acid, $modification_label_value,
      $modification_label_name) = split(/\t/, $sf);
    $self->{ptm_types}{$modification_label_name} = 1;
    if ($position ne "0") {
      $self->{sequence_sites}{ join("\t", $protein, $position) } = 1;
    }
    if ($pw->MolType($protein) eq $complex_type) {
##    don't print
    } else {
      $self->SequenceFeature($protein, $position, $amino_acid,
        $modification_label_value, $modification_label_name);
    }
  }
  $self->PTMtypes();
  $self->SequenceSites();
}

######################################################################
sub EvidenceId {
  my ($self, $code) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};

  return "pid_ev_$code";
}

######################################################################
sub EvidenceCodes {
  my ($self) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};

  for my $ec (keys %{ $self->{evidence_codes} }) {
    $self->PrLine(RDF_ID("bp:openControlledVocabulary",
        $self->EvidenceId($ec)));
    Indent();
    my $name = $evidence_codes{$ec};
    $self->PrLine(String("bp:TERM", $name));
    if (defined $evidence_code_xrefs{$ec}) {
      for my $db (keys %{ $evidence_code_xrefs{$ec} }) {
        my $id = $evidence_code_xrefs{$ec}{$db};
        my $compound_id = $db . "_$id";
        $self->PrLine(RDF_resource("bp:XREF", $compound_id));
        if ($db eq "ECO") {
          $self->{unification_xrefs}{$compound_id} =
              join("\t", $db, $id);
        } else {
## openControlledVocabulary cannot have a relationshipXref
##          $self->{relationship_xrefs}{$compound_id} =
##              join("\t", $db, $id);
          $self->{unification_xrefs}{$compound_id} =
              join("\t", $db, $id);
        }
      }
    }
    Exdent();
    $self->PrLine(CloseTag("bp:openControlledVocabulary"));
  }
}

######################################################################
sub Locations {
  my ($self) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};

  for my $loc (keys %{ $self->{locations} }) {
    my ($dummy, $name) = $lv->LabelValueToString($loc);
    $self->PrLine(RDF_ID("bp:openControlledVocabulary",
        $self->LocationId($name)));
    Indent();
    $self->PrLine(String("bp:TERM", $name));

    my @go = @{ $lv->GOTermsFor($loc) };
    for my $go (@go) {
      $go =~ s/GO://;
      $self->PrLine(RDF_resource("bp:XREF", "GO_$go"));
      $self->{unification_xrefs}{"GO_$go"} =
          join("\t", "GO", $go, "Gene Ontology cellular component");
    }
    Exdent();
    $self->PrLine(CloseTag("bp:openControlledVocabulary"));
  }
}

######################################################################
sub RelationshipXrefs {
  my ($self) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};

  for my $id (keys %{ $self->{relationship_xrefs} }) {
    my $xref = $self->{relationship_xrefs}{$id};
    my ($db, $entry, $relationship_type) = split(/\t/, $xref);
    $self->PrLine(relationshipXref($id, $id_type_lookup{$db}, $entry,
        $relationship_type));
  }
}

######################################################################
sub UnificationXrefs {
  my ($self) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};

  for my $id (keys %{ $self->{unification_xrefs} }) {
    my $xref = $self->{unification_xrefs}{$id};
    my ($db, $entry) = split(/\t/, $xref);
    $self->PrLine(unificationXref($id, $id_type_lookup{$db}, $entry));
  }
}

######################################################################
sub PublicationXrefs {
  my ($self) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};

  for my $id (keys %{ $self->{publication_xrefs} }) {
    my $xref = $self->{publication_xrefs}{$id};
    my ($db, $entry) = split(/\t/, $xref);
    $self->PrLine(publicationXref($id, $id_type_lookup{$db}, $entry));
  }
}

######################################################################
sub SequenceFeatureId {
  my ($protein_id, $position, $amino_acid, $modification_label_value,
      $modification_label_name) = @_;

  if ($position eq "") {
    $position = "0";
  }
  return join("_", "pid", "f",
      $protein_id, "$position$amino_acid$modification_label_value");
}

######################################################################
sub SequenceSiteId {
  my ($protein, $position) = @_;

  if ($position eq "") {
    $position = "0";
  }
  return join("_", "pid", "s", $protein, $position);
}

######################################################################
sub SequenceIntervalId {
  my ($protein, $begin, $end) = @_;

  if ($begin eq "") {
    $begin = "null";
  }
  if ($end eq "") {
    $end = "null";
  }
  return join("_", "pid", "s", $protein, $begin, $end);
}

######################################################################
sub SequenceFeature {
  my ($self, $protein, $position, $amino_acid, $modification_label_value,
      $modification_label_name) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};

  $self->PrLine(RDF_ID("bp:sequenceFeature", SequenceFeatureId($protein,
      $position, $amino_acid, $modification_label_value,
      $modification_label_name)));
  Indent();
  $self->PrLine(RDF_resource("bp:FEATURE-TYPE", $modification_label_name));
  if ($position ne "0") {
    $self->PrLine(RDF_resource("bp:FEATURE-LOCATION", SequenceSiteId($protein,
      $position)));
  }
  Exdent();
  $self->PrLine(CloseTag("bp:sequenceFeature"));

}

######################################################################
sub SequenceSite {
  my ($self, $protein, $position) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};

  my $status = "EQUAL";
  $self->PrLine(RDF_ID("bp:sequenceSite", SequenceSiteId($protein, $position)));
  Indent();
  $self->PrLine(String("bp:POSITION-STATUS", $status));
  $self->PrLine(Integer("bp:SEQUENCE-POSITION", $position));
  Exdent();
  $self->PrLine(CloseTag("bp:sequenceSite"));

}

######################################################################
sub SequenceInterval {
  my ($self, $id, $begin, $end) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};

  $self->PrLine(RDF_ID("bp:sequenceInterval", $id));
  Indent();
  $self->PrLine(Integer("bp:SEQUENCE-INTERVAL-BEGIN", $begin));
  $self->PrLine(Integer("bp:SEQUENCE-INTERVAL-END", $end));
  Exdent();
  $self->PrLine(CloseTag("bp:sequenceInterval"));

}

######################################################################
sub MoleculeInstances {
  my ($self) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};

  for my $id (keys %{ $self->{modifiedmols} }) {
    my ($molid, $labels, $ptms) = @{ $self->{modifiedmols}{$id} };
    $self->Molecule($molid, $labels, $ptms);
  }

  for my $id (keys %{ $self->{moleculeinstanceids} }) {
    my ($molid, $labels, $ptms, ) = @{ $self->{moleculeinstanceids}{$id} };
    my $moltype = $pw->MolType($molid);
    my $tag;
    if ($moltype eq $lv->StringToLabelValue("molecule-type", "complex")) {
      $tag = "bp:physicalEntityParticipant";
    } elsif ($moltype eq $lv->StringToLabelValue("molecule-type", "protein")) {
      $tag = "bp:sequenceParticipant";
    } elsif ($moltype eq $lv->StringToLabelValue("molecule-type", "rna")) {
      $tag = "bp:sequenceParticipant";
    } elsif ($moltype eq $lv->StringToLabelValue("molecule-type", "compound")) {
      $tag = "bp:physicalEntityParticipant";
    } else {
      $tag = "bp:sequenceParticipant";
    }
    $self->PrLine(RDF_ID($tag, $id));
    Indent();
    $self->PrLine(RDF_resource("bp:PHYSICAL-ENTITY",
        $self->PhysicalEntityId($molid, $labels, $ptms)));
    for my $lvid (@{ $labels }) {
      my ($label, $value) = $lv->LabelValueToString($lvid);
      if ($label eq "location") {
        $self->PrLine(RDF_resource("bp:CELLULAR-LOCATION",
            $self->LocationId($value)))
      } elsif ($label eq "activity-state") {
        my ($protein_id, $position, $amino_acid, $modification_label_value,
            $modification_label_name) = ($molid, "0", "X", $lvid, $value);
        if ($tag eq "bp:physicalEntityParticipant") {
          $self->PrLine("<bp:COMMENT rdf:datatype=\"http://www.w3.org/2001/XMLSchema#string\">");
          Indent();
          $self->PrLine(
              "bp:SEQUENCE-FEATURE-LIST " .
              SequenceFeatureId($protein_id, $position, $amino_acid,
              $modification_label_value, $modification_label_name)  .
              " ; " .
              "bp:FEATURE-TYPE " .
              $modification_label_name
            );
          Exdent();
          $self->PrLine("</bp:COMMENT>");
        } else {
          $self->PrLine(RDF_resource("bp:SEQUENCE-FEATURE-LIST",
            SequenceFeatureId($protein_id, $position, $amino_acid,
            $modification_label_value, $modification_label_name)));
        }
      }
    }
    for my $ptm (@{ $ptms }) {
      my ($protein_id, $position, $amino_acid, $modification_label_value,
          $modification_label_name) = @{ $ptm };
      $self->PrLine(RDF_resource("bp:SEQUENCE-FEATURE-LIST",
          SequenceFeatureId($protein_id, $position, $amino_acid,
          $modification_label_value, $modification_label_name)));
    }

    Exdent();
    $self->PrLine(CloseTag($tag));
  }
}

######################################################################
sub PrOWL {
  my ($self) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};

  $self->PrLine($BIOPAX_HEADER);
  $self->PrLine($HOMO_SAPIENS_DEF);
  $self->PrLine($PID_DATA_SOURCES);

  $self->PathwayAtomMap();
  $self->MoleculeList();
  $self->InteractionList();
  $self->SubnetList();
  $self->PathwayList();  

  $self->EvidenceCodes();
  $self->Locations();
  $self->SequenceFeatures();
  $self->UnificationXrefs();
  $self->RelationshipXrefs();
  $self->PublicationXrefs();
  $self->MoleculeInstances();

  $self->PrLine($BIOPAX_TRAILER);

}


######################################################################
1;
######################################################################
