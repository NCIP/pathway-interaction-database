#!/usr/local/bin/perl
package TableOutput;
require Exporter;
@ISA = qw(Exporter);
#@EXPORT = qw(
#);

use strict;
use Pathway;
use PWLabel;

######################################################################
sub new {
  my ($self, $pw, $lv, $organism, $sourceid,
      $pathway_id_map, $atom_id_map, $output_fh) = @_;
  my $x = {};
  $x->{pw} = $pw;
  $x->{lv} = $lv;
  $x->{organism} = $organism;
  $x->{sourceid} = $sourceid;
  $x->{pathwayidmap} = $pathway_id_map;
  $x->{atomidmap}    = $atom_id_map;
  $x->{output_fh} = $output_fh;
  return bless $x;
}

######################################################################
sub Lines {
  my ($self) = @_;
  return $self->{lines};
}

######################################################################
sub PrLine {
  my ($self, $x) = @_;
  my $fh = $self->{output_fh};
  if ($fh) {
    print $fh "$x\n";
  }
  push @{ $self->{lines} }, $x;
}

######################################################################
sub DoAll {
  my ($self) = @_;

  $self->DoPathways();
  $self->DoMols();
  $self->DoAtoms();
  $self->DoMaps();
}

######################################################################
sub DoMaps() {
  my ($self) = @_;

  my $sourceid = $self->{sourceid};
  my ($map);
  $map = $self->{pathwayidmap};

  for my $x (keys %{ $map }) {
    $self->PrLine("pathway_map\t$x\t$$map{$x}\t$sourceid");
  }
  $map = $self->{atomidmap};
  for my $x (keys %{ $map }) {
    $self->PrLine("atom_map\t$x\t$$map{$x}\t$sourceid");
  }
}

######################################################################
sub DoPathway {
  my ($self, $pid) = @_;

  my $pw = $self->{pw};
  my $source_id = $self->{sourceid};

  $self->PrLine(join("\t",
    "pathway",
    $pid,
    $pw->PathwayName($pid),
    $pw->PathwayOrg($pid),
##    $pw->PathwaySrcId($pid),
    $source_id,
    $pw->PathwayExId($pid)   
  ));

}

######################################################################
sub DoPathways {
  my ($self) = @_;

  my $pw = $self->{pw};
  my $highest_id = 0;

  for my $pid (@{ $pw->PathwayId }) {
    $self->DoPathway($pid);
    if ($pid > $highest_id) {
      $highest_id = $pid;
    }
  }
  $self->PrLine("highest_pathway_id\t$highest_id");
}

######################################################################
sub DoMol {
  my ($self, $mol_id) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};
  my $org = $self->{organism};

  my ($basic_mol_type, $id_type, $name_type, $i);

  $basic_mol_type = $lv->BasicMolType($pw->MolType($mol_id));
  $self->PrLine("mol\t$mol_id\t$basic_mol_type\t$org");
  if ($basic_mol_type eq "CX") {
    $self->DoComplex($mol_id);
  }

  my $mehash = $pw->MolExId($mol_id);
  for $id_type (keys %{ $mehash }) {
    for $i (keys %{ $$mehash{$id_type} }) {
      $self->PrLine("ext_mol_id\t$mol_id\t$i\t$id_type");
    }
  }

  my $mnhash = $pw->MolName($mol_id);
  for $name_type (keys %{ $mnhash }) {
    for $i (keys %{ $$mnhash{$name_type} }) {
      $self->PrLine("mol_name\t$mol_id\t$i\t$name_type");
    }
  }

}

######################################################################
sub DoMols {
  my ($self) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};

  for my $mol_id (@{ $pw->Mols() }) {
    $self->DoMol($mol_id);
  }
}

######################################################################
sub DoComplex {
  my ($self, $cx_id) = @_;

  my $pw = $self->{pw};

  my ($comp_id, $sub_mol_id, $lvid);

  for $comp_id (@{ $pw->Components($cx_id) }) {
    $sub_mol_id = $pw->ComponentMol($cx_id, $comp_id);
    $self->PrLine("component\t$cx_id\t$comp_id\t$sub_mol_id");
    for $lvid (@{ $pw->ComponentLabel($cx_id, $comp_id) }) {
      $self->PrLine("component_label\t$cx_id\t$comp_id\t$lvid");
    }
  }
}

######################################################################
sub DoAtoms {
  my ($self) = @_;

  my $pw = $self->{pw};
  my $highest_id = 0;

  for my $atom_id (@{ $pw->Atoms }) {
    $self->DoAtom($atom_id);
    if ($atom_id > $highest_id) {
      $highest_id = $atom_id;
    }
  }
  $self->PrLine("highest_atom_id\t$highest_id");
}

######################################################################
sub DoAtom {
  my ($self, $atom_id) = @_;

  my $pw = $self->{pw};
  my $org = $self->{organism};
  my $source_id = $self->{sourceid};

  my $atom_type = $pw->AtomType($atom_id);

  $self->PrLine("atom\t$atom_id\t$org\t$source_id");

  for my $lvid (@{ $pw->AtomLabel($atom_id) }) {
    $self->PrLine("atom_label\t$atom_id\t$lvid");
  }

  for my $lvid (@{ $pw->AtomCondition($atom_id) }) {
    $self->PrLine("atom_condition\t$atom_id\t$lvid");
  }

  for my $edge_id (@{ $pw->Edges($atom_id) }) {
    $self->DoEdge($atom_id, $edge_id);
  }

  for my $pathway_id (@{ $pw->AtomPathway($atom_id) }) {
    $self->PrLine("pathway_atom\t$pathway_id\t$atom_id");
  }

}

######################################################################
sub DoEdge {
  my ($self, $atom_id, $edge_id) = @_;

  my $pw = $self->{pw};

  my $mol_id    = $pw->EdgeMol($atom_id, $edge_id);
  my $edge_type = $pw->EdgeType($atom_id, $edge_id);

  $self->PrLine("edge\t$atom_id\t$edge_id\t$mol_id");

  for my $lvid (@{ $pw->EdgeLabel($atom_id, $edge_id) }) {
    $self->PrLine("edge_label\t$atom_id\t$edge_id\t$lvid");
  }

  for my $lvid (@{ $pw->MolLabel($atom_id, $edge_id) }) {
    $self->PrLine("mol_label\t$atom_id\t$edge_id\t$lvid");
  }
}

######################################################################
1;
######################################################################
