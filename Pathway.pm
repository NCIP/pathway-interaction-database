package Pathway;
require Exporter;
@ISA = qw(Exporter);
#@EXPORT = qw(
#  CreateIsA
#  CreateLabelSort
#);

use strict;
use PWLabel;
use Clan;

## pathway:   pid,(pname|porg|psrcid|pexid) -> String
## moltype:   molid             -> moltype        ## actually a label_value
## molexid:   molid, exidtype -> Set Of exid      ## external ids
## molname:   molid, nametype   -> Set Of name
## mollabel:  atomid, edgeid    -> Set Of label_value  ## i.e. mol instance
## moluse:    molid             -> Set Of atomid
## edge:      atomid, edgeid    -> molid
## edgetype:  atomid, edgeid    -> edgetype       ## actually a label_value
## edgelabel: atomid, edgeid    -> Set Of label_value
## atomtype:  atomid            -> atomtype       ## actually a label_value
## atomlabel: atomid            -> Set Of label_value
## component:      complexmolid, componentid -> molid
## r_component: molid           -> Set Of complex_molid
## componentlabel: complexmolid, componentid -> Set Of label_value
## atompathway:    atom_id                   -> Set Of pathway_id
## pathwaysubnet:  pathwayid                 -> Set Of subnetid
## subnetatom:     subnetid   -> Set Of atomid
## abstraction: atomid  -> pathway_id;
     ## for when an interaction stands for a connection to another pathway
## atomcondition: atom_id       -> Set Of label_value
     ## have to keep process-type [atomlabel] and process-condition
     ## separate
## macroprocess: atomid         -> label_value
## clans        -> List Of Clan
## family:    family_id -> Set Of member_id
## r_family:  member_id -> Set Of family_id
## partof:    whole_mol_id -> Set Of part_mol_id
## r_partof:  part_mol_id  -> Set Of whole_mol_id  (should be one-to-one)
## partbounds: part_mol_id -> start_site,end_site (comma-separated string)
## source2name: source_id -> source_name
## atom2source:  atom_id -> source_id
## pathway2source: pathway_id -> source_id
## edgeptm: atomid, edgeid -> Set Of (UniProt id, position, aa, modification)
## componentptm: complexmolid, componentid -> Set Of (UniProt id, position, aa, modification)
## memberlabel: family_id, member_id -> Set Of label_value
## memberptm: family_id, member_id -> Set Of (UniProt id, position, aa, modification)

## molinst: normal mol inst string -> locally unique mol inst id
## rmolinst: locally unique mol inst id -> normal mol inst str
## molinstidx: highest mol inst id
## molinstsubtype         ## reflexive
## genericlabelsubtype    ## reflexive

## map: for mapping one set of ids to another: molid -> molid

my $lv;       ## PWLabel

######################################################################
sub new {
  my ($self, $labeldb) = @_;
  $lv = $labeldb;
  return bless {};
}

######################################################################
sub r_numerically { $b <=> $a };

######################################################################
sub numerically { $a <=> $b };

######################################################################
sub AllUsedLabels {
  my ($self) = @_;

  my %uses;
  for my $atom (@{ $self->Atoms() }) {
    $uses{$self->AtomType($atom)} = 1;
    for my $lvid (@{ $self->AtomCondition($atom) }) {
      $uses{$lvid} = 1;
    }
    for my $lvid (@{ $self->AtomNegativeCondition($atom) }) {
      $uses{$lvid} = 1;
    }
    for my $lvid (@{ $self->AtomLabel($atom) }) {
      $uses{$lvid} = 1;
    }
    for my $edge (@{ $self->Edges($atom) }) {
      $uses{$self->EdgeType($atom, $edge)} = 1;
      for my $lvid (@{ $self->MolLabel($atom, $edge) }) {
        $uses{$lvid} = 1;
      }
      for my $lvid (@{ $self->EdgeLabel($atom, $edge) }) {
        $uses{$lvid} = 1;
      }
      for my $x (@{ $self->EdgePTM($atom, $edge) }) {
        my ($protein_id, $position, $amino_acid, $modification_label_value,
            $modification_label_name) = @{ $x };
        $uses{$modification_label_value} = 1;
      }   
    }
  }
  my $COMPLEX_TYPE = $lv->StringToLabelValue("molecule-type", "complex");
  for my $molid (@{ $self->Mols() }) {
    $uses{$self->MolType($molid)} = 1;
    if ($self->MolType($molid) eq $COMPLEX_TYPE) {
      for my $component (@{ $self->Components($molid) }) {
        for my $lvid (@{ $self->ComponentLabel($molid, $component) }) {
          $uses{$lvid} = 1;
        }
        for my $x (@{ $self->ComponentPTM($molid, $component) }) {
          my ($protein_id, $position, $amino_acid, $modification_label_value,
              $modification_label_name) = @{ $x };
          $uses{$modification_label_value} = 1;
        }  
      }
    }
  }
  return \%uses;
}

######################################################################
## Source
######################################################################
sub AddSource {
  my ($self, $source_id, $source_name) = @_;
  $self->{source2name}{$source_id} = $source_name;
}

######################################################################
sub SourceName {
  my ($self, $source_id) = @_;
  if (defined $self->{source2name}{$source_id}) {
    return $self->{source2name}{$source_id};
  } else {
    return "";
  }
}

######################################################################
sub AddAtomSource {
  my ($self, $atom_id, $source_id) = @_;
  $self->{atom2source}{$atom_id} = $source_id;
}

######################################################################
sub AtomSource {
  my ($self, $atom_id) = @_;
  if (defined $self->{atom2source}{$atom_id}) {
    return $self->{atom2source}{$atom_id};
  } else {
    return "";
  }
}

######################################################################
sub AddPathwaySource {
  my ($self, $pathway_id, $source_id) = @_;
  $self->{pathway2source}{$pathway_id} = $source_id;
}

######################################################################
sub PathwaySource {
  my ($self, $pathway_id) = @_;
  if (defined $self->{pathway2source}{$pathway_id}) {
    return $self->{pathway2source}{$pathway_id};
  } else {
    return "";
  }
}

######################################################################
## abstraction
######################################################################

sub AddAbstraction {
  my ($self, $atom_id, $pid, $pathway_name, $ext_pathway_id) = @_;

  if ($pid) {
    $self->{abstraction}{pathway_id}{$atom_id} = $pid;
  }
  if ($ext_pathway_id) {
    if ($pid) {
      $self->SetPathwayExId($pid, $ext_pathway_id);
    }
    $self->{abstraction}{ext_pathway_id}{$atom_id} = $ext_pathway_id;
  }
  if ($pathway_name) {
    if ($pid) {
      $self->SetPathwayName($pid, $pathway_name);
    }
    $self->{abstraction}{pathway_name}{$atom_id} = $pathway_name;
  }
}

sub AbstractionId {
  my ($self, $atom_id) = @_;
  if (defined $self->{abstraction}{pathway_id}{$atom_id}) {
    return $self->{abstraction}{pathway_id}{$atom_id}; 
  } else {
    return undef;
  }
}

sub AbstractionExtId {
  my ($self, $atom_id) = @_;
  if (defined $self->{abstraction}{ext_pathway_id}{$atom_id}) {
    return $self->{abstraction}{ext_pathway_id}{$atom_id}; 
  } else {
    return undef;
  }
}

sub AbstractionName {
  my ($self, $atom_id) = @_;
  if (defined $self->{abstraction}{pathway_name}{$atom_id}) {
    return $self->{abstraction}{pathway_name}{$atom_id}; 
  } else {
    return undef;
  }
}

######################################################################
## Pathway
######################################################################
sub AddPathwayId {
  my ($self, $pid) = @_;

  $self->{pathway}{$pid}{pname}  = "";
  $self->{pathway}{$pid}{pexid}  = "";
  $self->{pathway}{$pid}{porg}   = "";
  $self->{pathway}{$pid}{psrcid} = "";
  $self->{pathway}{$pid}{pname}  = "";
}

sub AddPathwayReferences {
  my ($self, $pid, $ref) = @_;
  return $self->{pathwayreference}{$pid}{$ref} = 1;
}

sub PathwayReferences {
  my ($self, $pid) = @_;
  return [ keys %{ $self->{pathwayreference}{$pid} } ]
}

sub SetPathwayName {
  my ($self, $pid, $name) = @_;
  $self->{pathway}{$pid}{pname} = $name;
}

sub SetPathwayExId {
  my ($self, $pid, $ext_id) = @_;
  $self->{pathway}{$pid}{pexid} = $ext_id;
}

sub SetPathwayOrg {
  my ($self, $pid, $org) = @_;
  $self->{pathway}{$pid}{porg} = $org;
}

sub SetPathwaySrcId {
  my ($self, $pid, $srcid) = @_;
  $self->{pathway}{$pid}{psrcid} = $srcid;
}

sub SetIsSubnet {
  my ($self, $pid, $is_subnet) = @_;
  if ($is_subnet == 1 || $is_subnet eq "Y") {
    $self->{pathway}{$pid}{issubnet} = 1;
  } else {
    $self->{pathway}{$pid}{issubnet} = 0;
  }
}

sub PathwayName {
  my ($self, $pid) = @_;
  if (defined $self->{pathway}{$pid}) {
    return $self->{pathway}{$pid}{pname};
  } else {
    return "";
  }
}

sub PathwayId {
  my ($self) = @_;
  return [ keys %{ $self->{pathway} } ];
}

sub PathwayExId {
  my ($self, $pid) = @_;
  if (defined $self->{pathway}{$pid}) {
    return $self->{pathway}{$pid}{pexid};
  } else {
    return "";
  }
}

sub PathwayOrg {
  my ($self, $pid) = @_;
  if (defined $self->{pathway}{$pid}) {
    return $self->{pathway}{$pid}{porg};
  } else {
    return "";
  }
}

sub PathwaySrcId {
  my ($self, $pid) = @_;
  if (defined $self->{pathway}{$pid}) {
    return $self->{pathway}{$pid}{psrcid};
  } else {
    return "";
  }
}

sub PathwayIsSubnet {
  my ($self, $pid) = @_;
  if (defined $self->{pathway}{$pid}{issubnet} &&
    $self->{pathway}{$pid}{issubnet}) {
    return 1;
  } else {
    return 0;
  }
}


######################################################################
## Pathway layers
######################################################################
sub SetLayer {
  my ($self, $atom_id, $layer_id) = @_;
  $self->{layer}{$atom_id}{$layer_id} = 1;
}

sub LayersOf {
  my ($self, $atom_id) = @_;
  if (!defined $self->{layer}{$atom_id}) {
    return [];
  } else {
    return [ keys %{ $self->{layer}{$atom_id} } ];
  }
}

######################################################################
## Curators and Reviewers
######################################################################

sub AddCurator {
  my ($self, $pathway, $curator) = @_;
  $self->{pathwaycurator}{$pathway}{$curator} = 1;
}

sub Curators {
  my ($self, $pathway) = @_;
  if (defined $self->{pathwaycurator}{$pathway}) {
    return [ keys %{ $self->{pathwaycurator}{$pathway} } ];
  } else {
    return [ ];
  }
}

sub AddReviewer {
  my ($self, $pathway, $reviewer) = @_;
  $self->{pathwayreviewer}{$pathway}{$reviewer} = 1;
}

sub Reviewers {
  my ($self, $pathway) = @_;
  if (defined $self->{pathwayreviewer}{$pathway}) {
    return [ keys %{ $self->{pathwayreviewer}{$pathway} } ];
  } else {
    return [ ];
  }
}

######################################################################
## For atom collections from multiple pathways
######################################################################

sub AddAtomPathway {
  my ($self, $atom_id, $pathway_id) = @_;
  $self->{atompathway}{$atom_id}{$pathway_id} = 1;
}

sub AtomPathway {
  my ($self, $atom_id) = @_;
  return [ keys %{ $self->{atompathway}{$atom_id} } ];
}

######################################################################
## Subnets
######################################################################
sub AddPathwaySubnet {
  my ($self, $pathway, $subnet) = @_;

  $self->{pathwaysubnet}{$pathway}{$subnet} = 1;
}

######################################################################
sub AddSubnetAtom {
  my ($self, $subnet, $atom) = @_;

  $self->{subnetatom}{$subnet}{$atom} = 1;
}

######################################################################
sub Subnets {
  my ($self) = @_;
  if (defined $self->{subnetatom}) {
    return [ keys %{ $self->{subnetatom} } ];
  } else {
    return [ ];
  }
}

######################################################################
sub SubnetsInPathway {
  my ($self, $pathway) = @_;

  if (defined $self->{pathwaysubnet}{$pathway}) {
    return [ keys %{ $self->{pathwaysubnet}{$pathway} } ];
  } else {
    return [ ];
  }
}

######################################################################
sub AtomsInSubnet {
  my ($self, $subnet) = @_;

  if (defined $self->{subnetatom}{$subnet}) {
    return [ keys %{ $self->{subnetatom}{$subnet} } ];
  } else {
    return [ ];
  }
}

######################################################################
sub CollapseSubnets {
  my ($self) = @_;

  ## order of collapse does not matter

  my (%atomsubnet, %edge_types);

  my $SUBNET_TYPE  = $lv->StringToLabelValue("process-type", "subnet");
  my $ACTIVITY_STATE  = $lv->StringToLabelValue("activity-state",
      "activity-state");
  my $LOCATION  = $lv->StringToLabelValue("location", "location");

  ## first invert subnetatom

  for my $subnet (@{ $self->Subnets() }) {
    for my $atom (@{ $self->AtomsInSubnet($subnet) }) {
      $atomsubnet{$atom}{$subnet} = 1;
    }
  }
  
  ## Find molinsts that are not local to a single subnet;
  ## specify corresonding mol-subnet edges.
  ## Must traverse every edge of every atom in every
  ## subnet, but can quit an inner cycle when one non-local
  ## atom is found for a given molinst

  for my $subnet (@{ $self->Subnets() }) {
    for my $atom (@{ $self->AtomsInSubnet($subnet) }) {
      for my $edge (@{ $self->Edges($atom) }) {
        my $molinst = $self->MolInstId($atom, $edge);
        my ($mol, $labels) =
            DecodeMolInstString($self->{rmolinst}{$molinst});
        my $found = 0;
        for my $other_atom (@{ $self->MolUse($mol) }) {
          if (defined $atomsubnet{$other_atom}{$subnet}) {
            next;
          }
          for my $other_edge (@{ $self->Edges($other_atom) }) {
            if ($self->MolInstId($other_atom, $other_edge) eq $molinst) {
              $found++;
              $edge_types{$subnet}{$molinst}
                  {$self->EdgeType($atom, $edge)} = 1;
              last;
            }
          }
          if ($found) {
            last;
          }
        }
      }
    }
  }

  ## Now, construct the new, collapsed-subnet atoms
  for my $subnet (@{ $self->Subnets() }) {
    my $atom = "subnet_$subnet";
    $self->AddAtom($atom, $SUBNET_TYPE);
    if (! defined $edge_types{$subnet}) {
      next;
    }
    my $edge = 0;
    for my $molinst (keys %{ $edge_types{$subnet} }) {
      for my $edge_type (keys %{ $edge_types{$subnet}{$molinst} }) {
        $edge++;
        my ($mol, $labels) =
            DecodeMolInstString($self->{rmolinst}{$molinst});
        $self->AddEdge($atom, $edge, $edge_type, $mol);
        for my $lvid (@{ $labels }) {
          if ($lv->IsA($lvid, $ACTIVITY_STATE) ||
              $lv->IsA($lvid, $LOCATION)) {
            $self->AddMolLabel($atom, $edge, $lvid);
          }
        }
      }
    }
  }

  ## And delete the atoms in the collapsed subnets
  for my $subnet (@{ $self->Subnets() }) {
    $self->PruneAtoms($self->AtomsInSubnet($subnet));
  }
}

######################################################################
sub ValidateSubnets {
  my ($self) = @_;

  ## For all subnets defined in the model, check to see that all of
  ## the atoms referenced atoms are in the model

  for my $subnet (@{ $self->Subnets() }) {
    my $ok = 1;
    for my $atom (@{ $self->AtomsInSubnet($subnet) }) {
      if (! defined $self->{atomtype}{$atom}) {
        $ok = 0;
        ## can kill right here; loop works off a copy of the list
        delete $self->{subnetatom}{$subnet};
        last;
      }
    }
  }

  ## Also, compare all surviving subnets. Delete any subnet S1 that is
  ## subsumed by another surviving subnet S2.

  my %temp;
  for my $subnet (@{ $self->Subnets() }) {
    push @{ $temp{scalar(@{ $self->AtomsInSubnet($subnet) })} }, $subnet;
  }
  my (@subnets, @kill);
  for my $size (sort r_numerically keys %temp) {
    for my $subnet (@{ $temp{$size} }) {
      push @subnets, $subnet;
    }
  }
  for (my $s1 = 0; $s1 < @subnets - 1; $s1++) {
    if ($kill[$s1]) {
      next;
    }
    for (my $s2 = $s1 + 1; $s2 < @subnets; $s2++) {
      if ($kill[$s2]) {
        next;
      }
      my $inter = IntersectSetsAsHashes($self->{subnetatom}{$subnets[$s1]},
        $self->{subnetatom}{$subnets[$s2]});
      my $inter_size  = scalar(keys %{ $inter });
      if ($inter_size && $inter_size ==
          scalar(keys %{ $self->{subnetatom}{$subnets[$s2]} })) {
        $kill[$s2] = 1;
      }
    }
  }
  for (my $s = 0; $s < @subnets; $s++) {
    if ($kill[$s]) {
      delete $self->{subnetatom}{$subnets[$s]};
    }
  }
}

######################################################################
sub IntersectSetsAsHashes {
  my ($s1, $s2) = @_;
  my %s3;
  for my $s (keys %{ $s1 }) {
    if (defined $$s2{$s}) {
      $s3{$s} = 1;
    }
  }
  return \%s3;
}


######################################################################
## Post-Translational Modifications ptm
## for now: a ptm expression is simply a conjunction of terms, nothing
## more complicated.
######################################################################

sub PTMList {
  my ($list) = @_;

  my (@ptms, %ptm_order);

  for my $ptm (keys %{$list }) {
    my $x = [ split(",", $ptm) ];
    my ($protein_id, $position, $amino_acid, $modification_label_value,
        $modification_label_name) = @{ $x };
    push @{ $ptm_order{$protein_id}{$position} }, $x;
  }
  for my $protein_id (sort keys %ptm_order) {
    for my $position (sort numerically keys
        %{ $ptm_order{$protein_id} }) {
      push @ptms, @{ $ptm_order{$protein_id}{$position} };
    }
  }
  return \@ptms;
}

sub AddEdgePTM {
  my ($self, $atomid, $edgeid, $protein_id, $position, $amino_acid,
      $modification_label_value, $modification_label_name) = @_;

  my $ptm_string = join(",", $protein_id, $position, $amino_acid,
      $modification_label_value, $modification_label_name);
  $self->{edgeptm}{$atomid}{$edgeid}{$ptm_string} = 1;
}

sub EdgePTM {
  my ($self, $atomid, $edgeid) = @_;

  if (defined $self->{edgeptm}{$atomid}) {
    if (defined $self->{edgeptm}{$atomid}{$edgeid}) {
      return PTMList($self->{edgeptm}{$atomid}{$edgeid});
    } else {
      return [];
    }
  } else {
    return [];
  }
}

sub AddComponentPTM {
  my ($self, $cxid, $compid, $protein_id, $position, $amino_acid,
      $modification_label_value, $modification_label_name) = @_;

  my $ptm_string = join(",", $protein_id, $position, $amino_acid,
      $modification_label_value, $modification_label_name);
  $self->{componentptm}{$cxid}{$compid}{$ptm_string} = 1;
}

sub ComponentPTM {
  my ($self, $cxid, $compid) = @_;

  if (defined $self->{componentptm}{$cxid}) {
    if (defined $self->{componentptm}{$cxid}{$compid}) {
      return PTMList($self->{componentptm}{$cxid}{$compid});
    } else {
      return [];
    }
  } else {
    return [];
  }
}

sub NormalPTMString {
  my ($ptmset) = @_;

  my @ptms;
  for my $x (@{ $ptmset }) {
    my ($protein_id, $position, $amino_acid, $modification_label_value,
        $modification_label_name) = @{ $x };
    push @ptms, "$position$amino_acid$modification_label_value";
  }
  return join(",", @ptms);
}

######################################################################
## Molecule
######################################################################
sub AddMol {
  my ($self, $molid, $moltype) = @_;
  $self->{moltype}{$molid} = $moltype;
}

sub AddMolExId {
  my ($self, $molid, $idtype, $id) = @_;
  $self->{molexid}{$molid}{$idtype}{$id} = 1;
}

sub AddMolName {
  my ($self, $molid, $nametype, $name) = @_;
  $self->{molname}{$molid}{$nametype}{$name} = 1;
}

sub AddMolLabel {
  my ($self, $atomid, $edgeid, $label_value) = @_;
  $self->{mollabel}{$atomid}{$edgeid}{$label_value} = 1;
}

sub Mols {
  my ($self) = @_;
  return [ keys %{ $self->{moltype} } ];
}

sub UsedMols {
  my ($self) = @_;

  my $COMPLEX_TYPE = $lv->StringToLabelValue("molecule-type", "complex");
  my %temp;
  for my $mol (keys %{ $self->{moluse} }) {
    $self->ComplexDescendents($mol, \%temp);
    $temp{$mol} = 1;
    for my $x (@{ $self->FamilyParent($mol) }) {
      $temp{$x} = 1;
    }
    for my $x (@{ $self->FamilyChildren($mol) }) {
      $temp{$x} = 1;
    }
    for my $x (@{ $self->MolPart($mol) }) {
      $temp{$x} = 1;
    }
    for my $x (@{ $self->MolWhole($mol) }) {
      $temp{$x} = 1;
    }
  }
  my $more = 1;
  while ($more) {
    $more = 0;
    for my $x (keys %temp) {
      for my $y (
          @{ $self->FamilyParent($x) },
          @{ $self->FamilyChildren($x) },
          @{ $self->MolPart($x) },
          @{ $self->MolWhole($x) }
        ) {
        if (! defined $temp{$y}) {
          $more++;
          $temp{$y} = 1;
        }
      }
    }
  }
  return [ keys %temp ];
}

sub MolUse {
  my ($self, $molid) = @_;
  if (defined $self->{moluse}{$molid}) {
    return [ keys %{ $self->{moluse}{$molid} } ];
  } else {
    return [];
  }
}

sub MolType {
  my ($self, $molid) = @_;
  return $self->{moltype}{$molid};
}
sub MolExId {
  my ($self, $molid) = @_;
  return $self->{molexid}{$molid};
}
sub MolName {
  my ($self, $molid) = @_;
  return $self->{molname}{$molid};
}
sub MolLabel {
  my ($self, $atomid, $edgeid) = @_;
  return [ keys %{ $self->{mollabel}{$atomid}{$edgeid} } ];
}

sub PickMolName {
  my ($self, $molid) = @_;

  my ($label_sort, $moltype) = $lv->LabelValueToString($self->MolType($molid));

  if ($moltype eq "protein" || $moltype eq "rna") {
    my ($official, $preferred, $alias);
    my $h = $self->MolName($molid);
    if (defined $$h{"PF"}) {
      for my $name (keys %{ $$h{"PF"} }) {
        return $name;
      }
    } elsif (defined $$h{"OF"}) {
      for my $name (keys %{ $$h{"OF"} }) {
        return $name;
      }
    }
  }

  ##
  ## Arbitrary!!! pick the shortest name or the longest ext id
  ##

  my ($n0, $n1, $n2, $name, $id);
  $n0 = $self->MolName($molid);
  for $n1 (keys %{ $n0 }) {
    for $n2 (keys %{ $$n0{$n1} }) {
      if ((not defined $name) || (length($n2) < length($name))) {
        $name = $n2;
      }
    }
  }
  $n0 = $self->MolExId($molid);
  for $n1 (keys %{ $n0 }) {
    for $n2 (keys %{ $$n0{$n1} }) {
      if ($n1 eq "EC") {
        $id = "$n1:$n2";
#        $name = $id;
      } elsif ((not defined $id) || (length("$n1:$n2") > length($id))) {
        $id = "$n1:$n2";
      }
    }
  }
  if (defined $name) {
    return $name;
  } elsif (defined $id) {
    return $id;
  } else {
#    return "???"
    return ""
  }
}

######################################################################
sub AddMolMap {
  my ($self, $old, $new) = @_;

  if ($old ne $new) {
    $self->{map}{$old} = $new;
  }
}

######################################################################
sub MolMap {
  my ($self, $old) = @_;

  if (defined $self->{map}{$old}) {
    return $self->{map}{$old};
  } else {
    return undef;
  }
}

######################################################################
sub RemapMols {
  my ($self) = @_;

  my $pb = $self;

  if ((! defined $pb->{map}) ||
      (keys %{ $pb->{map} } == 0) ) {
    return;
  }

  my $map = $pb->{map};

  ##
  ##  moltype, molexid, molname, moluse, partbounds
  ##
  for my $HASH (
      "moltype",
      "molexid",
      "molname",
      "moluse",
      "componentlabel",
      "partbounds"
      ) {
    my %temp;
    for my $mol (keys %{ $pb->{$HASH} }) {
      if (defined $$map{$mol}) {
        $temp{$$map{$mol}} = $pb->{$HASH}{$mol};
      } else {
        $temp{$mol} = $pb->{$HASH}{$mol};
      }
    }
    $pb->{$HASH} = \%temp;
  }


  ##
  ## edge
  ##
  for my $HASH (
      "edge"
      ) {
    my %temp;
    for my $x1 (keys %{ $pb->{$HASH} }) {
      for my $x2 (keys %{ $pb->{$HASH}{$x1} }) {
        my $mol = $pb->{$HASH}{$x1}{$x2};
        if (defined $$map{$mol}) {
          $temp{$x1}{$x2} = $$map{$mol};
        } else {
          $temp{$x1}{$x2} = $mol;
        }
      }
    }
    $pb->{$HASH} = \%temp;
  }

  ##
  ## component
  ##
  for my $HASH (
      "component",
      "r_component"
      ) {
    my %temp;
    for my $mol1 (keys %{ $pb->{$HASH} }) {
      my $m1 = $mol1;
      if (defined $$map{$mol1}) {
        $m1 = $$map{$mol1};
      }
      for my $x1 (keys %{ $pb->{$HASH}{$mol1} }) {
        my $mol2 = $pb->{$HASH}{$mol1}{$x1};
        my $m2 = $mol2;
        if (defined $$map{$mol2}) {
          $m2 = $$map{$mol2};
        }
        $temp{$m1}{$x1} = $m2;
      }
    }
    $pb->{$HASH} = \%temp;
  }

  ##
  ## family, r_family, partof, r_partof
  ##
  for my $HASH (
      "family",
      "r_family",
      "partof",
      "r_partof"
      ) {
    my %temp;
    for my $mol1 (keys %{ $pb->{$HASH} }) {
      my $m1 = $mol1;
      if (defined $$map{$mol1}) {
        $m1 = $$map{$mol1};
      }
      for my $mol2 (keys %{ $self->{$HASH}{$mol1} }) {
        my $mol2 = $pb->{$HASH}{$mol1}{$mol1};
        my $m2 = $mol2;
        if (defined $$map{$mol2}) {
          $m2 = $$map{$mol2};
        }
        $temp{$m1}{$m2} = $pb->{$HASH}{$mol1}{$mol2};
      }
    }
    $pb->{$HASH} = \%temp;
  }

}

######################################################################
# Converting interaction to string
######################################################################

sub StringToAtom {
  my ($self, $str) = @_;

  my ($a, $e, @edges, @parsed_edges, $atom, $atom_label, $atom_condition,
      $edge, $edge_label, $mol_label, $molid, $edge_type, $atom_negative_condition);
  my ($name, $value);

  ($a, @edges) = split("::", $str);
  if ($a =~ /^(\d+):\[([0-9,]*)\]\[([0-9,]*)\]$/) {
    ($atom, $atom_label, $atom_condition) = ($1, $2, $3);
  } else {
    print STDERR "bad atom string header: $a\n";
    return 0;
  }
  for $e (@edges) {
    if ($e eq "") {
      next;
    }
    if ($e =~ /^(\d+):(\d+)\[([0-9,]*)\]\[([0-9,]*)\]$/) {
      ($edge, $molid, $edge_label, $mol_label) = ($1, $2, $3, $4);
      push @parsed_edges, [ $edge, $molid, $edge_label, $mol_label ];
    } else {
      print STDERR "bad edge string: $e\n";
      return 0;
    }
  }
  for my $x (split(",", $atom_label)) {
    ($name, $value) = $lv->LabelValueToString($x);
    if ($name eq "process-type") {
      $self->AddAtom($atom, $x);
    }
    $self->AddAtomLabel($atom, $x);
  }
  for my $x (split(",", $atom_condition)) {
    $self->AddAtomCondition($atom, $x);
  }
  for (@parsed_edges) {
    ($edge, $molid, $edge_label, $mol_label) = @{ $_ };
    for my $x (split(",", $edge_label)) {
      ($name, $value) = $lv->LabelValueToString($x);
      if ($name eq "edge-type") {
        $self->AddEdge($atom, $edge, $edge_type, $molid);
      }
      $self->AddEdgeLabel($atom, $edge, $x);
    }
    for my $x (split(",", $mol_label)) {
      $self->AddMolLabel($atom, $edge, $x);
    }
  }

  return 1;
}

sub LabelListToString {
  my ($self, $list) = @_;

# label_list: [ comma-separated list of int ]
  if ($list && @{ $list } > 0) {
    return "[" . join(",", sort @{ $list }) . "]";
  }
  return "[]";
}

sub EdgeToString {
  my ($self, $atom, $edge) = @_;

# edge: edgeid:molid:edge_label_list:mol_label_list
  my $molid = $self->EdgeMol($atom, $edge);
  my $mol_label = $self->LabelListToString($self->MolLabel($atom, $edge));
  my $edge_label = $self->LabelListToString($self->EdgeLabel($atom, $edge));
  return "$edge:$molid" . $edge_label . $mol_label;

}

sub AtomToString {
  my ($self, $atom) = @_;

# atom: atomid::atom_label_list::atom_condition_list::edge1::edge2::...
  my @edges;
  for my $edge (@{ $self->Edges($atom) }) {
    push @edges, $self->EdgeToString($atom, $edge);
  }
  
  my $atom_label = $self->LabelListToString($self->AtomLabel($atom));
  my $atom_condition = $self->LabelListToString($self->AtomCondition($atom));
  my $atom_negative_condition = $self->LabelListToString($self->AtomNegativeCondition($atom));
  return "$atom:$atom_label" . $atom_condition . "::" .
    join("::", @edges);

}


######################################################################
## Generic label string
######################################################################
sub NormalLabelString {
  my ($labellist) = @_;
  return join(",", sort @{ $labellist });
}

sub BuildGenericLabelCache {
  my ($self) = @_;

  my (%strings, @strings, $i, $j, @si, @sj);
  my ($atom, $edge);
  for $atom (@{ $self->Atoms }) {
    my $x = NormalLabelString($self->AtomLabel($atom));
    $strings{NormalLabelString($self->AtomLabel($atom))} = 1;
  }
  ##
  ## Make the relation reflexive
  ##
  for (keys %strings) {
    $self->{genericlabelsubtype}{"$_:$_"} = 1;
  }
  @strings = keys %strings;
  for ($i = 0; $i < @strings-1; $i++) {
    for ($j = $i + 1; $j < @strings; $j++) {
      @si = split ",", $strings[$i];
      @sj = split ",", $strings[$j];
      if (LabelSetIsSubType(\@si, \@sj)) {
        $self->{genericlabelsubtype}{"$strings[$i]:$strings[$j]"} = 1;
      }
      if (LabelSetIsSubType(\@sj, \@si)) {
        $self->{genericlabelsubtype}{"$strings[$j]:$strings[$j]i]"} = 1;
      }
    }
  }

  undef %strings;
  for $atom (@{ $self->Atoms }) {
    for $edge (@{ $self->Edges($atom) }) {
      $strings{NormalLabelString($self->EdgeLabel($atom, $edge))} = 1;
    }
  }
  ##
  ## Make the relation reflexive
  ##
  for (keys %strings) {
    $self->{genericlabelsubtype}{"$_:$_"} = 1;
  }
  @strings = keys %strings;
  for ($i = 0; $i < @strings-1; $i++) {
    $self->{genericlabelsubtype}{"$strings[$i]:$strings[$i]"} = 1;
    for ($j = $i + 1; $j < @strings; $j++) {
      @si = split ",", $strings[$i];
      @sj = split ",", $strings[$j];
      if (LabelSetIsSubType(\@si, \@sj)) {
        $self->{genericlabelsubtype}{"$strings[$i]:$strings[$j]"} = 1;
      }
      if (LabelSetIsSubType(\@sj, \@si)) {
        $self->{genericlabelsubtype}{"$strings[$j]:$strings[$j]i]"} = 1;
      }
    }
  }
}

######################################################################
sub CollapseMols {
  my ($self) = @_;

  undef $self->{mollabel};
}
######################################################################
## Molecule instance
## molid X Seq Of label value id: defines a molecule instance
## Need to collapse these.
######################################################################

sub MolInstId {
  my ($self, $atom, $edge) = @_;
  my $s = $self->NormalMolInstString($atom, $edge);
  if (defined $self->{molinst}{$s}) {
    return $self->{molinst}{$s};
  } else {
    print STDERR "request for nonexistent mol inst id: $s\n";
    return 0;
  }
}

sub MolInstToString {
  my ($molid, $labelset, $ptmset) = @_;
  return join(":", $molid, NormalLabelString($labelset),
     NormalPTMString($ptmset));
}

sub MolInstIdToString {
  my ($self, $molinstid) = @_;
  if (defined $self->{rmolinst}{$molinstid}) {
    return $self->{rmolinst}{$molinstid};
  } else {
    print STDERR "request for unknown mol inst id $molinstid\n";
    return "";
  }
}

sub NormalMolInstString {
  my ($self, $atomid, $edgeid) = @_;

  return MolInstToString($self->EdgeMol($atomid, $edgeid),
      $self->MolLabel($atomid, $edgeid),
      $self->EdgePTM($atomid, $edgeid));
}

sub DecodeMolInstString {
  my ($s) = @_;
  my ($m, $labels, $ptms) = split ":", $s;
  my @labels = split ",", $labels;
  my @ptms   = split ",", $ptms;
  return ($m, \@labels, \@ptms);
}

sub AddMolInstToCache {
  my ($self, $value) = @_;
  if (defined $self->{molinst}{$value}) {
    return $self->{molinst}{$value};
  } else {
    $self->{molinstidx}++;
    $self->{molinst}{$value} = $self->{molinstidx};
    $self->{rmolinst}{$self->{molinstidx}} = $value;
    return $self->{molinstidx};
  }
}

sub CachedMolInstSubtype {
  my ($self, $a, $b) = @_;
  if (defined $self->{molinstsubtype}{"$a,$b"}) {
    return 1;
  } else {
    return 0;
  }
}

sub CheckForSubtype {
  my ($self, $molinstset) = @_;
  
  my @set = keys %{ $molinstset };
  my ($i, $j, $si, $sj);
  for ($i = 0; $i < @set-1; $i++) {
    $si = $self->MolInstIdToString($set[$i]);
    for ($j = $i + 1; $j < @set; $j++) {
      $sj = $self->MolInstIdToString($set[$j]);
#print STDERR "si = $si; sj = $sj\n";
      if ($self->MolInstIsSubType($si, $sj)) {
        $self->{molinstsubtype}{"$set[$i],$set[$j]"} = 1;
#print STDERR "found forward subtype $set[$i] -> $set[$j]\n";
      }
      if ($self->MolInstIsSubType($sj, $si)) {
        $self->{molinstsubtype}{"$set[$j],$set[$i]"} = 1;
#print STDERR "found reverse subtype $set[$j] -> $set[$i]\n";
      }
    }
  }
}

sub BuildMolInstCache {
  my ($self) = @_;
  my ($atom, $edge, $val, @strings);

  my (%temp, $mi);

  for $atom (@{ $self->Atoms }) {
    for $edge (@{ $self->Edges($atom) } ) {
      $val = $self->NormalMolInstString($atom, $edge);
##print STDERR "BuildMolInstCache: $val\n";
      $mi = $self->AddMolInstToCache($val);
      $temp{$self->EdgeMol($atom, $edge)}{$mi} = 1;
    }
  }

  ##
  ## Make the relation reflexive
  ##
  for (values %{ $self->{molinst} }) {
    $self->{molinstsubtype}{"$_,$_"} = 1;
  }

  for (keys %temp) {
    $self->CheckForSubtype($temp{$_});
  }

##
## for now, don't build the subtype relation -- two slow
##

##  @strings = keys %{ $self->{molinst} };
##
##  my ($i, $j);
##  for ($i = 0; $i < @strings-1; $i++) {
##    for ($j = $i + 1; $j < @strings; $j++) {
##      if ($self->MolInstIsSubType($strings[$i], $strings[$j])) {
##        $self->{molinstsubtype}{
##            "$self->{molinst}{$strings[$i]},$self->{molinst}{$strings[$j]}"
##            } = 1;
##      }
##      if ($self->MolInstIsSubType($strings[$j], $strings[$i])) {
##        $self->{molinstsubtype}{
##            "$self->{molinst}{$strings[$j]},$self->{molinst}{$strings[$i]}"
##            } = 1;
##      }
##    }
##  }
}

sub MolInstIsSubType {
  my ($self, $s1, $s2) = @_;
  # inputs are normal strings for molecule instances

  my ($m1, $m2, $lab1, $lab2, $p1);
  if ($s1 eq $s2) {
    return 1;
  }
  ($m1, $lab1) = DecodeMolInstString($s1);
  ($m2, $lab2) = DecodeMolInstString($s2);
  $p1 = $self->FamilyAncestors($m1);
  if ($m1 != $m2 && (not defined $$p1{$m2})) {
    return 0;
  }
  if (LabelSetIsSubType($lab1, $lab2)) {
    return 1;
  } else {
    return 0;
  }
}

######################################################################
sub Clans {
  my ($self) = @_;
  return $self->{clans};
}

######################################################################
sub BuildClanList {
  my ($self) = @_;

  ## dumb, but...
  ## this repeats processing from BuildMolInstCache
  ## but we don't want to build the initial clan list
  ## until after dropping duplicate atoms and detection of
  ## dup atoms depends on the molinstcache

  my ($atom, $edge, $val);

  my ($mi);

  for $atom (@{ $self->Atoms }) {
    my @molinstlist;
    my @mp;
    for $edge (@{ $self->Edges($atom) } ) {
      $val = $self->NormalMolInstString($atom, $edge);
      $mi = $self->AddMolInstToCache($val);
      push @molinstlist, $mi;
    }
    if ($self->IsMacroProcess($atom)) {
      push @mp, $self->{macroprocess}{$atom};
    }
    for my $c (@{ $self->AtomCondition($atom) }) {
      push @mp, $c;
    }
    push @{ $self->{clans} }, new Clan([ $atom ], \@molinstlist, \@mp) ;
  }
  if (defined $self->{clans}) {
    $self->{clans} = Clan::PartitionClans($self->{clans});
  } else {
    $self->{clans} = [];
  }
}

######################################################################
sub ForceIntoOneClan {
  my ($self) = @_;

  ## if we really want to put everything in one dot output...

  my ($atom, $edge, $val);

  my @atoms;
  my @molinstlist;
  my @mp;

  my ($mi);

  for $atom (@{ $self->Atoms }) {
    for $edge (@{ $self->Edges($atom) } ) {
      $val = $self->NormalMolInstString($atom, $edge);
      $mi = $self->AddMolInstToCache($val);
      push @molinstlist, $mi;
    }
    if ($self->IsMacroProcess($atom)) {
      push @mp, $self->{macroprocess}{$atom};
    }
    for my $c (@{ $self->AtomCondition($atom) }) {
      push @mp, $c;
    }
    push @atoms, $atom;
  }
  push @{ $self->{clans} }, new Clan(\@atoms, \@molinstlist, \@mp) ;
}

######################################################################
## Edge
######################################################################
sub AddEdge {
  my ($self, $atomid, $edgeid, $edgetype, $molid) = @_;
  $self->{edge}{$atomid}{$edgeid} = $molid;
  $self->{moluse}{$molid}{$atomid} = 1;
  $self->{edgetype}{$atomid}{$edgeid} = $edgetype;
}

sub AddEdgeLabel {
  my ($self, $atomid, $edgeid, $label_value) = @_;
  $self->{edgelabel}{$atomid}{$edgeid}{$label_value} = 1;
}

sub Edges {
  my ($self, $atomid) = @_;
  return [ keys %{ $self->{edge}{$atomid} } ];
}
sub EdgeType {
  my ($self, $atomid, $edgeid) = @_;
  return $self->{edgetype}{$atomid}{$edgeid};
}
sub EdgeMol {
  my ($self, $atomid, $edgeid) = @_;
  return $self->{edge}{$atomid}{$edgeid};
}
sub EdgeLabel {
  my ($self, $atomid, $edgeid) = @_;
  return [ keys %{ $self->{edgelabel}{$atomid}{$edgeid} } ];
}

######################################################################
## Atom
######################################################################
sub AddAtom {
  my ($self, $atomid, $processtype) = @_;
  $self->{atomtype}{$atomid} = $processtype;
  $self->{atomlabel}{$atomid}{$processtype} = 1;
}

sub ResetAtomType {
  my ($self, $atomid, $processtype) = @_;
  delete $self->{atomlabel}{$atomid}{$self->AtomType($atomid)};
  $self->AddAtom($atomid, $processtype);
}

sub AddAtomLabel {
  my ($self, $atomid, $label_value) = @_;
  $self->{atomlabel}{$atomid}{$label_value} = 1;
}

sub AddAtomCondition {
  my ($self, $atomid, $label_value) = @_;
  $self->{atomcondition}{$atomid}{$label_value} = 1;
}

sub AddAtomNegativeCondition {
  my ($self, $atomid, $label_value) = @_;
  $self->{atomnegativecondition}{$atomid}{$label_value} = 1;
}

sub AddAtomEvidence {
  my ($self, $atomid, $evidence_code) = @_;
  $self->{atomevidence}{$atomid}{$evidence_code} = 1;
}

sub AddAtomReferences {
  my ($self, $atomid, $ref) = @_;
  $self->{atomreference}{$atomid}{$ref} = 1;
}

sub AddAtomNotes {
  my ($self, $atomid, $note) = @_;
  $self->{atomnote}{$atomid}{$note} = 1;
}

sub Atoms {
  my ($self) = @_;
  return [ keys %{ $self->{atomtype} } ];
}

sub AtomType {
  my ($self, $atomid) = @_;
  return $self->{atomtype}{$atomid};
}

sub HasAtom {
  my ($self, $atomid) = @_;
  if (defined $self->{atomtype}{$atomid}) {
    return 1;
  } else {
    return 0;
  }
}

sub AtomLabel {
  my ($self, $atomid) = @_;
  return [ keys %{ $self->{atomlabel}{$atomid} } ]
}

sub AtomCondition {
  my ($self, $atomid) = @_;
  return [ keys %{ $self->{atomcondition}{$atomid} } ]
}

sub AtomNegativeCondition {
  my ($self, $atomid) = @_;
  return [ keys %{ $self->{atomnegativecondition}{$atomid} } ]
}

sub AtomEvidence {
  my ($self, $atomid) = @_;
  if (defined $self->{atomevidence}{$atomid}) {
    return [ keys %{ $self->{atomevidence}{$atomid} } ]
  } else {
    return [];
  }
}

sub AtomReferences {
  my ($self, $atomid) = @_;
  return [ keys %{ $self->{atomreference}{$atomid} } ]
}

sub AtomNotes {
  my ($self, $atomid) = @_;
  return [ keys %{ $self->{atomnote}{$atomid} } ]
}

sub MacroProcesses {
  my ($self) = @_;
  return $self->{macroprocess}
}

sub IsMacroProcess {
  my ($self, $atomid) = @_;
  if (defined $self->{macroprocess}{$atomid}) {
    return 1;
  } else {
    return 0;
  }
}

######################################################################
## Complex
######################################################################
sub AddComponent {
  my ($self, $complexmolid, $componentid, $componentmolid) = @_;
  $self->{component}{$complexmolid}{$componentid} = $componentmolid;
  $self->{r_component}{$componentmolid}{$complexmolid} = 1;
}

sub AddComponentLabel {
  my ($self, $complexmolid, $componentid, $label_value) = @_;
  $self->{componentlabel}{$complexmolid}{$componentid}{$label_value} = 1;
}

sub ComponentLabel {
  my ($self, $complexmolid, $componentid) = @_;
  return [ keys %{ $self->{componentlabel}{$complexmolid}{$componentid} } ];
}

sub Components {
  my ($self, $complexid) = @_;
  return [ keys %{ $self->{component}{$complexid} } ];
}

sub ComponentMol {
  my ($self, $complexid, $seqid) = @_;
  return $self->{component}{$complexid}{$seqid};
}

######################################################################
## Molecule Parts
######################################################################

sub AddMolPart {
  my ($self, $part, $whole, $start_site, $end_site) = @_;

  $self->{partof}{$whole}{$part} = 1;
  $self->{r_partof}{$part}{$whole} = 1;
  $self->{partbounds}{$part} = join(",", $start_site, $end_site);
}

sub PartBounds {
  my ($self, $part) = @_;
  if (defined $self->{partbounds}{$part}) {
    return $self->{partbounds}{$part};
  } else {
    return "";
  }
}

sub MolPart {
  my ($self, $whole_mol_id) = @_;

  if (defined $self->{partof}{$whole_mol_id}) {
    return [ keys %{ $self->{partof}{$whole_mol_id} } ];
  } else {
    return [];
  }
}

sub MolWhole {
  my ($self, $part_mol_id) = @_;

  if (defined $self->{r_partof}{$part_mol_id}) {
    return [ keys %{ $self->{r_partof}{$part_mol_id} } ];
  } else {
    return [];
  }
}

######################################################################
## Molecule families
######################################################################

sub AddFamilyMember {
  my ($self, $child, $parent) = @_;

  $self->{family}{$child}{$parent} = 1;
  $self->{r_family}{$parent}{$child} = 1;
}

sub FamilyParentClosure {
  my ($self, $child, $accum) = @_;
  if (not defined $$accum{$child}) {
    $$accum{$child} = 1;
    for (keys %{ $self->{family}{$child} }) {
      $$accum{$_} = 1;
      $self->FamilyParentClosure($_, $accum);
    }
  }
}

sub FamilyChildren {
  my ($self, $parent) = @_;
  if (defined $self->{r_family}{$parent}) {
    return [ keys %{ $self->{r_family}{$parent} } ]
  } else {
    return []
  }
}

sub FamilyParent {
  my ($self, $child) = @_;
  if (defined $self->{family}{$child}) {
    return [keys %{ $self->{family}{$child} } ]
  } else {
    return [];
  }
}

sub FamilyAncestors {
  my ($self, $child) = @_;
  my %temp;
  $self->FamilyParentClosure($child, \%temp);
  return \%temp;
}

sub AddMemberLabel {
  my ($self, $family, $member, $lvid) = @_;
  $self->{memberlabel}{$family}{$member}{$lvid} = 1;
}

sub MemberLabel {
  my ($self, $family, $member) = @_;
  if (defined $self->{memberlabel}{$family}{$member}) {
    return [keys %{ $self->{memberlabel}{$family}{$member} } ];
  } else {
    return [];
  }
}

sub AddMemberPTM {
  my ($self, $family, $member, $protein_id, $position, $amino_acid,
      $modification_label_value, $modification_label_name) = @_;

  my $ptm_string = join(",", $protein_id, $position, $amino_acid,
      $modification_label_value, $modification_label_name);
  $self->{memberptm}{$family}{$member}{$ptm_string} = 1;
}

sub MemberPTM {
  my ($self, $family, $member) = @_;

  if (defined $self->{memberptm}{$family}) {
    if (defined $self->{memberptm}{$family}{$member}) {
      return PTMList($self->{memberptm}{$family}{$member});
    } else {
      return [];
    }
  } else {
    return [];
  }
}


######################################################################
## Prune
######################################################################

######################################################################
sub PruneMol {
  my ($self, $prunable) = @_;

  my ($atom, $edgeseq);
  my ($molid);

  for $atom (keys %{ $self->{edge} }) {
    for $edgeseq (keys %{ $self->{edge}{$atom} }) {
      if (defined $$prunable{$self->{edge}{$atom}{$edgeseq}}) {
        $molid = $self->{edge}{$atom}{$edgeseq};
        delete $self->{moluse}{$molid}{$atom};
        delete $self->{edge}{$atom}{$edgeseq};
        delete $self->{edgetype}{$atom}{$edgeseq};
        delete $self->{edgelabel}{$atom}{$edgeseq};
        delete $self->{edgeptm}{$atom}{$edgeseq};
      }
    }
  }

}

######################################################################
sub PruneAtoms {
  my ($self, $atoms, $replace) = @_;

  my ($edge, $pid);
  my (%atomsubnet);

  for my $subnet (@{ $self->Subnets() }) {
    for my $atom (@{ $self->AtomsInSubnet() }) {
      $atomsubnet{$atom}{$subnet} = 1;
    }
  }

#print STDERR "In PruneAtoms\n";
  for my $atomid (@{ $atoms }) {

#print STDERR "PruneAtoms: pruning atomid = $atomid\n";

    my $src = $self->AtomSource($atomid);
    my $nas = $self->NormalAtomString($atomid);

    if (defined $self->{atompathway}{$atomid}) {
      for $pid (keys %{ $self->{atompathway}{$atomid} }) {
        if (defined $$replace{$src}{$nas}) {
          $self->{atompathway}{$$replace{$src}{$nas}}{$pid} = 1;
        }
      }
      delete $self->{atompathway}{$atomid};
    }

    if (defined $atomsubnet{$atomid}) {
      for my $subnet (keys %{ $atomsubnet{$atomid} }) {
        if (defined $$replace{$src}{$nas}) {
          $self->{subnetatom}{$subnet}{$$replace{$src}{$nas}} = 1;
        }
        delete $self->{subnetatom}{$subnet}{$atomid};
      }
    }

    if (defined $self->{atomtype}{$atomid}) {
      delete $self->{atomtype}{$atomid};
    }
    if (defined $self->{atomlabel}{$atomid}) {
      delete $self->{atomlabel}{$atomid};
    }
    if (defined $self->{edge}{$atomid}) {
      for $edge (@{ $self->Edges($atomid) }) {
        delete $self->{moluse}{$self->{edge}{$atomid}{$edge}}{$atomid};
      }
      delete $self->{edge}{$atomid};
    }
    if (defined $self->{edgelabel}{$atomid}) {
      delete $self->{edgelabel}{$atomid};
    }
    if (defined $self->{edgeptm}{$atomid}) {
      delete $self->{edgeptm}{$atomid};
    }
    if (defined $self->{edgetype}{$atomid}) {
      delete $self->{edgetype}{$atomid};
    }
    if (defined $self->{mollabel}{$atomid}) {
      delete $self->{mollabel}{$atomid};
    }
    if (defined $self->{atomcondition}{$atomid}) {
      delete $self->{atomcondition}{$atomid};
    }
    if (defined $self->{atomevidence}{$atomid}) {
      for my $ev (@{ $self->AtomEvidence($atomid) }) {
        $self->{atomevidence}{$$replace{$src}{$nas}}{$ev} = 1;
      }
      delete $self->{atomevidence}{$atomid};
    }
    if (defined $self->{atomreference}{$atomid}) {
      for my $ref (@{ $self->AtomReferences($atomid) }) {
        $self->{atomreference}{$$replace{$src}{$nas}}{$ref} = 1;
      }
      delete $self->{atomreference}{$atomid};
    }
  }

  ##
  ## not really necessary to prune orphaned mols since they are
  ## not reachable in an ordinary traversal
  ##
#  PruneOrphanedMols($self);
}

######################################################################
sub PruneOrphanedMols {
  my ($self) = @_;
  my (%live);
  my ($atomid, $edgeid);
  for $atomid (keys %{ $self->{edge} }) {
    for $edgeid (keys %{ $self->{edge}{$atomid} }) {
      my $mol = $self->{edge}{$atomid}{$edgeid};
      $live{$mol} = 1;
      if (defined $self->{family}{$mol}) {
        for my $member (keys %{ $self->{family}{$mol} }) {
          $live{$member} = 1;
        }
      }
      if (defined $self->{r_partof}{$mol}) {
        for my $whole (keys %{ $self->{r_partof}{$mol} }) {
          $live{$whole} = 1;
        }
      }
    }
  }
  my $more = 1;
  while ($more) {
    $more = 0;
    for my $cxid (keys %{ $self->{component} }) {
      if (defined $live{$cxid}) {
        for my $cseq (keys %{ $self->{component}{$cxid} }) {
          if (! defined $live{$self->{component}{$cxid}{$cseq}} ) {
            $more = 1;
            my $mol = $self->{component}{$cxid}{$cseq};
            $live{$mol} = 1;
            if (defined $self->{family}{$mol}) {
              for my $member (keys %{ $self->{family}{$mol} }) {
                $live{$member} = 1;
              }
            }
            if (defined $self->{r_partof}{$mol}) {
              for my $whole (keys %{ $self->{r_partof}{$mol} }) {
                $live{$whole} = 1;
              }
            }
          }
        }
      }
    }
  }
  for my $molid (keys %{ $self->{moltype} }) {
    if (! defined $live{$molid}) {

#print STDERR "PruneOrphanedMols: pruning molid = $molid\n";

      delete $self->{moltype}{$molid};
      if (defined $self->{moluse}{$molid}) {
        delete $self->{moluse}{$molid};
      }
      if (defined $self->{molname}{$molid}) {
        delete $self->{molname}{$molid};
      }
      if (defined $self->{molextid}{$molid}) {
        delete $self->{molextid}{$molid};
      }
      if (defined $self->{family}{$molid}) {
        for my $member (keys %{ $self->{family}{$molid} }) {
          delete $self->{r_family}{$member}{$molid};
        }
        delete $self->{family}{$molid};
      }
      if (defined $self->{r_partof}{$molid}) {
        for my $whole (keys %{ $self->{r_partof}{$molid} }) {
          delete $self->{partof}{$whole}{$molid}
        }
        delete $self->{r_partof}{$molid};
      }
    }
  }
  for my $cxid (keys %{ $self->{component} }) {
    if (! defined $live{$cxid}) {
      delete $self->{component}{$cxid};
      delete $self->{componentlabel}{$cxid};
      delete $self->{componentptm}{$cxid};
    }
  }
}

######################################################################
sub PruneDuplicateAtoms {
  my ($self) = @_;

  my ($a, $s, %prunable, %atoms_seen, $src);
  for $a (sort @{ $self->Atoms() }) {
    $s = $self->NormalAtomString($a);
    $src = $self->AtomSource($a);
    if (defined $atoms_seen{$src}{$s}) {
      $prunable{$a} = $atoms_seen{$src}{$s};
#      print STDERR "atom $a duplicates $atoms_seen{$src}{$s}\n";
    } else {
      $atoms_seen{$src}{$s} = $a;
    }
  }

  PruneAtoms($self, [ keys %prunable ], \%atoms_seen);

}

######################################################################
sub IdentifyMacroProcesses {
  my ($self) = @_;

  my ($atom_id, $atom_type);
  my ($name, $value);

  my $macro_process_type =
      $lv->PWLabel::StringToLabelValue("process-type", "macroprocess");

  my $function_type =
      $lv->PWLabel::StringToLabelValue("function", "function");

  for $atom_id (@{ $self->Atoms() }) {
    $atom_type = $self->AtomType($atom_id);
    ##
    ## Problem: "transcription" is a GO BP term AND a basic atom type
    ##
    ($name, $value) = $lv->PWLabel::LabelValueToString($atom_type);
    if ($value eq "transcription") {
      next;
    }
    if ($lv->PWLabel::IsA($atom_type, $macro_process_type)) {
      $self->{macroprocess}{$atom_id} = $atom_type;
    } elsif ($lv->PWLabel::IsA($atom_type, $function_type)) {
      $self->{macroprocess}{$atom_id} = $atom_type;
    }
  }
}


######################################################################
## equivalencing complexes
######################################################################

sub Max {
  my $max = 0;
  for my $x (@_) {
    if ($x > $max) {
      $max = $x;
    }
  }
  return $max;
}

######################################################################
sub FlattenComplexes {
  my ($self) = @_;

  ## NOTE: if there were any labels on a complex component that was
  ## itself a complex, these labels will be lost -- no place to put them

  my $complex_type = $lv->StringToLabelValue("molecule-type", "complex");
  my (%child2parent, %parent2child);

  for my $cx (@{ $self->Mols() }) {
    if ($self->MolType($cx) eq $complex_type) {
      for my $comp (@{ $self->Components($cx) }) {
        my $compmol = $self->ComponentMol($cx, $comp);
        if ($self->MolType($compmol) eq $complex_type &&
            @{ $self->Components($compmol) } > 0) {  ## skip empties
          $child2parent{$compmol}{$cx} = 1;
          $parent2child{$cx}{$compmol} = 1;
        }
      }
    }
  }

  my $more = 1;
  while ($more) {
    $more = 0;
    for my $cx (keys %parent2child) {
      my @fix;
      for my $comp (@{ $self->Components($cx) }) {
        my $compmol = $self->ComponentMol($cx, $comp);
        if ($self->MolType($compmol) eq $complex_type &&
            @{ $self->Components($compmol) } > 0) {
          push @fix, $comp;
        }
      }
      for my $comp (@fix) {
        my $n = Max(@{ $self->Components($cx) });
        my $child = $self->ComponentMol($cx, $comp);
        for my $addcomp (@{ $self->Components($child) }) {
          my $addmol = $self->ComponentMol($child, $addcomp);
          $n++;
          $self->AddComponent($cx, $n, $addmol);
          if ($self->MolType($addmol) eq $complex_type &&
            @{ $self->Components($addmol) } > 0) {
            $more++;
          }
          for my $lvid (@{ $self->ComponentLabel($child, $addcomp) }) {
            $self->AddComponentLabel($cx, $n, $lvid);
          }
          for my $ptm (
              @{ $self->ComponentPTM($child, $addcomp) }) {
            $self->AddComponentPTM($cx, $n, @{ $ptm });
          }
        }
        ## delete from the parent the child complex that
        ## was just expanded
        delete $self->{component}{$cx}{$comp};
        delete $self->{componentlabel}{$cx}{$comp};
        delete $self->{componentptm}{$cx}{$comp};
      }
    }
  }
}

######################################################################
sub EqComplex {
  my ($pa, $pb, $a_cx, $b_cx) = @_;

  # check for same number of components
  if (@{ $pa->Components($a_cx) } != @{ $pb->Components($b_cx) }) {
    return 0;
  }

  if (@{ $pa->Components($a_cx) } == 0) {
    return NamesIntersect($pa, $pb, $a_cx, $b_cx);
  }

  # index (labels and all) the "standard" candidate
  my %a_comps;
  for my $comp (@{ $pa->Components($a_cx) }) {
    my $str = $pa->ComponentMol($a_cx, $comp) . ":" .
        join(",", sort @{ $pa->ComponentLabel($a_cx, $comp) },
        Pathway::NormalPTMString($pa->ComponentPTM($a_cx, $comp)));
    $a_comps{$str} = 1;
  }

  # now see if there is a match for each component of the "b" mol
  for my $comp (@{ $pb->Components($b_cx) }) {
    my $compmol = $pb->ComponentMol($b_cx, $comp);
    if (defined $pb->{map}{$compmol}) {
      $compmol = $pb->{map}{$compmol};
    }
    my $str = $compmol . ":" .
        join(",", sort @{ $pb->ComponentLabel($b_cx, $comp) },
        Pathway::NormalPTMString($pb->ComponentPTM($b_cx, $comp)));
    if (! defined $a_comps{$str}) {
      return 0;
    }
  }

  # index (labels and all) the "standard" candidate
  my %b_comps1;
  for my $comp (@{ $pb->Components($b_cx) }) {
    my $str = $pb->ComponentMol($b_cx, $comp) . ":" .
        join(",", sort @{ $pb->ComponentLabel($b_cx, $comp) },
        Pathway::NormalPTMString($pb->ComponentPTM($b_cx, $comp)));
    $b_comps1{$str} = 1;

  }

  # now see if there is a match for each component of the "a" mol
  for my $comp (@{ $pa->Components($a_cx) }) {
    my $str2 = $pa->ComponentMol($a_cx, $comp) . ":" .
        join(",", sort @{ $pa->ComponentLabel($a_cx, $comp) },
        Pathway::NormalPTMString($pa->ComponentPTM($a_cx, $comp)));

    my $compmol = $pa->ComponentMol($a_cx, $comp);
    if (defined $pa->{map}{$compmol}) {
      $compmol = $pa->{map}{$compmol};
    }
    my $str = $compmol . ":" .
        join(",", sort @{ $pa->ComponentLabel($a_cx, $comp) },
        Pathway::NormalPTMString($pa->ComponentPTM($a_cx, $comp)));
    if (! defined $b_comps1{$str}) {
      return 0;
    }
  }

  return 1;
}

######################################################################
sub NamesIntersect {
  my ($pa, $pb, $a_mol, $b_mol) = @_;

  my %a_names;

  my $hash = $pa->MolName($a_mol);
  for my $nametype (keys %{ $hash }) {
    for my $name (keys %{ $$hash{$nametype} }) {
      $a_names{$nametype . ":" . lc($name)} = 1;
    }
  }
  my $hash = $pb->MolName($b_mol);
  for my $nametype (keys %{ $hash }) {
    for my $name (keys %{ $$hash{$nametype} }) {
      if (defined $a_names{$nametype . ":" . lc($name)}) {
        return 1;
      }
    }
  }
  return 0;
}

######################################################################
sub ExIdsIntersect {
  my ($pa, $pb, $a_mol, $b_mol) = @_;

  my %a_ids;

  my $hash = $pa->MolExId($a_mol);
  for my $idtype (keys %{ $hash }) {
    for my $id (keys %{ $$hash{$idtype} }) {
      $a_ids{$idtype . ":" . lc($id)} = 1;
    }
  }
  my $hash = $pb->MolExId($b_mol);
  for my $idtype (keys %{ $hash }) {
    for my $id (keys %{ $$hash{$idtype} }) {
      if (defined $a_ids{$idtype . ":" . lc($id)}) {
        return 1;
      }
    }
  }
  return 0;
}

######################################################################
sub ExIdsClash {
  my ($pa, $pb, $a_mol, $b_mol) = @_;

  my $a_hash = $pa->MolExId($a_mol);
  my $b_hash = $pb->MolExId($b_mol);

    ## precedence of id types
  for my $type ("UP", "LL", "CH") {
    my $a = join(",", sort map { lc } keys %{ $$a_hash{$type} });
    my $b = join(",", sort map { lc } keys %{ $$b_hash{$type} });
    if (($a || $b) && ($a ne $b)) {
      return 1;
    }
  }
  return 0;
}

######################################################################
sub EqSimplex {
  my ($pa, $pb, $a_mol, $b_mol) = @_;

  if ($pa->MolType($a_mol) ne $pb->MolType($b_mol)) {
    return 0;
  }
  if (! PartWholeMatch($pa, $pb, $a_mol, $b_mol)) {
    return 0;
  }
  if (ExIdsClash($pa, $pb, $a_mol, $b_mol)) {
    return 0;
  }

  my $names_intersect = NamesIntersect($pa, $pb, $a_mol, $b_mol);
#  my $exids_intersect = ExIdsIntersect($pa, $pb, $a_mol, $b_mol);

  if ($names_intersect) {
    return 1;
  } else {
    ## this condition added to accommodate anonymous subunits of whole
    ## molecules produced by parsing Reactome BioPAX -- this may, in fact,
    ## not be needed any more since we are propating UP ids from whole
    ## to part.
    my $anames = $pa->MolName($a_mol);
    my $bnames = $pb->MolName($b_mol);
    if ((! defined $anames || keys %{ $anames } == 0) &&
        (! defined $bnames || keys %{ $bnames } == 0)) {
      return 1;
    }
  }

  return 0;
}

######################################################################
sub PartWholeMatch {
  my ($pa, $pb, $a_mol, $b_mol) = @_;

  ## if one of them is a part and the other is not, they
  ## do not match

  if (@{ $pa->MolWhole($a_mol) } != @{ $pb->MolWhole($b_mol) }) {
    return 0;
  }

  ## if they are both parts, then their bounds must match
  if (@{ $pa->MolWhole($a_mol) } > 0) {
    if ($pa->PartBounds($a_mol) ne $pb->PartBounds($b_mol)) {
      return 0
    }
  }

  return 1;
}

######################################################################
sub MergeComplexes {
  my ($self) = @_;

  my (@cxs, %replace);

  my $complex_type = $lv->StringToLabelValue("molecule-type", "complex");

  # true: simplex equivalence has already determined/replaced
  # true: complexes have been flattened
  my %complex_invert;
  for my $mol (sort numerically @{ $self->Mols() }) {
    if ($self->MolType($mol) eq $complex_type) {
      for my $comp (@{ $self->Components($mol) }) {
        $complex_invert{$self->ComponentMol($mol, $comp)}{$mol} = 1;
      }
    }
  }

  my %candidates;
  for my $link (keys %complex_invert) {
    my @mols = sort numerically keys %{ $complex_invert{$link} };
    for (my $i = 0; $i < @mols - 1; $i++) {
      for (my $j = $i+1; $j < @mols; $j++) {
        $candidates{$mols[$i]}{$mols[$j]} = 1;
      }
    }
  }
  for my $mol1 (sort numerically keys %candidates) {
    for my $mol2 (sort numerically keys %{ $candidates{$mol1} }) {
      ## true: $mol1 ne $mol2
      ## true: $mol1 < $mol2
      if (defined $replace{$mol2}) {
        next;
      }
      if (EqComplex($self, $self, $mol1, $mol2)) {
         $replace{$mol2} = $mol1;
      }
    }
  }

  for my $old (keys %replace) {
    $self->ReplaceMolecule($old, $replace{$old});
    $self->DeleteMolecule($old);
  }

}

######################################################################
sub FlattenAndMergeComplexes {
  my ($self) = @_;

  $self->FlattenComplexes();
  $self->MergeComplexes();
}

######################################################################
sub MergeSimplexes {
  my ($self) = @_;

  my (@mols, %replace);

  my $complex_type = $lv->StringToLabelValue("molecule-type", "complex");

  my %simplex_invert;
  for my $mol (@{ $self->Mols() }) {
    my $hit = 0;
    if ($self->MolType($mol) ne $complex_type) {
      my $hash = $self->MolName($mol);
      for my $type (keys %{ $hash }) {
        for my $name (keys %{ $$hash{$type} }) {
          $simplex_invert{lc($name)}{$mol} = 1;
          $hit++;
        }
      }
      my $hash = $self->MolExId($mol);
      for my $type (keys %{ $hash }) {
        for my $ext_id (keys %{ $$hash{$type} }) {
          $simplex_invert{lc($ext_id)}{$mol} = 1;
          $hit++;
        }
      }
      if ($hit == 0) {
        $simplex_invert{"NIL"}{$mol} = 1;
      }
    }
  }
  my %candidates;
  for my $link (keys %simplex_invert) {
    my @mols = sort numerically keys %{ $simplex_invert{$link} };
    for (my $i = 0; $i < @mols - 1; $i++) {
      for (my $j = $i+1; $j < @mols; $j++) {
        $candidates{$mols[$i]}{$mols[$j]} = 1;
      }
    }
  }
  for my $mol1 (sort numerically keys %candidates) {
    for my $mol2 (sort numerically keys %{ $candidates{$mol1} }) {
      ## true: $mol1 ne $mol2
      ## true: $mol1 < $mol2
      if (defined $replace{$mol2}) {
        next;
      }
      if (EqSimplex($self, $self, $mol1, $mol2)) {
        $replace{$mol2} = $mol1;
      }
    }
  }

  for my $old (keys %replace) {
    $self->ReplaceMolecule($old, $replace{$old});
    $self->DeleteMolecule($old);
  }

}

######################################################################
sub MergeMolecules {
  my ($self) = @_;

  $self->MergeSimplexes();
  $self->FlattenAndMergeComplexes();
}

######################################################################
sub ComplexDescendents {
  my ($self, $mol, $accum) = @_;

  if (defined $$accum{$mol}) {
    return;
  } else {
    $$accum{$mol} = 1;
  }
  my $COMPLEX_TYPE = $lv->StringToLabelValue("molecule-type", "complex");
  if ($self->MolType($mol) eq $COMPLEX_TYPE) {
    for my $comp (@{ $self->Components($mol) }) {
      my $compmol = $self->ComponentMol($mol, $comp);
      $self->ComplexDescendents($compmol, $accum);
    }
  }
}

######################################################################
## Equivalence relation
######################################################################

######################################################################
sub LabelSetIsSubType {
  my ($alabs, $blabs) = @_;

  my ($found, $al, $bl, $asort, $bsort);

  for $bl (@{ $blabs }) {
    $found = 0;
    $bsort = $lv->PWLabel::LabelSort($bl);
    $asort = 0;
    for $al (@{ $alabs }) {
      if ($lv->PWLabel::LabelSort($al) == $bsort) {
        $asort = $bsort;    ## That is, mol a DOES have a spec
                            ## for this label, so there better be a
                            ## real match on subtype
      }
#print STDERR "LabelSetIsSubType: $al, $bl\n";
      if ($lv->PWLabel::IsA($al, $bl)) {
        $found = 1;
        last;
      }
    }
#print STDERR "asort = $asort, bsort = $bsort\n";
    if (not $found) {
      return 0;
    }
    ## If b DID NOT have a spec for this label, then we
    ## consider cosnider a's spec to be a subtype
  }
  return 1;
}

######################################################################
###  other parts of equivalence AtomIsSubType, AtomsEq, MolsEq,
###  EdgeIsSubType, EdgeMolIsSubType removed when molinst was made
###  local to one pathway
######################################################################

######################################################################
sub NormalAtomString {
  my ($self, $atom) = @_;

  my ($atom_label, $mol_inst, $edge, $edge_label, $cond, @temp);
  $atom_label = NormalLabelString($self->AtomLabel($atom));
  $cond = NormalLabelString($self->AtomCondition($atom));
  for $edge (@{ $self->Edges($atom) }) {
    $mol_inst = $self->MolInstId($atom, $edge);
    $edge_label = NormalLabelString($self->EdgeLabel($atom, $edge));
    push @temp, $mol_inst . ":" . $edge_label;
  }
  if (defined $self->AbstractionExtId($atom)) {
    return ($atom_label . "::" . $self->AbstractionExtId($atom) .
        "::" . $cond . "::" . join(",", sort @temp));
  } else{
    return ($atom_label . "::" . $cond . "::" . join(",", sort @temp));
  }
}

######################################################################
sub SplitTranscription {
  my ($self) = @_;

  my $TRANSCRIPTION =
      $lv->StringToLabelValue("process-type", "transcription");
  for my $atom (@{ $self->Atoms }) {
    if ($self->AtomType($atom) ne $TRANSCRIPTION) {
      next;
    }
    my (@ins, @outs);
    for my $edge (@{ $self->Edges($atom) }) {
      if ($lv->IsA($self->EdgeType($atom, $edge),
          $lv->StringToLabelValue("edge-type", "outgoing-edge"))) {
        push @outs, $edge;
      } elsif ($lv->IsA($self->EdgeType($atom, $edge),
          $lv->StringToLabelValue("edge-type", "incoming-edge"))) {
	push @ins, $edge;
      }
    }
    if (@outs < 2) {
      next;
    }
    my $n_split = 0;
    for my $out_e (@outs) {
      $n_split++;
      my $new_atom = $atom . "_" . $n_split;
      $self->AddAtom($new_atom, $TRANSCRIPTION);
      $self->AddAtomSource($new_atom, $self->AtomSource($atom));
      for my $x (@{ $self->AtomCondition($atom) }) {
        $self->AddAtomCondition($new_atom, $x);
      }
      for my $x (@{ $self->AtomNegativeCondition($atom) }) {
        $self->AddAtomNegativeCondition($new_atom, $x);
      }
      for my $x (@{ $self->AtomLabel($atom) }) {
        $self->AddAtomLabel($new_atom, $x);
      }
      for my $e ($out_e, @ins) {
        $self->AddEdge($new_atom, $e, $self->EdgeType($atom, $e),
          $self->EdgeMol($atom, $e));
        for my $x (@{ $self->EdgeLabel($atom, $e) }) {
          $self->AddEdgeLabel($new_atom, $e, $x);
        }
        for my $x (@{ $self->MolLabel($atom, $e) }) {
          $self->AddMolLabel($new_atom, $e, $x);
        }
      }
      for my $x (@{ $self->AtomPathway($atom) }) {
        $self->AddAtomPathway($new_atom, $x);
      }
    }
    ##
    ## Clean up old stuff
    ##
    for my $x (@{ $self->Edges($atom) }) {
      delete $self->{moluse}{$self->EdgeMol($atom, $x)}{$atom}
    }
    delete $self->{atompathway}{$atom};
    delete $self->{atomlabel}{$atom};
    delete $self->{atomcondition}{$atom};
    delete $self->{atomtype}{$atom};
    delete $self->{atom2source}{$atom};
    delete $self->{edge}{$atom};
    delete $self->{edgetype}{$atom};
    delete $self->{edgelabel}{$atom};
    delete $self->{edgeptm}{$atom};
    delete $self->{mollabel}{$atom};
  }
}

######################################################################
sub ReplaceMolecule {
  my ($self, $old, $new) = @_;

  ## moltype should not be necessary
##  $self->{moltype}{$new} = $self->{moltype}{$old};

  ## copy all names from old to new
  for my $nametype (keys %{ $self->{molname}{$old} }) {
    for my $name (keys %{ $self->{molname}{$old}{$nametype} }) {
      $self->{molname}{$new}{$nametype}{$name} = 1;
    }
  }

  ## copy all ids from old to new
  for my $idtype (keys %{ $self->{molexid}{$old} }) {
    for my $id (keys %{ $self->{molexid}{$old}{$idtype} }) {
      $self->{molexid}{$new}{$idtype}{$id} = 1;
    }
  }

  ## add uses in atoms
  for my $atom (keys %{ $self->{moluse}{$old} }) {
    for my $edge (keys %{ $self->{edge}{$atom} }) {
      if ($self->{edge}{$atom}{$edge} eq $old) {
        $self->{edge}{$atom}{$edge} = $new;
      }
    }
    $self->{moluse}{$new}{$atom} = 1;
  }

  ## fix uses in complex

  for my $cx (keys %{ $self->{r_component}{$old} }) {
    for my $comp (keys %{ $self->{component}{$cx} }) {
      if ($self->{component}{$cx}{$comp} eq $old) {
         $self->{component}{$cx}{$comp} = $new;
      }
    }
    $self->{r_component}{$new}{$cx} = 1;
  }

  ## family

  for my $m (keys %{ $self->{family}{$old} }) {
    $self->{family}{$new}{$m} = 1;
    $self->{r_family}{$m}{$new} = 1;
  }

  for my $f (keys %{ $self->{r_family}{$old} }) {
    $self->{r_family}{$new}{$f} = 1;
    $self->{family}{$f}{$new} = 1;
  }

  for my $m (keys %{ $self->{memberlabel}{$old} }) {
    for my $lvid (keys %{ $self->{memberlabel}{$old}{$m} }) {
      $self->{memberlabel}{$new}{$m}{$lvid} = 1;
    }
  }
  for my $m (keys %{ $self->{memberptm}{$old} }) {
    for my $ptm (keys %{ $self->{memberptm}{$old}{$m} }) {
      $self->{memberptm}{$new}{$m}{$ptm} = 1;
    }
  }
  for my $f (keys %{ $self->{r_family}{$old} }) {
    for my $lvid (keys %{ $self->{memberlabel}{$f}{$old} }) {
      $self->{memberlabel}{$f}{$new}{$lvid} = 1;
    }
    for my $ptm (keys %{ $self->{memberptm}{$f}{$old} }) {
      $self->{memberptm}{$f}{$new}{$ptm} = 1;
    }
  }

  ## part
  ## should not have to reset part stuff
#  for my $p (keys %{ $self->{partof}{$old} }) {
#    $self->{partof}{$new}{$p} = 1;
#    $self->{r_partof}{$p}{$new} = 1;
#  }
#
#  for my $w (keys %{ $self->{r_partof}{$old} }) {
#    $self->{r_partof}{$new}{$w} = 1;
#    $self->{partof}{$w}{$new} = 1;
#  }
#  $self->{partbounds}{$new} = $self->{partbounds}{$old};

}

######################################################################
sub DeleteMolecule {
  my ($self, $mol) = @_;

  delete $self->{moltype}{$mol};
  delete $self->{molname}{$mol};
  delete $self->{molexid}{$mol};
  delete $self->{moluse}{$mol};

  ## complex

  delete $self->{component}{$mol};
  delete $self->{r_component}{$mol};
  delete $self->{componentlabel}{$mol};
  delete $self->{componentptm}{$mol};

  ## family

  for my $m (keys %{ $self->{family}{$mol} }) {
    delete $self->{r_family}{$m}{$mol};
    delete $self->{memberlabel}{$mol}{$m};
    delete $self->{memberptm}{$mol}{$m};
  }
  delete $self->{family}{$mol};

  for my $f (keys %{ $self->{r_family}{$mol} }) {
    delete $self->{family}{$f}{$mol};
    delete $self->{memberlabel}{$f}{$mol};
    delete $self->{memberptm}{$f}{$mol};
  }
  delete $self->{r_family}{$mol};

  ## part

  for my $p (keys %{ $self->{partof}{$mol} }) {
    delete $self->{r_partof}{$p}{$mol};
  }
  delete $self->{partof}{$mol};

  for my $w (keys %{ $self->{r_partof}{$mol} }) {
    delete $self->{partof}{$w}{$mol};
  }
  delete $self->{r_family}{$mol};

  if (defined $self->{partbounds}{$mol}) {
    delete $self->{partbounds}{$mol};
  }
}

######################################################################
sub DeleteMolLabelsOfKind {
  my ($self, $label_kind_string) = @_;

  for my $atom (keys %{ $self->{mollabel} }) {
    for my $edge (keys %{ $self->{mollabel}{$atom} }) {
      for my $value (keys %{ $self->{mollabel}{$atom}{$edge} }) {
        if ($lv->LabelString($lv->LabelSort($value)) eq
            $label_kind_string) {
          delete $self->{mollabel}{$atom}{$edge}{$value};
        }
      }
    }
  }
}

######################################################################
sub CloneMolecules {
  my ($self, $clone) = @_;

  for my $molid (keys %{ $self->{moltype} }) {
    $clone->{moltype}{$molid} = $self->{moltype}{$molid};
    for my $nametype (keys %{ $self->{molname}{$molid} }) {
      for my $name (keys %{ $self->{molname}{$molid}{$nametype} }) {
        $clone->{molname}{$molid}{$nametype}{$name} = 1;
      }
    }
    for my $idtype (keys %{ $self->{molexid}{$molid} }) {
      for my $id (keys %{ $self->{molexid}{$molid}{$idtype} }) {
        $clone->{molexid}{$molid}{$idtype}{$id} = 1;
      }
    }
  }
  for my $cxid (keys %{ $self->{component} }) {
    for my $cseq (keys %{ $self->{component}{$cxid} }) {
      $clone->{component}{$cxid}{$cseq} = $self->{component}{$cxid}{$cseq};
    }
    if (defined $self->{componentlabel}{$cxid}) {
      for my $cseq (keys %{ $self->{componentlabel}{$cxid} }) {
        for my $value (keys %{ $self->{componentlabel}{$cxid}{$cseq} }) {
          $clone->{componentlabel}{$cxid}{$cseq}{$value} = 1;
        }
        if (defined $self->{componentptm}{$cxid}{$cseq}) {
          for my $ptm_string
              (keys %{ $self->{componentptm}{$cxid}{$cseq} }) {
            $clone->{componentptm}{$cxid}{$cseq}{$ptm_string} = 1;
          }
        }
      }
    }
  }
  for my $family (keys %{ $self->{family} }) {
    for my $member (keys %{ $self->{family}{$family} }) {
      $clone->{family}{$family}{$member} = 1;
      $clone->{r_family}{$member}{$family} = 1;
    }
  }
  for my $whole (keys %{ $self->{partof} }) {
    for my $part (keys %{ $self->{partof}{$whole} }) {
      $clone->{partof}{$whole}{$part} = 1;
      $clone->{r_partof}{$part}{$whole} = 1;
    }
  }
}

######################################################################
sub CloneAtoms {
  my ($self, $clone) = @_;
  for my $atom (keys %{ $self->{atomtype} }) {
    $clone->{atomtype}{$atom} = $self->{atomtype}{$atom};
    for my $value (keys %{ $self->{atomlabel}{$atom} }) {
      $clone->{atomlabel}{$atom}{$value} = 1;
    }
    for my $value (keys %{ $self->{atomcondition}{$atom} }) {
      $clone->{atomcondition}{$atom}{$value} = 1;
    }
    if (defined $self->{atom2source}{$atom}) {
      $clone->{atom2source}{$atom} = $self->{atom2source}{$atom};
    }
    for my $edge (keys %{ $self->{edge}{$atom} }) {
      $clone->{edge}{$atom}{$edge} = $self->{edge}{$atom}{$edge};
      $clone->{moluse}{$self->{edge}{$atom}{$edge}}{$atom} = 1;
      $clone->{edgetype}{$atom}{$edge} = $self->{edgetype}{$atom}{$edge};
      if (defined $self->{edgelabel}{$atom}) {
        for my $value (keys %{ $self->{edgelabel}{$atom}{$edge} }) {
          $clone->{edgelabel}{$atom}{$edge}{$value} = 1;
        }
      }
      if (defined $self->{mollabel}{$atom}) {
        for my $value (keys %{ $self->{mollabel}{$atom}{$edge} }) {
          $clone->{mollabel}{$atom}{$edge}{$value} = 1;
        }
      }
      if (defined $self->{edgeptm}{$atom}{$edge}) {
        for my $ptm_string (keys %{ $self->{edgeptm}{$atom}{$edge} }) {
          $clone->{edgeptm}{$atom}{$edge}{$ptm_string} = 1;
        }
      }
    }
    if (defined $self->{abstraction}{pathway_id}{$atom}) {
      $clone->{abstraction}{pathway_id}{$atom} =
          $self->{abstraction}{pathway_id}{$atom};
    }
    if (defined $self->{abstraction}{ext_pathway_id}{$atom}) {
      $clone->{abstraction}{ext_pathway_id}{$atom} = 
          $self->{abstraction}{ext_pathway_id}{$atom};
    }
    if (defined $self->{abstraction}{pathway_name}{$atom}) {
      $clone->{abstraction}{pathway_name}{$atom} =
          $self->{abstraction}{pathway_name}{$atom};
    }
  }
}

######################################################################
sub ClonePathwaysInModel {
  my ($self, $clone) = @_;

  for my $pid (keys %{ $self->{pathway} }) {
    $clone->{pathway}{$pid}{pname}  = $self->{pathway}{$pid}{pname};
    $clone->{pathway}{$pid}{pexid}  = $self->{pathway}{$pid}{pexid};
    $clone->{pathway}{$pid}{porg}   = $self->{pathway}{$pid}{porg};
    $clone->{pathway}{$pid}{psrcid} = $self->{pathway}{$pid}{psrcid};
    $clone->{pathway}{$pid}{pname}  = $self->{pathway}{$pid}{pname};
    $clone->{pathway}{$pid}{issubnet}  = $self->{pathway}{$pid}{issubnet};
    $clone->{pathway2source}{$pid}  = $self->{pathway2source}{$pid};
    for my $atom (keys %{ $self->{atompathway} }) {
      for my $pid (keys %{ $self->{atompathway}{$atom} }) {
        $clone->{atompathway}{$atom}{$pid} = 1;
      }
    }
  }
}

######################################################################
sub CloneSource {
  my ($self, $clone) = @_;
  for my $src (keys %{ $self->{source2name} }) {
    $clone->{source2name}{$src} = $self->{source2name}{$src};
  }
}

######################################################################
sub ClonePathway {
  my ($self, $clone) = @_;
  ## copy self into clone
  $self->CloneSource($clone);
  $self->CloneMolecules($clone);
  $self->CloneAtoms($clone);
  $self->ClonePathwaysInModel($clone);
}

######################################################################
sub CopyMol {
  my ($self, $copy, $mol) = @_;

  if (defined $copy->{moltype}{$mol}) {
    return;
  };

  my $COMPLEX_TYPE = $lv->StringToLabelValue("molecule-type", "complex");

  $copy->AddMol($mol, $self->MolType($mol));
  for my $type (keys %{ $self->{molname}{$mol} }) {
    for my $name (keys %{ $self->{molname}{$mol}{$type} }) {
      $copy->AddMolName($mol, $type, $name);
    }
  }
  for my $type (keys %{ $self->{molexid}{$mol} }) {
    for my $id (keys %{ $self->{molexit}{$mol}{$type} }) {
      $copy->AddMolExId($mol, $type, $id);
    }
  }
  if ($self->MolType($mol) eq $COMPLEX_TYPE) {
    for my $comp (@{ $self->Components($mol) }) {
      my $compmol = $self->ComponentMol($mol, $comp);
      $copy->AddComponent($mol, $comp, $compmol);
      for my $lvid (@{ $self->ComponentLabel($mol, $comp) }) {
        $copy->AddComponentLabel($mol, $comp, $lvid);
      }
      for my $ptm (@{ $self->ComponentPTM($mol, $comp) }) {
        my ($protein_id, $position, $amino_acid, $modification_label_value,
            $modification_label_name) = split(",", $ptm);
        $copy->AddComponentPTM($mol, $comp, $protein_id, $position, $amino_acid,
            $modification_label_value, $modification_label_name);
      }
      $self->CopyMol($copy, $compmol);
    }
  }
  for my $member (keys %{ $self->{memberlabel}{$mol} }) {
    for my $lvid (keys %{ $self->{memberlabel}{$mol}{$member} }) {
      $copy->{memberlabel}{$mol}{$member} = $lvid;
    }
  }
  for my $member (keys %{ $self->{memberptm}{$mol} }) {
    for my $ptm_string (keys %{ $self->{memberptm}{$mol}{$member} }) {
      $self->{memberptm}{$mol}{$member}{$ptm_string} = 1;
    }
  }
}

######################################################################
sub CopyAtom {
  my ($self, $copy, $atom) = @_;

  $copy->AddAtom($atom, $self->AtomType($atom));
  $copy->AddAtomSource($atom, $self->AtomSource($atom));
  for my $lvid (@{ $self->AtomLabel($atom) }) {
    $copy->AddAtomLabel($atom, $lvid);
  }
  for my $lvid (@{ $self->AtomCondition($atom) }) {
    $copy->AddAtomCondition($atom, $lvid);
  }
  for my $lvid (@{ $self->AtomNegativeCondition($atom) }) {
    $copy->AddAtomNegativeCondition($atom, $lvid);
  }
  for my $ref (@{ $self->AtomReferences($atom) }) {
    $copy->AddAtomReferences($atom, $ref);
  }
  for my $note (@{ $self->AtomNotes($atom) }) {
    $copy->AddAtomNotes($atom, $note);
  }
  for my $edge (@{ $self->Edges($atom) }) {
    my $mol = $self->EdgeMol($atom, $edge);
    my $edge_type = $self->EdgeType($atom, $edge);
    $copy->AddEdge($atom, $edge, $self->EdgeType($atom, $edge), $mol);
    for my $lvid (@{ $self->EdgeLabel($atom, $edge) }) {
      $copy->AddEdgeLabel($atom, $edge, $lvid);
    }
    for my $lvid (@{ $self->MolLabel($atom, $edge) }) {
      $copy->AddMolLabel($atom, $edge, $lvid);
    }
    for my $ptm (@{ $self->EdgePTM($atom, $edge) }) {
      my ($protein_id, $position, $amino_acid, $modification_label_value,
          $modification_label_name) = split(",", $ptm);
      $copy->AddEdgePTM($atom, $edge, $protein_id, $position, $amino_acid,
          $modification_label_value, $modification_label_name);
    }
    $self->CopyMol($copy, $mol);
  }
  if (defined $self->AbstractionExtId($atom) ||
      defined $self->AbstractionId($atom) ||
      defined $self->AbstractionName($atom) ) {
    $copy->AddAbstraction($atom, $self->AbstractionId($atom),
        $self->AbstractionName($atom), $self->AbstractionExtId($atom));
    if (defined $self->AbstractionId($atom)) {
      my $p = $self->AbstractionId($atom);
      $copy->SetIsSubnet($p, $self->PathwayIsSubnet($p));
    }
  }
}

######################################################################
sub CopyPathway {
  my ($self, $copy, $pid) = @_;

  ## just copy all the source names; won't hurt
  $self->CloneSource($copy);

  ## copy pathway info
  $copy->AddPathwayId($pid);
  $copy->SetPathwayName($pid, $self->PathwayName($pid));
  $copy->SetPathwayExId($pid, $self->PathwayExId($pid));
  $copy->SetPathwayOrg($pid, $self->PathwayOrg($pid));
  $copy->SetPathwaySrcId($pid, $self->PathwaySrcId($pid));
  $copy->AddPathwaySource($pid, $self->PathwaySource($pid));
  $copy->SetIsSubnet($pid, $self->PathwayIsSubnet($pid));
  for my $ref (@{ $self->PathwayReferences($pid) }) {
    $copy->AddPathwayReferences($pid, $ref);
  }
  for my $reviewer (@{ $self->Reviewers($pid) }) {
    $copy->AddReviewer($pid, $reviewer);
  }
  for my $curator (@{ $self->Curators($pid) }) {
    $copy->AddCurator($pid, $curator);
  }
  for my $atom (@{ $self->Atoms() }) {
    for my $p (@{ $self->AtomPathway($atom) }) {
      if ($p eq $pid) {
        $copy->AddAtomPathway($atom, $pid);
        $self->CopyAtom($copy, $atom);
      }
    }
  }
#  for my $s (@{ $self->SubnetsInPathway($pid) }) {
#    $copy->AddPathwaySubnet($pid, $s);
#    for my $a (@{ $self->AtomsInSubnet($s) }) {
#      $copy->AddSubnetAtom($s, $a);
#    }
#  }
}

######################################################################
1;
