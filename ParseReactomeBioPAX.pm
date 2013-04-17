

# Copyright SRA International
#
# Distributed under the OSI-approved BSD 3-Clause License.
# See http://ncip.github.com/pathway-interaction-database/LICENSE.txt for details.


package ParseReactomeBioPAX;

require Exporter;
@ISA = qw(Exporter);
#@EXPORT = qw(
#  CreateIsA
#  CreateLabelSort
#);

use strict;
use XML::Parser;
use PWLabel;
use Pathway;
use IdMap;

use constant REACTOME_SOURCE_ID   => 7;
use constant REACTOME_SOURCE_NAME => "Reactome";

#### global unique numeric ids
#### to be used for pathways, interactions, molecules, ptms
my ($next_id, $next_pathway_id, $next_interaction_id, $next_molecule_id,
    $next_ptm_id);
my (%final_id);
my (%rdf_id2id, %id2rdf_id);

#### global lv, pw
my ($lv, $pw, $idm);

my (@elem_stack, @attr_stack, %char_value, %id_def);

my ($unification_db, $unification_id);
my ($xref_db, $xref_id);
my (@components, @pathway_components, @step_interactions, @next_steps);

my (%id2type, %id2xref, %reaction2participant, %complex2component,
    %network2component,
    %catalysis2controller, %catalysis2controlled,
    %catalysis2direction, %catalysis2control_type,
    %r_catalysis,
    %id2name, %id2short_name, %id2synonym,
    %participant2stoichiometry, %participant2cellular_location,
    %participant2physical_entity,
    %participant2db_id,     ## populated from COMMENT's: internal Reactome DB_ID
    %edge2participant,
    %sequence_participant2feature, %sequence_feature2type,
    %feature2location, %feature2residue,
    %unification, %relationship, %publication,
    %site2position, %site2position_status, %interval2begin, %interval2end);
my (%uniprot2entrez_gene);
my ($pathway2reaction);


### Shiva's mapping of our terms to CHEBI
# Acetylation: MOD:00394; acetyl group (CHEBI:22190)
# Biotinylation: MOD:00126; biotin (CHEBI:15956)??
# Dephosphoryaltion
# Farnesyaltion: MOD:00437; farnesyl group (CHEBI:24017)
# Glysoaminoglycan: MOD:00224????; glycosaminoglycans (CHEBI:18085)
# Geranylgeranylation: MOD:00441; geranylgeranyl group (CHEBI:24231)
# Glycosylation: MOD:00693; glycosyl glycosides (CHEBI:24407)???
# Methylation: MOD:00599; methyl group (CHEBI:32875)
# Myristoylation: MOD:00438; myristoyl group (CHEBI:25456)
# Oxidation: MOD:00412; 
# Palmitoylation: MOD:00440; palmitoyl group (CHEBI:25839)
# Phosphorylation: MOD:00696; phosphate group (CHEBI:32958)??
# poly (ADP-ribosyl)ation: CHEBI:16960
# Sumoylation
# Ubiquitination

### Reactome mapping to ChEBI
# _1_4_alpha_D_glucosyl_n_group : ChEBI_15444
# L_selenocysteinyl_group : ChEBI_30003
# palmitoyl_group : ChEBI_25839
# phosphate_group : ChEBI_35780
# carboxyl_group : ChEBI_23025
# glycogen_group : ChEBI_28087
# myristoyl_group : ChEBI_25456
# acetyl_group : ChEBI_22190
# hydroxyl_group : ChEBI_24706
# half_cystyl_group : ChEBI_30770
# methyl_group : ChEBI_25251

my %translate_location;

my %reactome_feature_types = (
  "_1_4_alpha_D_glucosyl_n_group" => "_1_4_alpha_D_glucosyl_n_group",  ## not PID standard
  "_1_6__alpha_glucosyl__1_4_alpha_D_glucosyl_n_group" => "_1_6__alpha_glucosyl__1_4_alpha_D_glucosyl_n_group",  ## not PID standard
  "acetyl_group"                  => "acetylation",
  "acyl_GPI_group"                => "acyl_GPI_group",  ## not PID standard
  "carboxyl_group"                => "carboxyl_group",  ## not PID standard
  "glycogen_group"                => "glycogen_group",  ## not PID standard
  "GPI_anchor_group"              => "GPI_anchor_group",  ## not PID standard
  "half_cystyl_group"             => "half_cystyl_group",  ## not PID standard
  "hydroxyl_group"                => "hydroxyl_group",  ## not PID standard
  "L_selenocysteinyl_group"       => "L_selenocysteinyl_group",  ## not PID standard
  "limit_dextrin_group"           => "limit_dextrin_group",  ## not PID standard
  "methyl_group"                  => "methyl_group",  ## not PID standard
  "myristoyl_group"               => "myristoylation",
  "oligo__1_4_alpha_D_glucosyl__group"     => "oligo__1_4_alpha_D_glucosyl__group",  ## not PID standard
  "Orthophosphate__ChEBI_18367_1" => "Orthophosphate__ChEBI_18367_1",  ## not PID standard
  "palmitoyl_group"               => "palmitoylation",
  "Pantetheine_4__phosphate__ChEBI_16858_" => "Pantetheine_4__phosphate__ChEBI_16858_",  ## not PID standard
# Shiva: may be ok to map to phosphorylation if we allow unspecified
# amino acids
#  "phosphate_group"               => "phosphate_group"
  "phosphate_group"               => "phosphorylation"
);

my %residue_type = (
  "L-Threonine [ChEBI:6308]"  => "T",
  "L-Serine [ChEBI:6301]"     => "S",
  "L-Tyrosine [ChEBI:17895]"  => "Y",
  "Glycine [ChEBI:15428]"     => "G",
  "L-Aspartate [ChEBI:17053]" => "D",
  "L-Cysteine [ChEBI:17561]"  => "C",
  "L-Glutamate [ChEBI:16015]" => "E",
  "asparagino group"          => "X"
   ##Per Shiva: [aspartate + glutamine + ATP <=> asparagine + glutamate + AMP + pyrophosphate] 
);

my $p;

######################################################################
sub MaxValue {
  my $max = 0;
  for my $v (@_) {
    if ($v > $max) {
      $max = $v;
    }
  }
  return $max;
}

######################################################################
sub InferConnections {
  my ($self) = @_;

  my (%pid2abstraction);
  my $pathway_type =
    $lv->StringToLabelValue("process-type", "pathway");
  for my $atom (@{ $pw->Atoms() }) {
    if ($pw->AtomType($atom) eq $pathway_type) {
      if (defined $pw->AbstractionId($atom)) {
        $pid2abstraction{$pw->AbstractionId($atom)}{$atom} = 1;
#print STDERR "added pathway " . $pw->AbstractionId($atom) . ", atom $atom\n";
      } else {
        print STDERR "no abstraction defined for atom $atom of pathway type\n";
      }
    }
  }

  my $agent_type =
      $lv->StringToLabelValue("edge-type", "agent");
  my $inhibitor_type =
      $lv->StringToLabelValue("edge-type", "inhibitor");
  my $input_type =
      $lv->StringToLabelValue("edge-type", "input");
  my $output_type =
      $lv->StringToLabelValue("edge-type", "output");
  my $pathway_type =
      $lv->StringToLabelValue("process-type", "pathway");

  my (%mi2e, %m_p_e, %parent, %children, %siblings, %already_added);

  for my $atom (@{ $pw->Atoms() }) {
    my @pids =  @{ $pw->AtomPathway($atom) };
    if ($pw->AtomType($atom) eq "$pathway_type") {
      for my $pid (@pids) {
        $children{$pid}{$pw->AbstractionId($atom)} = 1;
        $parent{$pw->AbstractionId($atom)}{$pid} = 1;
      }
    }
    for my $edge (@{ $pw->Edges($atom) }) {
      my $edgetype = $pw->EdgeType($atom, $edge);
      my $molinstid = $pw->MolInstId($atom, $edge);
      ## need to save this so that we can recover the full PTM info
      ## for the instance
      if (! defined $mi2e{$molinstid}) {
        $mi2e{$molinstid} = "$atom,$edge";
      }
      for my $pid (@pids) {
        $m_p_e{$molinstid}{$pid}{$edgetype} = 1;
      }
    }
  }
  for my $pid (keys %children) {
    my @pids = keys %{ $children{$pid} };
    for (my $i = 0; $i < @pids - 1; $i++) {
      for (my $j = $i+1; $j < @pids; $j++) {
        $siblings{$pids[$i]}{$pids[$j]} = 1;
        $siblings{$pids[$j]}{$pids[$i]} = 1;
      }
    }
  }

  my $more = 1;
  while ($more) {
    $more = 0;
    for my $molinstid (keys %m_p_e) {
      my ($a,$e) = split(",", $mi2e{$molinstid});
      my @pids = keys %{ $m_p_e{$molinstid} };
#print STDERR "pids = ". scalar(@pids) . "\n";
      my %hits;
      for (my $i = 0; $i < @pids - 1; $i++) {
        for (my $j = $i+1; $j < @pids; $j++) {
          if (defined $siblings{$pids[$i]}{$pids[$j]}) {
            $hits{$pids[$i]}++;
            $hits{$pids[$j]}++;
          } elsif (defined $parent{$pids[$i]}{$pids[$j]}) {
            $hits{$pids[$i]}++;
          } elsif (defined $parent{$pids[$j]}{$pids[$i]}) {
            $hits{$pids[$j]}++;
          }
        }
      }
      @pids = keys %hits;
#print STDERR "hits = ". scalar(@pids) . "\n";
      if (@pids > 1) {
        for my $pid (@pids) {
#print STDERR "looking at pid $pid\n";
          if (! defined $pid2abstraction{$pid}) {
            print STDERR "no abstraction defined for pathway $pid\n";
            next;
          }
          for my $abstraction (keys %{ $pid2abstraction{$pid} }) {
#print STDERR "looking at abstraction $abstraction\n";
            my $next_edge = MaxValue(@{ $pw->Edges($abstraction) });
            for my $edgetype (keys %{ $m_p_e{$molinstid}{$pid} }) {
              if (! defined $already_added{"$molinstid,$abstraction,$edgetype"}) {
                $already_added{"$molinstid,$abstraction,$edgetype"} = 1;
                for my $parent (keys %{ $parent{$pid} }) {
                  $m_p_e{$molinstid}{$parent}{$edgetype} = 1
                }
                $next_edge++;
                $pw->AddEdge($abstraction,
                    $next_edge, $edgetype, $pw->EdgeMol($a, $e));
#print STDERR "adding edge $next_edge, edgetype $edgetype to $abstraction\n";
                for my $lvid (@{ $pw->MolLabel($a, $e) }) {
                  $pw->AddMolLabel($abstraction, $next_edge, $lvid);
                }
                for my $ptm (@{ $pw->EdgePTM($a, $e) }) {
                  my ($protein_id, $position, $amino_acid,
                      $modification_label_value,
                      $modification_label_name) = @{ $ptm };
                  $pw->AddEdgePTM($abstraction,
                      $next_edge, $protein_id, $position, $amino_acid,
                      $modification_label_value, $modification_label_name);
                }
                $more++;
              }
            }
          }
        }
      }
    }
  }
}

######################################################################
sub DumpXrefs {
  for my $x (keys %id2xref) {
    for my $y (keys %{ $id2xref{$x} }) {
      print "$x : $y\n";
    }
  }
}

######################################################################
sub Close {
  my ($relation) = @_;

  my %copy;
  for my $x (keys %{ $relation }) {
    for my $y (keys %{ $$relation{$x} }) {
      $copy{$x}{$y} = 1;
    }
  }

  my $more = 1;
  while ($more) {
    $more = 0;
    for my $x (keys %copy) {
      for my $y (keys %{ $copy{$x} }) {
        if (defined $copy{$y}) {
          for my $z (keys %{ $copy{$y} }) {
            if (!defined $copy{$x}{$z}) {
              $copy{$x}{$z} = 1;
              $more++;
            }
          }
        }
      }
    }
  }
  return \%copy;
}

######################################################################
sub CollapseCatalysis {

  for my $c (keys %catalysis2controlled) {
    my $controlled   = $catalysis2controlled{$c};
    if ($id2type{$catalysis2controlled{$c}} eq "bp:catalysis") {
      $r_catalysis{$catalysis2controlled{$c}} = $c;
#      print STDERR "found catalysis $c that controls a " .
#          "$id2type{$catalysis2controlled{$c}}\n";
    }
  }

  for my $c (keys %catalysis2controlled) {
    my $controller   = $catalysis2controller{$c};
    my $controlled   = $catalysis2controlled{$c};
    my $control_type = $catalysis2control_type{$c};
    my $direction    = $catalysis2direction{$c};    ##!!!!!!! huh
    if ($id2type{$catalysis2controlled{$c}} eq "bp:biochemicalReaction") {
      $reaction2participant{$controlled}{$control_type}{$controller} = 1;
      for my $xref (keys %{ $id2xref{$c} }) {
        if (defined $publication{$xref}) {
          $id2xref{$controlled}{$xref} = 1;
# print STDERR "transferred $xref from $c to $controlled\n";
        }
      }
    }
    if (defined $r_catalysis{$c}) {
      my $c1 = $c;
      my $level = 2;
      while (defined $r_catalysis{$c1}) {
        $c1 = $r_catalysis{$c1};
        my $controller1   = $catalysis2controller{$c1};
        my $controlled1   = $catalysis2controlled{$c1};
        my $control_type1 = $catalysis2control_type{$c1} . "_$level";
        my $direction1    = $catalysis2direction{$c1};    ##!!!!!!! huh
        $reaction2participant{$controlled}{$control_type1}{$controller1} = 1;
        $level++;
        for my $xref1 (keys %{ $id2xref{$c1} }) {
          if (defined $publication{$xref1}) {
            $id2xref{$controlled}{$xref1} = 1;
# print STDERR "transferred $xref1 from $c1 to $controlled\n";
          }
        }
      }
    }
  }
}

######################################################################
sub GetRootsAndLeaves {
  my ($relation) = @_;

  ## close the whole relation

  my $copy = Close($relation);

  ## find roots in the relation

  my %roots;
  for my $x (keys %{ $relation }) {
    $roots{$x} = 1;
  }
  for my $x (keys %{ $relation }) {
    for my $y (keys %{ $$relation{$x} }) {
      if (defined $roots{$y}) {
        $roots{$y} = 0;
      }
    }
  }

  ## remove from the closed relation anything that involves something
  ## that is both not a root and not a leaf

  for my $x (keys %roots) {
    if ($roots{$x} == 0) {
      delete $$copy{$x};
    }
  }
  for my $x (keys %{ $copy }) {
    for my $y (keys %{ $$copy{$x}}) {
      if (defined $roots{$y}) {
        delete $$copy{$x}{$y};
      }
    }
  }

  return $copy;

}
				   
######################################################################
sub dump {
  my ($self, $fh) = @_;

  for my $id (keys %id2type) {
    print $fh "#type\t$id\t$id2type{$id}\n";
  }
  for my $id (keys %id2xref) {
    for my $xref (keys %{ $id2xref{$id} }) {
      if (defined $unification{$xref}) {
        if (defined $unification{$xref}{UniProt}) {
          print $fh "#xref\t$id\tUniProt\t$unification{$xref}{UniProt}\n";
        }
      }
    }
  }
  for my $reaction (keys %reaction2participant) {
    for my $role (keys %{ $reaction2participant{$reaction} }) {
      for my $participant (keys %{ $reaction2participant{$reaction}{$role}}) {
        print $fh "#participant\t$reaction\t$role\t$participant\n";
      }
    }
  }
  for my $complex (keys %complex2component) {
    for my $component (keys %{ $complex2component{$complex} }) {
      print $fh "#complex_component\t$complex\t$component\n";
    }
  }
  for my $network (keys %network2component) {
    for my $component (keys %{ $network2component{$network} }) {
# we'll skip this in favor of the pathway2interaction
#      print "#network_component\t$network\t$component\n";
    }
  }
  for my $pathway (keys %{ $pathway2reaction }) {
    for my $interaction (keys %{ $$pathway2reaction{$pathway} }) {
      print $fh "#pathway_interaction\t$pathway\t$interaction\n";
    }
  }

}

######################################################################
sub DeleteSpontaneousInteractions {

  my $transcription_type =
      $lv->StringToLabelValue("process-type", "transcription");
  my $reaction_type =
      $lv->StringToLabelValue("process-type", "reaction");
  my $modification_type =
      $lv->StringToLabelValue("process-type", "modification");
  my $translocation_type =
      $lv->StringToLabelValue("process-type", "translocation");
  my $pathway_type =
      $lv->StringToLabelValue("process-type", "pathway");
  my %of_interest = (
    $transcription_type => 1,
    $reaction_type => 1,
    $modification_type => 1,
    $translocation_type => 1
  );
  my $agent_type =
      $lv->StringToLabelValue("edge-type", "agent");
  my $inhibitor_type =
      $lv->StringToLabelValue("edge-type", "inhibitor");
  my $input_type =
      $lv->StringToLabelValue("edge-type", "input");
  my %incoming = (
    $agent_type => 1,
    $inhibitor_type => 1,
    $input_type => 1
  );

  my (%pathway2atom, %atom2pathway);
  my (%pathway2abstraction);

  for my $atom (@{ $pw->Atoms() }) {
    my $atomtype = $pw->AtomType($atom);
    for my $pid (@{ $pw->AtomPathway($atom) }) {
      $pathway2atom{$pid}{$atom} = 1;
      $atom2pathway{$atom}{$pid} = 1;
    }
    if ($atomtype eq $pathway_type) {
      if (defined $pw->AbstractionId($atom)) {
        $pathway2abstraction{$pw->AbstractionId($atom)}{$atom} = 1;
      }
    }
  }

  my %kill;
  for my $atom (@{ $pw->Atoms() }) {
    my $atomtype = $pw->AtomType($atom);
    if (defined $of_interest{$atomtype}) {
      my $incoming = 0;
      for my $edge (@{ $pw->Edges($atom) }) {
        if (defined $incoming{$pw->EdgeType($atom, $edge)}) {
          $incoming++;
          last;
        }
      }
      if (@{ $pw->AtomCondition($atom) } > 0) {
        $incoming++;
      }
      if (! $incoming) {
print STDERR "#deleting_spontaneous_interaction\t$atom\n";
        $pw->PruneAtoms([ $atom ], {});
        $kill{$atom} = 1;
      }
    }
  }

  ## now we have to propagate kills
  my %fresh_kill_atom;
  my %kill_pathway;
  my %p2a_count;

  for my $pid (keys %pathway2atom) {
    $p2a_count{$pid} = scalar(keys %{ $pathway2atom{$pid} });
  }

  for my $atom (keys %kill) {
    $fresh_kill_atom{$atom} = 1;
  }
  my $more = 1;

  while ($more) {

    $more = 0;
    undef %kill;
    for my $atom (keys %fresh_kill_atom) {
      $kill{$atom} = 1;
    }
    undef %fresh_kill_atom;

    for my $atom (keys %kill) {
      for my $pid (keys %{ $atom2pathway{$atom} }) {
        $p2a_count{$pid}--;
        if ($p2a_count{$pid} == 0) {    ## that is, it just became == 0
          $kill_pathway{$pid} = 1;
          for my $abs (keys %{ $pathway2abstraction{$pid} }) {
            $fresh_kill_atom{$abs} = 1;
            $more++;
          }
        }
      }
    }
    for my $atom (keys %fresh_kill_atom) {
      $pw->PruneAtoms([ $atom ], {});
print STDERR "#deleting_dangling_abstraction\t$atom\n";
    }
  }

  for my $pid (keys %kill_pathway) {
print STDERR "#deleting_empty_pathway\t$pid\n";
    delete $pw->{pathway2source}{$pid};
    delete $pw->{pathway}{$pid};
    delete $pw->{pathwayreference}{$pid};
    delete $pw->{pathwayreviewer}{$pid};
    delete $pw->{pathwaycurator}{$pid};
  }
}

######################################################################
sub DoGraphStuff {
  my ($self) = @_;

  $pathway2reaction = GetRootsAndLeaves(\%network2component);
  CollapseCatalysis();
  if (defined $self->{uniprot_f}) {
    DoUniProt($self->{uniprot_f});
  }
  FinishComplexes();
  FinishMols();
  FinishInteractions();
  FinishPathways();
  PropagateWholeProteinIDs();
  DeleteSpontaneousInteractions();
#  for my $i (keys %id2rdf_id) {
#    print STDERR "#PID-to-Reactome\t$i\t$id2rdf_id{$i}\n";
#  }
}

######################################################################
sub DoUniProt {
  my ($f) = @_;

  my ($uniprot, $entrez_gene);

  open(UNI, $f) or die "cannot open $f";
  while (<UNI>) {
    chop;
    ($uniprot, $entrez_gene) = split /\t/;
    $uniprot2entrez_gene{$uniprot}{$entrez_gene} = 1;
  }
  close UNI;
}

######################################################################
sub FlattenPathwaySteps {

  my $more = 1;
  while ($more) {
    $more = 0;
    my %add;
    my %remove;
    for my $p (keys %network2component) {
      for my $q (keys %{ $network2component{$p} }) {
        if ($id2type{$q} eq "bp:pathwayStep") {
          $more++;
          for my $r (keys %{ $network2component{$q} }) {
            $add{$p}{$r} = 1;
            $remove{$q}{$r} = 1;
          }
          $remove{$p}{$q} = 1;
        }
      }
    }
    for my $q (keys %remove) {
      for my $r (keys %{ $remove{$q} }) {
        delete $network2component{$q}{$r};
      }
    }
    for my $p (keys %add) {
      for my $r (keys %{ $add{$p} }) {
        $network2component{$p}{$r} = 1;
      }
    }
  }
  for my $p (keys %network2component) {
    if ($id2type{$p} eq "bp:pathwayStep") {
      delete $network2component{$p};
    }
  }  


  ## delete remaining pathways that do not have components

  my $more = 1;
  while ($more) {
    $more = 0;
    my %kill;
    for my $x (keys %network2component) {
      for my $y (keys %{ $network2component{$x}}) {
        if ($id2type{$x} eq "bp:pathway") {
          if (! defined $network2component{$y} &&
              $id2type{$y} eq "bp:pathway") {
            $more++;
            $kill{$x}{$y} = 1;
          }
        }
      }
    }
    if ($more) {
      for my $x (keys %kill) {
        for my $y (keys %{ $kill{$x} }) {
          delete $network2component{$x}{$y};
#print STDERR "deleting $x, $y\n";
        }
      }
      for my $x (keys %network2component) {
        if (keys %{ $network2component{$x} } == 0) {
#print STDERR "deleting $x\n";
          delete $network2component{$x};
        }
      }
    }
  }

}

######################################################################
sub FinishPathways {

  FlattenPathwaySteps();

## May 30, 2007: we should add another type of flattening. If p is pathway
## that is not a top-level pathway and if p contains only ONE component c
## (where c is either interaction or pathway), then raise c to be a component
## of parent(p) and delete p.

  for my $pathway (keys %id2type) {
    if ($id2type{$pathway} ne "bp:pathway") {
      next;
    }
    my $pathway_id = FinalId("pathway", LookUp($pathway));
    $pw->AddPathwayId($pathway_id);
    $pw->AddPathwaySource($pathway_id, REACTOME_SOURCE_ID);
    $pw->SetPathwayExId($pathway_id, $pathway);
    $pw->SetPathwayOrg($pathway_id, "Hs");
    $pw->SetIsSubnet($pathway_id, 0);       ## initially, set everything to not-subnet
    for my $name (keys %{ $id2name{$pathway} }) {
      $pw->SetPathwayName($pathway_id, $name);
      last;    ##### take only one name
    }
    for my $xref (keys %{ $id2xref{$pathway} }) {
      if (defined $publication{$xref}{"Pubmed"}) {
        $pw->AddPathwayReferences($pathway_id, $publication{$xref}{"Pubmed"});
      }
    }

  }

  for my $pathway (keys %network2component) {
    for my $pway_component (keys %{ $network2component{$pathway} }) {
      if ($id2type{$pway_component} eq "bp:pathway") {
        if (defined $network2component{$pway_component}) {
          $pw->SetIsSubnet(FinalId("pathway", LookUp($pway_component)), 1);
        } else {
print STDERR "#empty_subnet\t$pway_component\n";
        }
      }
    }
  }

  for my $pathway (keys %network2component) {
    my $pathway_id = FinalId("pathway", LookUp($pathway));
    for my $pway_component (keys %{ $network2component{$pathway} }) {
      my $pway_component_id = FinalId("interaction", LookUp($pway_component));
      if ($id2type{$pway_component} eq "bp:pathway") {
        $pw->AddAtom($pway_component_id,
            $lv->StringToLabelValue("process-type", "pathway"));
        $pw->AddAtomSource($pway_component_id, REACTOME_SOURCE_ID);
        my $abstracted_p = FinalId("pathway", LookUp($pway_component));
        $pw->AddAbstraction($pway_component_id, $abstracted_p,
            $pw->PathwayName($abstracted_p),
            $pw->PathwayExId($abstracted_p));
      }
      $pw->AddAtomPathway($pway_component_id, $pathway_id);
    }
  }

}

######################################################################
sub FinishMols {

  for my $id (keys %id2type) {
    my $id_type = $id2type{$id};
    if ($id_type eq "bp:protein" || $id_type eq "bp:rna") {
      if (defined $id2xref{$id}) {
        for my $xref (keys %{ $id2xref{$id} }) {
          if (defined $unification{$xref}{UniProt}) {
            my $uniprot = $unification{$xref}{UniProt};
            $pw->AddMolExId(FinalId("molecule", LookUp($id)), "UP", $uniprot);
          }
        }
      }
    }
    if ($id_type eq "bp:smallMolecule") {
      if (defined $id2xref{$id}) {
        for my $xref (keys %{ $id2xref{$id} }) {
          if (defined $unification{$xref}{ChEBI}) {
            my $chebi = $unification{$xref}{ChEBI};
            $pw->AddMolExId(FinalId("molecule", LookUp($id)), "CH", $chebi);
          }
        }
      }
    }
    if ($id_type eq "bp:protein" || $id_type eq "bp:rna" ||
        $id_type eq "bp:complex" || $id_type eq "bp:smallMolecule" ||
        $id_type eq "bp:physicalEntity") {
      my (@preferred);
      for my $name (
          keys %{ $id2short_name{$id} },
          keys %{ $id2name{$id} },
          keys %{ $id2synonym{$id} }
        ) {
## lose the cellular location part of the name
        if ($name =~ / *\[[^\]]+\]$/) {
          print STDERR "dropping cellular location part of $name\n";
          $name =~ s/ *\[[^\]]+\]$//;
        }
        if (length($name) <= 200) {
          push @preferred, $name;
        } else {
          print STDERR "dropping name with length > 200: $name\n";
        }
      }
      if (@preferred) {
        for my $name (@preferred) {
          $pw->AddMolName(FinalId("molecule", LookUp($id)), "AS", $name);
        }
      } else {
        print STDERR "no name for $id\n";
      }
    }
  }
  AddHUGOSymbols();
}

######################################################################
sub PropagateWholeProteinIDs {

  for my $mol (@{ $pw->Mols() }) {
    my @parts = @{ $pw->MolPart($mol) };
    if (@parts) {
      my @ups;
      my $hash = $pw->MolExId($mol);
      my $name = $pw->PickMolName($mol);
      if (defined $$hash{"UP"}) {
        @ups = keys %{ $$hash{"UP"} };
        for my $part (@parts) {
          my $part_hash = $pw->MolExId($part);
          if (! defined $$part_hash{"UP"}) {
            for my $up (@ups) {
              $pw->AddMolExId($part, "UP", $up);
            }
          }
          my $name_hash = $pw->MolName($part);
          if (keys %{ $name_hash } == 0) {
            my $bounds = $pw->PartBounds($part);
            $bounds =~ s/,/_/;
            $pw->AddMolName($part, "PF", join("_", $name, "subunit",
                $bounds));
          }
        }
      }
    }
  }
}


######################################################################
sub DoFeatures {
  my ($what, $container_id, $sequence_id, $participant, $mol) = @_;

  ## container_id is interaction or complex
  ## sequence_id is edge_id or component_id

  my $physical_entity = $participant2physical_entity{$participant};
  for my $feature (keys %{ $sequence_participant2feature{$participant} }) {
    my $feature_type = $sequence_feature2type{$feature};
    if ($feature_type =~ /^LengthFeature/) {
      ## make a subunit
      my $fl = $feature2location{$feature};
      my $begin_site = $interval2begin{$fl};
      my $end_site   = $interval2end{$fl};
      my $begin_pos  = $site2position{$begin_site};
      my $end_pos    = $site2position{$end_site};
      my $whole = $mol;    ## they don't distinguish
      my $part = FinalId("molecule", LookUp(join("_",
          $mol, "subunit", $begin_pos, $end_pos)));
      ##  can't hurt to add it more than once, right?
      $pw->AddMol($part, $lv->BasicMolTypeCode("PR"));
      $pw->AddMolPart($part, $whole, $begin_pos, $end_pos);
      ## now replace the whole with the part
      $mol = $part;
    } else {
      my ($modification_label_name, $modification_label_value);
      my ($protein, $position, $amino_acid);
      $protein = "";
      $position = "0";
      $amino_acid = "X";
      ## try to nail down the protein
      for my $xref (keys %{ $id2xref{$physical_entity} }) {
        if (defined $unification{$xref}{UniProt}) {
          $protein = $unification{$xref}{UniProt};
          last;
        }
      }
      ## try to nail down the residue
      if (defined $feature2residue{$feature}) {
        my $residue = $feature2residue{$feature};
        if (defined $residue_type{$residue}) {
          $amino_acid = $residue_type{$residue};
        }
      }
      ## try to nail down the position
      if (defined $feature2location{$feature}) {
        my $site = $feature2location{$feature};
        if (defined $site2position{$site}) {
          $position = $site2position{$site};
        }
      }
      if (defined $reactome_feature_types{$feature_type}) {
        $modification_label_name = $reactome_feature_types{$feature_type};
      }
      if (!defined $modification_label_name ||
          $modification_label_name eq "") {
print STDERR "unknown_feature_type\t$feature_type for feature $feature\n";
        $modification_label_name = "ptm";
      }
      $modification_label_value = $lv->StringToLabelValue("ptm",
          $modification_label_name);
      if ($what eq "edge") {
        $pw->AddEdgePTM($container_id, $sequence_id, $protein,
            $position, $amino_acid,
            $modification_label_value, $modification_label_name);
      } elsif ($what eq "component") {
        $pw->AddComponentPTM($container_id, $sequence_id, $protein,
            $position, $amino_acid,
            $modification_label_value, $modification_label_name);
      } else {
print STDERR "DoFeatures: bad 'what' = $what\n";
      }
    }
  }
  return $mol;
}

######################################################################
sub FinishComplexes {

## Reactome components do not have location

  for my $complex (keys %complex2component) {
    my $complex_id = FinalId("molecule", LookUp($complex));
    my $comp;
    for my $component (keys %{ $complex2component{$complex} }) {
      my $stoichiometry;
      if (defined $participant2stoichiometry{$component}) {
        $stoichiometry = $participant2stoichiometry{$component};
      } else {
        $stoichiometry = 1;
      }
      if ($stoichiometry > 8) {
print STDERR "reducing stoichiometry from $stoichiometry to 2 for " .
    " component $component of complex $complex\n";
        $stoichiometry = 2;
      }
      while ($stoichiometry) {
        $comp++;
        $stoichiometry--;
        my $component_mol = FinalId("molecule",
            LookUp($participant2physical_entity{$component}));
        $component_mol = DoFeatures("component", $complex_id, $comp,
            $component, $component_mol);
        $pw->AddComponent($complex_id, $comp, $component_mol);
      }
    }
  }
  
}

######################################################################
sub CheckInteractionType {
  my ($id) = @_;

  ## check atoms if modification type to see if they should actually
  ## be declared to be transloactions

  my $atomtype_lvid = $pw->AtomType($id);
  my ($dummy, $atomtype) = $lv->LabelValueToString($atomtype_lvid);
  if ($atomtype ne "modification") {
    return;
  }

  my $translocation_type =
      $lv->StringToLabelValue("process-type", "translocation");

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
      my $ptms = $pw->EdgePTM($id, $edge);
      if (@{ $ptms }) {
        $ptms{$edge} = Pathway::NormalPTMString($ptms);
      }
    }    

  }

  if ($pr_in != $pr_out || $cx_in != $cx_out) {
    return;
  }

  ## if an input mol matches an output mol but there are 
  ## changes in ptms or state labels, then declare biochemicalConversion

  my ($location_change, $state_change);
  for my $mol (keys %mol2edge) {
    my @temp = keys %{ $mol2edge{$mol} };
    if (@temp == 1) {
      $no_match{$mol2type{$mol}}++;
    }
    for (my $i = 0; $i < @temp - 1; $i++) {
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
  if (! $state_change && $location_change) {
    $pw->ResetAtomType($id, $translocation_type);
  }
}

######################################################################
sub FinishInteractions {

##!! need to add mol labels and edge ptms

  for my $reaction (keys %reaction2participant) {
    if ($id2type{$reaction} eq "bp:catalysis") {
      next;
    }
    my $edge;
    my $atomid = FinalId("interaction", LookUp($reaction));
    for my $xref (keys %{ $id2xref{$reaction} }) {
      if (defined $publication{$xref}{"Pubmed"}) {
        $pw->AddAtomReferences($atomid, $publication{$xref}{"Pubmed"});
      }
    }
    for my $type (keys %{ $reaction2participant{$reaction} }) {
      for my $participant (keys
          %{ $reaction2participant{$reaction}{$type} }) {
        my $mol = FinalId("molecule",
            LookUp($participant2physical_entity{$participant}));
        my $stoichiometry;
        if (defined $participant2stoichiometry{$participant}) {
          $stoichiometry = $participant2stoichiometry{$participant};
        } else {
          $stoichiometry = 1;
        }
        if ($stoichiometry > 8) {
print STDERR "reducing stoichiometry from $stoichiometry to 2 for " .
    " component $participant of interaction $reaction\n";
          $stoichiometry = 2;
        }
        while ($stoichiometry) {
          $edge++;
          $stoichiometry--;
          my $edge_type;
          if ($type eq "bp:LEFT") {
            $edge_type = $lv->StringToLabelValue("edge-type", "input");
          } elsif ($type eq "bp:RIGHT") {
            $edge_type = $lv->StringToLabelValue("edge-type", "output");
          } elsif ($type eq "ACTIVATION" ||
                   $type eq "ACTIVATION_2" ||
                   $type eq "ACTIVATION-ALLOSTERIC_2") {
            $edge_type = $lv->StringToLabelValue("edge-type", "agent");
          } elsif ($type eq "INHIBITION" ||
                   $type eq "INHIBITION-COMPETITIVE" ||
                   $type eq "INHIBITION-COMPETITIVE_2") {
            $edge_type = $lv->StringToLabelValue("edge-type", "inhibitor");   
          } else {
            print STDERR "#unknown_edge_type\t$type\n";
          }
          if (defined $participant2cellular_location{$participant}) {
            my $location = $participant2cellular_location{$participant};
            if (defined $id2xref{$location}) {
              for my $xref (keys %{ $id2xref{$location} }) {
                if (defined $unification{$xref}{"GO"}) {
                  my $go_id = $unification{$xref}{"GO"};
                  my $lvid = $lv->StringToLabelValue("location", "GO:$go_id");
                  if ($lvid) {
                    if (defined $translate_location{$lvid}) {
                      $lvid = $translate_location{$lvid};
                    }
                    $pw->AddMolLabel($atomid, $edge, $lvid);
                  } else {
print STDERR "can't find label value for GO id \"$go_id\"\n";
                  }
                } else {
print STDERR "can't find GO unification for term \"$xref\"\n";
                }
              }
            } else {
print STDERR "can't find xref for term \"$location\"\n";
            }
          }
          $mol = DoFeatures("edge", $atomid, $edge, $participant, $mol);
          $pw->AddEdge($atomid, $edge, $edge_type, $mol);
          $edge2participant{$atomid}{$edge} = $participant;
        }
      }

    }
    CheckInteractionType($atomid);
  }
}

######################################################################
sub AddHUGOSymbols {

  use DBI;
  use constant MAX_LONG_LEN       => 16384;
  use constant ORACLE_LIST_LIMIT  => 500;

  my ($db_inst, $db_user, $db_pass, $schema) =
     ("cgprod", "web", "readonly", "cgap");
  my ($sql, $stm);
  my ($i, $list);
  my $protein_type = $lv->StringToLabelValue("molecule-type", "protein");
  my (%eg, %up);
  my ($eg, $up, $sym);

  for my $mol (@{ $pw->Mols }) {
    if ($pw->MolType($mol) eq $protein_type) {
      my $hash = $pw->MolExId($mol);
      for my $idtype (keys %{ $hash }) {
        if ($idtype eq "LL") {
          for my $eg (keys %{ $$hash{$idtype} }) {
            $eg{$eg}{$mol} = 1;
          }
        } elsif ($idtype eq "UP") {
          for my $up (keys %{ $$hash{$idtype} }) {
            $up{$up}{$mol} = 1;
          }
        }
      }
    }
  }

  my $db = DBI->connect("DBI:Oracle:" . $db_inst, $db_user, $db_pass);
  if (not $db or $db->err()) {
    print STDERR "Cannot connect to " . $db_user . "@" . $db_inst . "\n";
    die;
  }

  my @list = keys %up;
  for($i = 0; $i < @list; $i += ORACLE_LIST_LIMIT) {
 
    if(($i + ORACLE_LIST_LIMIT - 1) < @list) {
      $list = "'" . join("','", @list[$i..$i+ORACLE_LIST_LIMIT-1]) . "'";
    }
    else {
      $list = "'" . join("','", @list[$i..@list-1]) . "'";
    }
    $sql = qq!
select
  p.sp_id_or_secondary,
  g.symbol
from
  $schema.sp_primary p,
  $schema.ll2sp s,
  $schema.ll_gene g
where
      s.organism = 'Hs'
  and g.organism = 'Hs'
  and g.symbol_type = 'OF'
  and g.ll_id = s.ll_id
  and s.sp_primary = p.sp_primary
  and p.id_or_accession = 'a'
  and p.sp_id_or_secondary in ($list)
!;
    $stm = $db->prepare($sql);
    if(not $stm) {
      print STDERR "$sql\n";
      print STDERR "$DBI::errstr\n";
      print STDERR "prepare call failed\n";
      die;
    }
    if(!$stm->execute()) {
      print STDERR "$sql\n";
      print STDERR "$DBI::errstr\n";
      print STDERR "execute call failed\n";
      die;
    }
    while (($up, $sym) = $stm->fetchrow_array()) {
      for my $mol (keys %{ $up{$up} }) {
        if (@{ $pw->MolWhole($mol) }) {
          $pw->AddMolName($mol, "PF", "$sym subunit");
        } else {
          $pw->AddMolName($mol, "OF", $sym);
          if (defined $pw->{molname}{$mol}{AS}{$sym}) {
            delete $pw->{molname}{$mol}{AS}{$sym};
          }
        }
      }
    }
  }

  my @list = keys %eg;
  for($i = 0; $i < @list; $i += ORACLE_LIST_LIMIT) {
 
    if(($i + ORACLE_LIST_LIMIT - 1) < @list) {
      $list = "'" . join("','", @list[$i..$i+ORACLE_LIST_LIMIT-1]) . "'";
    }
    else {
      $list = "'" . join("','", @list[$i..@list-1]) . "'";
    }
    $sql = qq!
select  
  g.ll_id,
  g.symbol
from
  $schema.ll_gene g
where
      g.organism = 'Hs'
  and symbol_type = 'OF'
  and g.ll_id in ($list)
!;
    $stm = $db->prepare($sql);
    if(not $stm) {
      print STDERR "$sql\n";
      print STDERR "$DBI::errstr\n";
      print STDERR "prepare call failed\n";
      die;
    }
    if(!$stm->execute()) {
      print STDERR "$sql\n";
      print STDERR "$DBI::errstr\n";
      print STDERR "execute call failed\n";
      die;
    }
    while (($eg, $sym) = $stm->fetchrow_array()) {
      for my $mol (keys %{ $eg{$eg} }) {
        if (@{ $pw->MolWhole($mol) }) {
          $pw->AddMolName($mol, "PF", "$sym subunit");
        } else {
          $pw->AddMolName($mol, "OF", $sym);
          if (defined $pw->{molname}{$mol}{AS}{$sym}) {
            delete $pw->{molname}{$mol}{AS}{$sym};
          }
        }
      }
    }
  }

  $db->disconnect();

}

######################################################################
sub parse {
  my ($self, $fh) = @_;
  $p->parse($fh);
  $self->DoGraphStuff();

  ## Painful: Check to see if we can distinguish instances that
  ## have the same physical entity and the same explicit features
  ## and locations but differ, nevertheless, in internal Reactome DB_ID.
  ## On finding these, use dummy state labels to distinguish them.
  ## Note: we are only doing this on interaction edges, not on components
  ## of complexes.

  my @state_labels = (
    $lv->StringToLabelValue("activity-state", "state_1"),
    $lv->StringToLabelValue("activity-state", "state_2"),
    $lv->StringToLabelValue("activity-state", "state_3"),
    $lv->StringToLabelValue("activity-state", "state_4"),
    $lv->StringToLabelValue("activity-state", "state_5")
  );
  my %inst2db_id;
  for my $atom (@{ $pw->Atoms() }) {
    for my $edge (@{ $pw->Edges($atom) }) {
      my $s = $pw->NormalMolInstString($atom, $edge);
      $inst2db_id{$s}{$participant2db_id{$edge2participant{$atom}{$edge}}}
          {"$atom,$edge"} = 1;
    }
  }
  for my $i (keys %inst2db_id) {
    if (keys %{ $inst2db_id{$i} } > 1) {
      my $state_n = -1;
      for my $db_id (keys %{ $inst2db_id{$i} }) {
        $state_n++;
        if ($state_n < @state_labels) {
          for my $ae (keys %{ $inst2db_id{$i}{$db_id} }) {
            my ($atom, $edge) = split(",", $ae);
            $pw->AddMolLabel($atom, $edge, $state_labels[$state_n]);
          }
        } else {
          print STDERR "not enough distinct states for instance $i\n";
        }
      }
    }
  }

  $pw->BuildMolInstCache();
  $self->InferConnections();
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

    if ($x eq "rdf:ID") {
      $id_def{$elem} = $$attr_top{"rdf:ID"};
      $id2type{$$attr_top{"rdf:ID"}} = $elem;
    } elsif ($x eq "rdf:resource") {
      $$attr_top{"rdf:resource"} =~ s/^#//;
    }

  }

}

######################################################################
sub handle_end {
  my ($p, $elem) = @_;

  my ($attr_top, $elem_top) = (attr_top(), elem_top());

  if ($elem eq "bp:ORGANISM") {

  } elsif ($elem eq "bp:pathway") {

    my $pathway = $attr_stack[0]{"rdf:ID"};
    for my $pway_component (@pathway_components) {
      $network2component{$pathway}{$pway_component} = 1;
    }
    undef @pathway_components;

  } elsif ($elem eq "bp:AVAILABILITY") {

  } elsif ($elem eq "bp:NAME") {

    $id2name{$attr_stack[1]{"rdf:ID"}}{$char_value{$elem}} = 1;

  } elsif ($elem eq "bp:SHORT-NAME") {

    $id2short_name{$attr_stack[1]{"rdf:ID"}}{$char_value{$elem}} = 1;

  } elsif ($elem eq "bp:SYNONYMS") {

    $id2synonym{$attr_stack[1]{"rdf:ID"}}{$char_value{$elem}} = 1;

  } elsif ($elem eq "bp:PATHWAY-COMPONENTS") {

    push @pathway_components, $attr_stack[0]{"rdf:resource"};

  } elsif ($elem eq "bp:pathwayStep") {

    my $pathway_step = $attr_stack[0]{"rdf:ID"};
    for my $interaction (@step_interactions) {
      $network2component{$pathway_step}{$interaction} = 1;
    }

    undef @step_interactions;
    undef @next_steps;

  } elsif ($elem eq "bp:STEP-INTERACTIONS") {

    push @step_interactions, $attr_stack[0]{"rdf:resource"};

  } elsif ($elem eq "bp:NEXT-STEP") {

    push @next_steps, $attr_stack[0]{"rdf:resource"};

  } elsif ($elem eq "bp:transport") {

  } elsif ($elem eq "bp:biochemicalReaction") {

    my $edge;
    my $atomid = FinalId("interaction", LookUp($attr_stack[0]{"rdf:ID"}));
    $pw->AddAtom($atomid,
        $lv->StringToLabelValue("process-type", "modification"));
    $pw->AddAtomSource($atomid, REACTOME_SOURCE_ID);

  } elsif ($elem eq "bp:LEFT") {

    my $left = $attr_stack[0]{"rdf:resource"};
    my $reaction = $attr_stack[1]{"rdf:ID"};
    $reaction2participant{$reaction}{$elem}{$left} = 1;

  } elsif ($elem eq "bp:RIGHT") {

    my $right = $attr_stack[0]{"rdf:resource"};
    my $reaction = $attr_stack[1]{"rdf:ID"};
    $reaction2participant{$reaction}{$elem}{$right} = 1;

  } elsif ($elem eq "bp:EVIDENCE") {

  } elsif ($elem eq "bp:evidence") {

  } elsif ($elem eq "bp:EVIDENCE-CODE") {

  } elsif ($elem eq "bp:sequenceParticipant" ||
           $elem eq "bp:physicalEntityParticipant") {

    ## Looks like sequenceParticipant and physicalEntityParticipant
    ## can each have a PHYSICAL-ENTITY and a CELLULAR-LOCATION.
    ## Also, participants are defined outside of reactions and are
    ## referenced inside reactions.

  } elsif ($elem eq "bp:SEQUENCE") {

    print STDERR "$elem not expected\n";

  } elsif ($elem eq "bp:SEQUENCE-FEATURE-LIST") {

    my $sequence_participant = $attr_stack[1]{"rdf:ID"};
    $sequence_participant2feature{$sequence_participant}
        {$attr_stack[0]{"rdf:resource"}} = 1;

  } elsif ($elem eq "bp:sequenceFeature") {

  } elsif ($elem eq "bp:FEATURE-LOCATION") {

    $feature2location{$attr_stack[1]{"rdf:ID"}} =
        $attr_stack[0]{"rdf:resource"};

  } elsif ($elem eq "bp:FEATURE-TYPE") {

    my $sequence_feature = $attr_stack[1]{"rdf:ID"};
    $sequence_feature2type{$sequence_feature} = $attr_stack[0]{"rdf:resource"};

#print STDERR "added feature " . $attr_stack[0]{"rdf:resource"} . "\n";

  } elsif ($elem eq "bp:sequenceSite") {

  } elsif ($elem eq "bp:sequenceInterval") {

  } elsif ($elem eq "bp:POSITION-STATUS") {

    $site2position_status{$attr_stack[1]{id}} = $char_value{$elem};

  } elsif ($elem eq "bp:SEQUENCE-POSITION") {

    $site2position{$attr_stack[1]{"rdf:ID"}} = $char_value{$elem};

  } elsif ($elem eq "bp:SEQUENCE-INTERVAL-BEGIN") {

    $interval2begin{$attr_stack[1]{"rdf:ID"}} =
        $attr_stack[0]{"rdf:resource"};

  } elsif ($elem eq "bp:SEQUENCE-INTERVAL-END") {

    $interval2end{$attr_stack[1]{"rdf:ID"}} =
        $attr_stack[0]{"rdf:resource"};

  } elsif ($elem eq "bp:PARTICIPANTS") {

    print STDERR "$elem not expected\n";

  } elsif ($elem eq "bp:PHYSICAL-ENTITY") {

    my $participant = $attr_stack[1]{"rdf:ID"};
    my $physical_entity = $attr_stack[0]{"rdf:resource"};
    $participant2physical_entity{$participant} = $physical_entity;

  } elsif ($elem eq "bp:CELLULAR-LOCATION") {

    my $participant = $attr_stack[1]{"rdf:ID"};
    my $cellular_location = $attr_stack[0]{"rdf:resource"};
    $participant2cellular_location{$participant} = $cellular_location;

  } elsif ($elem eq "bp:STOICHIOMETRIC-COEFFICIENT") {

    my $participant = $attr_stack[1]{"rdf:ID"};
    my $stoichiometry = $char_value{$elem};
    $participant2stoichiometry{$participant} = $stoichiometry;

  } elsif ($elem eq "bp:INTERACTION-TYPE") {

  } elsif ($elem eq "bp:physicalInteraction") {

  } elsif ($elem eq "bp:modulation") {

  } elsif ($elem eq "bp:catalysis") {

  } elsif ($elem eq "bp:EC-NUMBER") {

  } elsif ($elem eq "bp:CONTROLLER") {

    my $catalysis = $attr_stack[1]{"rdf:ID"};
    my $controller = $attr_stack[0]{"rdf:resource"};
    $catalysis2controller{$catalysis} = $controller;

  } elsif ($elem eq "bp:CONTROLLED") {

    my $catalysis = $attr_stack[1]{"rdf:ID"};
    my $controlled = $attr_stack[0]{"rdf:resource"};
    $catalysis2controlled{$catalysis} = $controlled;

  } elsif ($elem eq "bp:SPONTANEOUS") {

  } elsif ($elem eq "bp:DIRECTION") {

    my $direction = $char_value{$elem};
    my $catalysis = $attr_stack[1]{"rdf:ID"};
    $catalysis2direction{$catalysis} = $direction;

  } elsif ($elem eq "bp:CONTROL-TYPE") {

    my $control_type = $char_value{$elem};
    my $catalysis = $attr_stack[1]{"rdf:ID"};
    $catalysis2control_type{$catalysis} = $control_type;

  } elsif ($elem eq "bp:physicalEntity") {

    my $molid = FinalId("molecule", LookUp($attr_stack[0]{"rdf:ID"}));
    $pw->AddMol($molid,
        $lv->StringToLabelValue("molecule-type", "molecule-type"));

  } elsif ($elem eq "bp:rna") {

    my $molid = FinalId("molecule", LookUp($attr_stack[0]{"rdf:ID"}));
    $pw->AddMol($molid, $lv->BasicMolTypeCode("RN"));

  } elsif ($elem eq "bp:protein") {

    my $molid = FinalId("molecule", LookUp($attr_stack[0]{"rdf:ID"}));
    $pw->AddMol($molid, $lv->BasicMolTypeCode("PR"));

  } elsif ($elem eq "bp:smallMolecule") {

    my $molid = FinalId("molecule", LookUp($attr_stack[0]{"rdf:ID"}));
    $pw->AddMol($molid, $lv->BasicMolTypeCode("CM"));

  } elsif ($elem eq "bp:complex") {

    my $molid = FinalId("molecule", LookUp($attr_stack[0]{"rdf:ID"}));
    $pw->AddMol($molid, $lv->BasicMolTypeCode("CX"));

    for my $participant (@components) {
      $complex2component{$attr_stack[0]{"rdf:ID"}}{$participant} = 1;
    }
    undef @components;

  } elsif ($elem eq "bp:COMPONENTS") {

    push @components, $attr_stack[0]{"rdf:resource"};

  } elsif ($elem eq "bp:openControlledVocabulary") {

  } elsif ($elem eq "bp:XREF") {

    $id2xref{$attr_stack[1]{"rdf:ID"}}{$attr_stack[0]{"rdf:resource"}} = 1;

  } elsif ($elem eq "bp:unificationXref") {

     $unification{$attr_stack[0]{"rdf:ID"}}{$xref_db} = $xref_id;

     undef $xref_db;
     undef $xref_id;

  } elsif ($elem eq "bp:relationshipXref") {

     $relationship{$attr_stack[0]{"rdf:ID"}}{$xref_db} = $xref_id;

     undef $xref_db;
     undef $xref_id;

  } elsif ($elem eq "bp:RELATIONSHIP-TYPE") {


  } elsif ($elem eq "bp:publicationXref") {

     $publication{$attr_stack[0]{"rdf:ID"}}{$xref_db} = $xref_id;

     undef $xref_db;
     undef $xref_id;

  } elsif ($elem eq "bp:DB") {

    if ($elem_stack[1] eq "bp:unificationXref") {
      $xref_db = $char_value{$elem}
    } elsif ($elem_stack[1] eq "bp:relationshipXref") {
      $xref_db = $char_value{$elem}
    } elsif ($elem_stack[1] eq "bp:publicationXref") {
      $xref_db = $char_value{$elem}
    }    

  } elsif ($elem eq "bp:ID") {

    if ($elem_stack[1] eq "bp:unificationXref") {
      $xref_id = $char_value{$elem}
    } elsif ($elem_stack[1] eq "bp:relationshipXref") {
      $xref_id = $char_value{$elem}
    } elsif ($elem_stack[1] eq "bp:publicationXref") {
      $xref_id = $char_value{$elem}
    }    

  } elsif ($elem eq "bp:ID-VERSION") {

  } elsif ($elem eq "bp:YEAR") {

  } elsif ($elem eq "bp:TITLE") {

  } elsif ($elem eq "bp:AUTHORS") {

  } elsif ($elem eq "bp:SOURCE") {

  } elsif ($elem eq "bp:COMMENT") {

    if ($char_value{$elem} =~ /^(Authored|Reviewed|Edited):\s*(.+$)/) {
      my ($date);
      my ($type, $name) = ($1, $2);
      if ($name =~ /(.+),?\s*(\d\d\d\d-\d\d-\d\d)\s*$/) {
        ($name, $date) = ($1, $2);
      }
      $name =~ s/,\s+$//;
      if ($elem_stack[1] eq "bp:pathway") {
        my $pid = FinalId("pathway", LookUp($attr_stack[1]{"rdf:ID"}));
        if ($type eq "Authored" || $type eq "Edited") {
          my $found;
          for my $c (@{$pw->Curators($pid) }) {
            if (index($c, $name) >= 0) {
              $found++;
              last;
            } elsif (index($name, $c) >= 0) {
              delete $pw->{pathwaycurator}{$pid}{$c};
              $pw->AddCurator($pid, $name);
              $found++;
              last;
            }
          }
          if (! $found) {
            $pw->AddCurator($pid, $name);
          }
        } elsif ($type eq "Reviewed") {
          $pw->AddReviewer($pid, $name);
        }
        if ($date) {
        }
      }
    } elsif ($char_value{$elem} =~ /^Reactome DB_ID:\s*(\d+)/) {
       ## This is always in a sequenceParticipant or a physicalEntityParticipant
       my $db_id = $1;
       $participant2db_id{$attr_stack[1]{"rdf:ID"}} = $db_id;
#       print STDERR "##DB_ID\t$1\t" . $attr_stack[1]{"rdf:ID"} . "\n";
     } elsif ($char_value{$elem} =~ /^residue:\s*(.+)/) {
       my $residue = $1;
       $feature2residue{$attr_stack[1]{"rdf:ID"}} = $residue;
     }

  } elsif ($elem eq "bp:SYNONYMS") {
  } elsif ($elem eq "bp:TAXON-XREF") {
  } elsif ($elem eq "bp:TERM") {
  } elsif ($elem eq "bp:bioSource") {
  } elsif ($elem eq "bp:DATA-SOURCE") {
  } elsif ($elem eq "bp:dataSource") {
  } elsif ($elem eq "bp:URL") {

  } else {

    print STDERR "#ignoring_tag\t$elem\n";

  }

  shift @elem_stack;
  shift @attr_stack;

}

######################################################################
sub FinalId {
  my ($id_type, $tmp_id) = @_;

  if ($id_type eq "pathway") {
    if (defined $final_id{pathway}{$tmp_id}) {
      return $final_id{pathway}{$tmp_id};
    } else {
      $next_pathway_id++;
      $final_id{pathway}{$tmp_id} = $next_pathway_id;
      return $next_pathway_id;
    }
  } elsif ($id_type eq "interaction") {
    if (defined $final_id{interaction}{$tmp_id}) {
      return $final_id{interaction}{$tmp_id};
    } else {
      $next_interaction_id++;
      $final_id{interaction}{$tmp_id} = $next_interaction_id;
      return $next_interaction_id;
    }
  } elsif ($id_type eq "molecule") {
    if (defined $final_id{molecule}{$tmp_id}) {
      return $final_id{molecule}{$tmp_id};
    } else {
      $next_molecule_id++;
      $final_id{molecule}{$tmp_id} = $next_molecule_id;
      return $next_molecule_id;
    }
  } else {
    print STDERR "unrecognized id type $id_type " .
        "provided to FinalId for id = $tmp_id\n";
    return 0;
  }
}

######################################################################
sub LookUp {
  my ($rdf_id) = @_;

  if (defined $rdf_id2id{$rdf_id}) {
    return $rdf_id2id{$rdf_id};
  } else {
    $next_id++;
    $rdf_id2id{$rdf_id} = $next_id;
    $id2rdf_id{$next_id} = $rdf_id;
    return $next_id;
  }
}

######################################################################
sub handle_char {
  my ($p, $s) = @_;

  $char_value{elem_top()} .= $s;
}

######################################################################
sub new {
  my ($self, $labeldb, $pathway, $idmap) = @_;

  $lv  = $labeldb;
  $pw  = $pathway;
  $idm = $idmap;
  $pw->AddSource(REACTOME_SOURCE_ID, REACTOME_SOURCE_NAME);
  $p = new XML::Parser(Handlers => {
    Start   => \&handle_start,
    End     => \&handle_end,
    Char    => \&handle_char
  } );
  for my $id ($idm->PTMExprId(), $idm->AtomId(),
      $idm->PathwayId(), $idm->MolId()) {
    if ($id > $next_id) {
      $next_id = $id;
    }
  }
  if ($next_id > 0) {
    $next_id--;
  }
  $next_pathway_id     = $idm->PathwayId();
  $next_interaction_id = $idm->AtomId();
  $next_molecule_id    = $idm->MolId();
  $next_ptm_id         = $idm->PTMExprId();

  $translate_location{$lv->StringToLabelValue("location", "cytosol")} =
      $lv->StringToLabelValue("location", "cytoplasm");
  $translate_location{$lv->StringToLabelValue("location", "nucleoplasm")} =
      $lv->StringToLabelValue("location", "nucleus");
  $translate_location{$lv->StringToLabelValue("location", "mitochonrdial matrix")} =
      $lv->StringToLabelValue("location", "mitochondria");
  $translate_location{$lv->StringToLabelValue("location", "endoplasmic reticulum lumen")} =
      $lv->StringToLabelValue("location", "endoplasmic reticulum");
  $translate_location{$lv->StringToLabelValue("location", "peroxisomal matrix")} =
      $lv->StringToLabelValue("location", "peroxisome");

  my $x = {};
  return bless $x;
}

######################################################################
1;
