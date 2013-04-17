#!/usr/local/bin/perl


# Copyright SRA International
#
# Distributed under the OSI-approved BSD 3-Clause License.
# See http://ncip.github.com/pathway-interaction-database/LICENSE.txt for details.


package StudyNet;
require Exporter;
@ISA = qw(Exporter);
#@EXPORT = qw(
#);

use strict;
use Pathway;
use PWLabel;

my $protein_type;
my $rna_type;
my $compound_type;
my $complex_type;

my $macro_process_type;
my $transcription_type;
my $translocation_type;
my $reaction_type;
my $modification_type;
my $binding_type;        ###!!! old

my $input_type;
my $agent_type;
my $inhibitor_type;
my $output_type;

my %is_output;

use constant NO_VAL => "-1";

######################################################################
sub new {
  my ($self, $pw, $lv) = @_;

  my $x = {};

  if (! defined $pw || ! defined $lv) {
    return bless $x;
  }

  $x->{pw} = $pw;
  $x->{lv} = $lv;

  $protein_type =
      $lv->StringToLabelValue("molecule-type", "protein");
  $compound_type =
      $lv->StringToLabelValue("molecule-type", "compound");
  $rna_type =
      $lv->StringToLabelValue("molecule-type", "rna");
  $complex_type =
      $lv->StringToLabelValue("molecule-type", "complex");

  $macro_process_type =
      $lv->StringToLabelValue("process-type", "macroprocess");
  $transcription_type =
      $lv->StringToLabelValue("process-type", "transcription");
  $translocation_type =
      $lv->StringToLabelValue("process-type", "translocation");
  $reaction_type =
      $lv->StringToLabelValue("process-type", "reaction");
  $modification_type =
      $lv->StringToLabelValue("process-type", "modification");
  $binding_type =    ###!!!!! old type
      $lv->StringToLabelValue("process-type", "binding");

  $input_type =
      $lv->StringToLabelValue("edge-type", "input");
  $agent_type =
      $lv->StringToLabelValue("edge-type", "agent");
  $inhibitor_type =
      $lv->StringToLabelValue("edge-type", "inhibitor");
  $output_type =
      $lv->StringToLabelValue("edge-type", "output");

  return bless $x;
}


######################################################################
sub and_3val {
  my (@vector) = @_;
  my $default = 1;
  for my $v (@vector) {
    if ($v eq "0" || $v eq "") {
      return 0;
    } elsif ($v eq NO_VAL) {
      $default = NO_VAL;
    }
  }
  return $default;
}

######################################################################
sub or_3val {
  my (@vector) = @_;
  my $default = 0;
  for my $v (@vector) {
    if ($v eq "1") {
      return 1;
    } elsif ($v eq NO_VAL) {
      $default = NO_VAL;
    }
  }
  return $default;
}

######################################################################
sub not_3val {
  my ($v) = @_;

  if ($v eq "1") {
    return 0;
  } elsif ($v eq "0" || $v eq "") {
    return 1;
  } else {
    return NO_VAL;
  }
}

######################################################################
sub CheckIfPresent {
  my ($self, $mol, $present_calls) = @_;

## Some rules for now:
## If compound, return true
## If protein or rna, if exists LL external id, then true iff
##    any external LL id for the mol is entered in %present_call
## If complex, return true iff CheckIfPresent for each component mol
##    (disregarding component labels)

  my $pw = $self->{pw};
  my $lv = $self->{lv};

  my $mol_type = $pw->MolType($mol);
  my $ext_ids = $pw->MolExId($mol);
  if ($mol_type eq $protein_type) {
    if (defined $$ext_ids{LL}) {
      my @temp;
      for my $ll (keys %{ $$ext_ids{LL} }) {
        if (defined $$present_calls{$ll}) {
          push @temp, $$present_calls{$ll};
        }
      }
      # check families; one hit is enough
      for my $child (@{ $pw->FamilyChildren($mol) }) {
        push @temp, $self->CheckIfPresent($child, $present_calls);
      }
      return or_3val(@temp);
    } else {
      return NO_VAL;
    }
  } elsif ($mol_type eq $rna_type) {
    if (defined $$ext_ids{LL}) {
      my @temp;
      for my $ll (keys %{ $$ext_ids{LL} }) {
        if (defined $$present_calls{$ll}) {
          push @temp, $$present_calls{$ll};
        }
      }
      return or_3val(@temp);
    } else {
      return NO_VAL;
    }
  } elsif ($mol_type eq $compound_type) {
    return 1;
  } elsif ($mol_type eq $complex_type) {
    my @temp;
    for my $comp (@{ $pw->Components($mol) }) {
      push @temp, $self->CheckIfPresent($pw->ComponentMol($mol, $comp),
          $present_calls);
    }
    return and_3val(@temp);
  }
}

######################################################################
sub FindActiveInteractions {
  my ($self, $atom_list, $present_calls, $active_calls, $ptms) = @_;

  # present_calls: Entrez Gene ID -> 0|1

  my $pw = $self->{pw};
  my $lv = $self->{lv};

  my ($more);
  my (%mol_inst_present, %mol_present, %atom_active,
      %molinst2mol,
      %input, %agent, %inhibitor, %output);
  my (@computable_atom);

  for my $mol (@{ $pw->UsedMols() }) {
    $mol_present{$mol} = $self->CheckIfPresent($mol, $present_calls);
  }

## rule for dealing with oscillation:
##   if the atom is ever active, then is is active

  for my $atom (@{ $atom_list } ) {

    if ($lv->IsA($pw->AtomType($atom), $macro_process_type)) {
      ## skip over macroprocess nodes
      next;
    }

    push @computable_atom, $atom;

    for my $edge (@{ $pw->Edges($atom) }) {

      my $mol = $pw->EdgeMol($atom, $edge);
      my $mol_inst = $pw->MolInstId($atom, $edge);
      my $edge_type = $pw->EdgeType($atom, $edge);
      $molinst2mol{$mol_inst} = $mol;

      ## if mol instance has no labels, then the mol instance
      ## is present if the mol is present

      if (scalar (@{ $pw->MolLabel($atom, $edge) }) == 0) {
        if (defined $mol_present{$mol}) {
          $mol_inst_present{$mol_inst} = $mol_present{$mol};
        }
      }

      if ($edge_type eq $input_type) {
        $input{$atom}{$mol_inst} = 1;
      } elsif ($edge_type eq $output_type) {
        $output{$atom}{$mol_inst} = 1;
      } elsif ($edge_type eq $agent_type) {
        $agent{$atom}{$mol_inst} = 1;
      } elsif ($edge_type eq $inhibitor_type) {
        $inhibitor{$atom}{$mol_inst} = 1;
      }

    }
  }

  my @temp;
  my $v;
  $more = 1;
  while ($more) {
    $more = 0;
    for my $atom (@computable_atom) {
      if ($atom_active{$atom} == 1) {
        next;
      }
      undef @temp;

      for my $mol_inst (keys %{ $input{$atom} }) {
        if (defined $mol_inst_present{$mol_inst}) {
          push @temp, $mol_inst_present{$mol_inst};
        } else {
          push @temp, NO_VAL;
        }
      }

      for my $mol_inst (keys %{ $agent{$atom} }) {
        if (defined $mol_inst_present{$mol_inst}) {
          push @temp, $mol_inst_present{$mol_inst};
        } else {
          push @temp, NO_VAL;
        }
      }

      for my $mol_inst (keys %{ $inhibitor{$atom} }) {
        if (defined $mol_inst_present{$mol_inst}) {
          push @temp, not_3val($mol_inst_present{$mol_inst});
        } else {
          push @temp, not_3val(NO_VAL);
        }
      }

      for my $mol_inst (keys %{ $output{$atom} }) {
        if (defined $molinst2mol{$mol_inst}) {
          if (defined $mol_present{$molinst2mol{$mol_inst}}) {
            push @temp, $mol_present{$molinst2mol{$mol_inst}};  ###!!!!!
          } else {
            push @temp, NO_VAL;
          }
        } else {
          push @temp, NO_VAL;
        }
      }

      $v = or_3val($atom_active{$atom}, and_3val(@temp));

      if (! defined $atom_active{$atom} || $v != $atom_active{$atom}) {
        $more = 1;
        $atom_active{$atom} = $v;
        for my $mol_inst (keys %{ $output{$atom} }) {
          $mol_inst_present{$mol_inst} =
              or_3val($mol_inst_present{$mol_inst}, $v);
        }
      }
    }
    last;
  }
  for my $atom (keys %atom_active) {
    $$active_calls{$atom} = $atom_active{$atom};
    for my $mol_inst (keys %{ $output{$atom} }) {
      my ($mol, $labels) = 
          Pathway::DecodeMolInstString(
          $pw->MolInstIdToString($mol_inst));
      if ($pw->MolType($mol) == $protein_type) {
        my (%labs, %ll_ids);
        for my $v (@{ $labels }) {
          my ($label, $name) = $lv->LabelValueToString($v);
          $labs{$name} = 1;
        }
        my $hash = $pw->MolExId($mol);
        if (defined $$hash{LL}) {
          for my $ll_id (keys %{ $$hash{LL} }) {
            $ll_ids{$ll_id} = 1;
          }
        }
        if (keys %labs) {
          $$ptms{join("\t",
              $pw->PickMolName($mol),
              join(",", keys %ll_ids),
              join(",", sort keys %labs)
              )} = 1;
        }
      }
    }
  }
}

######################################################################
sub FindRoots {
  my ($self) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};

  for my $atom (@{ $pw->Atoms() } ) {
    my (%ins, %outs);
    for my $edge (@{ $pw->Edges($atom) }) {
      my $mol = $pw->EdgeMol($atom, $edge);
      my $edge_type = $pw->EdgeType($atom, $edge);
      if ($edge_type eq $input_type) {
        $ins{$mol} = 1;
      } elsif ($edge_type eq $output_type) {
        $outs{$mol} = 1;
      }
    }
    my @ins  = keys %ins;
    my @outs = keys %outs;
    my ($mi, $mj, $inout);
    for (my $i = 0; $i < @outs; $i++) {
      $mi = $outs[$i];
      undef $mj;
      $inout = 0;
      for (my $j = 0; $j < @ins; $j++) {
        $mj = $ins[$j];
        if ($mi eq $mj) {
          $inout = 1;
          last;
        }
      }
      if (! $inout) {
        $is_output{$mi} = 1;
      }
    }
  }

  for my $mol (@{ $pw->Mols }) {
    if (! $is_output{$mol}) {
      my $mol_type = $lv->BasicMolType($pw->MolType($mol));
      print join("\t", "root", $mol, $mol_type,
          $pw->PickMolName($mol)) . "\n";

    }
  }

}

#require 'save_atom2atom_stuff.pm';
#require 'atom2atom_stuff.pm';

######################################################################

my %atom2pathway;
my %pathway2atom;
my %atom2atom;
my %prune_mols;
my @subgraph;

use constant MIN_SUBNET_SIZE => 3;

my @CCF = (    ## common co-factors
    'water',
    'atp',
    'nad',
    'nadh',
    'nadph',
    'nadp',
    'oxygen',
    'adp',
    'orthophosphate',
    'coa',
    'carbon dioxide',
    'diphosphate',
    'ammonia',
    'h+'
  );

######################################################################
sub DumpAtomAdjacencyMatrix {
  my ($lv, $pw) = @_;

  for my $p (keys %pathway2atom) {
    my $pname = $pw->PathwayName($p);
    for my $a1 (keys %{ $pathway2atom{$p} }) {
      for my $a2 (keys %{ $atom2atom{$a1} }) {
        if ($a1 eq $a2) {
          next;
        }
        for my $e1 (@{ $pw->Edges($a1) }) {
          my $m1 = $pw->EdgeMol($a1, $e1);
          my $mname1 = $pw->PickMolName($m1);
          if ($pw->MolType($m1) ne $protein_type) {
            next;
          }
          for my $e2 (@{ $pw->Edges($a2) }) {
            my $m2 = $pw->EdgeMol($a2, $e2);
            my $mname2 = $pw->PickMolName($m2);
            if ($pw->MolType($m2) ne $protein_type) {
              next;
            }
            print "adjacent\t$pname\t$a1\t$a2\t$m1\t$mname1\t$m2\t$mname2\n";
          }
        }
      }
    }
  }
}

######################################################################
sub SetUpAtomAdjacencyMatrix {
  my ($lv, $pw) = @_;

  my %mol2atom;
  my %atom2mol;

  CommonCoFactors($lv, $pw, \%prune_mols);

  for my $atom (@{ $pw->Atoms() }) {
    if (@{ $pw->AtomPathway($atom) } < 2) {
      next;
    }
    for my $pid (@{ $pw->AtomPathway($atom) }) {
      $atom2pathway{$atom}{$pid} = 1;
      $pathway2atom{$pid}{$atom} = 1;
    }    
    for my $edge (@{ $pw->Edges($atom) }) {
      my $molid = $pw->EdgeMol($atom, $edge);
      if (! defined $prune_mols{$molid}) {
        $atom2mol{$atom}{$molid} = 1;
        $mol2atom{$molid}{$atom} = 1;
      }
    }
  }
  for my $a1 (keys %atom2mol) {
    for my $m (keys %{ $atom2mol{$a1} }) {
      for my $a2 (keys %{ $mol2atom{$m} }) {
        if ($a1 ne $a2) {
          $atom2atom{$a1}{$a2} = 1;
          $atom2atom{$a2}{$a1} = 1;
        }
      }
    }
  }
}

######################################################################
sub AtomAdjacency {
  my ($self) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};
  SetUpAtomAdjacencyMatrix($lv, $pw);
  DumpAtomAdjacencyMatrix($lv, $pw);
}

######################################################################
sub XXX {
  my ($lv, $pw) = @_;

  my $n = 0;
  $subgraph[$n] = {};  ## (comma-sep atom list) -> (hash-set of pway)

  ## First, find all pairs of adjacent atoms that co-occur
  ## in multiple pathways. These are minimal shared subnets and
  ## may be the seeds of larger shared subnets

  for my $a1 (keys %atom2atom) {
    for my $a2 (keys %{ $atom2atom{$a1} }) {
      my $pway_set = IntersectSetsAsHashes(
          List2Hash([ @{ $pw->AtomPathway($a1) } ]),
          List2Hash([ @{ $pw->AtomPathway($a2) } ])
        );
      if (keys %{ $pway_set } > 1) {
        $subgraph[$n]{join(",", sort ($a1, $a2))} = $pway_set;
      }
    }
  }

  ## Now try to grow the subnets from seeds

  my $more = 1;
  while ($more) {
    $n++;
#    if ($n > 10) { last; }
    $subgraph[$n] = {};
    $more = 0;

    for my $atom_list (keys %{ $subgraph[$n-1] }) {
      my $pway_set = $subgraph[$n-1]{$atom_list};
      my $atom_set = List2Hash([ split(",", $atom_list) ]);
      for my $a1 (keys %{ $atom_set }) {
        for my $a2 (keys %{ $atom2atom{$a1} }) {
          if (defined $$atom_set{$a2}) {
            next;  ## a2 is already in this subnet; don't pursue
          }
          my $inter = IntersectSetsAsHashes($pway_set,
              List2Hash([ @{ $pw->AtomPathway($a2) } ]));
          if (keys %{ $inter } > 1) {
            $subgraph[$n]{join(",", sort ($a2, split(",", $atom_list)))}
                = $inter;
            $more = 1;
          }
        }
      }
    }

    my %temp;   ## (comma-sep pway list) -> (comma-sep atom list)
    for my $atom_list (keys %{ $subgraph[$n] }) {
      push @{ $temp{join(",",
          sort keys %{ $subgraph[$n]{$atom_list} })} }, $atom_list;
    }
    my %temp1;  ## (comma-sep pway list) -> (array of hash-sets of atoms)
    for my $pway_list (keys %temp) {
      $temp1{$pway_list} = ReFactor($temp{$pway_list});
    }

    ## Attempt to merge the newly expanded subnets

    ## Restore format: (comma-sep atom list) -> (hash-set of pway)
    $subgraph[$n] = {};
    for my $pway_list (keys %temp1) {
      for my $atom_hash (@{ $temp1{$pway_list } }) {
        my $alist = join(",", sort keys %{ $atom_hash });
        for my $p (split(",", $pway_list)) {
          $subgraph[$n]{$alist}{$p} = 1;
        }
      }
    }

  }

  my %final_subgraphs;
  for (my $i = 0; $i < @subgraph; $i++) {
    for my $atomlist (keys %{ $subgraph[$i] }) {
      if (defined $final_subgraphs{$atomlist}) {
        $final_subgraphs{$atomlist} = UnionSetsAsHashes(
            $final_subgraphs{$atomlist}, $subgraph[$i]{$atomlist});
      } else {
        $final_subgraphs{$atomlist} = $subgraph[$i]{$atomlist};
      }
    }
  }

  print "  <SubnetList>\n";
  my $id = 0;
  for my $atomlist (keys %final_subgraphs) {
    my @atoms = split(",", $atomlist);
    if (@atoms >= MIN_SUBNET_SIZE) {
      $id++;
      print "    <Subnet id=\"$id\">\n";
      print "      <SubnetComponentList>\n";
      for my $atom (@atoms) {
        print "        <SubnetComponent idref=\"$atom\" />\n";
      }
      print "      </SubnetComponentList>\n";
      print "      <SubnetUseList>\n";
      for my $pid (keys %{ $final_subgraphs{$atomlist} }) {
        print "        <SubnetUse pathway_idref=\"$pid\" />\n";
      }
      print "      </SubnetUseList>\n";
      print "    </Subnet>\n";
    }
  }
  print "  </SubnetList>\n";

}

######################################################################
sub ReFactor {
  my ($s ) = @_;

## $s is a list of items, each item being a comma-separated
##   string of elements
## reduce recursively by merging any pairs of list items that share at
##   least one element

  ## transform to a list of hash sets
  my @t;
  for my $i (@{ $s }) {
    my %u;
    for my $j (split(",", $i)) {
      $u{$j} = 1;
    }
    push @t, \%u;
  }
  my $more = 1;
  while ($more) {
    my @kill;
    for (my $i = 0; $i < @t; $i++) {
      $kill[$i] = 0;
    }
    $more = 0;
    for (my $i = 0; $i < @t - 1; $i++) {
      if (keys %{ $t[$i] } < 1) {
        next;    ## skip empty sets
      }
      for (my $j = $i+1; $j < @t; $j++) {
        if (keys %{ $t[$j] } < 1) {
          next;  ## skip empty sets
        }
        my $h = IntersectSetsAsHashes($t[$i], $t[$j]);
        if (keys %{ $h } > 0) {
#print STDERR "merging " . join(",", sort keys %{ $t[$i] }) .
#              " and " . join(",", sort keys %{ $t[$j] }) . "\n";
          $t[$i] = UnionSetsAsHashes($t[$i], $t[$j]);
#          $t[$j] = {};
          $kill[$j] = 1;   ## can't actually delete yet
          $more = 1;
        }
      }
    }
    ## copy the survivors
    my @zz;
    for (my $i = 0; $i < @kill; $i++) {
      if (! $kill[$i]) {
        push @zz, $t[$i];
      }
    }
    undef @t;
    for my $x (@zz) {
      push @t, $x;
    }
  }

  ## save non-empty sets
  my @v;
  for my $y (@t) {
    push @v, $y;
  }
  return \@v;
  
}

######################################################################
sub List2Hash {
  my ($list) = @_;
  my %hash;
  for my $x (@{ $list }) {
    $hash{$x} = 1;
  }
  return \%hash;
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
sub UnionSetsAsHashes {
  my ($s1, $s2) = @_;
  my %s3;
  for my $s (keys %{ $s1 }, keys %{ $s2 }) {
    $s3{$s} = 1;
  }
  return \%s3;
}

######################################################################
sub CommonCoFactors {
  my ($lv, $pw, $ccf) = @_;

  my $COMPOUND_TYPE = $lv->StringToLabelValue("molecule-type", "compound");
  for my $molid (@{ $pw->Mols() }) {
    if ($pw->MolType($molid) eq $COMPOUND_TYPE) {
      my $h = $pw->MolName($molid);
      for my $name_type (keys %{ $h }) {
        for my $c (@CCF) {
          if (defined $$h{$name_type}{$c}) {
            $$ccf{$molid} = 1;
          }
        }
      }
    }
  }
}

######################################################################
sub FactorIntoSubNets {
  my ($self) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};
  SetUpAtomAdjacencyMatrix($lv, $pw);
  XXX($lv, $pw);
}

######################################################################
sub oldFactorIntoSubNets {
  my ($self) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};

  my %temp1;
  my %temp2;

  ## build atom X pid
  for my $atom (@{ $pw->Atoms() }) {
    for my $pid (@{ $pw->AtomPathway($atom) }) {
      $temp1{$atom}{$pid} = 1;
    }
  }
  ## build { pid } X { atom }
  for my $atom (keys %temp1) {
    if (keys %{ $temp1{$atom} } > 1) {
      $temp2{join(",", sort keys %{ $temp1{$atom} })}{$atom} = 1;
    }
  }
  for my $pid_set (keys %temp2) {
    if (keys %{ $temp2{$pid_set} } > 1) {
      print "p-set X a-set\t$pid_set\t" . join(",",
          keys %{ $temp2{$pid_set} }) . "\n";
    }
  }
}

######################################################################
1;
