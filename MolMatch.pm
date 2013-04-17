

# Copyright SRA International
#
# Distributed under the OSI-approved BSD 3-Clause License.
# See http://ncip.github.com/pathway-interaction-database/LICENSE.txt for details.


package MolMatch;
require Exporter;
@ISA = qw(Exporter);
#@EXPORT = qw(
#  CreateIsA
#  CreateLabelSort
#);

use strict;
use PWLabel;
use Pathway;

## pathway pa is the "standard"
## pathway pb is what is to be adapted to the standard

## name_index: lc(name)->nametype->molid
## exid_index: lc(exid)->idtype->molid

######################################################################
sub new {
  my ($self, $lv, $pa, $pb) = @_;
  my $x = {};
  $x->{lv} = $lv;
  $x->{pa} = $pa;
  $x->{pb} = $pb;
  MakeIndexes($x);
  return bless $x;
}

######################################################################
sub MakeIndexes {
  my ($self) = @_;

  my $lv = $self->{lv};
  my $pa = $self->{pa};
  my $pb = $self->{pb};

  for my $mol (@{ $pa->Mols() }) {
    my $nhash = $pa->MolName($mol);
    for my $nametype (keys %{ $nhash }) {
      for my $name (keys %{ $$nhash{$nametype} }) {
        $self->{name_index}{$nametype . ":" . lc($name)}{$mol} = 1;
      }
    }
    my $ehash = $pa->MolExId($mol);
    for my $idtype (keys %{ $ehash }) {
      for my $id (keys %{ $$ehash{$idtype} }) {
        $self->{exid_index}{$idtype . ":" . lc($id)}{$mol} = 1;
      }
    }
  }
}

######################################################################
sub MatchComplexMols {
  my ($self) = @_;

  my $lv = $self->{lv};
  my $pa = $self->{pa};
  my $pb = $self->{pb};

  my (%comp_lookup, %a_empties, %b_empties);

  my $complex_type = $lv->StringToLabelValue("molecule-type", "complex");

  ## build reverse index on "standard": component mol -> complex mol
  for my $mol (@{ $pa->Mols() }) {
    if ($pa->MolType($mol) == $complex_type) {
      if (@{ $pa->Components($mol) } == 0) {
        $a_empties{$mol} = 1;
      } else {
        for my $comp (@{ $pa->Components($mol) }) {
          my $compmol = $pa->ComponentMol($mol, $comp);
          $comp_lookup{$compmol}{$mol} = 1;
        }
      }
    }
  }

  ## map the empties first; they need to look like simples
  for my $mol (@{ $pb->Mols() }) {
    if ($pb->MolType($mol) != $complex_type) {
      next;
    }
    if (@{ $pb->Components($mol) } == 0) {
      for my $a_empty (keys %a_empties) {
        if (Pathway::NamesIntersect($pa, $pb, $a_empty, $mol)) {
          $pb->AddMolMap($mol, $a_empty);
          last;    ## just take the first one
        }
      }
    }
  }

  for my $mol (@{ $pb->Mols() }) {
    if ($pb->MolType($mol) != $complex_type) {
      next;
    }
    my $found = 0;
    my $hit = "";
    if (@{ $pb->Components($mol) } > 0) {
      ## just the first component; then test all candidates that
      ## share this component
      my $comp = @{ $pb->Components($mol) }[0];
      my $compmol = $pb->ComponentMol($mol, $comp);
      if (defined $pb->MolMap($compmol)) {
        $compmol = $pb->MolMap($compmol);
      }
      if (defined $comp_lookup{$compmol}) {
        for my $candidate (keys %{ $comp_lookup{$compmol} }) {
          $found = Pathway::EqComplex($pa, $pb, $candidate, $mol);
          if ($found) {
            $hit = $candidate;
            last;
          }
        }
      }
    }
    if ($found) {
      $pb->AddMolMap($mol, $hit);
    }
  }
}

######################################################################
sub MatchSimpleMols {
  my ($self) = @_;

  my $lv = $self->{lv};
  my $pa = $self->{pa};
  my $pb = $self->{pb};

  my $complex_type = $lv->StringToLabelValue("molecule-type", "complex");
  for my $mol (@{ $pb->Mols() }) {
    if ($pb->MolType($mol) == $complex_type) {
      next;
    }
    my $found = 0;
    my $hit = "";
    my $ehash = $pb->MolExId($mol);
    for my $idtype (keys %{ $ehash }) {
      for my $id (keys %{ $$ehash{$idtype} }) {
        if (defined $self->{exid_index}{$idtype . ":" . lc($id)}) {
          for my $candidate_mol (
              keys %{ $self->{exid_index}{$idtype . ":" . lc($id)} }) {
            if (Pathway::EqSimplex($pa, $pb, $candidate_mol, $mol)) {
              if (! $found) {
                $hit = $candidate_mol;
                $found++;
              } elsif ($hit ne $candidate_mol) {
#print STDERR "MatchSimpleMols: $mol ex hit twice\n";
                $found++;
              }
            }
          }
        }
      }
    }
    if (! $found) {
      my $nhash = $pb->MolName($mol);
      for my $nametype (keys %{ $nhash }) {
        for my $name (keys %{ $$nhash{$nametype} }) {
          if (defined $self->{name_index}{$nametype . ":" . lc($name)}) {
            for my $candidate_mol (
                keys %{ $self->{name_index}{$nametype . ":" . lc($name)} }) {
              if (Pathway::EqSimplex($pa, $pb, $candidate_mol, $mol)) {
                if (! $found) {
                  $hit = $candidate_mol;
                  $found++;
                } elsif ($hit ne $candidate_mol) {
                  $found++;
#print STDERR "MatchSimpleMols: $mol name hit twice\n";
                }
              }
            }
          }
        }
      }
    }
    if ($found == 1) {
      $pb->AddMolMap($mol, $hit);
    }
  }
  
}

######################################################################
sub MatchMols {
  my ($self) = @_;

  my $lv = $self->{lv};
  my $pa = $self->{pa};
  my $pb = $self->{pb};

  ## create a {map} in $pb (which may be the same as $pa)

  $pa->FlattenComplexes();
  $pb->FlattenComplexes();
  $self->MatchSimpleMols();
  $self->MatchComplexMols();
}

######################################################################
sub r_numerically { $b <=> $a };

######################################################################
sub numerically { $a <=> $b };

######################################################################
1;

