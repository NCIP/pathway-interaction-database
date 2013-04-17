

# Copyright SRA International
#
# Distributed under the OSI-approved BSD 3-Clause License.
# See http://ncip.github.com/pathway-interaction-database/LICENSE.txt for details.


package ValueGenerator;
require Exporter;
@ISA = qw(Exporter);
#@EXPORT = qw(
#  CreateIsA
#  CreateLabelSort
#);

use strict;

use constant LO_VAL         => "1";
use constant HI_VAL         => "9";

######################################################################
sub new {
  my ($self, $undefined_atom_val, $undefined_mol_val,
      $atom_vals, $mol_vals, $transform_vector) = @_;

  my $x = {};
##  stored_x = m * y_in + b
##  y_out    = (stored_x - b ) / m
  if (defined $transform_vector) {
    my @tv = @{ $transform_vector };
    my ($y_lo, $y_hi) = ($tv[0], $tv[$#tv]);
    my $m = ($y_hi - $y_lo) / (HI_VAL - LO_VAL);
    my $b = $y_hi - ($m * HI_VAL);
    $x->{transform_m} = $m;
    $x->{transform_b} = $b;
#print STDERR "m = $m, b \ $b\n";
  }
  if (defined $x->{undefined_atomval}) {
    $x->{undefined_atomval} = LinearTransformIn($x, $undefined_atom_val);
  } else {
    $x->{undefined_atomval} = (HI_VAL + LO_VAL) / 2;
  }
  if (defined $x->{undefined_molval}) {
    $x->{undefined_molval}  = LinearTransformIn($x, $undefined_mol_val);
  } else {
    $x->{undefined_molval}  = (HI_VAL + LO_VAL) / 2;
  }
  if (defined $atom_vals) {
    for my $a (keys %{ $atom_vals }) {
      $x->{atomval}{$a} = LinearTransformIn($x, $$atom_vals{$a});
    }
  }
  if (defined $mol_vals) {
    for my $m (keys %{ $mol_vals }) {
      $x->{molval}{$m}{"\t"} = LinearTransformIn($x, $$mol_vals{$m}{"\t"});
    }
  }

  return bless $x;
}


######################################################################
sub LinearTransformIn {
  my ($self, $y) = @_;
  my ($x, $m, $b);
  if (defined $self->{transform_m}) {
    $m = $self->{transform_m};
    $b = $self->{transform_b};
    if ($m) {
      $x = ($y - $b) / $m;
      return $x;
    } else {
      return $b;
    }
  } else {
    return $y;
  }
}

######################################################################
sub LinearTransformOut {
  my ($self, $x) = @_;
  my ($y, $m, $b);
  if (defined $self->{transform_m}) {
    $m = $self->{transform_m};
    $b = $self->{transform_b};
    $y = $m * $x + $b;
    return $y;
  } else {
    return $x;
  }
}

######################################################################
sub SetAtomVal {
  my ($self, $atom, $value) = @_;
  $self->{atomval}{$atom} = $self->LinearTransformIn($value);
}

######################################################################
sub SetMolVal {
  my ($self, $dummy_molinst, $mol_id, $value) = @_;
  $self->{molval}{$mol_id}{"\t"} = $self->LinearTransformIn($value);
#print STDERR "set molval $mol_id to " . $self->{molval}{$mol_id}{"\t"} . "\n";
}

######################################################################
sub AtomVal {
  my ($self, $atom) = @_;

  if (defined $self->{atomval}{$atom}) {
    return $self->LinearTransformOut($self->{atomval}{$atom});
  } else {
    return $self->LinearTransformOut($self->{undefined_atomval});
  }
}

######################################################################
sub MolVal {
  my ($self, $dummy_molinst, $mol_id) = @_;
  if (defined $self->{molval}{$mol_id}{"\t"}) {
    return $self->LinearTransformOut($self->{molval}{$mol_id}{"\t"});
  } else {
    return $self->LinearTransformOut($self->{undefined_molval});
  }
}

######################################################################
sub MolHasValue {
  my ($self, $dummy_molinst, $mol_id) = @_;
  if (defined $self->{molval}{$mol_id}{"\t"}) {
    return 1;
  } else {
    return 0;
  }
}

######################################################################
sub AtomHasValue {
  my ($self, $atom) = @_;
  if (defined $self->{atomval}{$atom}) {
    return 1;
  } else {
    return 0;
  }
}

######################################################################
1;

