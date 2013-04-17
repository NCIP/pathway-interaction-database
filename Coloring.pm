

# Copyright SRA International
#
# Distributed under the OSI-approved BSD 3-Clause License.
# See http://ncip.github.com/pathway-interaction-database/LICENSE.txt for details.


package Coloring;
require Exporter;
@ISA = qw(Exporter);
#@EXPORT = qw(
#  CreateIsA
#  CreateLabelSort
#);

use strict;

######################################################################
sub new {
  my ($self, $breaks, $color_scale) = @_;
##  my ($self, $neutral_lo, $neutral_hi,
##      $breaks, $color_scale, $neutral) = @_;

  my $x = {};
#  $x->{neutral_lo}  = $neutral_lo;
#  $x->{neutral_hi}  = $neutral_hi;
  $x->{breaks}      = $breaks;
  $x->{color_scale} = $color_scale;
#  $x->{neutral}     = $neutral;

  return bless $x;
}


######################################################################
sub ColorScale {
  my ($self) = @_;
  return $self->{color_scale};
}

######################################################################
sub NumericBreaks {
  my ($self) = @_;
  return $self->{breaks};
}

######################################################################
sub Value2Color {
  my ($self, $v) = @_;

  my @color_scale = @{ $self->{color_scale} };
  my @breaks      = @{ $self->{breaks} };
#  my $neutral     = $self->{neutral};
#  my $neutral_lo  = $self->{neutral_lo};
#  my $neutral_hi  = $self->{neutral_hi};

#  if (($v <= $neutral_hi) && ($v >= $neutral_lo)) {
#    return $neutral;
#  }  

  for (my $i = 0; $i < @breaks; $i++) {
    if ($v <= $breaks[$i]) {
      return $color_scale[$i];
    }
  }
  return $color_scale[$#color_scale];
}

######################################################################
1;
