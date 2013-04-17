

# Copyright SRA International
#
# Distributed under the OSI-approved BSD 3-Clause License.
# See http://ncip.github.com/pathway-interaction-database/LICENSE.txt for details.


package IdMap;
require Exporter;
@ISA = qw(Exporter);
#@EXPORT = qw(
#  CreateIsA
#  CreateLabelSort
#);

use strict;

######################################################################
sub new {
  my ($self, $f) = @_;
  my $x = {};
  if (defined $f) {
    ReadIdMap($x, $f);
  }
  return bless $x;
}

######################################################################
sub PTMExprId {
  my ($self) = @_;
  if (defined $self->{start_ptm_id}) {
    return $self->{start_ptm_id};
  } else {
    return 0;
  }
}

######################################################################
sub AtomId {
  my ($self) = @_;
  if (defined $self->{start_atom_id}) {
    return $self->{start_atom_id};
  } else {
    return 0;
  }
}

######################################################################
sub MolId {
  my ($self) = @_;
  if (defined $self->{start_mol_id}) {
    return $self->{start_mol_id};
  } else {
    return 0;
  }
}

######################################################################
sub PathwayId {
  my ($self) = @_;
  if (defined $self->{start_pathway_id}) {
    return $self->{start_pathway_id};
  } else {
    return 0;
  }
}

######################################################################
sub SourceId {
  my ($self, $name) = @_;
  if (defined $self->{r_source_id}{$name}) {
    return $self->{r_source_id}{$name};
  } else {
    return 0;
  }
}

######################################################################
sub SourceName {
  my ($self, $id) = @_;
  if (defined $self->{source_id}{$id}) {
    return $self->{source_id}{$id};
  } else {
    return "";
  }
}

######################################################################
sub ReadIdMap {
  my ($self, $f) = @_;

  open(INF, $f) or die "cannot open $f";
  while (<INF>) {
    s/[\r\n]+//;
    s/^\s+//;
    s/\s+$//;
    if (/^$/ || /^#/) {
      next;
    }
    my  ($op, $id, $other) = split(/\t/, $_, 3);
    my $lcop = lc($op);
    if ($lcop eq "start_pathway_id") {
      $self->{$lcop} = $id;
    } elsif ($lcop eq "start_atom_id") {
      $self->{$lcop} = $id;
    } elsif ($lcop eq "start_mol_id") {
      $self->{$lcop} = $id;
    } elsif ($lcop eq "start_ptm_id") {
      $self->{$lcop} = $id;
    } elsif ($lcop eq "source_id") {
      $self->{source_id}{$id} = $other;
      $self->{r_source_id}{$other} = $id;
    }
  }
  close INF;
}

######################################################################
1;
######################################################################