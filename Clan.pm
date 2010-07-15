package Clan;
require Exporter;
@ISA = qw(Exporter);
#@EXPORT = qw(
#);

use strict;

######################################################################
sub new {
  my ($self, $atomset, $molinstset, $macropset) = @_;

  my $x = {};
  if (defined $atomset && $atomset ne "") {
    for my $a (@{ $atomset }) {
      $x->{atom}{$a} = 1;
    }
  }
  if (defined $molinstset && $molinstset ne "") {
    for my $m (@{ $molinstset }) {
      $x->{molinst}{$m} = 1;
    }
  }
  if (defined $macropset && $macropset ne "") {
    for my $p (@{ $macropset }) {
      $x->{macroprocess}{$p} = 1;
    }
  }
  return bless $x;
}

######################################################################
sub HasAtom {
  my ($self, $atom) = @_;
  if (defined $self->{atom}{$atom}) {
    return 1;
  } else {
    return 0;
  }
}

######################################################################
sub ClanAtoms {
  my ($self) = @_;
  return [ keys %{ $self->{atom} } ];
}

######################################################################
sub ClanSize {
  my ($self) = @_;

  return scalar(keys %{ $self->{atom} });
}

######################################################################
sub ClansTouch {
  my ($self, $x) = @_;

  for my $m (keys %{ $x->{molinst} }) {
    if (defined $self->{molinst}{$m}) {
      return 1;
    }
  }
  for my $p (keys %{ $x->{macroprocess} }) {
    if (defined $self->{macroprocess}{$p}) {
      return 1;
    }
  }
  return 0;
  
}

######################################################################
sub JoinClans {
  my ($self, $x) = @_;

  for my $a (keys %{ $x->{atom} }) {
    $self->{atom}{$a} = 1;
  }
  for my $m (keys %{ $x->{molinst} }) {
    $self->{molinst}{$m} = 1;
  }
  for my $p (keys %{ $x->{macroprocess} }) {
    $self->{macroprocess}{$p} = 1;
  }  
}

######################################################################
sub PartitionClans {
  my ($clanlist) = @_;

  my @clans = @{ $clanlist };
  my $n = @clans;
  my %killed;
  my $more;

  while (1) {
    $more = 0;
    for (my $i = 0; $i < $n - 1; $i++) {
      if ($killed{$i}) {
        next;
      }
      for (my $j = $i+1; $j < $n; $j++) {
        if ($killed{$j}) {
          next;
        }
        if ($clans[$i]->ClansTouch($clans[$j])) {
          $clans[$i]->JoinClans($clans[$j]);
          $more = 1;
          $killed{$j} = 1;
        }
      }
    }
    if (! $more) {
      last;
    }    
  }

  my @newclanlist;
  for (my $i = 0; $i < $n; $i++) {
    if (! $killed{$i}) {
      push @newclanlist, $clans[$i];
#print STDERR "[" . join(",", keys %{ $clans[$i]->{atom} }) . "], ";
#print STDERR "[" . join(",", keys %{ $clans[$i]->{molinst} }) . "]\n";
    }
  }

  return \@newclanlist;
}

