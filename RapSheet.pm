#!/usr/local/bin/perl


# Copyright SRA International
#
# Distributed under the OSI-approved BSD 3-Clause License.
# See http://ncip.github.com/pathway-interaction-database/LICENSE.txt for details.


package RapSheet;
require Exporter;
@ISA = qw(Exporter);
#@EXPORT = qw(
#);

use strict;
use Pathway;
use PWLabel;

my $transcription_type;
my $translocation_type;
my $reaction_type;
my $modification_type;
my $input_type;
my $agent_type;
my $inhibitor_type;
my $output_type;

my $complex_type;

my %error_types = (
  "PURE_INHIBITION"           => 1,
  "MODIFICATION_SPONTANEOUS"  => 1,
  "MODIFICATION_NO_OUTPUT"    => 1,
  "MODIFICATION_IDENTICAL_INPUT_OUTPUT" => 1,
  "REACTION_SPONTANEOUS"      => 1,
  "REACTION_NO_OUTPUT"        => 1,
  "REACTION_NO_INPUT"         => 1,
  "TRANSCRIPTION_NO_OUTPUT"   => 1,
  "TRANSCRIPTION_SPONTANEOUS" => 1,
  "TRANSCRIPTION_HAS_INPUT"   => 1,
  "TRANSCRIPTION_MULTIPLE_OUTPUT"       => 1,
  "TRANSLOCATION_NO_INPUT"    => 1,
  "TRANSLOCATION_NO_OUTPUT"   => 1,
  "TRANSLOCATION_UNSPECIFIED_INPUT"     => 1,
  "TRANSLOCATION_UNSPECIFIED_OUTPUT"    => 1,
  "TRANSLOCATION_NO_CHANGE"   => 1
);

######################################################################
sub new {
  my ($self, $pw, $lv) = @_;

  my $x = {};

  $x->{pw} = $pw;
  $x->{lv} = $lv;

  $transcription_type =
      $lv->StringToLabelValue("process-type", "transcription");
  $translocation_type =
      $lv->StringToLabelValue("process-type", "translocation");
  $reaction_type =
      $lv->StringToLabelValue("process-type", "reaction");
  $modification_type =
      $lv->StringToLabelValue("process-type", "modification");
  $input_type =
      $lv->StringToLabelValue("edge-type", "input");
  $agent_type =
      $lv->StringToLabelValue("edge-type", "agent");
  $inhibitor_type =
      $lv->StringToLabelValue("edge-type", "inhibitor");
  $output_type =
      $lv->StringToLabelValue("edge-type", "output");

  $complex_type = 
      $lv->StringToLabelValue("molecule-type", "complex");

  return bless $x;
}

######################################################################
sub RapErrors {
  my ($self) = @_;
  return $self->{errors};
}

######################################################################
sub AddError {
  my ($self, $atom, $errcode) = @_;
  if (! defined $error_types{$errcode}) {
    print STDERR "use of undefined error code $errcode\n";
  }
  $self->{errors}{$atom}{$errcode} = 1;
}

######################################################################
sub RapAtom {
  my ($self, $atom) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};

  my ($atom_type);
  my ($edge, $edge_type);
  my ($incoming, $outgoing);
  my ($agent, $inhibitor, $input, $output);
  my (%incoming, %outgoing, %input, %agent, %output, %inhibitor);
  my (%location, %activity, %ptm);
  my ($process_condition);
  my ($lvid, $label_name, $label_value);

  $atom_type = $pw->AtomType($atom);
  for $lvid (@{ $pw->AtomCondition($atom) }) {
    $process_condition++
  }

  for $edge (@{ $pw->Edges($atom) }) {
    $edge_type = $pw->EdgeType($atom, $edge);
    if ($edge_type == $agent_type) {
      $agent{$edge} = 1;
      $incoming{$edge} = 1;
    } elsif ($edge_type == $inhibitor_type) {
      $inhibitor{$edge} = 1;
      $incoming{$edge} = 1;
    } elsif ($edge_type == $input_type) {
      $input{$edge} = 1;
      $incoming{$edge} = 1;
    } elsif ($edge_type == $output_type) {
      $output{$edge} = 1;
      $outgoing{$edge} = 1;
    }
    for $lvid (@{ $pw->MolLabel($atom, $edge) }) {
      ($label_name, $label_value) = $lv->LabelValueToString($lvid);
      if ($label_name eq "location") {
        $location{$edge} = $lvid;
      } elsif ($label_name eq "activity-state") {
        $activity{$edge} = $lvid;
      }
    }
    $ptm{$edge} = Pathway::NormalPTMString($pw->EdgePTM($atom, $edge));
  }
  $agent     = scalar(keys %agent);
  $inhibitor = scalar(keys %inhibitor);
  $input     = scalar(keys %input);
  $incoming  = $agent + $inhibitor + $input;
  $output    = scalar(keys %output);
  $outgoing  = $output;

  ##
  ## transcription
  ##
  if ($atom_type == $transcription_type) {
    if ((! $agent) && (! $process_condition)) {
      $self->AddError($atom, "TRANSCRIPTION_SPONTANEOUS");
    } elsif (! $output) {
      $self->AddError($atom, "TRANSCRIPTION_NO_OUTPUT");
    } elsif ($input) {
      $self->AddError($atom, "TRANSCRIPTION_HAS_INPUT");
    }
    if ($output > 1) {
      $self->AddError($atom, "TRANSCRIPTION_MULTIPLE_OUTPUT");
    }
  ##
  ## translocation
  ##
  } elsif ($atom_type == $translocation_type) {
    if (! $input) {
      $self->AddError($atom, "TRANSLOCATION_NO_INPUT");
    } elsif (! $output) {
      $self->AddError($atom, "TRANSLOCATION_NO_OUTPUT");
    }
    for $edge (keys %input) {
      if (! defined $location{$edge}) {
        $self->AddError($atom, "TRANSLOCATION_UNSPECIFIED_INPUT");
      }
    }
    for $edge (keys %output) {
      if (! defined $location{$edge}) {
        $self->AddError($atom, "TRANSLOCATION_UNSPECIFIED_OUTPUT");
      } 
    }
    for $edge (keys %input) {
      for my $edge1 (keys %output) {
        if (defined $location{$edge} && defined $location{$edge1}) {
          if ($location{$edge} == $location{$edge1}) {
            $self->AddError($atom, "TRANSLOCATION_NO_CHANGE");
          }
        }
      }
    }

  ##
  ## modification
  ##
  } elsif ($atom_type == $modification_type) {
    if (! $output) {
      $self->AddError($atom, "MODIFICATION_NO_OUTPUT");
    } elsif ((! $agent) && (! $inhibitor) &&
        (! $process_condition) && (! $input)) {
      $self->AddError($atom, "MODIFICATION_SPONTANEOUS");
    } elsif ((! $agent) && (! $inhibitor) &&
        (! $process_condition) && ($input == 1)) {
      my @ei = keys %input;
      my $ei = $ei[0];
      if ($pw->MolType($pw->EdgeMol($atom, $ei)) ne $complex_type) {
        my $ok = 0;
        for my $eo (keys %output) {
          if ($pw->MolType($pw->EdgeMol($atom, $eo)) eq $complex_type) {
            $ok++;
          }
        }
        if (! $ok) {
          $self->AddError($atom, "MODIFICATION_SPONTANEOUS");
        }
      }
    }
    for my $edge (keys %input) {
      for my $edge1 (keys %output) {
        if ($pw->EdgeMol($atom, $edge) eq $pw->EdgeMol($atom, $edge1) &&
            $activity{$edge} eq $activity{$edge1} &&
            $ptm{$edge} eq $ptm{$edge1} &&
            $location{$edge} eq $location{$edge1}
            ) {
          $self->AddError($atom, "MODIFICATION_IDENTICAL_INPUT_OUTPUT");
        }
      }
    }


  ##
  ## reaction
  ##
  } elsif ($atom_type == $reaction_type) {
    if ((! $agent) && (! $inhibitor) &&
        (! $process_condition) && ($input < 2) && ($output == 1)) {
      $self->AddError($atom, "REACTION_SPONTANEOUS");
    } elsif (! $input) {
      $self->AddError($atom, "REACTION_NO_INPUT");
    } elsif (! $output) {
      $self->AddError($atom, "REACTION_NO_OUTPUT");
    }
  }

  if ($atom_type == $modification_type ||
      $atom_type == $reaction_type ||
      $atom_type == $translocation_type ||
      $atom_type == $transcription_type ) {
    if ($inhibitor == 1 && $agent == 0 && $input == 0 &&
        $process_condition == 0) {
      $self->AddError($atom, "PURE_INHIBITION");
    }    
  }
}

######################################################################
sub RapAll {
  my ($self) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};
  my $fh = $self->{output_fh};

  for my $atom (@{ $pw->Atoms() }) {
    $self->RapAtom($atom);
  }
}

######################################################################
1;
 
