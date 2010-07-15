#!/usr/local/bin/perl
package EquationOutput;
require Exporter;
@ISA = qw(Exporter);
#@EXPORT = qw(
#);

use strict;
use Pathway;
use PWLabel;

######################################################################
sub new {
  my ($self, $pw, $lv, $outfh) = @_;
  my $x = {};
  $x->{pw} = $pw;
  $x->{lv} = $lv;
  $x->{output_fh} = $outfh;
  return bless $x;
}

######################################################################
my $indent_level;
my $blanks =  "                                        ";
use constant INDENT_INCR => 2;

sub Indent {
  if ($indent_level*INDENT_INCR >= length($blanks)) {
    print STDERR "attempt to indent beyond limit\n";
  } else {
    $indent_level++;
  }
}

sub Exdent {
  if ($indent_level > 0) {
    $indent_level--;
  } else {
    print STDERR "attempt to exdent < 0\n";
  }
}

sub Lines {
  my ($self) = @_;
  return $self->{lines};
}

sub PrLine {
  my ($self, $x) = @_;
  my $fh = $self->{output_fh};
  if ($fh) {
    print $fh substr($blanks, 0, $indent_level*INDENT_INCR);
    print $fh "$x\n";
  }
  push @{ $self->{lines} },
      substr($blanks, 0, $indent_level*INDENT_INCR) . $x;
}

sub PrLineNoIndent {
  my ($self, $x) = @_;
  my $fh = $self->{output_fh};
  if ($fh) {
    print $fh "$x\n";
  }
  push @{ $self->{lines} }, $x;
}

######################################################################
sub numerically { $a <=> $b; }

######################################################################
sub MOLPAGE_URL {
  my ($molid) = @_;
#  return "http://pid.nci.nih.gov/search/MoleculePage?molid=$molid";
  return "/search/MoleculePage?molid=$molid";
}

######################################################################
sub ATOMPAGE_URL {
  my ($atomid) = @_;
#  return "http://pid.nci.nih.gov/search/InteractionPage?atomid=$atomid";
  return "/search/InteractionPage?atomid=$atomid";
}

######################################################################
sub EquationTerm {
  my ($self, $type, $mol, $labels, $ptms) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};

  my ($mol_name, $color, $condition, $macroprocess, @exids, $pathway);
  my (@locations, @modifications, @activity_state);
  my ($dummy);
  my (@lines);

  if ($type eq "condition") {
    $color = "green";
    for my $lvid (@{ $labels }) {
      ### only one label
      ($dummy, $condition) = $lv->LabelValueToString($lvid);
    }
  } elsif ($type eq "macroprocess") {
    $color = "black";
    for my $lvid (@{ $labels }) {
      ### only one label
      ($dummy, $macroprocess) = $lv->LabelValueToString($lvid);
    }
  } elsif ($type eq "pathway") {
    $color = "black";
    #### bad: this is NOT a label value id; it's an atom id
    for my $lvid (@{ $labels }) {
      $pathway = $pw->AbstractionName($lvid);
    }
  } else {
    if ($type eq "agent") {
      $color = "green";
    } elsif ($type eq "inhibitor") {
      $color = "red";
    } else {
      $color = "black";
    }
    $mol_name = $pw->PickMolName($mol);
    my $mehash = $pw->MolExId($mol);
    for my $idtype (keys %{ $mehash }) {
      for my $id (keys %{ $$mehash{$idtype} }) {
        push @exids, $idtype . ":" . $id;
      }
    }
    for my $lvid (@{ $labels }) {
      my ($label, $value) = $lv->LabelValueToString($lvid);
      if ($label eq "location") {
        push @locations, $value;
      } elsif ($label eq "activity-state") {
        push @activity_state, $value;
      }
    }
    for my $ptm (@{ $ptms }) {
      my ($protein_id, $position, $amino_acid, $label_value, $label_name)
          = @{ $ptm };
      push @modifications, "[ $amino_acid $position $label_name ]";
    }
  }

  push @lines, "<table border=0>";

  if ($type eq "condition") {
    push @lines, "<tr><td><font color=\"$color\">$condition</font></td></tr>";
    push @lines, "<tr><td><font color=\"green\">" .
        "<i>condition</i></font></td></tr>";
  } elsif ($type eq "macroprocess") {
    push @lines, "<tr><td><font color=\"$color\">$macroprocess</font></td></tr>";
  } elsif ($type eq "pathway") {
    push @lines, "<tr><td><font color=\"$color\">$pathway</font></td></tr>";
  } else {
    push @lines, "<tr><td><font color=\"$color\">" .
        "<a href=\"" . MOLPAGE_URL($mol) . "\">" .
        "$mol_name</a></font></td></tr>";
    if ($type eq "input" || $type eq "output") {
      push @lines, "<tr><td><font color=\"black\">" .
          "<i>$type</i></font></td></tr>";
    } elsif ($type eq "agent") {
      push @lines, "<tr><td><font color=\"green\">" .
          "<i>$type</i></font></td></tr>";
    } elsif ($type eq "inhibitor") {
      push @lines, "<tr><td><font color=\"red\">" .
          "<i>$type</i></font></td></tr>";
    }
  }
  for my $x (@exids, @locations, @activity_state, @modifications) {
    push @lines, "<tr><td>$x</td></tr>";
  }
  push @lines, "</table>";
  return join("\n", @lines);
}

######################################################################
sub Equation {
  my ($self, $atom) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};

  my (@conditions, @inhibitors, @agents, @inputs, @outputs);
  my (@row);

  push @row, "<table border=1>";
  push @row, "<tr valign=center>";

push @row, "";
  push @row, "<td>";
  push @row, "<table border=0>";
  push @row, "<tr><td><a href=\"" . ATOMPAGE_URL($atom) .
      "\">$atom</a></td></tr>";
  push @row, "</table>";
  push @row, "<td>";
push @row, "";

  for my $condition (@{ $pw->AtomCondition($atom) }) {
    my ($dummy, $value) = $lv->LabelValueToString($condition);
    push @conditions, [ "condition", undef, [ $value ], [ ] ];
  }
  for my $edge (@{ $pw->Edges($atom) }) {
    my ($dummy, $edge_type) =
      $lv->LabelValueToString($pw->EdgeType($atom, $edge));
    my $mol = $pw->EdgeMol($atom, $edge);
    my $mollabels = $pw->MolLabel($atom, $edge);
    my $ptms = $pw->EdgePTM($atom, $edge);
    if ($edge_type eq "agent") {
      push @agents, [ $edge_type, $mol, $mollabels, $ptms ];
    } elsif ($edge_type eq "inhibitor") {
      push @inhibitors, [ $edge_type, $mol, $mollabels, $ptms ];
    } elsif ($edge_type eq "input") {
      push @inputs, [ $edge_type, $mol, $mollabels, $ptms ];
    } elsif ($edge_type eq "output") {
      push @outputs, [ $edge_type, $mol, $mollabels, $ptms ];
    }
  }

  my $incoming = @conditions + @inhibitors + @agents + @inputs;
  my $outgoing = @outputs;

  for my $edge_data (@conditions, @inhibitors, @agents, @inputs) {
    push @row, "<td>";
    push @row, $self->EquationTerm(@{ $edge_data });
    push @row, "</td>";
push @row, "";
    $incoming--;
    if ($incoming) {
      push @row, "<td>";
      push @row, $self->InverseConnector("+");
      push @row, "</td>";
push @row, "";
    }
  }
  push @row, "<td>";
  push @row, $self->InverseConnector("->");
      push @row, "</td>";
push @row, "";
  if ($outgoing == 0) {
    my $atomtype = $pw->AtomType($atom);
    push @row, "<td>";
    if ($atomtype eq $lv->StringToLabelValue("process-type", "pathway")) {
      push @row, $self->EquationTerm("pathway", undef, [ $atom ], []);
    } else {
      push @row, $self->EquationTerm("macroprocess", undef, [ $atomtype ], []);
    }
    push @row, "</td>";
push @row, "";
  } else {
    for my $edge_data (@outputs) {
      push @row, "<td>";
      push @row, $self->EquationTerm(@{ $edge_data });
      push @row, "</td>";
push @row, "";
      $outgoing--;
      if ($outgoing) {
        push @row, "<td>";
        push @row, $self->InverseConnector("+");
        push @row, "</td>";
push @row, "";
      }
    }
  }

  push @row, "</tr>";
  push @row, "</table>";
  $self->PrLine(join("\n", @row));
}

######################################################################
sub Connector {
  my ($self, $connector) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};

  my @lines;
  push @lines, "<table border=0>";
  push @lines, "<tr><td>$connector</td><tr>";
  push @lines, "</table>";
  return join("\n", @lines);
}

######################################################################
sub InverseConnector {
  my ($self, $connector) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};

  my @lines;
  push @lines, "<table border=0>";
  push @lines, "<tr><td bgcolor=\"black\"><font color=\"white\">" .
      "$connector</font></td><tr>";
  push @lines, "</table>";
  return join("\n", @lines);
}

######################################################################
sub EquationOutput {
  my ($self) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};

  for my $atom (sort numerically @{ $pw->Atoms() }) {
    $self->PrLine("<hr>");
    $self->Equation($atom);
  }
  $self->PrLine("<hr>");
}

######################################################################
1;
######################################################################
