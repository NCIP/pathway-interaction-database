#!/usr/local/bin/perl
package HTMLOutput;
require Exporter;
@ISA = qw(Exporter);
#@EXPORT = qw(
#);

use strict;
use Pathway;
use PWLabel;

my $COLOR1 = "#3366cc";
my $COLOR2 = "#cccccc";

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
sub  numerically { $a <=> $b };

######################################################################
sub ENTREZ_GENE_URL {
  my ($id) = @_;
  return "http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=gene&" .
      "cmd=Retrieve&dopt=full_report&list_uids=$id";
}

######################################################################
sub UNIPROT_URL {
  my ($id) = @_;
  return "http://www.ebi.uniprot.org/entry/$id";
}

######################################################################
sub PUBMED_URL {
  my ($pubmed_id) = @_;
  return "http://www.ncbi.nlm.nih.gov:80/entrez/query.fcgi?".
      "cmd=Retrieve\&db=PubMed\&list_uids=$pubmed_id\&dopt=Abstract";
}

######################################################################
sub GO_URL {
  my ($id) = @_;
  return "http://www.godatabase.org/cgi-bin/amigo/go.cgi?action=query&" .
      "view=query&query=$id&search_constraint=terms";
}

######################################################################
sub CleanString {
  my ($s) = @_;

  while ($s =~ /([\"\'\200-\377])/) {
    my $c = $1;
    my $x = sprintf("&#x%x;", ord($c));
    $s =~ s/$c/$x/g;
  }
  return $s;
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
  my ($self, $line) = @_;

  my $fh = $self->{output_fh};
  if ($fh) {
    print $fh "$line\n";
  }
  push @{ $self->{lines} }, $line;
}

######################################################################
sub ExternalId {
  my ($self, $id_type, $id) = @_;

  if ($id_type eq "LL") {
    return "<a href=\"" . ENTREZ_GENE_URL($id) .
        "\" target=\"entrez\">LL:$id</a>"
  } elsif ($id_type eq "UP") {
    return "<a href=\"" . UNIPROT_URL($id) .
        "\" target=\"uniprot\">UP:$id</a>"
  } elsif ($id_type eq "GO") {
    return "<a href=\"" . GO_URL($id) .
     "\" target=\"go\">GO:$id</a>"
  } else {
    return "$id_type:$id"
  }
}

######################################################################
sub PrMolInst {
  my ($self, $mol, $mixed_labels, $mods) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};

  my ($label, $mol_type) = $lv->LabelValueToString($pw->MolType($mol));
  my $mt = $mol_type;

  if ($mol_type eq "protein") {
    my @whole  = @{ $pw->MolWhole($mol) };
    if (@whole) {
      $mt = "protein subunit [" . $pw->PartBounds($mol) . "]";
    }
  }

  my (@member_names);
  my @members = @{ $pw->FamilyChildren($mol) };
  my $not_a_subunit_family = 0;
  for my $member (@members) {
    if (@{ $pw->MolWhole($member) } == 0) {
      $not_a_subunit_family = 1;
    }
    push @member_names, $pw->PickMolName($member);
  }

  if (@members) {
    if ($mol_type eq "protein") {
      if ($not_a_subunit_family) {
        $mt = "protein family";
      } else {
        $mt = "protein subunit family";
      }
    } else {
      $mt = "$mol_type family";
    }
    $mt .= " [" . join(", ", @member_names) . "]";   
  }

  my $fontcolor;
  if ($mol_type eq "complex" && !@members && @{ $pw->Components($mol) }) {
    $fontcolor = $COLOR2;
  } else {
    $fontcolor = $COLOR1;
  }

  $self->PrLine("<tr>");

  ## name, type, xrefs

  $self->PrLine("<td width=\"275\">");
  $self->PrLine("<table name=molname>");
  $self->PrLine("<tr><td><font color=\"$fontcolor\">" . $pw->PickMolName($mol) .
      "</font></td></tr>");

  $self->PrLine("<tr><td><font color=\"$fontcolor\">" . $mt .
      "</font></td></tr>");

  my $mehash = $pw->MolExId($mol);
  for my $id_type (keys %{ $mehash }) {
    for my $i (keys %{ $$mehash{$id_type} }) {
      $self->PrLine("<tr><td><font color=\"$fontcolor\">" . 
          $self->ExternalId($id_type, $i) . "</font></td></tr>");
    }
  }

  $self->PrLine("</table>");
  $self->PrLine("</td>");

  ## labels

  $self->PrLine("<td width=\"150\">");
  if (@{ $mixed_labels }) {
    $self->PrLine("<table name=labels>");
    for my $value (@{ $mixed_labels }) {
      $self->PrLine("<tr><td><font color=\"$fontcolor\">$value</font></td></tr>");
    }
    $self->PrLine("</table>");
  } else {
    $self->PrLine("&nbsp;");
  }
  $self->PrLine("</td>");

  ## mods

  $self->PrLine("<td width=\"150\">");
  if (@{ $mods }) {
    $self->PrLine("<table name=mods>");
    for my $ptm (@{ $mods }) {
      my ($uniprot, $pos, $aa, $modification) = @{ $ptm };
      my ($label, $value) = $lv->LabelValueToString($modification);
      $self->PrLine("<tr><td><font color=\"$fontcolor\">" .
          $aa . "[" . $pos . "]" . "$value</font></td></tr>");
    }
    $self->PrLine("</table>");
  } else {
    $self->PrLine("&nbsp;");
  }
  $self->PrLine("</td>");

  $self->PrLine("<td width=\"125\">\&nbsp;</td>");

  ## close the row

  if ($mol_type eq "complex") {
    $self->PrLine("</font>");
  }
  $self->PrLine("</tr>");

  if ($mol_type eq "complex") {
    for my $comp (sort numerically @{ $pw->Components($mol) }) {
      my @mixed_labels;
      for my $lvid (@{ $pw->ComponentLabel($mol, $comp) }) {
        my ($label, $value) = $lv->LabelValueToString($lvid);
        push @mixed_labels, $value;
      }
      $self->PrMolInst($pw->ComponentMol($mol, $comp),
          \@mixed_labels,
          $pw->ComponentPTM($mol, $comp));
    }
  }
}
######################################################################
sub PrCondition {
  my ($self, $atom, $condition, $edge_add) = @_;
  my $pw = $self->{pw};
  my $lv = $self->{lv};

  if ($edge_add == 1) {
    $self->InteractionHeader();
  }

  my $edge = $condition;
  my ($label, $value) = $lv->LabelValueToString($condition);

  my $edge_type = "condition"; 
  $self->PrLine("<table name=edge>");
  $self->PrLine("<tr>");
  $self->PrLine("<td valign=center bgcolor=\"$COLOR1\" width=\"10%\">" .
      "<font color=\"white\"><b>");
  $edge_add++;
  $self->PrLine($atom . "." . $edge_add . "<br>$edge_type");
  $self->PrLine("</b></font><td>");

  $self->PrLine("</td>");

  $self->PrLine("<td width=\"90%\">");

  $self->PrLine("<table border=1 name=molinst>");
  $self->PrLine("<tr>");
  $self->PrLine("<td width=\"275\">");
 
  $self->PrLine("<table name=molname>");
  $self->PrLine("<tr><td><font color=\"#3366cccc\">$value</font></td></tr>");
  $self->PrLine("</table>"); 
  $self->PrLine("</td>");
  $self->PrLine("<td width=\"150\">");
  $self->PrLine("&nbsp;");
  $self->PrLine("</td>");
  $self->PrLine("<td width=\"150\">");
  $self->PrLine("&nbsp;");
  $self->PrLine("</td>");
  $self->PrLine("<td width=\"125\">");
  $self->PrLine("&nbsp;");
  $self->PrLine("</td></tr>");
  $self->PrLine("</table>");


  $self->PrLine("</td>");

  $self->PrLine("</tr>");
  $self->PrLine("</table>");

}
######################################################################
sub PrEdge {
  my ($self, $atom, $edge, $edge_add) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};

  if ($edge + $edge_add == 1) {
    $self->InteractionHeader();
  }

  my ($label, $edge_type) = $lv->LabelValueToString($pw->EdgeType($atom, $edge));

  $self->PrLine("<table name=edge>");
  $self->PrLine("<tr>");
  $self->PrLine("<td valign=center bgcolor=\"$COLOR1\" width=\"10%\">" .
      "<font color=\"white\"><b>");
  $edge_add = $edge_add + $edge;
  $self->PrLine($atom . "." . $edge_add . "<br>$edge_type");
  $self->PrLine("</b></font><td>");

  $self->PrLine("</td>");

  $self->PrLine("<td width=\"90%\">");

  my $mol = $pw->EdgeMol($atom, $edge);
  my @mixed_labels;
  for my $lvid (@{ $pw->MolLabel($atom, $edge) }, @{ $pw->EdgeLabel($atom, $edge) }) {
    my ($label, $value) = $lv->LabelValueToString($lvid);
    if ($label ne "edge-type") {
      push @mixed_labels, $value;
    }
  }

  $self->PrLine("<table border=1 name=molinst>");
  $self->PrMolInst($mol, \@mixed_labels, $pw->EdgePTM($atom, $edge));
  $self->PrLine("</table>");

  $self->PrLine("</td>");

  $self->PrLine("</tr>");
  $self->PrLine("</table>");
}

######################################################################
sub InteractionHeader {
  my ($self) = @_;

  $self->PrLine("<table name=edgeheader>");
  $self->PrLine("<tr>");
  $self->PrLine("<td valign=center bgcolor=\"$COLOR1\" width=\"10%\">" .
      "<font color=\"white\"><b>");
  $self->PrLine("\&nbsp;");
  $self->PrLine("</b></font><td>");
  $self->PrLine("</td>");
  $self->PrLine("<td width=\"90%\">");
  $self->PrLine("<table border=1 name=molinstheader>");
  $self->PrLine("<tr>" .
      "<td width=\"275\" align=center><b>Molecule</b></td>" .
      "<td width=\"150\" align=center><b>Location</b></td>" .
      "<td width=\"150\" align=center><b>State</b></td>" .
      "<td width=\"125\" align=center><b>Comments</b></td>" .
      "</tr>");
  $self->PrLine("</table>");
  $self->PrLine("</td>");
  $self->PrLine("</tr>");
  $self->PrLine("</table>");
}

######################################################################
sub PrEvidence {
  my ($self, $atom) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};

  my @temp;
  for my $e (@{ $pw->AtomEvidence($atom) }) {
    push @temp, $e;
  }
  if (@temp == 0) {
    return "<font color=\"$COLOR2\">no evidence code</font>";
  } else {
    return join(", ", sort @temp);
  }

}

######################################################################
sub PrReference {
  my ($self, $atom) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};

  my @temp;
  for my $r (@ { $pw->AtomReferences($atom) }) {
    push @temp, "<a href=\"" . PUBMED_URL($r) . "\" >" . 
        "$r" . "</a>";
  }
  if (@temp == 0) {
    return "<font color=\"$COLOR2\">no reference</font>";
  } else {
    return join(", ", @temp);
  }

}

######################################################################
sub PrAtom {
  my ($self, $atom) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};

  my ($dummy, $atomtype) = $lv->LabelValueToString($pw->AtomType($atom));

  $self->PrLine("<hr width=\"100%\" color=\"$COLOR2\">");

  $self->PrLine("<table border=0><tr><td>");
  $self->PrLine("<a name=\"interaction_$atom\">" .
      "<font color=\"$COLOR1\" size=4>" .
      "<b>Interaction $atom</b></font></a>");
  $self->PrLine("</td><td bgcolor=\"$COLOR2\">&nbsp;");
  $self->PrLine("</td><td>");
  $self->PrLine($atomtype);
  $self->PrLine("</td><td bgcolor=\"$COLOR2\">&nbsp;");
  $self->PrLine("</td><td>");
  $self->PrLine($self->PrEvidence($atom));
  $self->PrLine("</td><td bgcolor=\"$COLOR2\">&nbsp;");
  $self->PrLine("</td><td>");
  $self->PrLine($self->PrReference($atom));
  $self->PrLine("</td></table>");
  
  my $edge_add = 0; 
  for my $lvid (@{ $pw->AtomCondition($atom) }) {
    my ($label, $value) = $lv->LabelValueToString($lvid);
    $self->PrLine("<a name=\"condition_$atom" . "_" . "$lvid\"></a>");
    $self->PrCondition($atom, $lvid, $edge_add);
    $edge_add++;
  }
 
  for my $edge (sort numerically @{ $pw->Edges($atom) }) {
    $self->PrLine("<a name=\"edge_$atom" . "_" . "$edge\"></a>");
    $self->PrEdge($atom, $edge, $edge_add);
  }
}

######################################################################
sub PrHTML {
  my ($self) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};

  $self->PrLine("<html>");
  $self->PrLine("<head>");
  $self->PrLine("<STYLE TYPE=\"text/css\">");
  $self->PrLine("BODY {background-color: #ffffff; " .
      "font-family: Arial, Helvetica, sans-serif; font-size:9pt;}");
  $self->PrLine("</STYLE>");
  $self->PrLine("</head>");
  $self->PrLine("<body>");


  $self->PrLine("<p><font color=\"$COLOR1\" size=4>Go to interaction:</font>");
  $self->PrLine("<table border=0>");
  my $n = 0;
  for my $atom (sort numerically @{ $pw->Atoms() }) {
    if ($n % 20 == 0) {
      if ($n > 0) {
        $self->PrLine("</tr>");
      }
      $self->PrLine("<tr>");
    }
    $self->PrLine("<td bgcolor=\"$COLOR1\">" .
        "<a href=\"#interaction_$atom\"><font color=\"white\"><b>" .
        "$atom</b></font></a></td>");
    $n++;
  }
  $self->PrLine("</tr>");
  $self->PrLine("</table>");

  for my $atom (sort numerically @{ $pw->Atoms() }) {
    $self->PrAtom($atom);
  }

  $self->PrLine("</body>");
  $self->PrLine("<html>");

}


######################################################################
1;
######################################################################
