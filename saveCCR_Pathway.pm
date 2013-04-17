#!/usr/local/bin/perl


# Copyright SRA International
#
# Distributed under the OSI-approved BSD 3-Clause License.
# See http://ncip.github.com/pathway-interaction-database/LICENSE.txt for details.


use strict;
use Cache;

my $BASE;
use constant CACHE_ROOT       => "/share/content/CMAP/data/cache";
use constant PW_CACHE_PREFIX  => "PW";
use constant DOT_PGM          => "/share/content/CMAP/run/PW/dot";
my $GROUP_BY_LOCATION = 0;

my ($svg_height, $svg_width) = (700, 900);

my %element_types = (
  "protein"  => 1,
  "process"  => 1,
  "molecule" => 1,
  "complex"  => 1,
  "family"   => 1
);

my %attr_names = (
  "label"    => 1,
  "locusid"  => 1,
  "type"     => 1,
  "members"  => 1,
  "location" => 1,
  "color"    => 1
);

my %style;

$style{graph}{compound}{default} = "true";      ## allow edges to clusters
$style{graph}{color}{default}    = "lightgray"; ## fillcolor for clusters
$style{graph}{style}{default}    = "filled";    ## fill for clusters
$style{graph}{bgcolor}{default}  = "black";     ## fill for clusters

$style{node}{width}{default}     = "0";
$style{node}{style}{default}     = "filled";
$style{node}{fontcolor}{default} = "white";
$style{node}{fontname}{default}  = "Arial";
$style{node}{fillcolor}{default} = "lightgray";
$style{node}{fontsize}{default}  = "14";
$style{node}{shape}{default}     = "ellipse";
$style{node}{shape}{process}     = "box";

$style{edge}{color}{default}     = "white";
$style{edge}{arrowsize}{default} = "1.5";

my (@elements, %attr_values, %members, @edges, %rank, $pathway_name);
my (@elem_order);
my (@lines);

######################################################################
sub LL_URL {
  my ($locusid) = @_;

  return "http://cgap.nci.nih.gov/Genes/GeneInfo?ORG=Hs&amp;LLNO=$locusid";

}

######################################################################
sub HTML_Head {

  return
    qq!
<head>
<TITLE>CCR Pathway</TITLE>
<STYLE TYPE="text/css">
  BODY {background-color: #ffffff; font-family: Arial, Helvetica, sans-serif; font-size:10pt;}
  UL,OL,TH,TD,P,DD,DT,DL,BLOCKQOUTE,H1,H2,H3,H4,H5,H6
      {font-family: Arial, Helvetica, sans-serif; font-size:10pt; color:#336699}
  H1,H2{font-size:12pt}
</STYLE>
<script>
function spawn(url) {var w = window.open(url, "_blank")}
</script>
</head>
  !;
}

######################################################################
sub Draw {
  my ($pathf, $exprf) = @_;

  my $graphic_type = "svg";

  if ($pathf) {
    ParsePathway($pathf);
  }
  if ($exprf) {
    ParseExpression($exprf);
  }

  my $cache = new Cache(CACHE_ROOT, PW_CACHE_PREFIX);
  my ($cache_id, $filename) = $cache->MakeCacheFile();
  if ($cache_id != $CACHE_FAIL) {
    my $cmd = DOT_PGM .  " -T$graphic_type > $filename";
    open(OUTF, "| $cmd");
    WriteOutput();
    print OUTF join("\n", @lines) . "\n";
    close(OUTF);
    chmod 0666, $filename;
    return
      HTML_Head() .
      qq!
<body>
<h2>
$pathway_name
</h2>
<embed type="image/svg-xml" height=$svg_height width=$svg_width
src=CCR_Pathway.pl?cmd=image&iid=$cache_id&format=SVG>
</body>
      !;
  } else {
    return "cache failure\n";
  }


}

######################################################################
sub Color {
  my ($elem) = @_;

  if (defined $attr_values{$elem}{color}) {
    return $attr_values{$elem}{color};
  } else {
    return $style{fillcolor}{default};
  }

}

######################################################################
sub Line {
  my ($line) = @_;

  push @lines, $line;
}

######################################################################
sub Trim {
  my ($x) = @_;

  $x =~ s/"//g;
  $x =~ s/\s+/ /g;
  $x =~ s/^\s+//;
  $x =~ s/\s+$//;
  return $x;
}

######################################################################
sub Normalize {
  my ($x) = @_;

  return lc(Trim($x));
}

######################################################################
sub ParsePathway {
  my ($pathf, $exprf) = @_;

  my $state;

  while (<$pathf>) {
    s/\r//;
    s/\n//;
    my $line = Trim($_);

    if ($line =~ /^#/) {
      next;
    }
    if ($line =~ /^\s*$/) {
      next;
    }
    if ($line =~ /^pathway_name\s*=\s*(.+)/i) {
      $pathway_name = Trim($1);
      next;
    }
    if ($line =~ /^pathway_elements {/i) {
      $state = "pathway_elements";
      next;
    }
    if ($line =~ /^pathway {/i) {
      $state = "pathway";
      next;
    }
    if ($line =~ /^\}/) {
      if ($state ne "pathway_elements" &&
          $state ne "pathway") {
        print STDERR "syntax error: }\n";
      } else {
        undef $state;
      }
      next;
    }
    if ($state eq "pathway_elements") {
      if ($line =~ /\}/) {
        undef $state;
        next;
      }
      if ($line =~ /^(\S+)\s*\[([^\]]*)\]$/) {
        my ($elem, $attrs) = ($1, $2);
        push @elements, $elem;
        for my $a (split(/\s*;\s*/, $attrs)) {
          my ($name, $value) = split("=", $a);
          $name = Normalize($name);
          if (defined $attr_names{$name}) {
            if ($name ne "label") {
              $value = Normalize($value);
            }
            if ($name eq "type") {
              if (! defined $element_types{$value}) {
                print STDERR "syntax error: illegal type $value\n";
                next;
              }
            }
            if ($name eq "members") {
              for my $m (split(/\s*,\s*/, $value)) {
                push @{ $members{$elem} }, Trim($m);
              }
            } else {
              $attr_values{$elem}{$name} = $value;
            }
          } else {
            print STDERR "syntax error: illegal attribute name $name\n";
          }
        }
      } else {
        print STDERR "syntax error: $_\n";
      }
      next;
    }
    if ($state eq "pathway") {
      if ($line =~ /\}/) {
        undef $state;
        next;
      }
      my ($left, $op, $right) = split(/\s+/, $line);
      if ($op ne "->" && $op ne "-|") {
        print STDERR "syntax error: illegal operator $op\n";
      } else {
        push @edges, "$left\t$op\t$right";
      }
      next;
    }
  }
}

######################################################################
sub WriteOneComplex {
  my ($cx) = @_;

  my @members = (@{ $members{$cx} });
  my ($i);

  Line("subgraph cluster_$cx {");
  Line("style=filled;");
  Line("nodesep=0;");
  Line("edge [style=\"invis\", arrowhead=\"none\"];");
  for my $m (@members) {
    WriteOneElement($m);
  }
  for ($i = 0; $i < @members - 1; $i++) {
    Line("$members[$i] -> $members[$i+1];");
  }
  Line("};");
}

######################################################################
sub WriteOneElement {
  my ($elem) = @_;

  my ($type, @attrs);


  if (defined $attr_values{$elem}{type}) {
    $type = $attr_values{$elem}{type};
  } else {
    print STDERR "semantic error: no type for $elem\n";
#    return;
  }
  if (($type eq "complex" || $type eq "family")
      && ! defined $members{$elem}) {
    print STDERR "semantic error: $type $elem must have members\n";
    return;
  }
  if (($type ne "complex" && $type ne "family")
      && defined $members{$elem}) {
    print STDERR "semantic error: $type $elem cannot have members\n";
    return;
  }
  if (($type eq "complex" || $type eq "family")
      && defined $attr_values{locusid}) {
    print STDERR "semantic error: $type $elem cannot have locusid\n";
    return;
  }
  if ($type eq "complex") {
    WriteOneComplex($elem);
  } else {
    for my $a ("shape") {
      if (defined $style{node}{$a}{$type}) {
        push @attrs, "$a=\"$style{node}{$a}{$type}\"";
      }
    }
    push @attrs, "fillcolor=\"" . Color($elem) . "\"";
    if (defined $attr_values{$elem}{label}) {
      push @attrs, "label=\"$attr_values{$elem}{label}\"";
    }
    if (defined $attr_values{$elem}{locusid}) {
      push @attrs, "URL=\"javascript:spawn('" . LL_URL($attr_values{$elem}{locusid}) .
          "')\"";
    }
    Line("$elem [" . join(", ", @attrs) . "];");
  }
}

######################################################################
sub WriteElements {

  my (%rank, %elem_seen);

  for my $elem (keys %attr_values) {
    for my $name (keys %{ $attr_values{$elem} }) {
      if ($name eq "type") {
        if ($attr_values{$elem}{$name} eq "process") {
          push @{ $rank{process} }, $elem;
        }
      } elsif ($name eq "location") {
        if ($GROUP_BY_LOCATION) {
          push @{ $rank{$attr_values{$elem}{location}} }, $elem;
        }
      }
    }
  }

  for my $r (keys %rank) {
    if ($r eq "process") {
      Line("{ rank=same; ");
      for my $elem (@{ $rank{$r} }) {
        WriteOneElement($elem);
        $elem_seen{$elem} = 1;
      }
      Line("}");
    } else {
      Line("subgraph cluster_$r {");
      Line("style=\"\";");
      for my $elem (@{ $rank{$r} }) {
        WriteOneElement($elem);
        $elem_seen{$elem} = 1;
      }
      Line("}");
    }
  }

  for my $elem (@elements) {
    if (! defined $elem_seen{$elem}) {
      WriteOneElement($elem);
    }
  }

}

######################################################################
sub WriteOneEdge {
  my ($edge) = @_;

  my (@attrs);

  my ($left, $op, $right) = split(/\t/, $edge);
  if ($op eq "-|") {
    push @attrs, "arrowhead=\"tee\"";
  }
  if ($attr_values{$left}{type} eq "complex") {
    push @attrs, "ltail=\"cluster_$left\"";
    my $n = @{ $members{$left} } - 1;
    $left = $members{$left}[$n];
  }
  if ($attr_values{$right}{type} eq "complex") {
    push @attrs, "lhead=\"cluster_$right\"";
    my $n = @{ $members{$right} } - 1;
    $right = $members{$right}[$n];
  }
  Line("$left -> $right [" . join(", ", @attrs) . "];");
  
}

######################################################################
sub WriteEdges {
  for my $edge (@edges) {
    WriteOneEdge($edge);
  }
}

######################################################################
sub WriteOutput {

  Line("digraph G {");

  for my $x (keys %{ $style{graph} }) {
    if (defined $style{graph}{$x}{default}) {
      Line("$x=\"$style{graph}{$x}{default}\";");
    }
  }

  for my $y ("node", "edge") {
    for my $x (keys %{ $style{$y} }) {
      if (defined $style{$y}{$x}{default}) {
        Line("$y [$x=\"$style{$y}{$x}{default}\"];");
      }
    }
  }

  WriteElements();
  WriteEdges();
  Line("}");
}

######################################################################
sub GetPwImage_1 {
  my ($base, $id) = @_;

  $BASE = $base;

  my (@buf);
  my $i = 0;
  my $fn = CACHE_ROOT . "/" . PW_CACHE_PREFIX . ".$id";
  open(INF, $fn) or return "Failed to open cache file $fn";
  while (read INF, $buf[$i], 16384) {
    $i++;
  }
  close INF;
  return join("", @buf);
  
}

######################################################################
sub ConstructForm {

  return
    HTML_Head() .
    qq!
<body>
<h3>CCR Pathway Tool</h3>
Choose one of the following options:
<ul>
<li>Submit a pathway to be saved.
<li>Submit a pathway (with or without added expression data) to be drawn.
<li>Draw a saved pathway (with or without added expression data).
</ul>

<hr>
<h4>Submit a pathway to be saved</h4>
Not implemented.

<hr>
<h4>Submit a pathway to be drawn</h4>

<form enctype="multipart/form-data"
    action="CCR_Pathway.pl"
    method=POST>

<table border=1 cellspacing=1 cellpadding=4>
<tr>
  <td>Pathway file:</td><td><input type=file name=pathf></td><td>[required]</td>
</tr>
<tr>
  <td>Expression file:</td><td><input type=file name=exprf></td><td>[optional]</td>
</tr>
<tr>
  <td colspan=3 align=center><input type=submit name=cmd value=draw></td>
</tr>
</table>
</form>

<hr>
<h4>Draw a saved pathway</h4>
Not implemented.

<hr>

</body>
    !;
}

######################################################################
1;