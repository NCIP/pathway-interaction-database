#!/usr/local/bin/perl


# Copyright SRA International
#
# Distributed under the OSI-approved BSD 3-Clause License.
# See http://ncip.github.com/pathway-interaction-database/LICENSE.txt for details.


use CGI;
use strict;

BEGIN {
  my @path_elems = split("/", $0);
  pop @path_elems;
  push @INC, join("/", @path_elems);
}

if (-d "/app/oracle/product/dbhome/current") {
  $ENV{'ORACLE_HOME'} = "/app/oracle/product/dbhome/current";
}

use PathwayOld;

my ($db_user, $db_pass, $db_inst, $schema) =
    ("web", "readonly", "cgprod", "pid");

my $query     = new CGI;
my $gene      = $query->param("gene");
my $genes_a   = $query->param("genes_a");
my $genes_b   = $query->param("genes_b");
my $format    = $query->param("format");
my $pathway   = $query->param("pathway");
my $what      = $query->param("what");
my $source    = $query->param("source");
my $file_a    = $query->param("FILE_A");
my $file_b    = $query->param("FILE_B");

#my ($what, $genes, $format, $pathway, $source) = ("pathway", "207", "gif", 500055, "Reactome");
#my ($what, $genes_a, $format, $pathway, $source) = ("rank", "207,2208,2885,3122,3123,3159,3497,3561,3565,3566,3716,3718,6198,6464,6778", "text");
#my ($what, $genes_a, $format, $pathway, $source) = ("rank", lc("AKT1,BRCA1,BARD1"), "text");

my $db = DBI->connect("DBI:Oracle:" . $db_inst, $db_user, $db_pass);
if (not $db or $db->err()) {
  print STDERR "Cannot connect to " . $db_user . "@" . $db_inst . "\n";
  exit;
}

if ($what eq "rank") {
  if ($format eq "text") {
    print "Content-type: text/plain\n\n";
  } elsif ($format eq "html") {
    print "Content-type: text/html\n\n";
  } else {
    print "Content-type: text/plain\n\n";
    print "for what='rank', format is 'text' or 'html'\n";
    exit;
  } 
} elsif ($what eq "pathway") {
  if ($pathway !~ /^\d+$/) {
    print "Content-type: text/plain\n\n";
    print "for what='pathway', specify pathway=<numeric pathway id>\n";
    exit;
  }
  if ($format eq "svg") {
    print "Content-type: text/html\n\n";
  } elsif ($format eq "gif") {
    print "Content-type: text/html\n\n";
  } else {
    print "Content-type: text/plain\n\n";
    print "for what='pathway', format is 'svg' or 'gif'\n";
    exit;
  }
  if ($source eq "") {
    print "Content-type: text/plain\n\n";
    print "for what='pathway', specify source=<source name>\n";
    exit;
  }
} else {
  print "Content-type: text/html\n\n";
  my @lines;
  push @lines, "<html>";
  push @lines, HTML_HEAD();
  push @lines, "</head>";
  push @lines, "<body>";
  push @lines, "<h3>Rank Pathways By Gene Lists</h3>";
  push @lines, "<blockquote>";
  my $form = new CGI;
  push @lines, "<br>";
  push @lines, "<p>Upload one or two lists of Entrez Gene ids:";
  push @lines, "<br>";
  push @lines, $form->start_multipart_form(-method=>'POST',-action=>'PathwayGene.pl');
  push @lines, $form->hidden(-name=>'what',-value=>'rank');
  push @lines, "<blockquote>";
  push @lines, "Gene List A\&nbsp;";
  push @lines, $form->filefield(-name=>'FILE_A',-size=>30);
  push @lines, "<br><br>";
  push @lines, "Gene List B\&nbsp;";
  push @lines, $form->filefield(-name=>'FILE_B',-size=>30);
  push @lines, "<br>";
  push @lines, "<ul>";
  push @lines, "<li>Genes only in List A will be colored blue";
  push @lines, "<li>Genes only in List B will be colored red";
  push @lines, "<li>Genes in both lists (and complexes with components in both lists ) will be colored purple";
  push @lines, "</ul>";
  push @lines, "</blockquote>";
  push @lines, "Pick format:";
  push @lines, "<blockquote>";
  push @lines, $form->radio_group(-name=>'format', -linebreak=>'true',
      -values=>['html','text'], -labels=>{'html'=>'HTML','text'=>'text'},
      -default=>'html');
  push @lines, "</blockquote>";
  push @lines, $form->submit(-label=>'Rank Pathways');
  push @lines, "<br>";
  push @lines, $form->end_form;
  push @lines, "<br>";

  push @lines, "</blockquote>";
  print join("\n", @lines, "");
  exit;
}

my (%genes_a, %genes_b);
for my $g (split(",", $genes_a)) {
  $genes_a{$g} = 1;
}
if ($gene ne "") {
  $genes_a{$gene} = 1;
}
for my $g (split(",", $genes_b)) {
  $genes_b{$g} = 1;
}

my $g;
my $error;
if ($file_a) {
  while (<$file_a>) {
    s/\n//g;
    s/\r//g;
    s/^ +//;
    ($g) = split /\t/;
    $g = lc($g);
    if ($g eq "") {
      next;
    }
#    if ($g =~ /^\d+$/) {
      $genes_a{$g} = 1;
#    } else {
#      $error++;
#      print "File A has non-numeric (bad) Entrez Gene id: $g<br>\n";
#    }
  }
}
if ($file_b) {
  while (<$file_b>) {
    s/\n//g;
    s/\r//g;
    s/^ +//;
    ($g) = split /\t/;
    $g = lc($g);
    if ($g eq "") {
      next;
    }
#    if ($g =~ /^\d+$/) {
      $genes_b{$g} = 1;
#    } else {
#      $error++;
#      print "File B has non-numeric (bad) Entrez Gene id: $g<br>\n";
#    }
  }
}
if ($error) {
  exit;
}

if ($what eq "rank") {

  if ($format eq "html") {
    print HTML_HEAD();
    print "<h3>Pathways Ranked by Probability of Containing Genes from Query Lists</h3>";
  }
  print RankPathwaysByGeneHits($db, $schema, $format, \%genes_a, \%genes_b);

} elsif ($what eq "pathway") {

  my $mols_a = GetMolIdsForEntrezGene($db, $pathway, \%genes_a);
  my $mols_b = GetMolIdsForEntrezGene($db, $pathway, \%genes_b);
  print CreateGraphicFile($db, $schema, $pathway, $source, $mols_a, $mols_b, $format);
}

$db->disconnect();
