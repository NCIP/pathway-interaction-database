#!/usr/local/bin/perl

use DBI;
use FileHandle;
use strict;

BEGIN {
  my @path_elems = split("/", $0);
  pop @path_elems;
  push @INC, join("/", @path_elems);
}

use PWLabel;
use Parser;

my (
  $what,
  $inpf,
  $molmapf,
#  $pathway_id_base,
  $set_pathway_id,
  $atom_id_base,
  $organism,
  $source_id
) = @ARGV;


if (-d "/app/oracle/product/dbhome/current") {
  $ENV{'ORACLE_HOME'} = "/app/oracle/product/dbhome/current";
} elsif (-d "/app/oracle/product/8.1.7") {
  $ENV{'ORACLE_HOME'} = "/app/oracle/product/8.1.7";
} elsif (-d "/app/oracle/product/8.1.6") {
  $ENV{'ORACLE_HOME'} = "/app/oracle/product/8.1.6";
}

my $db = DBI->connect("DBI:Oracle:" . "lpgdev",
      "web", "readonly");
if (not $db or $db->err()) {
  print STDERR "#! Cannot connect to " . "web" . "@" .
      "lpgdev" . "\n";
  exit;
}

my $lv   = new PWLabel($db, "cgap");
my $pw   = new Pathway($lv);
my $xp   = new Parser($pw, $lv, $atom_id_base, $set_pathway_id);

$db->disconnect();

$xp->ReadMolMap($molmapf);

my $fh = new FileHandle;
open($fh, $inpf) or die "Cannot open input file $inpf";
if (! $xp->ReadFile($fh) ) {
  print join("", @{ $xp->ListErrors() });
  exit();
}
close $fh;
if (! $xp->Parse() ) { 
  print join("", @{ $xp->ListErrors() });
  exit();
}


$pw->BuildMolInstCache();
#print STDERR "return from BuildMolInstCache\n";

$pw->BuildGenericLabelCache();
#print STDERR "return from BuildGenericLabelCache\n";

$pw->PruneDuplicateAtoms();
#print STDERR "return from PruneDuplicateAtoms\n";

$pw->IdentifyMacroProcesses();
#print STDERR "return from IdentifyMacroProcesses\n";

if ($what eq "dot") {
  use DOToutput;
  my @lines;
  my $dot = new DOToutput($pw, $lv, "");
  $dot->ShowAtomIds(1);
  $dot->ShowSubTypeLines(1);
  $dot->DOTGraph();
  print join("\n", @lines, @{ $dot->Lines }) . "\n";
} elsif ($what eq "table") {
  use TableOutput;
  my @lines;
  my $pathway_id_map = $xp->PathwayIdMap();
  my $atom_id_map    = $xp->AtomIdMap();
  my $to = new TableOutput($pw, $lv, $organism, $source_id,
      $pathway_id_map, $atom_id_map, "");
  $to->DoAll();
  print join("\n", @lines, @{ $to->Lines }) . "\n";
} elsif ($what eq "lisp") {
  use LispOutput;
  my @lines;
  my $lisp = new LispOutput($pw, $lv);
  $lisp->PrPathway();
  print join("\n", @lines, @{ $lisp->Lines }) . "\n";
} elsif ($what eq "template") {
  use Templates;
  my @lines;
  my $pwt = new Templates($pw, $lv, "");
  $pwt->CollateAtomTypes($pw);
  $pwt->PrintTemplates($lv);
  return join("\n", @{ $pwt->Lines() }) . "\n";
} else {
  print STDERR "do nothing for what = $what\n";
}




