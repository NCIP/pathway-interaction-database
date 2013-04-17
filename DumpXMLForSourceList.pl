#!/usr/local/bin/perl


# Copyright SRA International
#
# Distributed under the OSI-approved BSD 3-Clause License.
# See http://ncip.github.com/pathway-interaction-database/LICENSE.txt for details.


BEGIN {
  my @path_elems = split("/", $0);
  pop @path_elems;
  push @INC, join("/", @path_elems);
}

use strict;
use Blocks;
use PWAppConfig;
use PWApp;
use Pathway;
use PathwayDB;
use PathParams;
use PWLabel;
use FileHandle;
use ServerSupport;
use URI::Escape;

my @ALL_EVIDENCE = (
  "NIL",
  "IAE",
  "IC",
  "IDA",
  "IFC",
  "IGI",
  "IMP",
  "IOS",
  "IPI",
  "RCA",
  "RGE",
  "TAS"
);

my @lines;

my ($s) = shift @ARGV;
if ($s ne "pid" && $s ne "pid2") {
  die "bad schema: $s";
} else {
  push @lines, "db_schema	$s";
  PWApp:SetSchema($s);
}

push @lines, "db_inst	cgprod";
push @lines, "db_user	web";
push @lines, "db_pass	readonly";
#push @lines, "db_schema	pid";
push @lines, "print	xml";
push @lines, "pathway_id	*";
for my $ic (@ALL_EVIDENCE) {
  push @lines, "evidence_code	$ic";
}
while (<>) {
  push @lines, "source_id	$_"; 
}
InitializeDatabase();
my $base         = "";
my $graphic_type = "text";
my $parm_string  = join("\n", @lines) . "\n";
my @response = split("\001", PrPath_1($base, $graphic_type, $parm_string));
print $response[2];

######################################################################
sub InitializeDatabase {

  InitializeLV();
  InitializeAtom2Atom();
}
