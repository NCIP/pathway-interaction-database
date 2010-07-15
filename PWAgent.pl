#!/usr/local/bin/perl

use strict;
use LWP::Simple;
use URI::Escape;

my ($what, $f) = @ARGV;

my $url = "http://cbiodev104.nci.nih.gov:8080/PW/DrawPathwayForAgent";

my @temp;
open(INP, $f) or die "cannot open $f";
while (<INP>) {
  chop;
  push @temp, $_;
}
close INP;

my $vars = "graphic_type=$what" .
    "&" .
    "params=" . uri_escape(join("\n", @temp) . "\n") ;

print get("$url?$vars");

 
