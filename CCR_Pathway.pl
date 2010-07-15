#!/usr/local/bin/perl

use strict;
use CGI;

BEGIN {
  my @path_elems = split("/", $0);
  pop @path_elems;
  push @INC, join("/", @path_elems);
}

use CCR_Pathway;

my $query            = new CGI;
my $cmd              = $query->param("cmd");
my $pathway_name     = $query->param("pathn");
my $pathway_rev      = $query->param("pathv");
my $pathway_data     = $query->param("pathf");
my $expression_data  = $query->param("exprf");
my $iid              = $query->param("iid");
my $format           = $query->param("format");

if ($cmd eq "") {
  print "Content-type: text/html\n\n";
  print ConstructForm();
} elsif ($cmd eq "get") {
  print "Content-type: text/plain\n\n";
} elsif ($cmd eq "draw") {
  print "Content-type: text/html\n\n";
  print Draw($pathway_data, $expression_data);
} elsif ($cmd eq "save") {
} elsif ($cmd eq "image") {
  print "Content-type: text/plain\n\n";
  print GetPwImage_1("", $iid);
} else {
  print "Unrecognized cmd $cmd\n";
  exit;
}
