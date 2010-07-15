#!/usr/local/bin/perl

BEGIN {
  my @path_elems = split("/", $0);
  pop @path_elems;
  push @INC, join("/", @path_elems);
}

use strict;
use Cache;
use CGI;
use PWAppConfig;

my $query  = new CGI;
my $id     = $query->param("GID");
my $format = $query->param("FORMAT");

my (@buf);
my $i = 0;
my $fn = CACHE_ROOT . "/" . PW_CACHE_PREFIX . ".$id";
if (!open(INF, $fn)) {
  print "Content-type: text/plain\n\n";
  print "Failed to open cache file $fn";
  exit;
}

while (read INF, $buf[$i], 16384) {
  $i++;
}
close INF;

if ($format eq "GIF") {
  $format = "JPG";
}

if ($format eq "SVG") {
  print "Content-type: image/svg-xml\n\n";
} elsif ($format eq "JPG") {
  print "Content-type: image/jpg\n\n";
}

print join("", @buf);
