#!/usr/local/bin/perl

BEGIN {

  my @path_elems = split("/", $0);
  pop @path_elems;
  push @INC, join("/", @path_elems);

}
######################################################################

use strict;
use Bayesian;
use CGI;

my $query = new CGI;
my $a = $query->param("a");
my $b = $query->param("b");
my $A = $query->param("A");
my $B = $query->param("B");
my $fold = $query->param("fold");

#my ($a, $b, $A, $B, $fold) = @ARGV;

print "Content-type: text/plain\n\n";

my $P;
if ($a/$A > $b/$B) {
  $P = sprintf "%.2e", 1 - Bayesian::Bayesian($fold, $a, $b, $A, $B);
}else {
  $P = sprintf "%.2e", 1 - Bayesian::Bayesian($fold, $b, $a, $B, $A);
}
print "a = $a, b  = $b, A = $A, B = $B, fold = $fold\n";

print "p = $P\n";

