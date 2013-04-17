#!/usr/local/bin/perl


# Copyright SRA International
#
# Distributed under the OSI-approved BSD 3-Clause License.
# See http://ncip.github.com/pathway-interaction-database/LICENSE.txt for details.


use strict;
use CGAPConfig;

BEGIN {
  my @path_elems = split("/", $0);
  pop @path_elems;
  push @INC, join("/", @path_elems);
  push @INC, join("/", "/app/oracle/product/10gClient");
  push @INC, join("/", "/app/oracle/product/10gClient/lib");
}

## use lib "/usr/lib64/perl5/site_perl/5.8.5"; 
use lib "/app/oracle/product/10gClient/lib"; 


print "Content-type: text/plain\n\n";

print "7777: I am here\n";

for my $env (keys %ENV) {
  my $cmd = "export $env";
  system($cmd);
  print "$env: $ENV{$env}\n";
}

use DBI;
my $db = DBI->connect("DBI:Oracle:" . DB_INSTANCE, DB_USER, DB_PASS);
print "8888: I am here\n";
if (not $db or $db->err()) {
  print "Cannot connect to " . DB_USER . "@" . DB_INSTANCE . "\n";
}
print "9999: I am here\n";
