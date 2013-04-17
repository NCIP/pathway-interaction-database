#!/usr/local/bin/perl


# Copyright SRA International
#
# Distributed under the OSI-approved BSD 3-Clause License.
# See http://ncip.github.com/pathway-interaction-database/LICENSE.txt for details.


use strict;
use DBI;

open (HELLO, "> /cgap/webcontent/PID/dev/run/hello.txt");
print HELLO "hello\n";
my $db_inst = 'cgprod';
my $db_user = 'pid';
my $db_pass = 'pid24prod';

my $db = DBI->connect("DBI:Oracle:" . $db_inst, $db_user, $db_pass);
# my $e = $db->err();
#print HELLO "dberr:  $e\n";
#if (not $db ) {
  print HELLO "Cannot connect to " . $db_user . "@" . $db_inst . "\n";
#  die;
#}
close(HELLO);

