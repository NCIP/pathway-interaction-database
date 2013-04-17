#!/usr/local/bin/perl


# Copyright SRA International
#
# Distributed under the OSI-approved BSD 3-Clause License.
# See http://ncip.github.com/pathway-interaction-database/LICENSE.txt for details.


use strict;
use DBI;

# open (HELLO, "> /share/content/PID/run/hello.txt");
print "hello\n";
# my $db_inst = 'cgprod';
my $db_inst = 'cgprod.nci.nih.gov';
my $db_user = 'pid';
my $db_pass = 'pid24prod';

my $db = DBI->connect("DBI:Oracle:" . $db_inst, $db_user, $db_pass);
my $rc  = $db->disconnect;

print "rc is $rc\n";
