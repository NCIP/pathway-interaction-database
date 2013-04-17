#!/usr/local/bin/perl


# Copyright SRA International
#
# Distributed under the OSI-approved BSD 3-Clause License.
# See http://ncip.github.com/pathway-interaction-database/LICENSE.txt for details.


my $src       = "/local/content/biopaths/PathwaysNature/files";
my $dest      = "/share/content/PID/data/transfer";
my $last_dest = "/share/content/PID/data/last_transfer";

print "Content-type: text/plain\n\n";

#system("ls -l $src/*.js");
#system("ls -l $dest");

system("cp -p $dest/*.js $last_dest");
system("cp -p $dest/Stats $last_dest");
system("rm $dest/*.js");
system("rm $dest/Stats");
system("cp -p $src/*.js $dest");
system("cp -p $src/Stats $dest");

print "Done.\n";
