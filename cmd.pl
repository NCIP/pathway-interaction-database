#!/usr/local/bin/perl

use CGI;

my $query = new CGI;

my @cmds = $query->param("cmd");

print "Content-type: text/plain\n\n";

for my $cmd (@cmds) {
  $cmd =~ s/[\r\n]+//;
  print "## $cmd\n\n";
  system($cmd);
}
exit;


#my $cmd = "ls -l /share/content/CMAP/run/biopaths_tmp";
#my $cmd = "cp /local/content/biopaths/cgi-bin/PathwaysNature/* /share/content/CMAP/run/biopaths_tmp";

#my $cmd = "ls -l /local/content/biopaths/PathwaysNature";

#my $cmd = "ls -l /local/content/biopaths/cgi-bin/PathwaysNature";
#my $cmd = "cat /local/content/biopaths/PathwaysNature/files/rhodopsinpathway2.js";

my $cmd = "chown schaefec /share/content/CMAP/run/biopaths_tmp/*";
#print "## $cmd\n\n";

#system $cmd;





