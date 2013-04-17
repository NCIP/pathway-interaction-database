#!/usr/local/bin/perl


# Copyright SRA International
#
# Distributed under the OSI-approved BSD 3-Clause License.
# See http://ncip.github.com/pathway-interaction-database/LICENSE.txt for details.


use strict;
use CGI;

BEGIN {
  my @path_elems = split("/", $0);
  pop @path_elems;
  push @INC, join("/", @path_elems);
}

my $src       = "/share/content/PID/data/biopaths/PathwaysNature/files";
my $dest      = "/share/content/PID/data/tmp";
my $rundir    = "/share/content/PID/run";

my $CLEAVED = "zCleaved_Definitions3.js";
my $FAMILY  = "zFamily_Definitions3.js";

my $query  = new CGI;

my $js_f   = $query->param("js");

my $pex    = $js_f;
$pex       =~ s/\d*\.js$//;
$pex       = lc($pex);

my $err_f       = "$pex." . time() . ".err";

my $rule_f      = "$pex.rule";
my $base_dot_f  = "$pex.dot";
my $dot_f       = "$base_dot_f.1";
my $svg_f       = "$pex.svg";
my $gif_f       = "$pex.gif";
my $map_f       = "$pex.gif.html";
my $html_f      = "$pex.html";
my $xml_f       = "$pex.xml";
my $tar_f       = "$pex.tar";
#my $zip_f       = "$pex.zip";
my $zip_f       = "$pex.tar";

#print "Content-type: text/plain\n\n";
print "Content-type: application/octet-stream\n\n";
#print "Content-type: application/zip\n\n";

system("cp $src/$js_f $dest/$js_f > $dest/$err_f 2>&1");
system('echo "" >> $dest/$js_f');
system("cat $src/$FAMILY  >> $dest/$js_f");
system('echo "" >> $dest/$js_f');
system("cat $src/$CLEAVED >> $dest/$js_f");
system('echo "" >> $dest/$js_f');

my $cmd = join(" ",
  "$rundir/Driver.pl",
  "-it js",
  "-if $dest/$js_f",
  "-ot xml",
  "-of $dest/$xml_f",
#  "-db_user web",
#  "-db_pass readonly",
#  "-db_inst cgprod",
#  "-schema pid",
  "-ontf $rundir/pid_ontology.xml",
  "-standalone",
  ">> $dest/$err_f 2>&1"
);

system("$cmd");

open(INF, "$dest/$xml_f");
open(OUT, ">$dest/$xml_f.1");
while (<INF>) {
  s/[\r\n]+//;
  print OUT "$_\r\n";
}
close OUT;
close INF;
system("mv $dest/$xml_f.1 $dest/$xml_f");

my $cmd = join(" ",
  "$rundir/Driver.pl",
  "-it xml",
  "-if $dest/$xml_f",
  "-ot rule",
  "-of $dest/$rule_f",
#  "-db_user web",
#  "-db_pass readonly",
#  "-db_inst cgprod",
#  "-schema pid",
  "-standalone",
  ">> $dest/$err_f 2>&1"
);

system("$cmd");

open(INF, "$dest/$rule_f");
open(OUT, ">$dest/$rule_f.1");
while (<INF>) {
  s/[\r\n]+//;
  print OUT "$_\r\n";
}
close OUT;
close INF;
system("mv $dest/$rule_f.1 $dest/$rule_f");

my $cmd = join(" ",
  "$rundir/Driver.pl",
  "-it xml",
  "-if $dest/$xml_f",
  "-ot dot",
  "-of $dest/$base_dot_f",
#  "-db_user web",
#  "-db_pass readonly",
#  "-db_inst cgprod",
#  "-schema pid",
  "-standalone",
  ">> $dest/$err_f 2>&1"
);

system("$cmd");

system("$rundir/dot -Tsvg $dest/$dot_f > $dest/$svg_f 2>>$dest/$err_f");
system("$rundir/dot -Tgif $dest/$dot_f  > $dest/$gif_f 2>>$dest/$err_f");
system("$rundir/dot -Tcmap $dest/$dot_f > $dest/$map_f 2>>$dest/$err_f");

open(OUT, ">$dest/$map_f.1");
print OUT qq!
<html>
<head>
</head>
<body>
<img src="$pex.gif" usemap="#map">
<map name="map">
!;

open(INF, "$dest/$map_f");
while (<INF>) {
  print OUT;
}
close INF;

print OUT qq!
</map>
</body>
</html>
!;

system("mv $dest/$map_f.1 $dest/$map_f");

my $cmd = join(" ",
  "$rundir/Driver.pl",
  "-it xml",
  "-if \"$dest/$xml_f\"",
  "-ot html",
  "-of \"$dest/$html_f\"",
#  "-db_user web",
#  "-db_pass readonly",
#  "-db_inst cgprod",
#  "-schema pid",
  "-standalone",
  ">> $dest/$err_f 2>&1"
);

system("$cmd");

open(INF, "$dest/$err_f");
open(OUT, ">$dest/$err_f.1");
while (<INF>) {
  s/[\r\n]+//;
  print OUT "$_\r\n";
}
close OUT;
close INF;
system("mv $dest/$err_f.1 $dest/$err_f");

## -j option to junk the path

#system("zip -j $dest/$zip_f $dest/$html_f $dest/$svg_f " .
#    "$dest/gif_f $dest/map_f $dest/$err_f");
system("tar -cf $dest/$zip_f " .
    "-C $dest $svg_f " .
    "-C $dest $html_f " .
    "-C $dest $map_f " .
    "-C $dest $gif_f " .
    "-C $dest $xml_f " .
    "-C $dest $rule_f " .
    "-C $dest $err_f");

for my $f ($js_f, $err_f, $dot_f, $svg_f, $html_f, $gif_f, $map_f,
    $xml_f, $rule_f) {
  system("chmod 777 $dest/$f > /dev/null 2>&1");
#  system("rm -f $dest/$f > /dev/null 2>&1");
}

for my $f ($zip_f) {
  system("chmod 777 $dest/$f > /dev/null 2>&1");
}

system("cat $dest/$zip_f");

for my $f ($zip_f) {
#  system("rm -f $dest/$f > /dev/null 2>&1");
}



