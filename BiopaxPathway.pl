#!/usr/local/bin/perl


# Copyright SRA International
#
# Distributed under the OSI-approved BSD 3-Clause License.
# See http://ncip.github.com/pathway-interaction-database/LICENSE.txt for details.


BEGIN {
  my @path_elems = split("/", $0);
  pop @path_elems;
  push @INC, join("/", @path_elems);
}

use strict;
use DBI;
use CGI;

if (-d "/app/oracle/product/dbhome/current") {
  $ENV{'ORACLE_HOME'} = "/app/oracle/product/dbhome/current";
} elsif (-d "/app/oracle/product/8.1.7") {
  $ENV{'ORACLE_HOME'} = "/app/oracle/product/8.1.7";
} elsif (-d "/app/oracle/product/8.1.6") {
  $ENV{'ORACLE_HOME'} = "/app/oracle/product/8.1.6";
}

my $PREDEFINED_PATH = "/share/content/PID/data/predefined";

my ($db_inst, $db_user, $db_pass, $schema) = (
    "cgprod", "web", "readonly", "pid"
);

my $query = new CGI;

my $word = $query->param("word");
my $cmd  = $query->param("cmd");
my $pid  = $query->param("pid");

$word = lc($word);
$cmd  = lc($cmd);

my $db = DBI->connect("DBI:Oracle:" . $db_inst,
      $db_user, $db_pass);
if (not $db or $db->err()) {
  print "Content-type: text/plain\n\n";
  print STDERR "Cannot connect to database\n";
  exit;
}

print "Content-type: text/plain\n\n";

if ($cmd eq "search") {
  Search($word);
} else {
  Retrieve($pid);
}

######################################################################
sub Retrieve {

  system("cat $PREDEFINED_PATH/$pid.bpx");

}

######################################################################
sub Search {

  my ($sql, $stm);
  my ($source, $pname, $ext_id, $pid, $pname1);

  $sql = "select " .
      "s.source_name, p.pathway_name, p.ext_pathway_id, p.pathway_id " .
      "from $schema.pw_source s, $schema.pw_pathway p " .
      "where s.source_id = p.pathway_source_id " .
      "order by s.source_name, p.pathway_name";

  $stm = $db->prepare($sql);
  if(not $stm) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "prepare call failed\n";
    exit;
  }
  if(!$stm->execute()) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "execute call failed\n";
    exit;
  }
  while (($source, $pname, $ext_id, $pid) =
      $stm->fetchrow_array()) {
    $pname1 = lc($pname);
    $source =~ s/NATURE/NCI-Nature Curated/;
    if ($pname1 =~ /$word/ || $ext_id =~ /$word/) {
      print join("\t",
          $source,
          $pname,
          $ext_id,
          $pid,
      ) . "\n";
    }
  }
  $db->disconnect();
}
