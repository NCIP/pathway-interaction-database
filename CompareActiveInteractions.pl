#!/usr/local/bin/perl

BEGIN {
  my @path_elems = split("/", $0);
  pop @path_elems;
  push @INC, join("/", @path_elems);
}

if (-d "/app/oracle/product/dbhome/current") {
  $ENV{'ORACLE_HOME'} = "/app/oracle/product/dbhome/current";
} elsif (-d "/app/oracle/product/8.1.7") {
  $ENV{'ORACLE_HOME'} = "/app/oracle/product/8.1.7";
} elsif (-d "/app/oracle/product/8.1.6") {
  $ENV{'ORACLE_HOME'} = "/app/oracle/product/8.1.6";
}

use strict;
use Pathway;
use PathwayDB;
use PWLabel;
use StudyNet;
use ParsePathwayXML;
use FileHandle;
use CGI;

my (
  $present_fA,
  $present_fB,
  $format,
  $what
) = @ARGV;

## $format = svg | gif | xml
## $what   = A | B | AB

my (
  $db_inst,
  $db_user,
  $db_pass,
  $schema
) = ("cgprod", "web", "readonly", "cgap");

my $query     = new CGI;
my $present_fA = $query->param("fA");
my $present_fB = $query->param("fB");
my $what       = $query->param("what");
my $format     = $query->param("format");

use constant NO_VAL               => "-1";
use constant MAX_DOT_INTERACTIONS => 500;
use constant CACHE_ROOT           => "/share/content/CMAP/data/cache";
use constant PW_CACHE_PREFIX      => "PW";
use constant DOT_PGM              => "/share/content/CMAP/run/PW/dot";
use constant XML_F                => "/share/content/CMAP/data/BioCarta.xml";

my $xml_f = XML_F;

if ($what ne "A" && $what ne "B" && $what ne "AB") {
  print "Content-type: text/plain\n\n";
  print STDERR "bad value for parameter 'what'\n";
  exit;
}
if ($format ne "xml" && $format ne "svg" && $format ne "gif") {
  print "Content-type: text/plain\n\n";
  print STDERR "bad value for parameter 'format'\n";
  exit;
}

my $lv = InitializeLV($db_inst, $db_user, $db_pass, $schema);
my $pw = ReadPathwayXML($lv, $xml_f);
$pw->DeleteMolLabelsOfKind("location");
$pw->BuildMolInstCache();
DoComparison();

######################################################################
sub r_numerically { $b <=> $a };

######################################################################
sub DoComparison {

  my $active_callsA = DoOneSample($present_fA);
  my $active_callsB = DoOneSample($present_fB);
  my (%active_calls, @prune);

  if ($what eq "A") {         ## in A - B
    for my $a (keys %{ $active_callsA }) {
      if ($$active_callsA{$a} eq "1" && $$active_callsB{$a} ne "1") {
        $active_calls{$a} = 1;
      }
    }
  } elsif ($what eq "B") {    ## in B - A
    for my $a (keys %{ $active_callsB }) {
      if ($$active_callsB{$a} eq "1" && $$active_callsA{$a} ne "1") {
        $active_calls{$a} = 1;
      }
    }
  } elsif ($what eq "AB") {   ## in A * B
    for my $a (keys %{ $active_callsA }) {
      if ($$active_callsA{$a} eq "1" && $$active_callsB{$a} eq "1") {
        $active_calls{$a} = 1;
      }
    }
  } else {
  }

  for my $a (@{ $pw->Atoms() }) {
    if (not defined $active_calls{$a}) {
      push @prune, $a;
    }
  }

  $pw->PruneAtoms(\@prune, {});

  if ($format eq "xml") {

    use XMLOutput;

    my @lines;
    my $xml = new XMLOutput($pw, $lv);
    $xml->PrXML();

    print "Content-type: text/plain\n\n";
    print join ("\n", @{ $xml->Lines() }) . "\n";

  } elsif ($format eq "svg" || $format eq "gif") {
 
    use Cache;
    use Clan;
    use DOToutput;

    my @lines;

    my ($size, $clan, @clans, %sizeof, @cache_ids);

    print "Content-type: text/html\n\n";
    print "<html>\n";
    print "<body>\n";

    $pw->BuildClanList();
    for $clan (@{ $pw->Clans }) {
      $size = $clan->ClanSize();
      push @{ $sizeof{$size} }, $clan;
    }
    for $size (sort r_numerically keys %sizeof) {
      for $clan (@{ $sizeof{$size} }) {
        push @clans, $clan;
      }
    }
    for $clan (@clans) {
      $size = $clan->ClanSize();
      if ($size > MAX_DOT_INTERACTIONS) {
        print "Number of interactions $size exceeds maximum " .
            MAX_DOT_INTERACTIONS . " in a single subgraph.\n";
        print "Try XML format\n";
        exit;
      }
      my $cache = new Cache(CACHE_ROOT, PW_CACHE_PREFIX);
      my ($cache_id, $filename) = $cache->MakeCacheFile();
      if ($cache_id != $CACHE_FAIL) {
        my $cmd = DOT_PGM .  " -T$format > $filename";
        open(OUTF, "| $cmd");
        my $dot = new DOToutput($pw, $lv, *OUTF);
        $dot->DOTGraph($clan);
        close(OUTF);
        chmod 0666, $filename;
        push @cache_ids, $cache_id;
      } else {
        print "cache failure: " . CACHE_ROOT . "\n";
      }
    }
    my $n;
    for my $id (@cache_ids) {
      $n++;
      if ($format eq "svg") {
        print "<p>Subgraph $n\n";
        print "<p>\n";
        print "<embed type=\"image/svg-xml\" " .
            "height=600 width=800 " .
            "src=\"/cmapcgi/PW/GetPwImage.pl?" .
            "GID=$id&FORMAT=SVG\">\n";
      } else {
        print "<p>Subgraph $n\n";
        print "<p>\n";
        print "<img src=\"/cmapcgi/PW/GetPwImage.pl?" .
            "GID=$id&FORMAT=GIF\">\n";
      }
    }
  }

}

######################################################################
sub DoOneSample {
  my ($present_f) = @_;

  my %present_calls;

  ReadPresentCalls($present_f, \%present_calls);
  my $sn = new StudyNet($pw, $lv);
  return IterateOverPathways($pw, $sn, \%present_calls);
}

######################################################################
sub IterateOverPathways {
  my ($pw, $sn, $present_calls) = @_;

  my %pathway2atom;
  for my $atom (@{ $pw->Atoms() }) {
    for my $pid (@{ $pw->AtomPathway($atom) }) {
      push @{ $pathway2atom{$pid} }, $atom;
    }
  }

##  for my $lid (keys %{ $present_calls }) {
    my %active_calls;
    my %ptms;
    $sn->FindActiveInteractions($pw->Atoms(),
##        $present_calls{$lid}, \%active_calls, \%ptms);
        $present_calls, \%active_calls, \%ptms);
    for my $pid (keys %pathway2atom) {
      my $pexid = $pw->PathwayExId($pid);
      for my $atom (@{ $pathway2atom{$pid} }) {
        if (defined $active_calls{$atom}) {
#          print "$lid\t$pexid\t$atom\t$active_calls{$atom}\n";
        }
      }
      for my $ptm (keys %ptms) {
#        print "#ptm\t$lid\t$pexid\t$ptm\n";
      }
    }
##  }
  return \%active_calls;
}

######################################################################
sub ReadPathwayXML {
  my ($lv, $xml_f) = @_;

  my $fh = new FileHandle;
  open($fh, $xml_f) or die "cannot open $xml_f";
  my $pw   = new Pathway($lv);
  my $parser = new ParsePathwayXML($lv, $pw);
  $parser->parse($fh);
  close $fh;
  return $pw;
}

######################################################################
sub InitializeLV {
  my ($db_inst, $db_user, $db_pass, $schema) = @_;

  my $db = DBI->connect("DBI:Oracle:" . $db_inst,
        $db_user, $db_pass);
  if (not $db or $db->err()) {
    print STDERR "Cannot connect to " . $db_user . "@" .
        $db_pass . "\n";
    exit;
  }
  my $lv   = new PWLabel($db, $schema);
  $db->disconnect();
  return $lv;
}

######################################################################
sub ReadPresentCalls {
  my ($f, $present_calls) = @_;

#  my ($lid, $gene_id, $call);
  my ($lid, $gene_id, $call);

##  open(INF, $f) or die "cannot open $f";
##  while (<INF>) {
  while (<$f>) {
    chop;
#    ($lid, $gene_id, $call) = split /\t/;
    ($gene_id, $call) = split /\t/;
    if ($call eq "P") {
#      $$present_calls{$lid}{$gene_id} = 1;
      $$present_calls{$gene_id} = 1;
    } elsif ($call eq "A") {
#      if (! defined $$present_calls{$lid}{$gene_id}) {
#        $$present_calls{$lid}{$gene_id} = 0;
#      } elsif ($$present_calls{$lid}{$gene_id} eq "P") {
#      } elsif ($$present_calls{$lid}{$gene_id} eq NO_VAL) {
#      }
      if (! defined $$present_calls{$gene_id}) {
        $$present_calls{$gene_id} = 0;
      } elsif ($$present_calls{$gene_id} eq "P") {
      } elsif ($$present_calls{$gene_id} eq NO_VAL) {
      }
    } elsif ($call eq "M") {
#      if (! defined $$present_calls{$lid}{$gene_id}) {
#        $$present_calls{$lid}{$gene_id} = -1;
#      } elsif ($$present_calls{$lid}{$gene_id} eq "P") {
#      } elsif ($$present_calls{$lid}{$gene_id} eq NO_VAL) {
#      }
      if (! defined $$present_calls{$gene_id}) {
        $$present_calls{$gene_id} = -1;
      } elsif ($$present_calls{$gene_id} eq "P") {
      } elsif ($$present_calls{$gene_id} eq NO_VAL) {
      }
    }
  }
  close $f;
##  close INF;
}
