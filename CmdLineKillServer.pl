#!/usr/local/bin/perl

use strict;
use FileHandle;
use Socket;

BEGIN {
  my @path_elems = split("/", $0);
  pop @path_elems;
  push @INC, join("/", @path_elems);
}

use Blocks;

######################################################################
sub CatchPipe {
  print STDERR "Caught PIPE signal\n";
}
$SIG{PIPE} = \&CatchPipe;

######################################################################

my ($host, $port)  = @ARGV;

my $proto          = getprotobyname('tcp');
my $request = "KillServer()";
my $fh = new FileHandle;

if ($host eq "") {
  print "host missing\n";
  exit;
}
if ($port eq "") {
  print "port missing\n";
  exit;
}
if ($host !~ /\./) {
  $host .= ".nci.nih.gov";
}

print "Reseting host = $host, port = $port\n";

my $iaddr          = gethostbyname($host);
my $sin            = sockaddr_in($port, $iaddr);

if( !socket($fh, PF_INET, SOCK_STREAM, $proto) ) {
  print "Cannot open socket to $host:$port\n";
  exit;
}
if( !connect($fh, $sin) ) {
  print "Cannot connect to $host:$port, $!\n";
  exit;
}
if( !SendBlocks($fh, \$request) ) {
  print "SendBlocks failed for $host:$port\n";
  exit;
}
close($fh);

print "done\n";
