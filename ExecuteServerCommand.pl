#!/usr/local/bin/perl


# Copyright SRA International
#
# Distributed under the OSI-approved BSD 3-Clause License.
# See http://ncip.github.com/pathway-interaction-database/LICENSE.txt for details.


use strict;
use FileHandle;
use CGI;
use Socket;

BEGIN {
  my @path_elems = split("/", $0);
  pop @path_elems;
  push @INC, join("/", @path_elems);
}

use Blocks;

my $query = new CGI;
my $what  = $query->param("what");
my $host  = $query->param("host");
my $port  = $query->param("port");
my $cmd   = $query->param("cmd");

if ($what eq "exec") {
  DoExec();
} else {
  print "Content-type: text/html\n\n";
  DoForm();
}

######################################################################
sub CatchPipe {
  print STDERR "Caught PIPE signal\n";
}
$SIG{PIPE} = \&CatchPipe;

######################################################################
sub DoForm {

print qq!
<head>
</head>
<form action=\"ExecuteServerCommand.pl\" method=POST>
<input type=hidden name=what value=exec>
host: <input type=text size=30 name=host><br>
port: <input type=text size=6 name=port><br>
command: <input type=text size=60 name=cmd><br>
<input type=submit>
</form>
<body>
</body>
!;

}

######################################################################
sub DoExec {

  my $proto   = getprotobyname('tcp');
  my $request = $cmd;
  my $fh      = new FileHandle;

  if ($host eq "") {
    print "Content-type: text/html\n\n";
    print "host missing\n";
    DoForm();
    exit;
  }
  if ($port eq "") {
    print "Content-type: text/html\n\n";
    print "port missing\n";
    DoForm();
    exit;
  }
  if ($host !~ /\./) {
    $host .= ".nci.nih.gov";
  }

  my $iaddr          = gethostbyname($host);
  my $sin            = sockaddr_in($port, $iaddr);

  if( !socket($fh, PF_INET, SOCK_STREAM, $proto) ) {
    print "Content-type: text/html\n\n";
    print "Cannot open socket to $host:$port\n";
    DoForm();
    exit;
  }
  if( !connect($fh, $sin) ) {
    print "Content-type: text/html\n\n";
    print "Cannot connect to $host:$port, $!\n";
    DoForm();
    exit;
  }
  if( !SendBlocks($fh, \$request) ) {
    print "Content-type: text/html\n\n";
    print "SendBlocks failed for $host:$port\n";
    DoForm();
    exit;
  }
  my $response;
  RecvBlocks($fh, \$response);
  close($fh);
  print "Content-type: text/plain\n\n";
  print $response;

}
