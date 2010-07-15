#!/usr/local/bin/perl

use strict;
use constant HOURS_TO_LIVE => 5;
use constant HOURS_TO_LIVE_FOR_GENOMICS => 148;
use constant TIME_TO_LIVE  => HOURS_TO_LIVE * 60 * 60;
use constant TIME_TO_LIVE_FOR_GENOMICS  => HOURS_TO_LIVE_FOR_GENOMICS * 60 * 60;

my (
    $CACHE_PATH,
    $PREFIXES
) = @ARGV;


my ($dev, $ino, $mode, $nlink, $uid, $gid, $rdev, $size,
    $atime, $mtime, $ctime, $blksize, $blocks);

my ($filename, %prefixes, $prefix, $suffix); 
my $now = time();

for $prefix (split ",", $PREFIXES) {
  $prefixes{$prefix} = 1;
}

opendir (CACHE_DIR, $CACHE_PATH) || die "can not open $CACHE_PATH";  
while( $filename = readdir(CACHE_DIR) ) {
  if( !($filename =~ /^\./) ) {
    if( $filename =~ /^[a-zA-Z]+\.\d+$/) {
      ($prefix, $suffix) = split /\./, $filename;
      if (defined $prefixes{$prefix}){
        $filename = $CACHE_PATH . "/" . $filename;
        ($dev, $ino, $mode, $nlink, $uid, $gid, $rdev, $size,
            $atime, $mtime, $ctime, $blksize, $blocks) = stat ($filename);
        if ( $prefix eq "GENOMICS") {
          if ($now - $atime > TIME_TO_LIVE_FOR_GENOMICS) {
            print "deleting file: $filename\n";
            unlink($filename);
          }
        }
        else {
          if ($now - $atime > TIME_TO_LIVE) {
            print "deleting file: $filename\n";
            unlink($filename);
          }
        }
      }
    }
  }
}
closedir CACHEDIR;

