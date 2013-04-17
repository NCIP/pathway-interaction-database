

# Copyright SRA International
#
# Distributed under the OSI-approved BSD 3-Clause License.
# See http://ncip.github.com/pathway-interaction-database/LICENSE.txt for details.


package Cache;
require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(MakeCacheFile, FindCacheFile, $CACHE_FAIL);

my $CACHE_FAIL = 0;

use constant MAX_CACHE_ID   => 500;
use constant MAX_TRIES      => 30;
use constant FIRST_CACHE_ID => 1;

######################################################################
sub new {
  my ($class, $cache_path, $cache_prefix) = @_;
  my $self = {};

  $self->{CACHE_PATH} = $cache_path;
  $self->{CACHE_PREFIX} = $cache_prefix;

  $self->{CACHE_PATH} =~ s/\/$//;

  if (not (-e $self->{CACHE_PATH})) {
    die "cache path $self->{CACHE_PATH} does not exist";
  }
  if (not (-w $self->{CACHE_PATH})) {
    die "cache path $self->{CACHE_PATH} is not writable";
  }

  my $status = "not_ok";
  my $sequence_file  = "$self->{CACHE_PATH}/$self->{CACHE_PREFIX}.txt";
  my $temp_file1     = "$self->{CACHE_PATH}/$self->{CACHE_PREFIX}.tmp1";
  my $temp_file2     = "$self->{CACHE_PATH}/$self->{CACHE_PREFIX}.tmp2";
  for (my $tries = 1; $tries <= MAX_TRIES; $tries++) {
    if (-e $sequence_file) {
      $status = "ok";
      last;
    } else {
      sleep 1;
    }
  }
  if ($status eq "ok") {
    return bless $self;
  } else {
    if (-e $temp_file1) {
      unlink $temp_file1;
    }
    if (-e $temp_file2) {
      unlink $temp_file2;
    }
    print STDERR "Unlinked both files\n";
    if (open(OUT, ">$sequence_file")) {
      print OUT "1";
      close OUT;
      chmod 0666, $sequence_file;
      return bless $self;
    } else {
      die "Could not create PW.txt";
    }
 
  }


}

######################################################################
sub MakeCacheFile {
  my ($self) = @_;

  my $cache_id = GetCacheId($self);
  if ($cache_id == $CACHE_FAIL) {
    return ($CACHE_FAIL, "");
  } else {
    my $fname = "$self->{CACHE_PATH}/$self->{CACHE_PREFIX}.$cache_id";
    return ($cache_id, $fname);
  }

}

######################################################################
sub FindCacheFile {
  my ($self, $cache_id) = @_;

  if ($cache_id == $CACHE_FAIL) {
    return "";
  } elsif( -e "$self->{CACHE_PATH}/$self->{CACHE_PREFIX}.$cache_id" ) {
    return "$self->{CACHE_PATH}/$self->{CACHE_PREFIX}.$cache_id";
  } else {
    return "";
  } 
}

######################################################################
sub GetCacheId {

  my ($self) = @_;
  my $sequence_file  = "$self->{CACHE_PATH}/$self->{CACHE_PREFIX}.txt";
  my $temp_file1     = "$self->{CACHE_PATH}/$self->{CACHE_PREFIX}.tmp1";
  my $temp_file2     = "$self->{CACHE_PATH}/$self->{CACHE_PREFIX}.tmp2";

  my $id = $CACHE_FAIL;

  for (my $tries = 1; $tries <= MAX_TRIES; $tries++) {
    if (-e $sequence_file) {
      if (rename $sequence_file, $temp_file1) {
        if (open(IN, $temp_file1)) {
          while (<IN>) {
            $id = $_;
          }
          close $temp_file1;
          if ($id >= MAX_CACHE_ID) {
            $id = FIRST_CACHE_ID;
          } else {
            $id = $id + 1;
          }
          if (open(OUT, ">$temp_file2")) {
            print OUT $id;
            close OUT;
            unlink $temp_file1;
            rename $temp_file2, $sequence_file;
            chmod 0666, $sequence_file;
          } else {
            $id = $CACHE_FAIL;
            print STDERR "Cannot open $temp_file2\n";
            rename $temp_file1, $sequence_file;
            chmod 0666, $sequence_file;
          }
          last;
        } else {
          print STDERR "Cannot open $temp_file1\n";
          rename $temp_file1, $sequence_file;
          chmod 0666, $sequence_file;
          last;
        }
      } else {
        print STDERR "Cannot rename $sequence_file to $temp_file1\n";
        last;
      }
    } else {
      sleep 1;
    }
  }

  return $id;

}

######################################################################
1;
