######################################################################
# Scan.pm
######################################################################


# Copyright SRA International
#
# Distributed under the OSI-approved BSD 3-Clause License.
# See http://ncip.github.com/pathway-interaction-database/LICENSE.txt for details.


use Exporter ();
@ISA = qw(Exporter);
@EXPORT = qw (
  Scan
);

######################################################################
sub Scan {
  my $ret_val = 0; 
  for my $input (@_) {
    my $ret = exam($input); 
    $ret_val = $ret_val + $ret;
  } 
  return $ret_val;
}

######################################################################
sub exam {
  my ($inp) = @_; 
  my $input = $inp;
  my $ret_val = 0;
  ## if( ($input =~ /javascript/i) or ($input =~ /<script>.+<\/script>/i) ) {
  if( ($input =~ /javascript/i) or ($input =~ /\<script\>/i) or ($input =~ /\<\/script\>/i) or ($input =~ /vbscript/i) ) {
    ## print "<br><b><center>Error in input: $input</b>!</center>";
    print "<br><b><center>Error in input</b>!</center>";
    return 1;
  }
    
  if( $input =~ /=/ ) {
    my @tmp = split "=", $input;
    for (my $i=0; $i<@tmp; $i+2) {
      $tmp[$i] =~ s/\s+$//;
      my @left = split /\s+/, $tmp[$i]; 
      $tmp[$i+1] =~ s/^\s+//;
      my @right = split /\s+/, $tmp[$i+1]; 
      my $left_index = @left;
      if( $left[$left_index-1] eq $right[0] ) {
        print "<br><b><center>Error in input</b>!</center>";
        ## print "<br><b><center>Error in input: $input</b>!</center>";
        return 1;
      }
    }
  }
  if( $input =~ /\|\|/ ) {
    print "<br><b><center>Error in input</b>!</center>";
    ## print "<br><b><center>Error in input: $input</b>!</center>";
    return 1 ;
  }
  ## if( ($input =~ /\'\-\-/) or ($input =~ /\'\s+\-\-/) ) {
  if( $input =~ /\-\-/ ) {
    print "<br><b><center>Error in input</b>!</center>";
    ## print "<br><b><center>Error in input: $input</b>!</center>";
    return 1;
  }
  return $ret_val;
}

######################################################################

1;

