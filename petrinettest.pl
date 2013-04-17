#!/usr/local/bin/perl


# Copyright SRA International
#
# Distributed under the OSI-approved BSD 3-Clause License.
# See http://ncip.github.com/pathway-interaction-database/LICENSE.txt for details.


##
## More options to consider:
##   - save the input config file and allow retrieval
##   - allow uploading a config file
##   - limit to atoms of a given data source (maybe also with NOT)
##   - allow retrieval by template
##   - allow retrieval by degree (on incoming and/or outgoing edge) of
##     separation (but don't follow pruned mols)
##

BEGIN {
  my @path_elems = split("/", $0);
  pop @path_elems;
  push @INC, join("/", @path_elems);
}

use strict;
use Blocks;
use PWAppConfig;
use PWApp;
use Pathway;
use PathwayDB;
use PathParams;
use PWLabel;
use FileHandle;
use ServerSupport;
use URI::Escape;

######################################################################
sub InitializeDatabase {

  InitializeLV();
}


######################################################################
sub GetPwImage {
  my ($base, $id) = @_;

  return GetPwImage_1 ($base, $id);  
}

######################################################################
sub PrPath {
  my ($base, $graphic_type, $parm_string) = @_;

  $parm_string = uri_unescape($parm_string);
  return PrPath_1 ($base, $graphic_type, $parm_string);
}

######################################################################
sub MoleculePage {
  my ($base, $molid) = @_;

  return MoleculePage_1($base, $molid);
}

######################################################################
sub ListAllPathways {
  my ($base) = @_;

  return ListAllPathways_1($base);
}

######################################################################
sub AtomPage {
  my ($base, $atomid) = @_;

  return AtomPage_1($base, $atomid);
}

######################################################################

######################################################################
#
# main
#

SetProgramName($0);

SetSafe(
    "PrPath",
    "GetPwImage",
    "MoleculePage",
    "ListAllPathways",
    "AtomPage"
);

SetForkable(
    "PrPath",
    "GetPwImage",
    "MoleculePage",
    "ListAllPathways",
    "AtomPage"
);

InitializeDatabase();

my $value = PrPath_1("", "text", join("\n",
  "db_inst	cgprod",
  "db_user	web",
  "db_pass	readonly",
  "db_schema	cgap",
  "print	petrinet",
  "pathway_name	g1/s check point",
  "value	9	mol_name	tgfb1",
  ""
	    ));

my @tmp = split /\001/, $value;
print $tmp[2];
