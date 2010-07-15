#!/usr/local/bin/perl

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
  InitializeComplex2Component();
  InitializeAtom2Atom();
}


######################################################################
sub GetPwImage {
  my ($base, $id) = @_;

  return GetPwImage_1 ($base, $id);  
}

######################################################################
sub MakeCommandFile {
  my ($base, $graphic_type, $parm_string) = @_;

  $parm_string = uri_unescape($parm_string);
  return MakeCommandFile_1 ($base, $graphic_type, $parm_string);
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
sub MoleculeInstance {
  my ($base, $inst) = @_;

  return MoleculeInstance_1($base, $inst);
}

######################################################################
sub PathwayPage {
  my ($base, $pathid) = @_;

  return PathwayPage_1($base, $pathid);
}

######################################################################
sub ListAllPathways {
  my ($base) = @_;

  return ListAllPathways_1($base);
}

######################################################################
sub MoleculeHeader {
  my ($base, $molecule, $source_id, $format) = @_;

  return MoleculeHeader_1($base, $molecule, $source_id, $format);
}

######################################################################
sub PathwayHeader {
  my ($base, $pid, $format) = @_;

  return PathwayHeader_1($base, $pid, $format);
}

######################################################################
sub SearchPathwayKeywords {
  my ($base, $word, $format) = @_;

  return SearchPathwayKeywords_1($base, $word, $format);
}

######################################################################
sub Comment {
  my ($base, $name, $email, $organization, $atomid, $molname, $molid,
      $pathname, $pathid, $comment) = @_;

  return Comment_1($base, $name, $email, $organization,
      $atomid, $molname, $molid, $pathname, $pathid, $comment);
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
    "KillServer",
    "ResetServer",
    "PrPath",
    "GetPwImage",
    "MoleculePage",
    "MoleculeInstance",
    "PathwayPage",
    "ListAllPathways",
    "AtomPage",
    "MakeCommandFile",
    "SearchPathwayKeywords",
    "PathwayHeader",
    "MoleculeHeader",
    "Comment"
);

SetForkable(
    "PrPath",
    "GetPwImage",
    "MoleculePage",
    "MoleculeInstance",
    "PathwayPage",
    "ListAllPathways",
    "AtomPage",
    "MakeCommandFile",
    "SearchPathwayKeywords",
    "PathwayHeader",
    "MoleculeHeader"
###    "Comment"          ##!!! do not make this forkable
);

InitializeDatabase();
# StartServer(PW_SERVER_PORT, "PWAppServer");

#print PrPath("", "text", join("\n",
#  "db_inst	lpgprod",
#  "db_user	web",
#  "db_pass	readonly",
#  "db_schema	cgap2",
#  "print	dot",
#  "mol_name	brca1",
#  "value	9	mol_name	brca1"
#	    ));
exit();

#print PrPath_1("","text", join("\n",
#"db_user\tweb",
#"db_inst\tlpgprod",
#"db_pass\treadonly",
#"db_schema\tcgap2",
#"print\tdot",
#"mol_id\t20035") . "\n");

exit();

