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
  # unshift @INC, ".";
  my @path_elems = split("/", $0);
  pop @path_elems;
  unshift @INC, join("/", @path_elems);
}

if (-d "/app/oracle/product/dbhome/current") {
  $ENV{'ORACLE_HOME'} = "/app/oracle/product/dbhome/current";
}


use strict;
use Blocks;
use PWAppConfig;
use PWApp;
use Scan;
use Pathway;
use PathwayDB;
use PathwayXML;
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
  my @lines;
  my $scancheck = Scan($id);
  if ($scancheck > 0) {
    push @lines, "Bad data in parameters\n";
    return join("\n", @lines) . "\n";
  }

 
  return GetPwImage_1 ($base, $id);  
}

######################################################################
sub MakeCommandFile {
  my ($base, $graphic_type, $parm_string) = @_;

  my @lines;
  my $scancheck = Scan($base) + Scan($graphic_type) + Scan($parm_string);
  if ($scancheck > 0) {
    push @lines, "Bad data in parameters\n";
    return join("\n", @lines) . "\n";
  }

  $parm_string = uri_unescape($parm_string);
  return MakeCommandFile_1 ($base, $graphic_type, $parm_string);
}

######################################################################
sub PrPath {
  my ($base, $graphic_type, $parm_string) = @_;
  my @lines;
  my $scancheck = Scan($base) + Scan($graphic_type) + Scan($parm_string);
  if ($scancheck > 0) {
    push @lines, "Bad data in parameters\n";
    return join("\n", @lines) . "\n";
  }

  $parm_string = uri_unescape($parm_string);
  return PrPath_1 ($base, $graphic_type, $parm_string);
}

######################################################################
sub PrBatch {
  my ($base, $graphic_type, $parm_string) = @_;
  my @lines;
  my $scancheck = Scan($base) + Scan($graphic_type) + Scan($parm_string);
  if ($scancheck > 0) {
    push @lines, "Bad data in parameters\n";
    return join("\n", @lines) . "\n";
  }

  $parm_string = uri_unescape($parm_string);
  return PrBatch_1 ($base, $graphic_type, $parm_string);
}

######################################################################
sub MoleculePage {
  my ($base, $molid) = @_;
  my @lines;
  my $scancheck = Scan($base) + Scan($molid);
  if ($scancheck > 0) {
    push @lines, "Bad data in parameters\n";
    return join("\n", @lines) . "\n";
  }

  return MoleculePage_1($base, $molid);
}

######################################################################
sub MoleculeInstance {
  my ($base, $inst) = @_;

  my @lines;
  my $scancheck = Scan($base) + Scan($inst);
  if ($scancheck > 0) {
    push @lines, "Bad data in parameters\n";
    return join("\n", @lines) . "\n";
  }

  return MoleculeInstance_1($base, $inst);
}

######################################################################
sub PathwayPage {
  my ($base, $pathid) = @_;
  my @lines;
  my $scancheck = Scan($base) + Scan($pathid);
  if ($scancheck > 0) {
    push @lines, "Bad data in parameters\n";
    return join("\n", @lines) . "\n";
  }

  return PathwayPage_1($base, $pathid);
}

######################################################################
sub ListAllPathways {
  my ($base) = @_;
  my @lines;
  my $scancheck = Scan($base);
  if ($scancheck > 0) {
    push @lines, "Bad data in parameters\n";
    return join("\n", @lines) . "\n";
  }

  return ListAllPathways_1($base);
}

######################################################################
sub SearchHeader {
  my ($base, $atomid, $format) = @_;
  my @lines;
  my $scancheck = Scan($base) + Scan($atomid) + Scan($format);
  if ($scancheck > 0) {
    push @lines, "Bad data in parameters\n";
    return join("\n", @lines) . "\n";
  }

  return SearchHeader_1($base, $atomid, $format);
}

######################################################################
sub InteractionHeader {
  my ($base, $atomid, $format) = @_;
  my @lines;
  my $scancheck = Scan($base) + Scan($atomid) + Scan($format);
  if ($scancheck > 0) {
    push @lines, "Bad data in parameters\n";
    return join("\n", @lines) . "\n";
  }

  return InteractionHeader_1($base, $atomid, $format);
}

######################################################################
sub MoleculeHeader {
  my ($base, $molecule, $source_id, $format) = @_;
  my @lines;
  my $scancheck = Scan($base) + Scan($molecule) + Scan($source_id) + Scan($format);
  if ($scancheck > 0) {
    push @lines, "Bad data in parameters\n";
    return join("\n", @lines) . "\n";
  }

  return MoleculeHeader_1($base, $molecule, $source_id, $format);
}

######################################################################
sub IntermediatePage {
  my ($base, $molecule) = @_;

  my @lines;
  my $scancheck = Scan($base) + Scan($molecule);
  if ($scancheck > 0) {
    push @lines, "Bad data in parameters\n";
    return join("\n", @lines) . "\n";
  }

  return IntermediatePage_1($base, $molecule);
}

######################################################################
sub AdvancedHeader {
  my ($base, $molecule, $source_id, $evidence_code, $format) = @_;
  my @lines;
  my $scancheck = Scan($base) + Scan($molecule) + Scan($source_id) + Scan($evidence_code) + Scan($format);
  if ($scancheck > 0) {
    push @lines, "Bad data in parameters\n";
    return join("\n", @lines) . "\n";
  }

  return AdvancedHeader_1($base, $molecule, $source_id, $evidence_code, $format);
}

######################################################################
sub NetworkHeader {
  my ($base, $molecule, $source_id, $format) = @_;
  my @lines;
  my $scancheck = Scan($base) + Scan($source_id) + Scan($format);
  if ($scancheck > 0) {
    push @lines, "Bad data in parameters\n";
    return join("\n", @lines) . "\n";
  }

  return NetworkHeader_1($base, $molecule, $source_id, $format);
}
######################################################################
sub NetworkMolListHeader {
  my ($base, $molecule, $source) = @_;
  my @lines;
  my $scancheck = Scan($base) + Scan($molecule) + Scan($source);
  if ($scancheck > 0) {
    push @lines, "Bad data in parameters\n";
    return join("\n", @lines) . "\n";
  }

  return NetworkMolListHeader_1($base, $molecule, $source);
}

######################################################################
sub BatchHeader {
  my ($base, $molecule, $source_id, $format) = @_;
  my @lines;
  my $scancheck = Scan($base) + Scan($molecule) + Scan($source_id) + Scan($format);
  if ($scancheck > 0) {
    push @lines, "Bad data in parameters\n";
    return join("\n", @lines) . "\n";
  }

  return BatchHeader_1($base, $molecule, $source_id, $format);
}
######################################################################
sub ConnectHeader {
  my ($base, $molecule, $source_id, $format) = @_;
  my @lines;
  my $scancheck = Scan($base) + Scan($molecule) + Scan($source_id) + Scan($format);
  if ($scancheck > 0) {
    push @lines, "Bad data in parameters\n";
    return join("\n", @lines) . "\n";
  }

  return ConnectHeader_1($base, $molecule, $source_id, $format);
}

######################################################################
sub PathwayHeader {
  my ($base, $pid, $genes_a, $genes_b, $format) = @_;
  my @lines;
  my $scancheck = Scan($base) + Scan($pid) + Scan($format);
  if ($scancheck > 0) {
    push @lines, "Bad data in parameters\n";
    return join("\n", @lines) . "\n";
  }

  return PathwayHeader_1($base, $pid, $genes_a, $genes_b, $format);
}

######################################################################
sub CitationHeader {
  my ($base, $pid, $format) = @_;

  my @lines;
  my $scancheck = Scan($base) + Scan($pid) + Scan($format);
  if ($scancheck > 0) {
    push @lines, "Bad data in parameters\n";
    return join("\n", @lines) . "\n";
  }
 return CitationHeader_1($base, $pid, $format);
}

######################################################################
sub MolListHeader {
  my ($base, $pid, $format) = @_;
  my @lines;
  my $scancheck = Scan($base) + Scan($pid) + Scan($format);
  if ($scancheck > 0) {
    push @lines, "Bad data in parameters\n";
    return join("\n", @lines) . "\n";
  }

  return MolListHeader_1($base, $pid, $format);
}

######################################################################
sub SearchPathwayKeywords {
  my ($base, $word, $format) = @_;
  my @lines;
  my $scancheck = Scan($base) + Scan($word) + Scan($format);
  if ($scancheck > 0) {
    push @lines, "Bad data in parameters\n";
    return join("\n", @lines) . "\n";
  }

  return SearchPathwayKeywords_1($base, $word, $format);
}
######################################################################
sub SearchMoleculeKeywords {
  my ($base, $gene_id) = @_;
  my @lines;
  my $scancheck = Scan($base) + Scan($gene_id);
  if ($scancheck > 0) {
    push @lines, "Bad data in parameters\n";
    return join("\n", @lines) . "\n";
  }

  return SearchMoleculeKeywords_1($base, $gene_id);
}

######################################################################
sub Comment {
  my ($base, $name, $email, $organization, $atomid, $molname, $molid,
      $pathname, $pathid, $comment) = @_;
  my @lines;
  my $scancheck = Scan($base) + Scan($name) + Scan($email) + Scan($organization) + Scan($atomid) + Scan($molname) + Scan($molid) + Scan($pathname) + Scan($pathid) + Scan($comment);
  if ($scancheck > 0) {
    push @lines, "Bad data in parameters\n";
    return join("\n", @lines) . "\n";
  }

  return Comment_1($base, $name, $email, $organization,
      $atomid, $molname, $molid, $pathname, $pathid, $comment);
}

######################################################################
sub AtomPage {
  my ($base, $atomid) = @_;
  my @lines;
  my $scancheck = Scan($base) + Scan($atomid);
  if ($scancheck > 0) {
    push @lines, "Bad data in parameters\n";
    return join("\n", @lines) . "\n";
  }

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
    "PrBatch",
    "ServerInfo",
    "GetPwImage",
    "MoleculePage",
    "MoleculeInstance",
    "PathwayPage",
    "ListAllPathways",
    "AtomPage",
    "MakeCommandFile",
    "SearchPathwayKeywords",
    "SearchMoleculeKeywords",
    "PathwayHeader",
    "CitationHeader",
    "MolListHeader",
    "NetworkMolListHeader",
    "MoleculeHeader",
    "IntermediatePage",
    "AdvancedHeader",
    "NetworkHeader",
    "BatchHeader",
    "ConnectHeader",
    "SearchHeader",
    "InteractionHeader",
    "Comment"
);

SetForkable(
    "PrPath",
    "PrBatch",
    "GetPwImage",
    "MoleculePage",
    "MoleculeInstance",
    "PathwayPage",
    "ListAllPathways",
    "AtomPage",
    "MakeCommandFile",
    "SearchPathwayKeywords",
    "SearchMoleculeKeywords",
    "PathwayHeader",
    "CitationHeader",
    "MolListHeader",
    "MoleculeHeader",
    "NetworkMolListHeader",
    "IntermediatePage",
    "AdvancedHeader",
    "NetworkHeader",
    "BatchHeader",
    "ConnectHeader",
    "SearchHeader",
    "InteractionHeader"
###    "Comment"          ##!!! do not make this forkable
);

InitializeDatabase();
StartServer(PW_SERVER_PORT, "PWAppServer");

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

