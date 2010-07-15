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
use PWLabel;
use FileHandle;
use IdMap;

my (
  $db_inst,
  $db_user,
  $db_pass,
  $schema,
  $input_type,
  $input_file,
  $output_type,
  $output_file,
  $ontology_file,
  $id_file,
  $moldefs_file,
  $standalone,
  $atomids,
  $subtypelines
);

my %output_types = (
  "sql"    => 1,
  "table"  => 1,
  "xml"    => 1,
  "dot"    => 1,
  "biopax" => 1,
  "html"   => 1
);
my %input_types = (
  "js"     => 1,
  "xml"    => 1,
  "biopax" => 1
);
my ($lv, $pw, $idm);

##
## Set a timer
##
my $query_timeout = 500;       ## max length (in seconds) for a request
                               ## to run before exiting/resetting server
# $SIG{ALRM} = \&CatchAlarm;
alarm $query_timeout;         ## Set timer
##
##
##

ReadOptions();
print STDOUT "After ReadOptions\n";
##
## Read (if was given) a file that tells where to start numbering
## new entities (pathways, atoms, molecules, ptms)
##

$idm = new IdMap($id_file);
print STDOUT "After IdMap\n";

##
## Read the ontology from an xml input or from the database
##

if ($input_type eq "xml") {
  $lv = new PWLabel(undef, undef);
} elsif (defined $ontology_file ||
    (defined $db_user && defined $db_pass && defined $db_inst &&
      defined $schema)) {
  $lv = ReadOntology($ontology_file, $db_user, $db_pass, $db_inst, $schema);
} else {
  die "no ontology input specified";
}

##
## Read existing mol defs from an xml input or from the database
##

my $moldefs = new Pathway($lv);
if ($moldefs_file) {
  ReadMolDefs($lv, $moldefs_file, $moldefs);
} elsif (defined $db_user && defined $db_pass && defined $db_inst &&
      defined $schema) {
  use PathwayDB;
  use DBI;
  my $db = DBI->connect("DBI:Oracle:" . $db_inst, $db_user, $db_pass);
  if (not $db or $db->err()) {
    print STDERR "Cannot connect to " . $db_user . "@" . $db_inst . "\n";
    exit;
  }
  my $pdb = new PathwayDB($db, $schema, "", "", $moldefs, $lv);
  $pdb->GetAllMolecules();
  $db->disconnect();
}

##
## Read the input to be translated/processed
##
print STDOUT "Before ReadInput\n";
$pw = ReadInput($input_type, $input_file, $lv);
print STDOUT "After ReadInput\n";
exit(0);
##
## And produce the translation/analysis
##
print STDOUT "Before WriteOutput\n";
WriteOutput($output_type, $output_file, $lv, $pw);

######################################################################
sub ReadMolDefs {
  my ($lv, $f, $moldefs) = @_;

  use ParsePathwayXML;

  my $fh = new FileHandle;
  open($fh, $f) or die "cannot open $f";
  my $parser = new ParsePathwayXML($lv, $moldefs);
  $parser->parse($fh);
  close $fh;
  $moldefs->MergeMolecules();
}

######################################################################
sub ReadInput {
  my ($input_type, $f, $lv) = @_;

  my $pw;
  if ($input_type eq "xml") {
    $pw = ReadPathwayXML($lv, $f);
    $pw->MergeMolecules();
  } elsif ($input_type eq "js") {
    use ParseJS;
    use MolMatch;
    $pw = new Pathway($lv);
    my $js = new ParseJS($lv, $pw, $idm);
    print STDOUT "After new ParseJS\n";
    $js->ParseJS($f);
    print STDOUT "After ParseJS\n"; 
    $pw->MergeMolecules();
    print STDOUT "After MergeMolecules\n";
    my $mm = new MolMatch($lv, $moldefs, $pw);
    $mm->MatchMols();
    $pw->RemapMols();
  } elsif ($input_type eq "biopax") {
  }

  $pw->BuildMolInstCache();
  $pw->BuildGenericLabelCache();
  if (0) {
    $pw->PruneDuplicateAtoms();
  }
  $pw->IdentifyMacroProcesses();
  $pw->ValidateSubnets();
  $pw->CollapseSubnets();
  return $pw;
}

######################################################################
sub WriteOutput {
  my ($output_type, $output_file, $lv, $pw) = @_;

  if ($output_type eq "html") {

    use HTMLOutput;
    open(OUTF, ">$output_file") or die "cannot open $output_file";
    my $xml = new HTMLOutput($pw, $lv, *OUTF);
    $xml->PrHTML();
    close OUTF;

  } elsif ($output_type eq "xml") {

    use XMLOutput;
    open(OUTF, ">$output_file") or die "cannot open $output_file";
    my $xml = new XMLOutput($pw, $lv, *OUTF);
    $xml->PrXML();
    close OUTF;

  } elsif ($output_type eq "biopax") {

    use BioPAXOutput;
    open(OUTF, ">$output_file") or die "cannot open $output_file";
    my $xml = new BioPAXOutput($pw, $lv, *OUTF);
    $xml->PrOWL();
    close OUTF;

  } elsif ($output_type eq "table") { 

    use SQLOutput;
    open(OUTF, ">$output_file") or die "cannot open $output_file";
    my $sql = new SQLOutput($pw, $lv, *OUTF, "table");
    $sql->SetPTMExprId($idm->PTMExprId());
    my %list;
    for my $molid (@{ $pw->Mols() }) {
      if ($idm->MolId() > 0 && $molid < $idm->MolId()) {
        $list{$molid} = 1;
      }
    }
    $sql->SetExistingDefs(\%list);
    $sql->PrSQL();
    close OUTF;

  } elsif ($output_type eq "sql") {

    use SQLOutput;
    open(OUTF, ">$output_file") or die "cannot open $output_file";
    my $sql = new SQLOutput($pw, $lv, *OUTF, "sql");
    $sql->PrSQL();
    close OUTF;

  } elsif ($output_type eq "dot") { 
    use Cache;
    use Clan;
    use DOToutput;

    my @lines;
    my ($size, $clan, @clans, %sizeof, @cache_ids);

    if ($standalone) {
      ## Why? There are some pathways that have non-overlapping
      ## graphs showing processes involving different members of
      ## a molecule-family. We could try to box them (or some such
      ## graphic trick), but that might force a closer relation than
      ## we want; for example, all forms (ptm) and locations of each
      ## family member would go in the same box.
      $pw->ForceIntoOneClan();
    } else {
      $pw->BuildClanList();
    }
    for $clan (@{ $pw->Clans }) {
      $size = $clan->ClanSize();
      push @{ $sizeof{$size} }, $clan;
    }
    for $size (sort r_numerically keys %sizeof) {
      for $clan (@{ $sizeof{$size} }) {
        push @clans, $clan;
      }
    }
    my $n = 0;
    for $clan (@clans) {
      $n++;
      open(OUTF, ">$output_file.$n") or die "cannot open $output_file.$n";
      my $dot = new DOToutput($pw, $lv, *OUTF);
      if ($atomids) {
        $dot->ShowAtomIds(1);
      }
      if ($standalone) {
        $dot->SetStandalone(1);
        $dot->ShowAtomIds(1);
        $dot->ShowSubTypeLines(1);
        my $html_fn = $output_file;
        $html_fn =~ s/.*\///;
        $html_fn =~ s/\.dot$//;
        $html_fn .= ".html";
        $dot->SetHTMLFileName($html_fn);
      }
      if ($subtypelines) {
        $dot->ShowSubTypeLines(1);
      }
      $dot->DOTGraph($clan);
      close(OUTF);
    }
  }
}

######################################################################
sub r_numerically { $b <=> $a };

######################################################################
sub ReadOntology {
  my ($ontology_f, $db_user, $db_pass, $db_inst, $schema) = @_;

  my $lv;
  if ($ontology_f) {
    use ParsePathwayXML;
    my $fh = new FileHandle;
    $lv = new PWLabel(undef, undef);
    open($fh, $ontology_f) or die "cannot open $ontology_f";
    my $pw0 = new Pathway($lv);
    my $parser = new ParsePathwayXML($lv, $pw0);
    $parser->parse($fh);
    close $fh;
  } else {
    my $db = DBI->connect("DBI:Oracle:" . $db_inst, $db_user, $db_pass);
    if (not $db or $db->err()) {
      print STDERR "Cannot connect to " . $db_user . "@" . $db_pass . "\n";
      exit;
    }
    $lv   = new PWLabel($db, $schema);
    $db->disconnect();
  }
  return $lv;
}

######################################################################
sub ReadPathwayXML {
  my ($lv, $xml_f) = @_;

  use ParsePathwayXML;
  my $fh = new FileHandle;
  open($fh, $xml_f) or die "cannot open $xml_f";
  my $pw   = new Pathway($lv);
  my $parser = new ParsePathwayXML($lv, $pw);
  $parser->parse($fh);
  close $fh;
  return $pw;
}

######################################################################
sub ReadOptions {

  use Getopt::Long;

  if (@ARGV == 0) {
    print join("\n\t",
      "options:",
      "it: input type (xml, js, biopax)",
      "ot: output type (xml, biopax, dot, table, sql)",
      "if: input file",
      "of: output file",
      "ontf: ontology file (if not providing db params)",
      "idf: file with starting ids (pathway, atom, mol, ptm) and source ids",
      "moldef: xml pathway file with existing molecule defs " .
          "(into which new stuff is to be mapped)",
      "db_user",
      "db_pass",
      "db_inst",
      "schema",
      "standalone: standalone dot (local urls)",
      "atomids: show atomids in dot",
      "subtypelines: show subtypelines in dot"
    ) . "\n";
    exit;
  }

  GetOptions (
    "it:s"             => \$input_type,
    "ot:s"             => \$output_type,
    "if:s"             => \$input_file,
    "of:s"             => \$output_file,
    "ontf:s"           => \$ontology_file,
    "idf:s"            => \$id_file,
    "moldef:s"         => \$moldefs_file,
    "db_inst:s"        => \$db_inst,
    "db_user:s"        => \$db_user,
    "db_pass:s"        => \$db_pass,
    "schema:s"         => \$schema,
    "standalone"       => \$standalone,
    "atomids"          => \$atomids,
    "subtypelines"     => \$subtypelines,
  ) or die "exiting";

  if (! defined $input_types{$input_type}) {
    die "input type must be one of " .
        join(", ", keys %input_types) . "\n";
  }
  if (! defined $output_types{$output_type}) {
    die "output type must be one of " .
        join(", ", keys %output_types) . "\n";
  }
  if ($input_file eq "") {
    die "specify input file";
  }
  if ($output_file eq "") {
    die "specify output file";
  }
  if ($input_file eq $output_file) {
    die "input and output file names are identical";
  }
  if ($input_type ne "xml" && $ontology_file eq "" && $db_user eq "") {
    die "input type must be xml or else must specify ontology file " .
        "or database parameters";
  }
}

######################################################################
