#!/usr/local/bin/perl


# Copyright SRA International
#
# Distributed under the OSI-approved BSD 3-Clause License.
# See http://ncip.github.com/pathway-interaction-database/LICENSE.txt for details.


use strict;
use CGI;
use DBI;
use Cache;

BEGIN {
  my @path_elems = split("/", $0);
  pop @path_elems;
  push @INC, join("/", @path_elems);
}

if (-d "/app/oracle/product/dbhome/current") {
  $ENV{'ORACLE_HOME'} = "/app/oracle/product/dbhome/current";
}

use constant CACHE_ROOT      => "/share/content/PID/data/cache" ;
use constant PW_CACHE_PREFIX => "PW" ;
#use constant DOT_PGM         => "/usr/local/bin/dot";
use constant DOT_PGM         => "/share/content/PID/run/dot";
use constant PREDEF_DIR      => "/share/content/PID/data/predefined";

my $color        = "blue";
my ($db_inst, $db_user, $db_pass, $schema) =
   ("cgprod", "web", "readonly", "pid");

my $query  = new CGI;
#print "Content-type: text/plain\n\n";
print "Content-type: text/html\n\n";

my $pathway_id     = $query->param("pid");
my $entrez_gene_id = $query->param("egid");
my $source         = $query->param("source");
my $format         = $query->param("format");

#my ($pathway_id, $entrez_gene_id, $source) = ("200018", "5781", "NATURE");

my (%mol2color);

my $db = DBI->connect("DBI:Oracle:" . $db_inst, $db_user, $db_pass);
if (not $db or $db->err()) {
  print STDERR "Cannot connect to " . $db_user . "@" . $db_inst . "\n";
  exit;
}

my $mols = GetMolIdsForEntrezGene($db, $pathway_id, $entrez_gene_id);

$db->disconnect();

CreateGraphicFile($pathway_id, $source, $format);

######################################################################
sub CreateGraphicFile {
  my ($pid, $source) = @_;

  my %source_map = (
    "NATURE"   => "nature",
    "BioCarta" => "biocarta",
    "Reactome" => "reactome"
  );

  my $map_file;

  my $inf = PREDEF_DIR . "/" . "$source_map{$source}/$pid.dot";  

  open(INF, $inf) or die "cannot open $inf";
  my $cache = new Cache(CACHE_ROOT, PW_CACHE_PREFIX);
  my ($cache_id, $filename) = $cache->MakeCacheFile();
  if ($cache_id != $CACHE_FAIL) {
    open(OUTF, ">$filename.dot");
    while (<INF>) {
      s/[\r\n]+//;
      if (/^ *\/\//) {
        next;
      }
      if (/^ *M_(\d+)/ && ! / -> /) {
        my $mol = $1;
        if (defined $mol2color{$mol}) {
          my $color = $mol2color{$mol};
          if (/fontcolor="([^"]+)"/) {
            s/fontcolor="[^"]+"/fontcolor="$color"/;
          } else {
            s/\];/ fontcolor="$color"];/;
          }
        }
      }
      print OUTF "$_\n";
    }
    close(OUTF);
    chmod 0666, "$filename.dot";
    my $cmd = DOT_PGM .  " $filename.dot -T$format > $filename";
    system($cmd);
    chmod 0666, $filename;
    if ($format eq "gif") {
      my $cmd = DOT_PGM .  " $filename.dot -Tcmap > $filename.map";
      system($cmd);
      chmod 0666, "$filename.map";
      $map_file = "$filename.map"
    }
    unlink("$filename.dot");
  }
  print "<html><head></head><body>\n";
  if ($format eq "svg") {
    print "<embed type=\"image/svg-xml\" height=800 width=1000 " .
       "src=\"http://pid.nci.nih.gov/PathwayGraphic?GID=$cache_id&FORMAT=SVG\">";
  } elsif ($format eq "gif") {
    print "<img " .
       "src=\"http://pid.nci.nih.gov/PathwayGraphic?GID=$cache_id&FORMAT=GIF\">";
  }
}

######################################################################
sub GetMolIdsForEntrezGene {
  my ($db, $pathway_id, $entrez_gene_id) = @_;

  if ($entrez_gene_id =~ /^\s*$/) {
    return [];
  }

  my $sql = qq!/*
                   (mol_id_1, mol_id_2, relation):
identity:          (mol, mol, 'i')
complex component: (component, complex, 'c')
family member:     (family, member, 'm')
protein subunit:   (whole, part, 's')

Have to allow for cases where a component is a family and cases
where a family is a complex with a component that is a family.

When collecting mols for an interaction, find components in a complex
and find members of a family, but do not find families with a member or
complexes with a component.
*/

select unique
  mm_outer_family.mol_id_1
from
  pid.pw_pathway_atom pa,
  pid.pw_edge e,
  pid.pw_mol_mol mm_outer_family,
  pid.pw_mol_mol mm_inner_family,
  pid.pw_mol_mol mm_complex,
  pid.pw_mol_srch s,
  cgap.ll_gene g
where
      mm_outer_family.mol_id_2 = mm_complex.mol_id_2
  and mm_complex.mol_id_1 = mm_inner_family.mol_id_1
  and e.mol_id = mm_outer_family.mol_id_1
  and mm_complex.relation in ('s','c','i')
  and mm_outer_family.relation in ('s','m','i')
  and mm_inner_family.relation in ('s','m','i')
  and e.atom_id = pa.atom_id
  and pa.pathway_id = $pathway_id
  and mm_inner_family.mol_id_2 = s.mol_id
  and s.map_name = '$entrez_gene_id'
  and to_char(g.ll_id) = s.map_name
!;

  my $stm = $db->prepare($sql);
  if(not $stm) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "prepare call failed\n";
    exit;
  }
  if(!$stm->execute()) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "execute call failed\n";
    exit;
  }
  my ($mol, @tmp);
  while (($mol) = $stm->fetchrow_array()) {
    $mol2color{$mol} = $color;
    push @tmp, $mol;
  }
  return \@tmp;
}
