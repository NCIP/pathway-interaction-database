#!/usr/local/bin/perl


# Copyright SRA International
#
# Distributed under the OSI-approved BSD 3-Clause License.
# See http://ncip.github.com/pathway-interaction-database/LICENSE.txt for details.


use strict;
use CGI;
use DBI;

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

my ($db_inst, $db_user, $db_pass, $schema) =
   ("cgprod", "web", "readonly", "pid");

my (@tag, @sql);

push @tag, "NCI-Nature Curated: Number of pathways: ";
push @sql, "select count(distinct pathway_id) from $schema.pw_pathway where pathway_source_id = 5 and subnet = 'N' and pathway_id > 199999";

push @tag, "NCI-Nature Curated: Number of subnets: ";
push @sql, "select count(*) from $schema.pw_pathway where pathway_source_id = 5 and subnet = 'Y'";

push @tag, "NCI-Nature Curated: Number of interactions: ";
push @sql, "select count(*) from $schema.pw_atom where atom_source_id = 5" ;

push @tag, "NCI-Nature Curated: Total number of proteins: ";
push @sql, "select count(*) from $schema.pw_mol where basic_mol_type = 'PR' and mol_id >= 200000 and mol_id < 500000" ;

push @tag, "NCI-Nature Curated: Total number of compounds: ";
push @sql, "select count(*) from $schema.pw_mol where basic_mol_type = 'CM' and mol_id >= 200000 and mol_id < 500000" ;

push @tag, "NCI-Nature Curated: Total number of complexes: ";
push @sql, "select count(*) from $schema.pw_mol where basic_mol_type = 'CX' and mol_id >= 200000 and mol_id < 500000" ;

push @tag, "NCI-Nature Curated: Number of proteins used simply: ";
push @sql, "select count(unique m.mol_id) from $schema.pw_edge e, $schema.pw_mol m, $schema.pw_atom a where atom_source_id = 5 and a.atom_id = e.atom_id and e.mol_id = m.mol_id and m.basic_mol_type = 'PR'" ;

push @tag, "NCI-Nature Curated: Number of proteins used in complexes: ";
push @sql, "select count(unique m.mol_id) from $schema.pw_complex_component c, $schema.pw_mol m where m.mol_id >= 200000 and mol_id < 500000 and m.basic_mol_type = 'PR' and m.mol_id = c.component_mol_id" ;

push @tag, "NCI-Nature Curated: Number of compounds used simply: ";
push @sql, "select count(unique m.mol_id) from $schema.pw_edge e, $schema.pw_mol m, $schema.pw_atom a where atom_source_id = 5 and a.atom_id = e.atom_id and e.mol_id = m.mol_id and m.basic_mol_type = 'CM'" ;

push @tag, "NCI-Nature Curated: Number of compounds used in complexes: ";
push @sql, "select count(unique m.mol_id) from $schema.pw_complex_component c, $schema.pw_mol m where m.mol_id >= 200000 and mol_id < 500000 and m.basic_mol_type = 'CM' and m.mol_id = c.component_mol_id" ;

push @tag, "NCI-Nature Curated: Number of unique PUBMED ids: ";
push @sql, "select distinct count(distinct a.pmid) from $schema.pw_pubmed a, $schema.pw_references b, $schema.pw_pathway_atom c where c.pathway_id between 200000 and 299999 and a.pmid = b.pmid and b.atom_id = c.atom_id" ;

push @tag, "NCI-Nature Curated: Number of total PUBMED ids: ";
push @sql, "select count(a.pmid) from $schema.pw_pubmed a, $schema.pw_references b, $schema.pw_pathway_atom c where c.pathway_id between 200000 and 299999 and a.pmid = b.pmid and b.atom_id = c.atom_id" ;


######################################################################

push @tag, "BioCarta Imported: Number of pathways: ";
push @sql, "select count(*) from $schema.pw_pathway where pathway_source_id in (2,3)";

push @tag, "BioCarta Imported: Number of interactions: ";
push @sql, "select count(*) from $schema.pw_atom where atom_source_id in (2,3)" ;

push @tag, "BioCarta Imported: Total number of proteins: ";
push @sql, "select count(*) from $schema.pw_mol where basic_mol_type = 'PR' and mol_id < 200000" ;

push @tag, "BioCarta Imported: Total number of compounds: ";
push @sql, "select count(*) from $schema.pw_mol where basic_mol_type = 'CM' and mol_id < 200000" ;

push @tag, "BioCarta Imported: Total number of complexes: ";
push @sql, "select count(*) from $schema.pw_mol where basic_mol_type = 'CX' and mol_id < 200000" ;

push @tag, "BioCarta Imported: Number of proteins used simply: ";
push @sql, "select count(unique m.mol_id) from $schema.pw_edge e, $schema.pw_mol m, $schema.pw_atom a where atom_source_id in (2,3) and a.atom_id = e.atom_id and e.mol_id = m.mol_id and m.basic_mol_type = 'PR'" ;

push @tag, "BioCarta Imported: Number of proteins used in complexes: ";
push @sql, "select count(unique m.mol_id) from $schema.pw_complex_component c, $schema.pw_mol m where m.mol_id < 200000 and m.basic_mol_type = 'PR' and m.mol_id = c.component_mol_id" ;

push @tag, "BioCarta Imported: Number of compounds used simply: ";
push @sql, "select count(unique m.mol_id) from $schema.pw_edge e, $schema.pw_mol m, $schema.pw_atom a where atom_source_id in (2,3) and a.atom_id = e.atom_id and e.mol_id = m.mol_id and m.basic_mol_type = 'CM'" ;

push @tag, "BioCarta Imported: Number of compounds used in complexes: ";
push @sql, "select count(unique m.mol_id) from $schema.pw_complex_component c, $schema.pw_mol m where m.mol_id < 200000 and m.basic_mol_type = 'CM' and m.mol_id = c.component_mol_id" ;

######################################################################

push @tag, "Reactome Imported: Number of pathways: ";
push @sql, "select count(*) from $schema.pw_pathway where pathway_source_id = 7 and subnet = 'N'";

push @tag, "Reactome Imported: Number of subnets: ";
push @sql, "select count(*) from $schema.pw_pathway where pathway_source_id = 7 and subnet = 'Y'";

push @tag, "Reactome Imported: Number of interactions: ";
push @sql, "select count(*) from $schema.pw_atom where atom_source_id = 7" ;

push @tag, "Reactome Imported: Total number of proteins: ";
push @sql, "select count(*) from $schema.pw_mol where basic_mol_type = 'PR' and mol_id >= 500000" ;

push @tag, "Reactome Imported: Total number of compounds: ";
push @sql, "select count(*) from $schema.pw_mol where basic_mol_type = 'CM' and mol_id >= 500000" ;

push @tag, "Reactome Imported: Total number of complexes: ";
push @sql, "select count(*) from $schema.pw_mol where basic_mol_type = 'CX' and mol_id >= 500000" ;

push @tag, "Reactome Imported: Number of proteins used simply: ";
push @sql, "select count(unique m.mol_id) from $schema.pw_edge e, $schema.pw_mol m, $schema.pw_atom a where atom_source_id = 7 and a.atom_id = e.atom_id and e.mol_id = m.mol_id and m.basic_mol_type = 'PR'" ;

push @tag, "Reactome Imported: Number of proteins used in complexes: ";
push @sql, "select count(unique m.mol_id) from $schema.pw_complex_component c, $schema.pw_mol m where m.mol_id >= 500000 and m.basic_mol_type = 'PR' and m.mol_id = c.component_mol_id" ;

push @tag, "Reactome Imported: Number of compounds used simply: ";
push @sql, "select count(unique m.mol_id) from $schema.pw_edge e, $schema.pw_mol m, $schema.pw_atom a where atom_source_id = 7 and a.atom_id = e.atom_id and e.mol_id = m.mol_id and m.basic_mol_type = 'CM'" ;

push @tag, "Reactome Imported: Number of compounds used in complexes: ";
push @sql, "select count(unique m.mol_id) from $schema.pw_complex_component c, $schema.pw_mol m where m.mol_id >= 500000 and m.basic_mol_type = 'CM' and m.mol_id = c.component_mol_id" ;

my $query  = new CGI;
print "Content-type: text/plain\n\n";

my $db = DBI->connect("DBI:Oracle:" . $db_inst, $db_user, $db_pass);
if (not $db or $db->err()) {
  print STDERR "Cannot connect to " . $db_user . "@" . $db_inst . "\n";
  exit;
}

my (@lines);
push @lines, "<pre>";
for (my $i = 0; $i < @sql; $i++) {
  push @lines, $tag[$i] . " " . DoSql($db, $sql[$i]);
}
push @lines, "</pre>";

$db->disconnect();

print join("\n", @lines) . "\n";

######################################################################
sub DoSql {
  my ($db, $sql) = @_;

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
  my ($n) = $stm->fetchrow_array();
  $stm->finish();
  return $n;
}


