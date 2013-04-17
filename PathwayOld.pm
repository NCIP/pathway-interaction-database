#!/usr/local/bin/perl


# Copyright SRA International
#
# Distributed under the OSI-approved BSD 3-Clause License.
# See http://ncip.github.com/pathway-interaction-database/LICENSE.txt for details.


use DBI;
use Cache;
use strict;

use constant ORACLE_LIST_LIMIT => 500;
use constant MAX_LONG_LEN       => 16384;

use constant CACHE_ROOT      => "/share/content/PID/data/cache" ;
use constant PW_CACHE_PREFIX => "PW" ;
#use constant DOT_PGM         => "/usr/local/bin/dot";      ## on dev
use constant DOT_PGM         => "/share/content/PID/run/dot"; ## on stage, prod
use constant PREDEF_DIR      => "/share/content/PID/data/predefined";

my $PID_SITE = "http://pid.nci.nih.gov";

my $color_a        = "blue";
my $color_b        = "red";
my $color_ab       = "purple";

######################################################################
sub r_numerically { $b <=> $a; }

######################################################################
sub numerically { $a <=> $b; }

######################################################################
sub HTML_HEAD {
  return qq!
<head>
<STYLE TYPE=\"text/css\">
BODY {font-family: Arial}
</STYLE>
</head>
!;
}

######################################################################
sub GetPathwayHitsForEntrezGenes {
  my ($db, $schema, $genes, $pd, $gd) = @_;

  my $list;
  my @genes = keys %{ $genes };

  my ($source_name, $pathway_id, $ext_pathway_id, $pathway_name,
        $entrez_gene_id);

  for(my $i = 0; $i < @genes; $i += ORACLE_LIST_LIMIT) {
 
    if(($i + ORACLE_LIST_LIMIT - 1) < @genes) {
      $list = join(",", @genes[$i..$i+ORACLE_LIST_LIMIT-1]);
    }
    else {
      $list = join(",", @genes[$i..@genes-1]);
    }
    my $sql = qq!
select unique
  c.source_name,
  p.pathway_id,
  p.ext_pathway_id,
  p.pathway_name,
  pg.entrez_gene_id
from
  $schema.pw_pathway_gene pg,
  $schema.pw_pathway p,
  $schema.pw_source c
where
      p.pathway_id = pg.pathway_id
  and c.source_id = p.pathway_source_id
  and pg.entrez_gene_id in ($list)
!;

    my $stm = $db->prepare($sql);
    if(not $stm) {
      print STDERR "$sql\n";
      print STDERR "$DBI::errstr\n";
      die "prepare call failed\n";
    }
    if(!$stm->execute()) {
      print STDERR "$sql\n";
      print STDERR "$DBI::errstr\n";
      die "execute call failed\n";
    }
    while (($source_name, $pathway_id, $ext_pathway_id, $pathway_name,
        $entrez_gene_id) = $stm->fetchrow_array()) {
      $pd->{gene}{$pathway_id}{$entrez_gene_id} = 1;
      $pd->{source}{$pathway_id}     = $source_name;
      $pd->{ext_id}{$pathway_id}     = $ext_pathway_id;
      $pd->{name}{$pathway_id}       = $pathway_name;
      $gd->{pathway}{$entrez_gene_id}{$pathway_id} = 1;

    }
  }
}

######################################################################
sub GetTotalGenesInPathways {
  my ($db, $schema) = @_;

  my $sql = qq!
select
  count (unique pg.entrez_gene_id)
from
  $schema.pw_pathway_gene pg
!;
  my $stm = $db->prepare($sql);
  if(not $stm) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    die "prepare call failed\n";
  }
  if(!$stm->execute()) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    die "execute call failed\n";
  }
  my ($total_genes_in_pathways) = $stm->fetchrow_array();
  $stm->finish();
  return $total_genes_in_pathways;
}

######################################################################
sub MakeSymbolList {
  my ($list, $gd) = @_;

  my @symbols;
  for my $g (@{ $list }) {
    if (defined $gd->{symbol}{$g}{OF}) {
      push @symbols, $gd->{symbol}{$g}{OF};
    } else {
      for my $t (keys %{ $gd->{symbol}{$g} }) {
        push @symbols, $gd->{symbol}{$g}{$t};
        last; ##  just take one
      }
    }
  }
  return \@symbols;
}

######################################################################
sub CountGenesInPathways {
  my ($db, $schema, $pd) = @_;

  my $list;
  my ($pathway_id, $count);
  my @pids = keys %{ $pd->{gene} };

  for(my $i = 0; $i < @pids; $i += ORACLE_LIST_LIMIT) {
 
    if(($i + ORACLE_LIST_LIMIT - 1) < @pids) {
      $list = join(",", @pids[$i..$i+ORACLE_LIST_LIMIT-1]);
    }
    else {
      $list = join(",", @pids[$i..@pids-1]);
    }
    my $sql = qq!
select unique
  pg.pathway_id,
  count (unique pg.entrez_gene_id)
from
  $schema.pw_pathway_gene pg
where
      pg.pathway_id in ($list)
group by
  pg.pathway_id
!;
    my $stm = $db->prepare($sql);
    if(not $stm) {
      print STDERR "$sql\n";
      print STDERR "$DBI::errstr\n";
      die "prepare call failed\n";
    }
    if(!$stm->execute()) {
      print STDERR "$sql\n";
      print STDERR "$DBI::errstr\n";
      die "execute call failed\n";
    }
    while (($pathway_id, $count) = $stm->fetchrow_array()) {
      $pd->{genecount}{$pathway_id} = $count;
    }
  }
}

######################################################################
sub RankPathwaysByGeneHits {
  my ($db, $schema, $source, $format, $terms_a, $terms_b) = @_;

  my %terms;
  for my $t (keys %{ $terms_a }) {
    $terms{$t} = 1;
  }
  for my $t (keys %{ $terms_b }) {
    $terms{$t} = 1;
  }
  my @terms = keys %terms;
  if (@terms == 0) {
    return;
  }

  my ($sql, $stm, $list);
  my ($source_name, $pathway_id, $ext_pathway_id, $pathway_name,
      $entrez_gene_id, $count, $symbol, $symbol_type);

  my %pd;   ## pathway data
  my %gd;   ## gene data

  my $total_genes_in_pathways = GetTotalGenesInPathways($db, $schema);

  my ($genes_a, $genes_b) = GetEntrezGeneIds($db, $schema, $terms_a, $terms_b, \%gd);
  GetPathwayHitsForEntrezGenes($db, $schema, $genes_a, \%pd, \%gd);
  GetPathwayHitsForEntrezGenes($db, $schema, $genes_b, \%pd, \%gd);

  my @genes = keys %{ $gd{pathway} };
  my $n_query = scalar(@genes);

  CountGenesInPathways($db, $schema, \%pd);

  my (%tmp, %pval_cache);
  for my $pid (keys %{ $pd{gene} }) {
    my (@alist, @blist);
    my $hits_target = 0;
    for my $g (keys %{ $pd{gene}{$pid} }) {
      if (defined $$genes_a{$g}) {
        push @alist, $g;
      }
      if (defined $$genes_b{$g}) {
        push @blist, $g;
      }
    }
    for my $g (@genes) {
      if (defined $pd{gene}{$pid}{$g}) {
        $hits_target++;
      }
    }

    my $pval;
    my $good = $pd{genecount}{$pid};
    if (defined $pval_cache{"$good,$hits_target"}) {
      $pval = $pval_cache{"$good,$hits_target"};
    } else {
      $pval = hypergeom(
          $good,
          $total_genes_in_pathways - $good,
          $n_query,
          $hits_target
        );
      $pval = sprintf("%.2e", $pval);
      $pval_cache{"$good,$hits_target"} = $pval;
    }

## remove any possible duplicate gene ids (from multiple symbols)
    my %cull;
    for my $g (@alist) {
      $cull{$g} = 1;
    }
    undef @alist;
    @alist = keys %cull;
    my %cull;
    for my $g (@blist) {
      $cull{$g} = 1;
    }
    undef @blist;
    @blist = keys %cull;
## make symbol lists

    my @a_syms = @{ MakeSymbolList(\@alist, \%gd) };    
    my @b_syms = @{ MakeSymbolList(\@blist, \%gd) };    

    if ($format eq "text") {
      push @{ $tmp{$pval} }, join("\t",
        $pval,
        $good,
        $total_genes_in_pathways - $good,
        $n_query,
        $hits_target,
        join(",", @alist),
        join(",", @a_syms),
        join(",", @blist),
        join(",", @b_syms),
        $pid,
        $pd{source_name}{$pid},
        $pd{ext_id}{$pid},
        $pd{name}{$pid}
      );
    } elsif ($format eq "html") {
      my $a_syms;
      if (@a_syms == 0) {
        $a_syms = "&nbsp;";
      } else {
        $a_syms = join(", ", sort @a_syms);
      }
      my $b_syms;
      if (@b_syms == 0) {
        $b_syms = "&nbsp;";
      } else {
        $b_syms = join(", ", sort @b_syms);
      }
      push @{ $tmp{$pval} }, join("\t",
        "<tr>",
        "<td>$pval</td>",
        "<td>$a_syms</td>",
        "<td>$b_syms</td>",
        "<td><a href=\"pathway_landing.shtml?" .
          "what=graphic\&" .
          "pathway_id=$pid\&" .
          "gif=on\&" .
          "source=$pd{source}{$pid}\&" .
          "genes_a=" . join(",", @alist) . "\&" .
          "genes_b=" . join(",", @blist) .
          "\" target=_>$pd{name}{$pid}</a></td>",
        "</tr>"
      );
    }
  }
  my @rows;
  push @rows, "<h1 class=\"pagetitle\">Batch query results for $source </h1>";

  for my $pval (sort numerically keys %tmp) {
    for my $r (@{ $tmp{$pval} }) {
      push @rows, $r;
    }
  }
  if ($format eq "html") {
    unshift @rows, "<table border=1 cellspacing=1 cellpadding=4>\n" .
        "<tr><td>P-value</td><td>Genes (A)</td><td>Genes (B)</td><td>Pathway Name</td></tr>";   
    push @rows, "</table>";
  } elsif ($format eq "text") {
    unshift @rows, join("\t",
      "P-value",
      "N genes in pathway",
      "N genes NOT in pathway",
      "N genes in Query",
      "N Hits",
      "Gene ids in A list",
      "Gene symbols in A list",
      "Gene ids in B list",
      "Gene symbols in B list",
      "Pathway id",
      "Pathway source",
      "Pathway ext id",
      "Pathway name");
  }
  return join("\n", @rows, "");
}

#####################################################################
sub GetEntrezGeneIds {
  my ($db, $schema, $terms_a, $terms_b, $gd) = @_;

  ## limit the search to Entrez Gene ids that are in some pathway

  my %terms;
  my (%genes_a, %genes_b);

  for my $t (keys %{ $terms_a }, keys %{ $terms_b }) {
    $terms{$t} = 1;
  }

  my (@numeric_terms, @non_numeric_terms);

  for my $t (keys %terms) {
    if ($t =~ /^\d+$/) {
      push @numeric_terms, $t;
    } else {
      push @non_numeric_terms, $t;
    }
  }

  ## do numeric

  my $list;
  for(my $i = 0; $i < @numeric_terms; $i += ORACLE_LIST_LIMIT) {
 
    if(($i + ORACLE_LIST_LIMIT - 1) < @numeric_terms) {
      $list = join(",", @numeric_terms[$i..$i+ORACLE_LIST_LIMIT-1]);
    }
    else {
      $list = join(",", @numeric_terms[$i..@numeric_terms-1]);
    }
    my $sql = qq!
select unique
  pg.entrez_gene_id,
  lower(e.symbol),
  e.symbol,
  e.symbol_type
from
  $schema.pw_pathway_gene pg,
  cgap.ll_gene e
where
      pg.entrez_gene_id = e.ll_id
  and e.organism = 'Hs'
  and e.ll_id in ($list)
!;
    my $stm = $db->prepare($sql);
    if(not $stm) {
      print STDERR "$sql\n";
      print STDERR "$DBI::errstr\n";
      die "prepare call failed\n";
    }
    if(!$stm->execute()) {
      print STDERR "$sql\n";
      print STDERR "$DBI::errstr\n";
      die "execute call failed\n";
    }
    my ($entrez_gene_id, $lc_symbol, $symbol, $symbol_type);
    while (($entrez_gene_id, $lc_symbol, $symbol, $symbol_type) =
        $stm->fetchrow_array()) {
      $gd->{symbol}{$entrez_gene_id}{$symbol_type} = $symbol;
      if (defined $$terms_a{$entrez_gene_id}) {
        $genes_a{$entrez_gene_id} = 1;
      }
      if (defined $$terms_b{$entrez_gene_id}) {
        $genes_b{$entrez_gene_id} = 1;
      }
    }
  }

  ## do non_numeric, gene symbols

  my $list;
  for(my $i = 0; $i < @non_numeric_terms; $i += ORACLE_LIST_LIMIT) {
 
    if(($i + ORACLE_LIST_LIMIT - 1) < @non_numeric_terms) {
      $list = lc("'" . join("','", @non_numeric_terms[$i..$i+ORACLE_LIST_LIMIT-1]) . "'");
    }
    else {
      $list = lc("'" . join("','", @non_numeric_terms[$i..@non_numeric_terms-1]) . "'");
    }
    my $sql = qq!
select unique
  pg.entrez_gene_id,
  lower(e.symbol),
  e.symbol,
  e.symbol_type
from
  $schema.pw_pathway_gene pg,
  cgap.ll_gene e
where
      pg.entrez_gene_id = e.ll_id
  and pg.pathway_id < 500000
  and e.organism = 'Hs'
  and lower(e.symbol) in ($list)
!;
    my $stm = $db->prepare($sql);
    if(not $stm) {
      print STDERR "$sql\n";
      print STDERR "$DBI::errstr\n";
      die "prepare call failed\n";
    }
    if(!$stm->execute()) {
      print STDERR "$sql\n";
      print STDERR "$DBI::errstr\n";
      die "execute call failed\n";
    }
    my ($entrez_gene_id, $lc_symbol, $symbol, $symbol_type);
    while (($entrez_gene_id, $lc_symbol, $symbol, $symbol_type) =
        $stm->fetchrow_array()) {
      $gd->{symbol}{$entrez_gene_id}{$symbol_type} = $symbol;
      if (defined $$terms_a{$lc_symbol}) {
        $genes_a{$entrez_gene_id} = 1;
      }
      if (defined $$terms_b{$lc_symbol}) {
        $genes_b{$entrez_gene_id} = 1;
      }
    }
  }

  ## do non_numeric, UniProt ids (we still have significantly more
  ## UP-to-EG associations in ll2sp than in ll2acc)

  my $list;
  for(my $i = 0; $i < @non_numeric_terms; $i += ORACLE_LIST_LIMIT) {
 
    if(($i + ORACLE_LIST_LIMIT - 1) < @non_numeric_terms) {
      $list = lc("'" . join("','", @non_numeric_terms[$i..$i+ORACLE_LIST_LIMIT-1]) . "'");
    }
    else {
      $list = lc("'" . join("','", @non_numeric_terms[$i..@non_numeric_terms-1]) . "'");
    }
    my $sql = qq!
select unique
  pg.entrez_gene_id,
  lower(e.symbol),
  e.symbol,
  e.symbol_type
from
  $schema.pw_pathway_gene pg,
  cgap.ll_gene e,
  cgap.ll2sp s,
  cgap.sp_primary p
where
      pg.entrez_gene_id = e.ll_id
  and pg.pathway_id < 500000
  and e.organism = 'Hs'
  and e.ll_id = s.ll_id
  and s.sp_primary = p.sp_primary
  and lower(p.sp_id_or_secondary) in ($list)
!;
    my $stm = $db->prepare($sql);
    if(not $stm) {
      print STDERR "$sql\n";
      print STDERR "$DBI::errstr\n";
      die "prepare call failed\n";
    }
    if(!$stm->execute()) {
      print STDERR "$sql\n";
      print STDERR "$DBI::errstr\n";
      die "execute call failed\n";
    }
    my ($entrez_gene_id, $lc_symbol, $symbol, $symbol_type);
    while (($entrez_gene_id, $lc_symbol, $symbol, $symbol_type) =
        $stm->fetchrow_array()) {
      $gd->{symbol}{$entrez_gene_id}{$symbol_type} = $symbol;
      if (defined $$terms_a{$lc_symbol}) {
        $genes_a{$entrez_gene_id} = 1;
      }
      if (defined $$terms_b{$lc_symbol}) {
        $genes_b{$entrez_gene_id} = 1;
      }
    }
  }

  return (\%genes_a, \%genes_b);
}

#####################################################################
sub oldGetEntrezGeneIds {
  my ($db, $schema, $terms_a, $terms_b, $gd) = @_;

  ## limit the search to Entrez Gene ids that are in some pathway

  my %terms;
  my (%genes_a, %genes_b);

  for my $t (keys %{ $terms_a }, keys %{ $terms_b }) {
    $terms{$t} = 1;
  }
  my @terms = keys %terms;

  my $list;
  for(my $i = 0; $i < @terms; $i += ORACLE_LIST_LIMIT) {
 
    if(($i + ORACLE_LIST_LIMIT - 1) < @terms) {
      $list = lc("'" . join("','", @terms[$i..$i+ORACLE_LIST_LIMIT-1]) . "'");
    }
    else {
      $list = lc("'" . join("','", @terms[$i..@terms-1]) . "'");
    }
    my $sql = qq!
select unique
  m1.map_name,
  pg.entrez_gene_id,
  e.symbol,
  e.symbol_type
from
  $schema.pw_pathway_gene pg,
  cgap.ll_gene e,
  $schema.pw_mol_srch m1,
  $schema.pw_mol_srch m2
where
      pg.entrez_gene_id = e.ll_id
  and pg.pathway_id < 500000
  and e.organism = 'Hs'
  and to_char(pg.entrez_gene_id) = m2.map_name
  and m1.mol_id not in (select family_mol_id from $schema.pw_family_member)
  and m1.mol_id = m2.mol_id
  and m1.map_name in ($list)
!;

    my $stm = $db->prepare($sql);
    if(not $stm) {
      print STDERR "$sql\n";
      print STDERR "$DBI::errstr\n";
      die "prepare call failed\n";
    }
    if(!$stm->execute()) {
      print STDERR "$sql\n";
      print STDERR "$DBI::errstr\n";
      die "execute call failed\n";
    }
    my ($term, $entrez_gene_id, $symbol, $symbol_type);
    while (($term, $entrez_gene_id, $symbol, $symbol_type) =
        $stm->fetchrow_array()) {
      $gd->{symbol}{$entrez_gene_id}{$symbol_type} = $symbol;
      if (defined $$terms_a{$term}) {
        $genes_a{$entrez_gene_id} = 1;
      }
      if (defined $$terms_b{$term}) {
        $genes_b{$entrez_gene_id} = 1;
      }
    }
  }
  return (\%genes_a, \%genes_b);
}

#####################################################################
sub CreateGraphicFile {
  my ($db, $schema, $pid, $source, $mols_a, $mols_b, $format) = @_;

  my %source_map = (
    "NATURE"   => "nature",
    "BioCarta" => "biocarta",
    "Reactome" => "reactome"
  );

  my %nice_source = (
    "NATURE"   => "NCI-Nature Curated",
    "BioCarta" => "BioCarta",
    "Reactome" => "Reactome"
  );

  my $sql = qq!
select
  p.pathway_name
from
  $schema.pw_pathway p
where
      p.pathway_id = $pid
!;

  my $stm = $db->prepare($sql);
  if(not $stm) {
    print STDERR "$sql\n";
    exit;
  }
  if(!$stm->execute()) {
    print STDERR "$sql\n";
    exit;
  }

  my ($pathway_name) = $stm->fetchrow_array();
  $stm->finish();

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
        if (defined $$mols_a{$mol}) {
          if (defined $$mols_b{$mol}) {
            if (/fontcolor="([^"]+)"/) {
              s/fontcolor="[^"]+"/fontcolor="$color_ab"/;
            } else {
              s/\];/ fontcolor="$color_ab"];/;
            }
          } else {
            if (/fontcolor="([^"]+)"/) {
              s/fontcolor="[^"]+"/fontcolor="$color_a"/;
            } else {
              s/\];/ fontcolor="$color_a"];/;
            }
          }
        } elsif (defined $$mols_b{$mol}) {
          if (/fontcolor="([^"]+)"/) {
            s/fontcolor="[^"]+"/fontcolor="$color_b"/;
          } else {
            s/\];/ fontcolor="$color_b"];/;
          }
        }
      }
      print OUTF "$_\n";
    }
    close(OUTF);
    chmod 0666, "$filename.dot";
    if ($format ne 'svg') {
      $format = 'gif';
    }
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
  my @tmp;
  push @tmp, "<html>";
  push @tmp, HTML_HEAD();
  push @tmp, "<body>";
  push @tmp, "<p>";

  if ($format eq "svg") {
    push @tmp, 
        "<embed type=\"image/svg-xml\" height=800 width=1000 " .
         "src=\"$PID_SITE/PathwayGraphic?" .
         "GID=$cache_id&FORMAT=SVG\">";
  } elsif ($format eq "gif") {
    push @tmp, "<img " .
        "src=\"$PID_SITE/PathwayGraphic?" .
        "GID=$cache_id&FORMAT=GIF\" usemap=#map_$cache_id>";
    push @tmp, "<map name=map_$cache_id>";
    open(MAP, $map_file) or die "cannot open $map_file";
    while (<MAP>) {
      chomp;
      push @tmp, "$_";
    } 
    close(MAP);
 
    push @tmp, "</map>";
    # Need to add map information here
  }
  return join("\n", @tmp, "");
}

######################################################################
sub GetMolIdsForEntrezGene {
  my ($db, $pathway_id, $genes) = @_;

  my @genes = keys %{ $genes };
  if (@genes == 0) {
    return {};
  }

  my ($mol, %mols, $list);

  for(my $i = 0; $i < @genes; $i += ORACLE_LIST_LIMIT) {
 
    if(($i + ORACLE_LIST_LIMIT - 1) < @genes) {
      $list = join("','", @genes[$i..$i+ORACLE_LIST_LIMIT-1]);
    }
    else {
      $list = join("','", @genes[$i..@genes-1]);
    }
    $list = "'" . $list . "'";

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
  and s.map_name in ($list)
  and to_char(g.ll_id) = s.map_name
!;

    my $stm = $db->prepare($sql);
    if(not $stm) {
      print STDERR "$sql\n";
      exit;
    }
    if(!$stm->execute()) {
      print STDERR "$sql\n";
      exit;
    }
    while (($mol) = $stm->fetchrow_array()) {
      $mols{$mol} = 1;
    }
  }
  return \%mols;
}


######################################################################
sub logfact {
   return gammln(shift(@_) + 1.0);
}

sub hypergeom {
   # There are m "bad" and n "good" balls in an urn.
   # Pick N of them. The probability of i or more successful selections:
   # (m!n!N!(m+n-N)!)/(i!(n-i)!(m+i-N)!(N-i)!(m+n)!)
   my ($n, $m, $N, $i) = @_;

   my $loghyp1 = logfact($m)+logfact($n)+logfact($N)+logfact($m+$n-$N);
   my $loghyp2 = logfact($i)+logfact($n-$i)+logfact($m+$i-$N)+logfact(+$N-$i)+logfact($m+$n);
   return exp($loghyp1 - $loghyp2);
}

sub gammln {
  my $xx = shift;
  my @cof = (76.18009172947146, -86.50532032941677,
             24.01409824083091, -1.231739572450155,
             0.12086509738661e-2, -0.5395239384953e-5);
  my $y = my $x = $xx;
  my $tmp = $x + 5.5;
  $tmp -= ($x + .5) * log($tmp);
  my $ser = 1.000000000190015;
  for my $j (0..5) {
     $ser += $cof[$j]/++$y;
  }
  -$tmp + log(2.5066282746310005*$ser/$x);
}

######################################################################
1;
######################################################################
