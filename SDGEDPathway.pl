#!/usr/local/bin/perl



# Copyright SRA International
#
# Distributed under the OSI-approved BSD 3-Clause License.
# See http://ncip.github.com/pathway-interaction-database/LICENSE.txt for details.



######################################################################
# 
#
######################################################################

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
use DBI;
use Bayesian;
use CGI;

# constant determining whether to use Oracle or flat file
# as data source
use constant MANY_SAGE_LIBRARIES => 20;
use constant MINIMUM_SEQUENCES => 40000;
use constant ORACLE_LIST_LIMIT => 500;

my $CGAP_SCHEMA = "cgap2";
use constant SAGEFREQ => "/share/content/CGAP/data/sagefreq.dat";
use constant DB_USER         => "web";
use constant DB_PASS         => "readonly";
use constant DB_INSTANCE     => "lpgprod";

my @CCF = (    ## common co-factors
    'water',
    'atp',
    'nad',
    'nadh',
    'nadph',
    'nadp',
    'oxygen',
    'adp',
    'orthophosphate',
    'coa',
    'carbon dioxide',
    'diphosphate',
    'ammonia',
    'h+'
  );

#my (
#  $org,
#  $factor,
#  $pvalue,
#  $tissue_a, $histology_a,
#  $tissue_b, $histology_b  
#) = @ARGV;

my (
  $org,
  $factor,
  $pvalue,
  $tissue_a, $histology_a,
  $tissue_b, $histology_b  
);

my $query = new CGI;
$org         = $query->param("ORG");
$factor      = $query->param("FACTOR");
$pvalue      = $query->param("PVAL");
$tissue_a    = $query->param("TISSUE_A");
$tissue_b    = $query->param("TISSUE_B");
$histology_a = $query->param("HIST_A");
$histology_b = $query->param("HIST_B");

print "Content-type: text/plain\n\n";
print SDGEDPathway_1($org, $factor, $pvalue,
      $tissue_a, $histology_a, $tissue_b, $histology_b);


######################################################################
## GXS variables

my (%setA, %setB, %seqA, %seqB, %libA, %libB, %tallied_clus);

######################################################################
sub Init {
  undef %setA;
  undef %setB;
  undef %seqA;
  undef %seqB;
  undef %libA;
  undef %libB;
  undef %tallied_clus;
}

######################################################################
sub MayBeDifferent {
  my ($factor, $a, $b, $A, $B) = @_;

  my ($big, $small, $a_ratio, $b_ratio);
  $a_ratio = $a/$A;
  $b_ratio = $b/$B;

  if ($a_ratio == $b_ratio) {
    return 0;
  } elsif ($factor == 1) {
    return 1;
  } else {
    if ($a_ratio > $b_ratio) {
      $big   = $a_ratio;
      $small = $b_ratio;
    } else {
      $big   = $b_ratio;
      $small = $a_ratio;
    }
    if ($big > $factor * $small) {
      return 1;
    } else {
      return 0;
    }
  }
}

######################################################################
sub OddsRatio {
  my ($G_A, $G_B, $TotalA, $TotalB) = @_;
  ## G_A: number of ESTs in Set A that hit gene G
  ## TotalA: number of ESTs in Set A
  ## BarG_A: number of ESTS in Set A that do not hit gene G

  my $BarG_A = $TotalA - $G_A;
  my $BarG_B = $TotalB - $G_B;

  if ($G_A == 0) {
    return 0;
  } elsif ($G_B == 0) {
    return "NaN";
  } else {
    return sprintf("%.2f", ($G_A * $BarG_B) / ($G_B * $BarG_A));
  }
  
}

######################################################################
sub numerically { $a <=> $b };
sub r_numerically { $b <=> $a };

######################################################################
# SAGE
######################################################################

######################################################################
sub GetSymsOfCIDs {
  my ($db, $cids_hit, $cid2sym, $cid2title, $cid2loc) = @_;

  my ($sql, $stm);
  my ($cluster_number, $gene, $title, $loc);
  my ($i, $list);
  my @cids = keys %{ $cids_hit };

  if (@cids == 0) {
    return;
  }

  for($i = 0; $i < @cids; $i += ORACLE_LIST_LIMIT) {
 
    if(($i + ORACLE_LIST_LIMIT - 1) < @cids) {
      $list = join(",", @cids[$i..$i+ORACLE_LIST_LIMIT-1]);
    }
    else {
      $list = join(",", @cids[$i..@cids-1]);
    }

    $sql = "select cluster_number, gene, description, locuslink " . 
        "from $CGAP_SCHEMA.hs_cluster " .
        "where cluster_number in ($list)";

    $stm = $db->prepare($sql);
    if (not $stm) {
      print STDERR "$sql\n";
      print STDERR "$DBI::errstr\n";
      print STDERR "prepare call failed\n";
      return;
    }
    if (not $stm->execute()) {
      print STDERR "$sql\n";
      print STDERR "$DBI::errstr\n";
      print STDERR "execute call failed\n";
      return;
    }
    while (($cluster_number, $gene, $title, $loc) =
        $stm->fetchrow_array()) {
      if ($gene) {
        $$cid2sym{$cluster_number} = $gene;
      }
      if ($title) {
        $$cid2title{$cluster_number} = $title;
      }
      if ($loc) {
        $$cid2loc{$cluster_number} = $loc;
      }
    }
  }
}

######################################################################
sub GetBestMapsOfTags {
  my ($db, $tags, $tag2cid, $cid2title, $cid2sym, $cid2loc) = @_;

  my ($sql, $stm);
  my ($i, $list, $tag, $cluster_number, $loc,
      %tags_seen, @missing_cids, %cids_hit);

  for($i = 0; $i < @{ $tags }; $i += ORACLE_LIST_LIMIT) {
    if(($i + ORACLE_LIST_LIMIT - 1) < @{ $tags }) {
      $list = join("','", @{ $tags }[$i..$i+ORACLE_LIST_LIMIT-1]);
    }
    else {
      $list = join("','", @{ $tags }[$i..@{ $tags }-1]);
    }

    $sql = "select s.tag, s.cluster_number " .
        "from $CGAP_SCHEMA.sagebest_tag2clu s " .
        "where tag in ('$list')";

    $stm = $db->prepare($sql);
    if(not $stm) {
      print STDERR "$sql\n";
      print STDERR "$DBI::errstr\n";
      print STDERR "prepare call failed\n";
    }
    if(!$stm->execute()) {
      print STDERR "$sql\n";
      print STDERR "$DBI::errstr\n";
      print STDERR "execute call failed\n";
    }

    while(($tag, $cluster_number) = $stm->fetchrow_array()) {
      $tags_seen{$tag} = 1;
      $$tag2cid{$tag} = $cluster_number;
      $cids_hit{$cluster_number} = 1;
    }
  }

  GetSymsOfCIDs($db, \%cids_hit, $cid2sym, $cid2title, $cid2loc);

  for $tag (@{ $tags }) {
    if (not defined $tags_seen{$tag}) {
      push @missing_cids, $tag;
    }
  }

}

######################################################################
sub SummarizeTagFreqs {
  my ($db, $lid_list, $tagfreqs, $totalfreqs,
      $tags_hit) = @_;

  my ($sql, $stm);
  my ($tag, $tagfreq, $libcount);

  $sql = "select " .
      "tag, sum(frequency), count(sage_library_id) ".
      "from $CGAP_SCHEMA.sagefreq " .
      "where sage_library_id in ($lid_list) " .
      "group by tag";

  $stm = $db->prepare($sql);
  if(not $stm) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "prepare call failed\n";
  }
  if(!$stm->execute()) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "execute call failed\n";
  }
  while(($tag, $tagfreq, $libcount) = $stm->fetchrow_array()) {
    $$tags_hit{$tag} = 1;
    $$tagfreqs{$tag} = $tagfreq;
    $$totalfreqs += $tagfreq;
  }

}

######################################################################
sub Readsagefreqdat {
  my ($setA, $tagfreqsA, $totalfreqsA,
      $setB, $tagfreqsB, $totalfreqsB,
      $tags_hit) = @_;

  my ($sql, $stm);
  my ($tag, $tagfreq);

  open (SAGEFREQIN, SAGEFREQ)  or die "Cannot open " . SAGEFREQ;
  while(<SAGEFREQIN>) {
    chop;
    my ($tag, $lid, $freq) = split "\t", $_;
    if ($$setA{$lid} == 1) {
      $$tags_hit{$tag} = 1;
      $$tagfreqsA{$tag} += $freq;
      $$totalfreqsA += $freq;
    }
    if ($$setB{$lid} == 1) {
      $$tags_hit{$tag} = 1;
      $$tagfreqsB{$tag} += $freq;
      $$totalfreqsB += $freq;
    }
  }  
  close(SAGEFREQIN);

}

######################################################################
sub  GetLLsUsedInPW {
  my ($db, $used_lls) = @_;

  my ($stm, $sql);
  my ($loc);

  my $CGAP_SCHEMA = "cgap2";

  $sql = qq!
select
  e.ext_mol_id
from
  $CGAP_SCHEMA.pw_ext_mol_id e
where
  e.id_type = 'LL'
  !;
  $stm = $db->prepare($sql);
  if (not $stm) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "prepare call failed\n";
    return;
  }
  if (not $stm->execute()) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "execute call failed\n";
    return;
  }
  while (($loc) = $stm->fetchrow_array()) {
    $$used_lls{$loc} = 1;
  }
}

######################################################################
sub SAGEComputeDifferences {
  my ($db, $factor, $pvalue, $tags_hit,
        $total_seqsA, $total_seqsB,
        $freqsA,     $freqsB, 
        $result
    ) = @_;

  my %exists;
 
  if ($factor < 1) {
    $factor = 1;
  }

##  if ($pvalue =~ /^e/) {
##      $pvalue = "1" . $pvalue;
##  } 

  my ($tag);
  my ($odds_ratio, $P);
  my ($G_A, $G_B);
  my (%order, %NaN);
  my (@tags_hit, %tag2cid, %tag2title, %cid2sym, %cid2loc,
      %tagdat, %used_lls);
  my ($cid, $sym);

  my (@HI, @LO);

  for $tag (keys %{ $tags_hit }) {

    $G_A = $$freqsA{$tag}; $G_A or $G_A = 0;
    $G_B = $$freqsB{$tag}; $G_B or $G_B = 0;

    if (not MayBeDifferent($factor,
        $G_A, $G_B, $total_seqsA, $total_seqsB)) {
      next;
    }

    if ($G_A or $G_B) {
      $odds_ratio =
          OddsRatio($G_A, $G_B, $total_seqsA, $total_seqsB);

      if( defined $exists{"$G_A,$G_B"} ) {
        $P = $exists{"$G_A,$G_B"};
      }
      else {
 
        if ($G_A/$total_seqsA > $G_B/$total_seqsB) {
          $P = 1 - Bayesian::Bayesian($factor,
                   $G_A, $G_B, $total_seqsA, $total_seqsB);
          $exists{"$G_A,$G_B"} = $P;
        } else {
          $P = 1 - Bayesian::Bayesian($factor,
                   $G_B, $G_A, $total_seqsB, $total_seqsA);
          $exists{"$G_A,$G_B"} = $P;
        }
 
      }
      if ($P <= $pvalue) {
        push @tags_hit, $tag;
        $P = sprintf "%.2f", $P;
        $tagdat{$tag} =
            "A = $G_A/$total_seqsA; B = $G_B/$total_seqsB; P = $P";
        if ($odds_ratio eq "NaN") {
          push @HI, $tag;
        } elsif ($odds_ratio > 1) {
          push @HI, $tag;
        } else {
          push @LO, $tag;
        }
      }
    }
  }

  GetBestMapsOfTags($db, \@tags_hit, \%tag2cid, \%tag2title, \%cid2sym,
      \%cid2loc);

  GetLLsUsedInPW($db, \%used_lls);

  my ($loc);

  if (@HI > 0 || @LO > 0) {
    for my $c (@CCF) {
      push @{ $result }, "prune_mol_name\t$c";
    }
  }

  for $tag (@HI) {
    if (defined $tag2cid{$tag}) {
      $cid = $tag2cid{$tag};
      if (defined $cid2loc{$cid}) {
        $loc = $cid2loc{$cid};
        if (defined $used_lls{$loc}) {
          push @{ $result }, "#$tag\tHs.$cid\t$cid2sym{$cid}\t" .
            "LL $loc\t$tagdat{$tag}";
          push @{ $result }, "mol_ext_id\tLL\t$loc";
          push @{ $result }, "value\t9\tmol_ext_id\tLL\t$loc";
          push @{ $result }, "connect_mol_ext_id\tLL\t$loc";
        } else {
          push @{ $result }, "#$tag\tHs.$cid\t$cid2sym{$cid}\t" .
            "LL $loc\t$tagdat{$tag}";
          push @{ $result }, "#NO PATHWAY FOR: mol_ext_id\tLL\t$loc";
          push @{ $result }, "#NO PATHWAY FOR: value\t9\tmol_ext_id\tLL\t$loc";
        }
      }
    }
  }

  for $tag (@LO) {
    if (defined $tag2cid{$tag}) {
      $cid = $tag2cid{$tag};
      if (defined $cid2loc{$cid}) {
        $loc = $cid2loc{$cid};
        if (defined $used_lls{$loc}) {
          push @{ $result }, "#$tag\tHs.$cid\t$cid2sym{$cid}\t" .
            "LL $loc\t$tagdat{$tag}";
          push @{ $result }, "mol_ext_id\tLL\t$loc";
          push @{ $result }, "value\t1\tmol_ext_id\tLL\t$loc";
          push @{ $result }, "connect_mol_ext_id\tLL\t$loc";
        } else {
          push @{ $result }, "#$tag\tHs.$cid\t$cid2sym{$cid}\t" .
            "LL $loc\t$tagdat{$tag}";
          push @{ $result }, "#NO PATHWAY FOR: mol_ext_id\tLL\t$loc";
          push @{ $result }, "#NO PATHWAY FOR: value\t1\tmol_ext_id\tLL\t$loc";
        }
      }
    }
  }

}

######################################################################
sub ComputeSDGED_1 {

  my ($org, $factor, $pvalue, $sA, $sB) = @_;

  my ($total_seqsA, $total_seqsB);
  my (@order);

  Init();

  for (@{ $sA }) {
    $setA{$_} = 1 ;
  }
  my $setA = join(",", keys %setA);
  my $total_libsA = scalar(keys %setA);

  for (@{ $sB }) {
    $setB{$_} = 1 ;
    if (defined $setA{$_}) {
      push @order, "Library id $_ is in Pool A and in Pool B";
      return \@order;
    }
  }
  my $setB = join(",", keys %setB);
  my $total_libsB = scalar(keys %setB);

  if ($setA eq "") {
    push @order, "No libraries in Pool A";
    return \@order;
  }
  if ($setB eq "") {
    push @order, "No libraries in Pool B";
    return \@order;
  }

  my ($db);
  my (%freqsA, %freqsB, %tags_hit);

  $db = DBI->connect("DBI:Oracle:" . DB_INSTANCE, DB_USER, DB_PASS);
  if (not $db or $db->err()) {
    push @order, "Cannot connect to " . DB_USER . "@" . DB_INSTANCE . "\n";
    push @order, "$DBI::errstr\n";
    return \@order;
  }

  if ($total_libsA > MANY_SAGE_LIBRARIES ||
      $total_libsB > MANY_SAGE_LIBRARIES) {
    Readsagefreqdat(\%setA, \%freqsA, \$total_seqsA,
        \%setB, \%freqsB, \$total_seqsB, \%tags_hit);
  } else {
    SummarizeTagFreqs($db, $setA, \%freqsA, \$total_seqsA,
        \%tags_hit);
    SummarizeTagFreqs($db, $setB, \%freqsB, \$total_seqsB,
        \%tags_hit);
  }

  if ($total_seqsA == 0) {
    push @order, "No sequences in A";
    return \@order;
  }
  if ($total_seqsB == 0) {
    push @order, "No sequences in B";
    return \@order;
  }

  if ($total_seqsA < MINIMUM_SEQUENCES) {
    push @order, "$total_seqsA sequences in Pool A below minimum " .
        MINIMUM_SEQUENCES;
    return \@order;
  }

  if ($total_seqsB < MINIMUM_SEQUENCES) {
    push @order, "$total_seqsB sequences in Pool B below minimum " .
        MINIMUM_SEQUENCES;
    return \@order;
  }

  SAGEComputeDifferences($db, $factor, $pvalue, \%tags_hit,
      $total_seqsA, $total_seqsB,
      \%freqsA,     \%freqsB, 
      \@order
  );

  $db->disconnect();

  return \@order;

}

######################################################################
sub SDGEDLibrarySelect_1 {
  my (
      $tissue_a,
      $histology_a,
      $tissue_b,
      $histology_b,
      $set_a,
      $set_b
  ) = @_;

  my ($db, $sql, $stm);

  $db = DBI->connect("DBI:Oracle:" . DB_INSTANCE, DB_USER, DB_PASS);
  if (not $db or $db->err()) {
    print STDERR "Cannot connect to " . DB_USER . "@" . DB_INSTANCE . "\n";
    exit();
  }

  SAGESelectLibrarySet($db, $set_a, $tissue_a, $histology_a);
  SAGESelectLibrarySet($db, $set_b, $tissue_b, $histology_b);

  $db->disconnect();

}

######################################################################
sub SAGESelectLibrarySet {
  ## This does a "liberal" selection. That is, if tissue=liver
  ## is specified, it will select any library that has a keyword
  ## k such that IsKindOf($k, 'liver'), even if the library is also
  ## keyworded for tissues that are not 'liver').

  my ($db, $items_ref, $tissue, $hist, $good_size, $good_qual) = @_;

  my ($choices);
  my ($sql, $stm);

  $sql =
      "select distinct i.sage_library_id " .
      "from $CGAP_SCHEMA.sagelibinfo i ";

  if ($tissue) {
    $choices = "('" . join("', '", split(",", $tissue)). "')";
    $sql .= " where " .
            "exists (select k1.sage_library_id " .
        "from $CGAP_SCHEMA.sagekeywords k1 " .
        "where k1.sage_library_id = i.sage_library_id and " .
        "k1.keyword in $choices)"
  }

  if ($hist) {
    $choices = "('" . join("', '", split(",", $hist)). "')";
    $sql .= " and exists (select k4.sage_library_id " .
        "from $CGAP_SCHEMA.sagekeywords k4 " .
        "where k4.sage_library_id = i.sage_library_id and " .
        "k4.keyword in $choices)"
  }

  my ($lid);

  $stm = $db->prepare($sql);
  if(not $stm) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "prepare call failed\n";
  }
  if(!$stm->execute()) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "execute call failed\n";
  }
  while (($lid) = $stm->fetchrow_array()) {
    push @{ $items_ref }, $lid;
  }      

}

######################################################################
sub SDGEDPathway_1 {
  my ($org, $factor, $pvalue,
      $tissue_a, $histology_a, $tissue_b, $histology_b) = @_;

  my (@set_a, @set_b);

  SDGEDLibrarySelect_1(
        $tissue_a,
        $histology_a,
        $tissue_b,
        $histology_b,
        \@set_a,
        \@set_b
    );

  my $order = ComputeSDGED_1($org, $factor, $pvalue,
      \@set_a, \@set_b);

  return join("\n", @{ $order }) . "\n";
  
}

######################################################################
1;
######################################################################

