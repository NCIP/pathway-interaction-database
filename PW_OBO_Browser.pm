#!/usr/local/bin/perl


# Copyright SRA International
#
# Distributed under the OSI-approved BSD 3-Clause License.
# See http://ncip.github.com/pathway-interaction-database/LICENSE.txt for details.


use strict;
use PIDConfig;
use DBI;
use URI::Escape;

if (-d "/app/oracle/product/10gClient") {
  $ENV{'ORACLE_HOME'} = "/app/oracle/product/10gClient"
} elsif (-d "/app/oracle/product/dbhome/current") {
  $ENV{'ORACLE_HOME'} = "/app/oracle/product/dbhome/current";
} elsif (-d "/app/oracle/product/8.1.7") {
  $ENV{'ORACLE_HOME'} = "/app/oracle/product/8.1.7";
} elsif (-d "/app/oracle/product/8.1.6") {
  $ENV{'ORACLE_HOME'} = "/app/oracle/product/8.1.6";
}

my $BASE;

use constant MAX_ROWS_PER_FETCH => 1000;
use constant MAX_LONG_LEN       => 16384;
use constant ORACLE_LIST_LIMIT  => 500;

##
## GO stuff
##

my $GO_ROOT = "1";
my %GO_CLASS_ID = (
  "Signaling Pathways" => "3",
  "Regulatory Pathways" => "4",
  "3"=>"SP",
  "4"=>"RP",
  "SP"=>"Signaling Pathways",
  "RP"=>"Regulatory pathways"
);

######################################################################
sub numerically { $a <=> $b; }

######################################################################
sub GO_IMAGE_TAG {
  my ($image_name) = @_;
  my $tmp = $image_name;
  $tmp =~ s/.jpg$//; 
  return "<image src=\"$BASE/images/$image_name\" " .
      "width=15 height=15 alt=\"$tmp\" border=0>";
}

######################################################################
sub GO_GENE_URL {
  my ($go_id, $org) = @_;
  return "$BASE/Genes/GoGeneQuery?PAGE=1&ORG=$org&" .
      "GOID=$go_id";
}

######################################################################
sub GetGeneCounts {
  my ($db, $id2name, $counts) = @_;

  my ($sql, $stm);
  my ($go_id, $organism, $count);
  my $list = join(",", keys %{ $id2name });

  my $list = 0;
  $sql = qq!
select
  p.pathway_id,
  p.organism,
  5 
from
  $CGAP_SCHEMA.pw_category p
where
      p.pathway_id in ($list)
  !;
  $stm = $db->prepare($sql);
  if (!$stm) {
    ## print STDERR "sql: $sql\n";
    ## print STDERR "$DBI::errstr\n";
    print "<br><b><center>GetGenecountsError in input</b>!</center>";
    $db->disconnect();
    return "";
  } 
  if (!$stm->execute()) {
    ## print STDERR "$sql\n";
    print "execute call failed\n";
    $db->disconnect();
    return "";
  }
  while (($go_id, $organism, $count) = $stm->fetchrow_array()) {
    if ($count) {
      $$counts{$go_id}{$organism} = $count;
    }
  }

}

######################################################################
sub AddAncestors {
  my ($db, $focal_node, $nodes) = @_;

  my ($sql, $stm);
  my ($go_ancestor_id, $go_class);

  $$nodes{$focal_node} = 1;
  my $list = join(",", keys %{ $nodes });

  $sql = qq!
select 
  pm.pathway_ancestor_id, pathway_class 
from
  $CGAP_SCHEMA.pw_pathway_mapping pm,
  $CGAP_SCHEMA.pw_category p
where
      pm.pathway_ancestor_id in ($list)
  and pm.pathway_id = p.pathway_id
  !;

  # print "SQL:  $sql\n";
  $stm = $db->prepare($sql);
  if (!$stm) {
    ## print STDERR "sql: $sql\n";
    ## print STDERR "$DBI::errstr\n";
    print "<br><b><center>AddAncestorsError in input</b>!</center>";
    $db->disconnect();
    return "";
  } 
  if (!$stm->execute()) {
    ## print STDERR "$sql\n";
    print "execute call failed\n";
    $db->disconnect();
    return "";
  }

  while (($go_ancestor_id, $go_class) = $stm->fetchrow_array()) {
   # if ($go_ancestor_id == 4) {
   #   $go_class = 'RP';
   # }
   $$nodes{$go_ancestor_id} = 1;
    
    # print "$go_ancestor_id $go_class <p>";
    ## ancestor table does not include the three top-level ancestors
    ## so add them now, if appropriate.
    ## Here, add Signaling Pathways and Regulatory Pathways 
    if (defined $GO_CLASS_ID{$go_class}) {
      $$nodes{$GO_CLASS_ID{$GO_CLASS_ID{$go_class}}} = 1; 
    }
  }
}

######################################################################
sub CloseNodes {
  my ($db, $focal_node, $nodes) = @_;

  my ($sql, $stm);
  my ($go_id);

  delete $$nodes{$focal_node};

  if (defined $GO_CLASS_ID{$focal_node}) {
    $sql = qq!
select
  pm.pathway_id 
from
  $CGAP_SCHEMA.pw_pathway_mapping pm
where
      pm.pathway_ancestor_id = $focal_node
  or pm.pathway_ancestor_id in (
    select
      p.pathway_ancestor_id
    from
      $CGAP_SCHEMA.pw_pathway_mapping p
    where
      p.pathway_ancestor_id = $focal_node
  )
    !;
  } else {
    $sql = qq!
select
  p.pathway_id 
from
  $CGAP_SCHEMA.pw_pathway_mapping p
where
      p.pathway_ancestor_id = $focal_node
    !;
  }


  # print "SQL:  $sql\n"; 
  $stm = $db->prepare($sql);
  if (!$stm) {
    print "<br><b><center>CloseNodesError in input</b>!</center>";
    $db->disconnect();
    return "";
  } 
  if (!$stm->execute()) {
    ## print STDERR "$sql\n";
    print "execute call failed\n";
    $db->disconnect();
    return "";
  }
  while (($go_id) = $stm->fetchrow_array()) {
    delete $$nodes{$go_id};
  }

}

######################################################################
sub Descend {
  my ($level, $node, $kids_name, $name2id, $lines) = @_;

  my ($name, $parent_type);
  # print "In Descend\n";
  my $nodecount = 0;
  # print "$level $node $kids_name $name2id<P> $lines\n"; 
  if ($level > 30) {
    return;
  }
  for $name (sort keys %{ $$kids_name{$node} }) {
    # print "<P>In Descend for"; 
    $parent_type = $$kids_name{$node}{$name};
    my $pathway_id = $$name2id{$name};
    my $class_type = "categ";
    if ($pathway_id < 5) {
      $class_type = "categ top-node";
    }
    if ($pathway_id > 199999) {
      push @{ $lines }, "$level\t$parent_type\t<a class=\"pathway\" href=http://pid.nci.nih.gov/search/pathway_landing.shtml?pathway_id=$pathway_id&source=NATURE&what=graphic&jpg=on>$name</a>\t$$name2id{$name}";
    } else {
      push @{ $lines }, "$level\t$parent_type\t<span class=\"$class_type\">$name</span>\t$$name2id{$name}";

    } 
    Descend($level+1, $name, $kids_name, $name2id, $lines);
  }
  return;
}

######################################################################
sub GetAllParents {
  my ($db, $all_parents) = @_;

  my ($sql, $stm);
  my ($go_parent_id);

  $sql = qq!
select unique
  pm.pathway_parent
from
  $CGAP_SCHEMA.pw_pathway_parent pm
  !;


  # print "SQL:  $sql\n";
  $stm = $db->prepare($sql);
  if (!$stm) {
    ## print STDERR "sql: $sql\n";
    ## print STDERR "$DBI::errstr\n";
    print "<br><b><center>GetAllParentsError in input</b>!</center>";
    $db->disconnect();
    return "";
  } 
  if (!$stm->execute()) {
    print "execute call failed\n";
    $db->disconnect();
    return "";
  }
  while (($go_parent_id) = $stm->fetchrow_array()) {
    $$all_parents{$go_parent_id} = 1;
  }
}

######################################################################
sub PIDTree {
  my ($db, $nodes, $kids_id, $kids_name, $id2name, $name2id, $lines) = @_;

  my ($sql, $stm);
  my ($go_name, $go_class, $go_id, $go_parent_name,
      $go_parent_id, $parent_type);
  my $list = join(",", keys %{ $nodes });

  if (substr($list,0,1) eq ',') {
    $list = substr($list,1);
  }

  $sql = qq!
select
  p2.pathway_name,
  pm.pathway_class,
  p2.pathway_id,
  p1.pathway_name,
  pm.pathway_parent,
  'P'
from
  $CGAP_SCHEMA.pw_pathway_parent pm,
  $CGAP_SCHEMA.pw_category p1,
  $CGAP_SCHEMA.pw_category p2 
where
      pm.pathway_parent in ($list)
  and p1.pathway_id = pm.pathway_parent
  and p2.pathway_id = pm.pathway_id
  !;

  $stm = $db->prepare($sql);
  # print "SQL:  $sql\n";
  if (!$stm) {
    print "<br><b><center>GotreeError in input</b>!</center>";
    $db->disconnect();
    return "";
  } 
  if (!$stm->execute()) {
    ## print STDERR "$sql\n";
    print "execute call failed\n";
    $db->disconnect();
    return "";
  }
  while (($go_name, $go_class, $go_id, $go_parent_name,
      $go_parent_id, $parent_type) = $stm->fetchrow_array()) {

    #if (defined $GO_OBSOLETE{$go_id}) {
    #  next;
    #}
    $$kids_id{$go_parent_id}{$go_id}       = $parent_type;
    $$kids_name{$go_parent_name}{$go_name} = $parent_type;
    $$id2name{$go_id}   = $go_name;
    $$name2id{$go_name} = $go_id;
    $$id2name{$go_parent_id}   = $go_parent_name;
    $$name2id{$go_parent_name} = $go_parent_id;
  } 

  Descend(0, "", $kids_name, $name2id, $lines);
}

######################################################################
sub FormatOneLine {
  my ($what, $target, $focal_node, $nodes, $all_parents, $counts, $line,
      $lines_out) = @_;

  my ($level, $type, $name, $id) = split(/\t/, $line);
  my ($open_close_url, $count);
  if (defined $$nodes{$id}) {
    ##
    ## id is currently "open"
    ##
    $open_close_url = "<a href=\"javascript:" .
        "document.bf.NODE.value='$id';" .
        "document.bf.CMD.value='close';" .
        "document.bf.submit()\">" .
        GO_IMAGE_TAG("minus.gif") . "</a>";
  } elsif (defined $$all_parents{$id}) {
    ##
    ## id is currently "closed" but has kids
    ##
    $open_close_url = "<a href=\"javascript:" .
        "document.bf.NODE.value='$id';" .
        "document.bf.CMD.value='open';" .
        "document.bf.submit()\">" .
        GO_IMAGE_TAG("plus.gif") . "</a>";
  } else {
    $open_close_url = GO_IMAGE_TAG($type eq "P" ? "partof.gif" : "isa.gif");
  }
  if ($$counts{$id}) {
    if ($what eq "PROTEIN") {
      $count = "&nbsp;<a href=\"" . GO_PROTEIN_URL($name) .
          "\" target=$target>\[$$counts{$id}\]<a>";
    } elsif ($what eq "GENE") {
      my $h = $$counts{$id};
      for my $org ("Hs", "Mm") {
        my $c = $$h{$org};
        if ($c) {
          $count .= "&nbsp;<a href=\"" . GO_GENE_URL($id, $org) .
            "\" target=$target>\[$org:$c\]</a>";
        }
      }
    }
  }
  if ($id eq $focal_node) {
    $name = "<font color=green>$name</font>";
  }
  my $indent;
  for (my $i = 1; $i <= $level; $i++) {
#    $indent .= "&nbsp;&nbsp;&nbsp;&nbsp;"
    $indent .= GO_IMAGE_TAG("blank.gif");
  }
  push @{ $lines_out },
      $indent . $open_close_url . "&nbsp;" . $name . $count . "<br>";
  if ($id eq $focal_node) {
    return 1;
  } else {
    return 0;
  }
}

######################################################################
sub FormatLines {
  my ($what, $url, $target, $focal_node, $nodes, $all_parents, $counts, $lines,
      $lines_out) = @_;

  push @{ $lines_out },
      "<form name=bf method=post action=\"$url\">";
  push @{ $lines_out },
      "<input type=hidden name=CMD>";
  push @{ $lines_out },
      "<input type=hidden name=NODE>";
  for my $n (keys %{ $nodes }) {
    if ($n) {
      push @{ $lines_out },
          "<input type=hidden name=GOIDS value=\"$n\">";
    }
  }

  my ($i, $is_focal, $focal_point);
  for my $line (@{ $lines }) {
    $i++;
    $is_focal =
        FormatOneLine($what, $target, $focal_node, $nodes,
            $all_parents, $counts, $line, $lines_out);
    if ($is_focal) {
      $focal_point = $i;
    }
  }

  push @{ $lines_out }, "</form>\n";

  return $focal_point;
}

######################################################################
sub GOBrowser_1 {
  my ($base, $cmd, $url, $target,
      $focal_node, $context_node_list) = @_;

  $BASE = $base;

  my (%nodes, @lines, %all_parents, %counts,
      %kids_id, %kids_name, %id2name, %name2id);

  my $dbtest = DB_INSTANCE . " " . DB_USER . " " . DB_PASS;
  my $db = DBI->connect("DBI:Oracle:" . DB_INSTANCE, DB_USER, DB_PASS);
  if (not $db or $db->err()) {
    ## print STDERR "Cannot connect to " . DB_USER . "@" . DB_INSTANCE . "\n";
    print STDERR "Cannot connect to database \n";
    ## print STDERR "$DBI::errstr\n";
    return;
  }
 
  $nodes{$GO_ROOT} = 1;
  for my $i (split(",", $context_node_list)) {
    $nodes{$i} = 1;
  }

  if ($cmd eq "expand") {
    my $sql = qq!
select
  pathway_id
from
  $CGAP_SCHEMA.pw_category
!;   
 
    my $stm = $db->prepare($sql);
    if (not $stm) {
      $db->disconnect();
      return "";
    }
    if (!$stm->execute()) {
      print "Execute failed\n";
      $db->disconnect();
      return;
    }

    while ((my $pathway_id) = $stm->fetchrow_array()) {
      $nodes{$pathway_id} = 1; 
    }
    $stm->finish();

  }

  if ($cmd eq "") {
    $cmd = "open";
  }
  if ($focal_node eq "") {
    $focal_node = $GO_ROOT;
  }


  if ($cmd eq "open") {
    $nodes{$focal_node} = 1;
    AddAncestors($db, $focal_node, \%nodes);
  } elsif ($cmd eq "close") {
    CloseNodes($db, $focal_node, \%nodes);
  } else {
    print STDERR "Illegal action: $cmd\n";
  }
 
  PIDTree($db, \%nodes, \%kids_id, \%kids_name, \%id2name, \%name2id, \@lines);
  GetAllParents($db, \%all_parents);
  GetGeneCounts($db, \%id2name, \%counts);
  $db->disconnect();

  my @lines_out;
  my $BANNER_ETC_LINES = 12;
  my $gene_or_prot = '';
  my $focal_line = FormatLines($gene_or_prot, $url, $target, $focal_node,
      \%nodes, \%all_parents, \%counts, \@lines, \@lines_out) + 
      $BANNER_ETC_LINES;
  my $PIXELS_PER_LINE = 16;
  my $y_coord = $focal_line * $PIXELS_PER_LINE;
  push @lines_out,
      qq!
<script>
  window.scrollTo(0,$y_coord);
</script>
      !;
  return join("\n", @lines_out) . "\n";

}

######################################################################
sub GetGOTotals {
  my ($db, $org, $go2cid, $total, $direct_total) = @_;

  my ($i, $list, $sql, $stm);
  my ($go_id, $tot, $direct_tot);
  my ($row, $rowcache);

  my @go_ids = keys %{ $go2cid };

  for($i = 0; $i < @go_ids; $i += ORACLE_LIST_LIMIT) {
    if(($i + ORACLE_LIST_LIMIT - 1) < @go_ids) {
      $list = join(",", @go_ids[$i..$i+ORACLE_LIST_LIMIT-1]);
    }
    else {
      $list = join(",", @go_ids[$i..@go_ids-1]);
    }
    $sql = qq!
select
  go_id,
  ug_count,
  direct_ug_count
from
  $CGAP_SCHEMA.go_count
where
      organism = '$org'
  and go_id in ($list)
    !;
    $stm = $db->prepare($sql);
    if (not $stm) {
      print "<br><b><center>GetGoTotals Error in input</b>!</center>";
      $db->disconnect();
      return "";
    }
    if (!$stm->execute()) {
      print "Execute failed\n";
      $db->disconnect();
      return;
    }
    while ($rowcache = $stm->fetchall_arrayref(undef, MAX_ROWS_PER_FETCH)) {
      for $row (@{ $rowcache }) {
        ($go_id, $tot, $direct_tot) = @{ $row };
        $$total{$go_id} = $tot;
        $$direct_total{$go_id} = $direct_tot;
      }
    }
  }
}

######################################################################
sub GetGONames {
  my ($db, $go2cid, $go2name, $go2class) = @_;

  my ($list, $i, $sql, $stm);
  my ($go_id, $go_name, $go_class);
  my ($rowcache, $row);

  my @go_ids = keys %{ $go2cid };

  for($i = 0; $i < @go_ids; $i += ORACLE_LIST_LIMIT) {
    if(($i + ORACLE_LIST_LIMIT - 1) < @go_ids) {
      $list = join(",", @go_ids[$i..$i+ORACLE_LIST_LIMIT-1]);
    }
    else {
      $list = join(",", @go_ids[$i..@go_ids-1]);
    }
    $sql = qq!
select
  p.pathway_id,
  p.pathway_name,
  pathway_class 
from
  $CGAP_SCHEMA.pw_pathway p, $CGAP_SCHEMA.pw_pathway_mapping pp

where
      p.pathway_id in ($list)
      and p.pathway_id = pp.pathway_id
    !;

    $stm = $db->prepare($sql);
    if (not $stm) {
      print "<br><b><center>GetGoNames:  Error in input</b>!</center>";
      $db->disconnect();
      return "";
    }
    if (!$stm->execute()) {
      print "Execute failed\n";
      $db->disconnect();
      return;
    }
    while ($rowcache = $stm->fetchall_arrayref(undef, MAX_ROWS_PER_FETCH)) {
      for $row (@{ $rowcache }) {
        ($go_id, $go_name, $go_class) = @{ $row };
        $$go2name{$go_id} = $go_name;
        $$go2class{$go_id} = $go_class;
      }
    }
  }
}

######################################################################
sub GetGOAnnotations {
  my ($db, $org, $cids, $go2cid) = @_;

  my ($list, $i, $sql, $stm);
  my ($go_id, $cid);
  my ($rowcache, $row);

  for($i = 0; $i < @{ $cids }; $i += ORACLE_LIST_LIMIT) {
    if(($i + ORACLE_LIST_LIMIT - 1) < @{ $cids }) {
      $list = join(",", @{ $cids }[$i..$i+ORACLE_LIST_LIMIT-1]);
    }
    else {
      $list = join(",", @{ $cids }[$i..@{ $cids } - 1]);
    }
    $sql = qq!
select
  a.go_ancestor_id,
  u.cluster_number
from
  $CGAP_SCHEMA.ll_go l,
  $CGAP_SCHEMA.gene2unigene u,
  $CGAP_SCHEMA.go_ancestor a
where
      l.ll_id = u.gene_id
  and a.go_id = l.go_id
  and u.cluster_number in ($list)
  and u.organism = '$org'
    !;

    $stm = $db->prepare($sql);
    if (not $stm) {
      print "<br><b><center>GetGOAnnotationsError in input</b>!</center>";
      $db->disconnect();
      return "";
    }
    if (!$stm->execute()) {
      print "Execute failed\n";
      $db->disconnect();
      return;
    }
    while ($rowcache = $stm->fetchall_arrayref(undef, MAX_ROWS_PER_FETCH)) {
      for $row (@{ $rowcache }) {
        ($go_id, $cid) = @{ $row };
        $$go2cid{$go_id}{$cid} = 1;
      }
    }
  }
}

######################################################################
1;







