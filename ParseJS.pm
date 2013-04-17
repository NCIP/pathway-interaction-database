

# Copyright SRA International
#
# Distributed under the OSI-approved BSD 3-Clause License.
# See http://ncip.github.com/pathway-interaction-database/LICENSE.txt for details.


package ParseJS;
require Exporter;
@ISA = qw(Exporter);
#@EXPORT = qw(
#  CreateIsA
#  CreateLabelSort
#);

BEGIN {
  my @path_elems = split("/", $0);
  pop @path_elems;
  push @INC, join("/", @path_elems);
}

use PWLabel;
use Pathway;
use FileHandle;
use strict;

use constant OP              => 0;
use constant NAME            => 1;
use constant EXT_IDS         => 2;
use constant ACTIVITY_LABELS => 3;
use constant LOC_FUNC_LABELS => 4;
use constant ROLE            => 5;
use constant COMPONENTS      => 6;
use constant MISC_INFO       => 7;

use constant PTMS            => COMPONENTS;
use constant REFNAME         => EXT_IDS;
use constant PMID            => NAME;
use constant NOTE            => NAME;

# MISC_INFO is ;-separated:
#   revsersible/irreversible
#   ? atom condition ?
#   moltype
#   macroprocess subtype

my %abbrev2mod = (
  "Ac" => "acetylation",
  "Bi" => "biotinylation",
  "De" => "dephosphorylation",
  "Fa" => "farnesylation",
  "Gag"=> "glycosaminoglycan",
  "Ge" => "geranylgeranylation",
  "Gl" => "glycosylation",
  "OH" => "hydroxylation", 
  "Me" => "methylation",
  "My" => "myristoylation",
  "Ox" => "oxidation",
  "Pa" => "palmitoylation",
  "Ph" => "phosphorylation",
  "pADPr" => "poly (ADP-ribosyl)ation",
  "Su" => "sumoylation",
  "Ub" => "ubiquitination"
);

my %interaction_types = (
  "Modification"      => "modification",
  "NoModification"    => "modification",
  "Transcription"     => "transcription",
  "NoTranscription"   => "transcription",
  "Translocation"     => "translocation",
  "NoTranslocation"   => "translocation",
  "MacroProcess"      => "macroprocess",
  "NoMacroProcess"    => "macroprocess",
  "Reaction"          => "reaction",
  "NoReaction"        => "reaction",
  "Pathway"           => "pathway",
  "NoPathway"         => "pathway",
  "Subnet"            => "subnet",
  "NoSubnet"          => "subnet"
);

my %evidence_ops = (
  "IAE" => 1,
  "IC"  => 1,
  "IDA" => 1,
  "IEA" => 1,
  "IEP" => 1,
  "IFC" => 1,
  "IGI" => 1,
  "IMP" => 1,
  "IOS" => 1,
  "IPI" => 1,
  "ISS" => 1,
  "NAS" => 1,
  "ND"  => 1,
  "NR"  => 1,
  "RCA" => 1,
  "RGE" => 1,
  "TAS" => 1
);

## @tmp_interactions[interaction_id]{interation_element}: TAB-separated
## string of filtered,original-text fields of Shownamed. Cleaned out with
## each new pathway that is parsed in.

## %interactions{pathway_id}{interaction_id}{edge_id}: array of
##   [ operator, molname, ext_ids, activity, location/function, role,
##     moltype, ptms_and_components, misc_info ]
## (where misc_info is for elements that represent the interaction
## itself; has: reversible, condition, moltype, macroprocess_subtype).
## Accumulated across parsed pathways.
## pathway_ids and interaction_ids are not yet mapped.

## %pathway2curator{pathway_id}[]: curator
## %interaction2type{pathway_id}[interaction_id]: interaction_type
## %interaction2abstraction{pathway_id}[interaction_id]: subnet or
##      pathway external id
## %interaction2condition{pathway_id}[interaction_id]{condition}: 1

## %interaction2evidence{pathway_id}{interaction_id}{evidence_op}: 1
## %interaction2reference{pathway_id}{interaction_id}[]: pubmedid
## %interaction2note{pathway_id}{interaction_id}[]: note text

## moldefs{pathway_id}[interaction_id]{edge}: array of
##   [ molname, ext_ids, activity, location/function, moltype,
##     components, ptms ]

## familydefs{pathway_id}[interaction_id]{edge}: array of
##   [ molname, ext_ids, activity, location/function, moltype,
##     components, ptms ]

## emptyfamilytype{pathway_id}[interaction_id]{edge}

## %mapped_interactions{pathway_id}[interaction_id]:
##    mapped (global) interaction_id
## %mapped_pathways{pathway_id}: mapped (global) pathway id
## %edge_mols{pathway_id}[interaction_id]{edge}: mapped (global) mol id
## %comp_mols{pathway_id}[interaction_id]{edge}{comp_seq}:
##    mapped (global) mol id
## %family_mols{pathway_id}[interaction_id]{edge}[]:
##    mapped (global) mol id
## %edge_types{pathway_id}[interaction_id]{edge}: edgetype

my ($lv, $pw, $idm);
my ($pathway_id, %pathway);
my (@tmp_interactions, %interactions, %moldefs, %familydefs,
    %mapped_interactions, %mapped_pathways,
    %emptyfamilytype,
    %edge_mols, %comp_mols, %edge_types);
my (%interaction2evidence, %interaction2reference, %interaction2note,
    %interaction2type, %interaction2condition, %interaction2abstraction,
    %pathway2curator);

######################################################################
sub new {
  my ($self, $lv0, $pw0, $idm0) = @_;
  $lv = $lv0;
  $pw = $pw0;
  $idm = $idm0;
  return bless {};
}

######################################################################
sub ParseJS {
  my ($self, $js_f) = @_;

  ReadJSInput($js_f);
  Transform();
  MapInteractions();
  MapMols();
  MapPathways();
  Convert();
}

######################################################################
sub numerically { $a <=> $b };

######################################################################
sub ReadJSInput {
  my ($f) = @_;
  open(INF, $f) or die "cannot open $f";
  while (<INF>) {
    s/[\r\n]+//;
    s/^\s+//;
    s/\s+$//;
    ## if file did not end properly with a \n, then
    ## there might not be a proper separation between the files
    if (/^\}./) {
print STDERR "## found one $_\n";
      ParseLine("}");
      s/^\}\s*//;
print STDERR "## reduced to $_\n";
    }
    ParseLine($_);
  }
  close INF;
}

######################################################################
sub NormalizeLabels {
  my ($label_list) = @_;

  my %labels;
  my @labels = split(";", $label_list);
  for my $x (@labels) {
    $x =~ s/^\s+//;
    $x =~ s/\s+$//;
    $x =~ s/\s+/ /g;
    $x = lc($x);
    $labels{$x} = 1;
  }
  return join(";", sort keys %labels);
}

######################################################################
sub NormalizePTMs {
  my ($expr) = @_;

  my %expr;
  for my $term (split("_", $expr)) {
    $expr{$term} = 1;
  }
  return join("_", sort keys %expr);
}

######################################################################
sub EqComplex {
  my ($components1, $components2) = @_;

  if ($components1 eq "" && $components2 eq "") {
    return 1;
  }

  my @components1 = @{ $components1 };
  my @components2 = @{ $components2 };

  if (@components1 != @components2) {
    return 0;
  }
  my $found;
  for (my $i = 0; $i < @components1; $i++) {
    my ($name1, $ext_ids1, $activity1, $location1, $moltype1, $comps1)
        = @{ $components1[$i] };
    my $m1 = [ $name1, $ext_ids1, $activity1, $location1, $moltype1,
        $comps1 ];
    my $ptms1 = $comps1;   ## there aren't any components of components
    for (my $j = 0; $j < @components2; $j++) {
      my ($name2, $ext_ids2, $activity2, $location2, $moltype2, $comps2)
          = @{ $components2[$j] };
      my $m2 = [ $name2, $ext_ids2, $activity2, $location2, $moltype2,
          $comps2 ];

      my $ptms2 = $comps2;   ## there aren't any components of components
      if (EqMolecule($m1, $m2) && $activity1 eq $activity2 &&
          $location1 eq $location2 && $ptms1 eq $ptms2) {
        $found++;
        last;
      }
    }
    if ($found < $i+1) {
      last;
    }
  }
  if ($found == @components1) {
    return 1;
  } else {
    return 0;
  }
}

######################################################################
sub MergeNames {
  my ($i, $j) = @_;

  my ($pex1, $inter1, $edge1, $seq1) = split(",", $i, 4);
  my ($pex2, $inter2, $edge2, $seq2) = split(",", $j, 4);
  my ($replaced_name, $replacing_name, $merged_name);

  if ($seq2) {
    $replaced_name = $moldefs{$pex2}[$inter2]{$edge2}[5][$seq2-1][0];
   } else {
    $replaced_name = $moldefs{$pex2}[$inter2]{$edge2}[0];
  }
  if ($seq1) {
    $replacing_name = $moldefs{$pex1}[$inter1]{$edge1}[5][$seq1-1][0];
  } else {
    $replacing_name = $moldefs{$pex1}[$inter1]{$edge1}[0];
  }

  my (%seen, @result);
  for my $n (split(";", $replacing_name), split(";", $replaced_name)) {
    if ($n && ! defined $seen{$n}) {
      push @result, $n;
      $seen{$n} = 1;
    }
  }
  $merged_name = join(";", @result);

  if ($seq1) {
    $moldefs{$pex1}[$inter1]{$edge1}[5][$seq1-1][0] = $merged_name;
  } else {
    $moldefs{$pex1}[$inter1]{$edge1}[0] = $merged_name;
  }
  # also replace the loser's name, since this may later be used in Convert
  if ($seq2) {
    $moldefs{$pex2}[$inter2]{$edge2}[5][$seq2-1][0] = $merged_name;
  } else {
    $moldefs{$pex2}[$inter2]{$edge2}[0] = $merged_name;
  }

}

######################################################################
sub NamesIntersect {
  my ($n1, $n2) = @_;

  my (%n1);
  for my $n (split(";", $n1)) {
    if ($n) {
      $n1{lc($n)} = 1;
    }
  }
  for my $n (split(";", $n2)) {
    if (defined $n1{lc($n)}) {
      return 1;
    }
  }
  return 0;
}

######################################################################
sub EqMolecule {
  my ($m1, $m2) = @_;

  my ($name1, $ext_ids1, $activity1, $location1, $moltype1, $components1)
      = @{ $m1 };
  my ($name2, $ext_ids2, $activity2, $location2, $moltype2, $components2)
      = @{ $m2 };
#print "comparing " . join(",", @{$m1}) . " and " . join(",", @{$m2}) . "\n";

  if ($moltype1 ne $moltype2) {
#    print STDERR "## types not equal: '$moltype1', '$moltype2'\n";
    if ($moltype1 ne "cleave" && $moltype1 ne "family") {
      return 0;
    }
  }
  if ($ext_ids1 ne $ext_ids2) {
#print "ext ids not equal\n";
    return 0;
  }
  if ($moltype1 eq "complex") {
    if ($components2 ne "" && ref($components2) ne "ARRAY") {
      my $reftype = ref($components2);
      # $components2 = @{ $components2 };

      print STDERR "INTERNAL ERROR: components2 is not an array: $components2 Reftype:  $reftype\n";
      return 0;
    }
    if ($components2 eq "" || @{ $components2 } == 0) {
      if (NamesIntersect($name1, $name2)) {
        return 1;
      } else {
        return 0;
      }
    } else {
      # return EqComplex($components1, $components2);
      return 0;
    }
  }
  if (! NamesIntersect($name1, $name2)) {
    return 0;
  }
  return 1;
}

######################################################################
sub MapInteractions {

  my $n = $idm->AtomId();
  for my $pid (keys %interactions) {
    for (my $inter = 1; $inter < @{ $interactions{$pid} }; $inter++) {
      $n++;
      $mapped_interactions{$pid}[$inter] = $n;
    }
  }
}

######################################################################
sub MapPathways {

  my $n = $idm->PathwayId();
  for my $pid (keys %interactions) {
    $n++;
    $mapped_pathways{$pid} = $n;
  }
  for my $pid (keys %interaction2abstraction) {
    for (my $i = 1; $i < @{ $interaction2abstraction{$pid} }; $i++) {
      my $ext_pathway_id = $interaction2abstraction{$pid}[$i];
      if (! defined $mapped_pathways{$ext_pathway_id}) {
        $n++;
        $mapped_pathways{$ext_pathway_id} = $n;
      }
    } 
  }
}

######################################################################
sub MapMols {

  use constant MOLTYPE_OFFSET    => 4;
  use constant COMPONENTS_OFFSET => 5;

  ## order: (1) families and whole/parts; (2) simple molecules;
  ## (3) "outer" complexes; (4) "inner" complexes

  my (%replace,
      @all_moldefs,   @regular,   @family_cleave,   @outer,   @inner,
      %r_all_moldefs, %r_regular, %r_family_cleave, %r_outer, %r_inner
     );
  my ($n,             $regularn,  $family_cleaven,  $outern,  $innern) =
     (0,              0,          0,                0,        0);

  for my $pid (keys %moldefs) {
    for (my $i = 0; $i < @{ $moldefs{$pid} }; $i++) {
      for my $j (sort numerically keys %{ $moldefs{$pid}[$i] }) {
        my $moltype = $moldefs{$pid}[$i]{$j}[MOLTYPE_OFFSET];
        my $comps   = $moldefs{$pid}[$i]{$j}[COMPONENTS_OFFSET];
        if ($moltype eq "family" || $moltype eq "cleave") {
          $family_cleave[$family_cleaven]   = $moldefs{$pid}[$i]{$j};
          $r_family_cleave{$family_cleaven} = "$pid,$i,$j";
          $family_cleaven++;
        } elsif ($moltype eq "complex") {
          $outer[$outern]   = $moldefs{$pid}[$i]{$j};
          $r_outer{$outern} = "$pid,$i,$j";
          $outern++;
        } else {
          $regular[$regularn]   = $moldefs{$pid}[$i]{$j};
          $r_regular{$regularn} = "$pid,$i,$j";
          $regularn++;
        }
        ## if a complex or family or cleave, make a separate entry for each component
        my $comp_seq = 0;
        if (ref($comps) eq "ARRAY") {
          for my $comp (@{ $comps }) {
            $comp_seq++;
            my $moltype2 = $$comp[MOLTYPE_OFFSET];
            if ($moltype2 eq "complex") {
              $inner[$innern]   = $comp;
              $r_inner{$innern} = "$pid,$i,$j,$comp_seq";
              $innern++;
            } else {
              $regular[$regularn]   = $comp;
              $r_regular{$regularn} = "$pid,$i,$j,$comp_seq";
              $regularn++;
            }
          }
        }
#        print join("\t", "##moldefs:", "$pid,$i,$j",
#            join(",", @{ $moldefs{$pid}[$i]{$j} }[0..4],
#              join(",", @{ $moldefs{$pid}[$i]{$j}[5] }))) . "\n";
      }
    }
  }

  for (my $i = 0; $i < @family_cleave; $i++) {
    $all_moldefs[$n]   = $family_cleave[$i];
    $r_all_moldefs{$n} = $r_family_cleave{$i};
    $n++;
  }
  for (my $i = 0; $i < @regular; $i++) {
    $all_moldefs[$n]   = $regular[$i];
    $r_all_moldefs{$n} = $r_regular{$i};
    $n++;
  }
  for (my $i = 0; $i < @outer; $i++) {
    $all_moldefs[$n]   = $outer[$i];
    $r_all_moldefs{$n} = $r_outer{$i};
    $n++;
  }
  for (my $i = 0; $i < @inner; $i++) {
    $all_moldefs[$n]   = $inner[$i];
    $r_all_moldefs{$n} = $r_inner{$i};
    $n++;
  }
  for (my $i = 0; $i < @all_moldefs - 1; $i++) {
    if (defined $replace{$i}) {
      next;
    }

    for (my $j = $i+1; $j < @all_moldefs; $j++) {
      if (EqMolecule($all_moldefs[$i], $all_moldefs[$j])) {
        ## true: j > i
        ## if this is the first time that j has been found equivalent to a
        ## i, then mark j for replacement by i (guarantees that
        ## j will be replaced by the lowest possible i)
        if (! defined $replace{$j}) {

          $replace{$j} = $i;

          #
          # fill in empty family types, if necessary
          #
          my ($name1, $ext_ids1, $activity1, $location1, $moltype1,
              $components1) = @{ $all_moldefs[$i] };
          if ($moltype1 eq "family" &&
              ($components1 eq "" || @{ $components1 } == 0)) {
            my ($name2, $ext_ids2, $activity2, $location2, $moltype2,
                $components2) = @{ $all_moldefs[$j] };
            if ($moltype2 eq "protein" ||
                $moltype2 eq "compound" ||
                $moltype2 eq "complex" ||
                $moltype2 eq "rna") {
              my ($pex, $inter, $edge) = $r_all_moldefs{$i};
              if (defined $emptyfamilytype{$pex}[$inter]{$edge}) {
                print STDERR "## empty mol type for family $name1 already " .
                    "defined " .
                    "to be $emptyfamilytype{$pex}[$inter]{$edge}\n";
              } else {
                $emptyfamilytype{$pex}[$inter]{$edge} = $moltype2;
              }
            }
          }

        }
      }
    }
  }

  # merge names of replaced mols (j) into names of replacing mols (i)
  for my $j (keys %replace) {
    my $i = $replace{$j};
    MergeNames($r_all_moldefs{$i}, $r_all_moldefs{$j});
  }

  my $next_mol = $idm->MolId();
  my (%mol_remap);
  for (my $i = 0; $i < @all_moldefs; $i++) {
    if (! defined $replace{$i}) {
      if (! defined $mol_remap{$i}) {
        $next_mol++;
        $mol_remap{$i} = $next_mol;
      }
    }
  }

  for (my $i = 0; $i < @all_moldefs; $i++) {
    my ($pid, $inter, $edge, $comp) = split(",", $r_all_moldefs{$i});
    if (defined $replace{$i}) {
      if ($comp) {
        $comp_mols{$pid}[$inter]{$edge}{$comp} = $mol_remap{$replace{$i}};
      } else {
        $edge_mols{$pid}[$inter]{$edge} = $mol_remap{$replace{$i}};
      }
    } else {
      if ($comp) {
        $comp_mols{$pid}[$inter]{$edge}{$comp} = $mol_remap{$i};
      } else {
        $edge_mols{$pid}[$inter]{$edge} = $mol_remap{$i};
      }
    }
  }
}

######################################################################
sub TransformNote {
  my ($i, $a) = @_;

  if ($$a[NOTE] ne "") {
    ## older js files did not have a ReferenceGif so used a NotesGif
    if ($$a[NOTE] =~ /^\s*PMID:\s*(\d+)\s*(.*)/) {
      my ($pmid, $refname) = ($1, $2);
      $$a[NOTE] = $pmid;
      $$a[REFNAME] = $refname;
      TransformReference($i, $a);
    } else {
      push @{ $interaction2note{$pathway_id}{$i} }, $$a[NOTE];
    }
  }
}

######################################################################
sub TransformEvidence {
  my ($i, $op) = @_;

  $interaction2evidence{$pathway_id}{$i}{$op} = 1;
}

######################################################################
sub TransformReference {
  my ($i, $a) = @_;

  if ($$a[PMID] ne "" || $$a[REFNAME] ne "") {
    my $pmid = $$a[PMID];
    $pmid =~ s/^\s+//;
    $pmid =~ s/\s+$//;
    my $name = $$a[REFNAME];
    $name =~ s/^\s+//;
    $name =~ s/\s+$//;
    $name =~ s/\s+/ /g;
    push @{ $interaction2reference{$pathway_id}{$i} }, $pmid;
  }
}

######################################################################
sub NormalizeExtIds {
  my ($e)= @_;

  if ($e eq "") {
    return $e;
  }

  my ($type, $id);

  if ($e =~ /^(O|P|Q)[\dA-Z]{5}(-\d+)?$/) {
    ($type, $id) = ("UP", $e);
  } elsif ($e =~ /^\d+$/) {
    ($type, $id) = ("LL", $e);
  } elsif ($e =~ /(CAS)(:?\s*)(.+)/) {
    ($type, $id) = ("CA", $3);
  } elsif ($e =~ /(LL|UP|EC|CA|KG|GO)(:?\s*)(.+)/) {
    ($type, $id) = ($1, $3);
  } else {
    ($type, $id) = ("ZZ", $e);
  }
  return "$type:$id";
}

######################################################################
sub TransformMolecule {
  my ($i, $edge, $a) = @_;

  if (defined $interaction2type{$pathway_id}[$i]) {
    $$a[ROLE] = "output";
  } else {
    $$a[ROLE] = "input";
  }
  $$edge++;
  $edge_types{$pathway_id}[$i]{$$edge} = $$a[ROLE];
  $interactions{$pathway_id}[$i]{$$edge} = $a
}

######################################################################
sub TransformInteractionType {
  my ($i, $j, $op, $edge, $a) = @_;

  my ($reversible, $condition, $moltype, $macroprocess_subtype) =
      split(";", $$a[MISC_INFO], 1000);
#print STDERR "## i = $i\n";
#print STDERR "## j = $j\n";
#print STDERR "## op = $op\n";
#print STDERR "## edge = $$edge\n";
#print STDERR "## reversible = $reversible\n";
#print STDERR "## condition = $condition\n";
#print STDERR "## moltype = $moltype\n";
#print STDERR "## macroprocess_subtype = $macroprocess_subtype\n";
  for my $x ($reversible, $condition, $moltype, $macroprocess_subtype) {
    $x =~ s/^\s+//;
    $x =~ s/^\s+$//;
    $x =~ s/\s+/ /g;
    $x = lc($x);
    $x =~ s/^go\s?:?\s?(\d+)$/GO:$1/;
  }
  if (defined $interaction2type{$pathway_id}[$i] &&
      $interaction2type{$pathway_id}[$i] ne $interaction_types{$op}) {
    print STDERR "## multiple interaction types for pathway $pathway_id, " .
        "interaction $i: $interaction2type{$pathway_id}[$i], $op\n";
  } else {
    if ($interaction_types{$op} eq "macroprocess") {
      $interaction2type{$pathway_id}[$i] = $macroprocess_subtype;
    } else {
      $interaction2type{$pathway_id}[$i] = $interaction_types{$op};
    }
    if ($interaction_types{$op} eq "pathway" ||
        $interaction_types{$op} eq "subnet") {
      $interaction2abstraction{$pathway_id}[$i] = $macroprocess_subtype;
    }
  }

  if ($condition ne "") {
    $interaction2condition{$pathway_id}[$i]{$condition} = 1;
  }

  if ($op =~ /^No/) {
    $$a[ROLE] = "inhibitor";
  } elsif ($interaction_types{$op} eq "pathway" || 
      $interaction_types{$op} eq "subnet") {
    $$a[ROLE] = "input";
  } else {
    $$a[ROLE] = "agent";
  }

  if ($moltype eq "") {
##    print STDERR "## empty moltype $pathway_id, $i, $$edge: $tmp_interactions[$i]{$j}\n";
  } elsif ($$a[NAME] eq "" && $$a[EXT_IDS] eq "") {
    print STDERR "## no mol name or ext_id $pathway_id, $i, $$edge: $tmp_interactions[$i]{$j}\n";
  } else {
    $$edge++;
    $edge_types{$pathway_id}[$i]{$$edge} = $$a[ROLE];
    $interactions{$pathway_id}[$i]{$$edge} = $a;
    $moldefs{$pathway_id}[$i]{$$edge} = [
        $$a[NAME],
        $$a[EXT_IDS],
        $$a[ACTIVITY_LABELS],
        $$a[LOC_FUNC_LABELS],
        lc($moltype),                     ##!! lc moltype
        TransformComponents($i, $j, $a)
      ];
  }
}

######################################################################
sub TransformComponents {
  my ($i, $j, $a) = @_;

  my @temp;
  my $comp_ptms = $$a[COMPONENTS];
  my ($components, $ptms0);
  my @comp_ptms = split(";", $comp_ptms, 1000);
#print STDERR join("\n", "comp_ptms:", @comp_ptms, "");
  # following complexity to work around some bad js, as a result of
  # which we cannot always expect the last ;-split-fragment being
  # the ptms
  if (@comp_ptms > 2) {
    if ($comp_ptms[$#comp_ptms-1] !~ /\&/) {
      print STDERR "illegal final ';' in comp_ptms: " . $$a[COMPONENTS] . "\n";
      $ptms0 = $comp_ptms[$#comp_ptms-1];
      $components = join(";", @comp_ptms[0..$#comp_ptms-2]);
    } else {
      $ptms0 = $comp_ptms[$#comp_ptms];
      $components = join(";", @comp_ptms[0..$#comp_ptms-1]);
    }
  } else {
    ($components, $ptms0) = ($comp_ptms[0], $comp_ptms[1]);
  }

  if ($components ne "") {
    my @fields = split(/\&/, $components, 1000);
    my $mod;
    if (@fields % 6 == 0) {
      $mod = 6;
    } elsif (@fields % 5 == 0) {
      $mod = 5;
    } else {
      print STDERR "# number fields neither mod 6 nor mod 5 in " .
          "$pathway_id, " . "$i, $j:\t$components\n";
    }
    while (@fields > 1) {
      my ($name, $ext_ids, $activity, $location, $moltype, $ptms);
      if ($mod == 6) {
        ($name, $ext_ids, $activity, $location, $moltype, $ptms) =
          splice(@fields, 0, 6);
      } elsif ($mod == 5) {
        ($name, $ext_ids, $activity, $location, $moltype) =
            splice(@fields, 0, 5);
      } else {
        last;
      }
      $ext_ids  = NormalizeExtIds($ext_ids);
      $activity = NormalizeLabels($activity);
      $location = NormalizeLabels($location);
      push @temp, [ $name, $ext_ids, $activity, $location, $moltype, $ptms ];
#        print STDERR "##comp name $pathway_id, $i, $j\t$name\n";
#        print STDERR "##comp ext_ids $pathway_id, $i, $j\t$ext_ids\n";
#        print STDERR "##comp activity $pathway_id, $i, $j\t$activity\n";
#        print STDERR "##comp location $pathway_id, $i, $j\t$location\n";
#        print STDERR "##comp moltype $pathway_id, $i, $j\t$moltype\n";
#        print STDERR "##comp ptms $pathway_id, $i, $j\t$ptms\n";
    }
  }
  return (\@temp, $ptms0);
}

######################################################################
sub TransformReactionElement {
  my ($i, $j, $edge) = @_;

  my @a = split("\t", $tmp_interactions[$i]{$j}, 1000);
  my $op = $a[OP];

  if (defined $evidence_ops{$op}) {
    TransformEvidence($i, $op);
  } elsif ($op eq "Reference") {
    TransformReference($i, \@a);
  } elsif ($op eq "Notes")  {
    TransformNote($i, \@a);
  } elsif ($op eq "Associates") {
    # Associates is a pure operator, does not define a participant
  } elsif (defined $interaction_types{$op}) {
    $a[EXT_IDS] = NormalizeExtIds($a[EXT_IDS]);
    TransformInteractionType($i, $j, $op, $edge, \@a);
  } elsif ($op eq "Complex" || $op eq "Protein" || $op eq "RNA" || 
        $op eq "Compound" || $op eq "Family" || $op eq "Cleave") {
    $a[EXT_IDS]         = NormalizeExtIds($a[EXT_IDS]);
    $a[ACTIVITY_LABELS] = NormalizeLabels($a[ACTIVITY_LABELS]);
    $a[LOC_FUNC_LABELS] = NormalizeLabels($a[LOC_FUNC_LABELS]);
    if ($a[NAME] eq "" && $a[EXT_IDS] eq "") {
      print STDERR "## no mol name or ext_id $pathway_id, $i, $j: $tmp_interactions[$i]{$j}\n";
    } else {
      TransformMolecule($i, $edge, \@a);
      my ($components, $ptms) = TransformComponents($i, $j, \@a);
      $moldefs{$pathway_id}[$i]{$$edge} = [
          $a[NAME],
          $a[EXT_IDS],
          $a[ACTIVITY_LABELS],
          $a[LOC_FUNC_LABELS],
          lc($op),                # moltype,
          $components,
          $ptms
        ];
      if ($op eq "Family" && @{ $components } == 0) {
        $emptyfamilytype{$pathway_id}[$i]{$$edge} = "";
      }
    }
  }
}

######################################################################
sub TransformOneInteraction {
  my ($i) = @_;

# It seems that every interaction has at least one of these:
# (No)Modification, (No)Transcription, (No)Translocation, (No)MacroProcess.
# Associates is a pure operator (like "+")

  my $edge = 0;

  for my $j (sort numerically keys %{ $tmp_interactions[$i] }) {
    TransformReactionElement($i, $j, \$edge);
  }
}

######################################################################
sub Transform {
  for (my $i = 0; $i < @tmp_interactions; $i++) {
    if ($tmp_interactions[$i]) {
      TransformOneInteraction($i);
    }
  }
}

######################################################################
sub ConvertPTM {
  my ($expr) = @_;

  my @temp;
  for my $term (split("_", $expr)) {
    my ($aa, $pos, $modabbrev);
    my ($modstring, $modvalue);
    if ($term =~ /([A-Z])(\d+)(\w+)/) {
      ($aa, $pos, $modabbrev) = ($1, $2, $3);
      if (defined $abbrev2mod{$modabbrev}) {
        $modstring = $abbrev2mod{$modabbrev};
        $modvalue = $lv->StringToLabelValue("ptm", $modstring);
        if ($modvalue eq 0) {
## ouch: some ptm terms are also GO process terms
          $modvalue = $lv->StringToLabelValue("process-type", $modstring);
        }
      } else {
        print STDERR "##bad mod abbrev in term: $term\n";
      }
      push @temp, [$pos, $aa, $modstring, $modvalue ];
# print STDERR "#convertPTM: $pos, $aa, $modstring, $modvalue\n";
    } else {
      print STDERRR "##bad ptm: $term\n";
    }
  }
  return \@temp;
}

######################################################################
sub ConvertCleaved {
  my ($pex, $i, $edge, $whole_mol, $parts) = @_;

  my $seq;

  for my $part (@{ $parts }) {

    $seq++;

    my ($name, $ext_ids, $activity_labels, $loc_func_labels, $moltype,
        $ptms) = @{ $part };

    my $part_mol = $comp_mols{$pex}[$i]{$edge}{$seq};

    my $moltype_lvid = $lv->StringToLabelValue("molecule-type", $moltype);
    if ($moltype_lvid == 0) {
      print STDERR "## no moltype for part $pex, interaction $i, edge $edge, " .
          "part $seq\n";
    }
    $pw->AddMol($part_mol, $moltype_lvid);
    my ($start_site, $end_site) = sort numerically split("-", $loc_func_labels);
    $start_site =~ s/\s+//g;
    $end_site   =~ s/\s+//g;
    $pw->AddMolPart($part_mol, $whole_mol, $start_site, $end_site);

    my ($pf_name, @as_names) = split(";", $name);
    if ($pf_name) {
      $pw->AddMolName($part_mol, "PF", $pf_name);
    }
    for my $as_name (@as_names) {
      if ($as_name) {
        $pw->AddMolName($part_mol, "AS", $as_name);
      }
    }

    for my $x (split(";", $ext_ids)) {
      if ($x) {
        my ($type, $id) = split(":", $x);
        $pw->AddMolExId($part_mol, $type, $id);
## sometimes it seems that the whole mol has not been provided with
## a protein id
        if ($type eq "UP") {
          $pw->AddMolExId($whole_mol, $type, $id);
        }
      }
    }
  }
}

######################################################################
sub ConvertFamily {
  my ($pex, $i, $edge, $family_mol, $members) = @_;

  if (@{ $members } == 0) {
    return $emptyfamilytype{$pex}[$i]{$edge};
  }

  my $seq;
  my $member_mol_type;

  for my $member (@{ $members }) {

    $seq++;
    my ($name, $ext_ids, $activity_labels, $loc_func_labels, $moltype,
        $ptms) = @{ $member };
    my $protein_id;
    if ($moltype eq "protein" || $moltype eq "compound" || $moltype eq "complex"
        || $moltype eq "rna") {
      if (defined $member_mol_type && $member_mol_type ne $moltype) {
        print STDERR "multiple member type for family $name\n";
      } else {
        $member_mol_type = $moltype;
      }
    }
    my $member_mol = $comp_mols{$pex}[$i]{$edge}{$seq};

    my $moltype_lvid = $lv->StringToLabelValue("molecule-type", $moltype);
    if ($moltype_lvid == 0) {
      print STDERR "# no moltype for member $pex, interaction $i, edge $edge " .
          "member $seq\n";
    }
    $pw->AddMol($member_mol, $moltype_lvid);

    $pw->AddFamilyMember($member_mol, $family_mol);

    my ($pf_name, @as_names) = split(";", $name);
    if ($pf_name) {
      $pw->AddMolName($member_mol, "PF", $pf_name);
    }
    for my $as_name (@as_names) {
      if ($as_name) {
        $pw->AddMolName($member_mol, "AS", $as_name);
      }
    }

    for my $x (split(";", $ext_ids)) {
      if ($x) {
        my ($type, $id) = split(":", $x);
        if ($type eq "UP") {
          $protein_id = $id;
        }
        $pw->AddMolExId($member_mol, $type, $id);
      }
    }

    for my $x (split(";", $activity_labels)) {
      if ($x) {
        my $lvid = $lv->StringToLabelValue("activity-state", $x);
        if ($lvid ne 0) {
          $pw->AddMemberLabel($family_mol, $member_mol, $lvid);
        }
      }
    }
    for my $x (@{ ConvertPTM($ptms) }) {
      my ($pos, $aa, $modstring, $modvalue) = @{ $x };
      $pw->AddMemberPTM($family_mol, $member_mol, $protein_id, $pos, $aa,
          $modvalue, $modstring);
    }

  }

  return $member_mol_type;
}

######################################################################
sub ConvertComplex {
  my ($pex, $i, $edge, $cxmol, $components) = @_;

  my $seq;

  for my $comp (@{ $components }) {

    $seq++;
    my ($name, $ext_ids, $activity_labels, $loc_func_labels, $moltype,
        $ptms) = @{ $comp };
    my $protein_id;
    my $compmol = $comp_mols{$pex}[$i]{$edge}{$seq};
    my $moltype_lvid = $lv->StringToLabelValue("molecule-type", $moltype);
    if ($moltype_lvid == 0) {
      print STDERR "# no moltype for $pex, interaction $i, edge $edge, " .
          "component $seq\n";
    }
    $pw->AddMol($compmol, $moltype_lvid);

    $pw->AddComponent($cxmol, $seq, $compmol);

    my ($pf_name, @as_names) = split(";", $name);
    if ($pf_name) {
      $pw->AddMolName($compmol, "PF", $pf_name);
    }
    for my $as_name (@as_names) {
      if ($as_name) {
        $pw->AddMolName($compmol, "AS", $as_name);
      }
    }

    for my $x (split(";", $ext_ids)) {
      if ($x) {
        my ($type, $id) = split(":", $x);
        if ($type eq "UP") {
          $protein_id = $id;
        }
#        print STDERR "##convertcx: $i, $edge, $seq, ext id: $x\n";
        $pw->AddMolExId($compmol, $type, $id);
      }
    }
    for my $x (split(";", $activity_labels)) {
      if ($x) {
#        print STDERR "##convertcx: $i, $edge, $seq, activity label: $x\n";
        my $lvid = $lv->StringToLabelValue("activity-state", $x);
        if ($lvid ne 0) {
          $pw->AddComponentLabel($cxmol, $seq, $lvid);
        }
      }
    }
    for my $x (split(";", $loc_func_labels)) {
      if ($x) {
#        print STDERR "##convertcx: $i, $edge, $seq, loc_func label: $x\n";
        my $loc = $lv->StringToLabelValue("location", $x);
        if ($loc ne 0) {
          $pw->AddComponentLabel($cxmol, $seq, $loc);
          next;
        }
        my $func = $lv->StringToLabelValue("function", $x);
        if ($func ne 0) {
          $pw->AddComponentLabel($cxmol, $seq, $func);
          next;
        }
        print STDERR "##unrecognized loc_func label in complex: " .
            "$pex, $i, $edge, $seq: $x\n";
      }
    }
# print STDERR "about  to call ptm from $i, cx = $cxmol, $seq\n";
    for my $x (@{ ConvertPTM($ptms) }) {
      my ($pos, $aa, $modstring, $modvalue) = @{ $x };
      $pw->AddComponentPTM($cxmol, $seq, $protein_id, $pos, $aa,
          $modvalue, $modstring);
    }
  }
}

######################################################################
sub Convert {

  my %complex_seen;

  for my $pex (keys %mapped_pathways) {

    ## if this pathway was only referenced, not defined, then don't
    ## process it
    if (!defined $pathway{$pex}{name}) {
      next;
    }

    my $SOURCE_NAME  = $pathway{$pex}{source};
    my $SOURCE_VALUE = $idm->SourceId($SOURCE_NAME);
    $pw->AddSource($SOURCE_VALUE, $SOURCE_NAME);

    my $pid = $mapped_pathways{$pex};

    $pw->AddPathwayId($pid);
    $pw->SetPathwayName($pid, $pathway{$pex}{name});
    $pw->SetPathwayExId($pid, $pex);
    $pw->SetPathwayOrg($pid, $pathway{$pex}{organism});
    $pw->SetPathwaySrcId($pid, $SOURCE_VALUE);
    $pw->AddPathwaySource($pid, $SOURCE_VALUE);

    $pw->SetIsSubnet($pid, $pathway{$pex}{issubnet});

    for my $curator (@{ $pathway2curator{$pex} }) {
#      $pw->AddCurator($pid, $curator);
# actually, these are reviewers
      $pw->AddReviewer($pid, $curator);
    }

    if (! defined $interactions{$pex}) {
      next;
    }

    for (my $i = 1; $i < @{ $interactions{$pex} }; $i++) {

      my $atom_id = $mapped_interactions{$pex}[$i];
      my $atom_type = $lv->StringToLabelValue("process-type",
          $interaction2type{$pex}[$i]);

      if ($atom_type == 0 && $interaction2type{$pex}[$i] ne "") {
        print STDERR "## interaction type $interaction2type{$pex}[$i] " .
            "not defined as process-type " .
            "in $pex, interaction $i; " .
            "will look for definition as function\n";
        $atom_type = $lv->StringToLabelValue("function",
            $interaction2type{$pex}[$i]);
      }

      if ($atom_type != 0) {

        if (defined $interaction2abstraction{$pex}[$i]) {
          my $ext_pathway_id = $interaction2abstraction{$pex}[$i];
          my ($pid, $pathway_name) = ("", $pathway{$ext_pathway_id}{name});
          $pw->AddAbstraction($atom_id, $pid, $pathway_name,
              $ext_pathway_id);
        }

        $pw->AddAtom($atom_id, $atom_type);
        $pw->AddAtomSource($atom_id, $SOURCE_VALUE);

        for my $code (keys %{ $interaction2evidence{$pex}{$i} }) {
          $pw->AddAtomEvidence($atom_id, $code);
        }

        for my $ref (@{ $interaction2reference{$pex}{$i} }) {
          $pw->AddAtomReferences($atom_id, $ref);
        }

        for my $note (@{ $interaction2note{$pex}{$i} }) {
          $pw->AddAtomNotes($atom_id, $note);
        }

        for my $condition (keys %{ $interaction2condition{$pex}[$i] }) {
          my $lvid = $lv->StringToLabelValue("process-type", $condition);
          if ($lvid == 0 && $condition ne "") {
            print STDERR "## condition $condition not defined as process-type " .
                "in $pex, interaction $i; " .
                "will look for definition as function\n";
            $lvid = $lv->StringToLabelValue("function", $condition);
            if ($lvid == 0 && $condition ne "") {
              print STDERR "## condition $condition not defined as process-type " .
                  "or as function in $pex, interaction $i\n";
            } else {
              $pw->AddAtomCondition($atom_id, $lvid);
            }
	  } else {
            $pw->AddAtomCondition($atom_id, $lvid);
          }
        }

        for my $edge (sort numerically keys %{ $edge_mols{$pex}[$i] }) {
          my $mol = $edge_mols{$pex}[$i]{$edge};
          my $edge_type = $lv->StringToLabelValue("edge-type",
              $edge_types{$pex}[$i]{$edge});
          $pw->AddEdge($atom_id, $edge, $edge_type, $mol);
        }
      } elsif ($interaction2type{$pex}[$i] ne "") {
        print STDERR "## process type " .
            "$interaction2type{$pex}[$i] not defined " .
                "in $pex, interaction $i\n";
      }

      ## use moldefs here because ext ids, etc. have been normalized
      for my $edge (sort numerically keys %{ $moldefs{$pex}[$i] }) {
        my $protein_id;
        my ($name, $ext_ids, $activity_labels, $loc_func_labels, $moltype,
            $components, $ptms) = @{ $moldefs{$pex}[$i]{$edge} };
        my $mol = $edge_mols{$pex}[$i]{$edge};

        if ($moltype eq "family") {
          my $member_type = ConvertFamily($pex, $i, $edge, $mol, $components);
          my $moltype_lvid = $lv->StringToLabelValue("molecule-type", $member_type);
          if ($moltype_lvid == 0) {
            print STDERR "## no moltype for family $pex, interaction $i, edge $edge\n";
          }
          $pw->AddMol($mol, $moltype_lvid);
        } elsif ($moltype eq "cleave") {
          my $member_type = ConvertCleaved($pex, $i, $edge, $mol, $components);
          my $moltype_lvid =
              $lv->StringToLabelValue("molecule-type", "protein");
          if ($moltype_lvid == 0) {
            print STDERR "# no moltype for cleave $pex, interaction $i, edge $edge\n";
          }
          $pw->AddMol($mol, $moltype_lvid);
        } else {
          if ($moltype eq "complex") {
            if (! defined $complex_seen{$mol}) {
              ConvertComplex($pex, $i, $edge, $mol, $components);
              $complex_seen{$mol} = 1;
            }
          }
          my $moltype_lvid = $lv->StringToLabelValue("molecule-type", $moltype);
          if ($moltype_lvid == 0) {
            print STDERR "# no moltype for $pex, interaction $i, edge $edge\n";
          }
          $pw->AddMol($mol, $moltype_lvid);
        }

        my ($pf_name, @as_names) = split(";", $name);
        if ($pf_name) {
          $pw->AddMolName($mol, "PF", $pf_name);
        }
        for my $as_name (@as_names) {
          if ($as_name) {
            $pw->AddMolName($mol, "AS", $as_name);
          }
        }

        for my $x (split(";", $ext_ids)) {
          if ($x) {
            my ($type, $id) = split(":", $x);
#            print STDERR "##convert: $i, $edge, ext id: $x\n";
            $pw->AddMolExId($mol, $type, $id);
            if ($type eq "UP") {
              $protein_id = $id;
            }
          }
        }

        if ($atom_type != 0) {
          for my $x (split(";", $activity_labels)) {
            if ($x) {
#              print STDERR "##convert: $i, $edge, activity label: $x\n";
              my $lvid = $lv->StringToLabelValue("activity-state", $x);
              if ($lvid ne 0) {
                $pw->AddMolLabel($atom_id, $edge, $lvid);
              }
            }
          }
          for my $x (split(";", $loc_func_labels)) {
            if ($x) {
#              print STDERR "##convert: $i, $edge, loc_func label: $x\n";
              my $loc = $lv->StringToLabelValue("location", $x);
              if ($loc ne 0) {
                $pw->AddMolLabel($atom_id, $edge, $loc);
                next;
              }
              my $func = $lv->StringToLabelValue("function", $x);
              if ($func ne 0) {
                $pw->AddEdgeLabel($atom_id, $edge, $func);
                next;
              }
              print STDERR "##unrecognized loc_func label: $pex, $i, $edge: $x\n";
            }
          }
          if ($moltype ne "complex") {
            if ($ptms) {
              for my $x (@{ ConvertPTM($ptms) }) {
                my ($pos, $aa, $modstring, $modvalue) = @{ $x };
                $pw->AddEdgePTM($atom_id, $edge, $protein_id, $pos, $aa,
                    $modvalue, $modstring);
              }
            }
          }
        }

      }

      $pw->AddAtomPathway($atom_id, $pid);
    }

  }
}

######################################################################
sub ParseLine {
  my ($line) = @_;

  if ($line =~ /^pathway\((.*)\);$/) {
    DoPathway($1);
  } elsif ($line =~ /^function\s+\w+\(\)\s+\{/) {
    if (defined $pathway_id) {
      Transform();
      undef @tmp_interactions;
    }
  } elsif ($line =~ /^curators\((.*)\);$/) {
    DoCurators($1);
  } elsif ($line =~ /^shownamed\((.*)\);$/) {
    DoShownamed($1);
  } else {
    print STDERR "## ignoring: $line\n";
  }
}

######################################################################
sub DoPathway {
  my ($args) = @_;
  my @args = split(",", $args, 1000);
  for my $a (@args) {
    $a =~ s/^\s*'?\s*//;
    $a =~ s/\s*'?\s*$//;
  }
  my ($pid, $name, $organism, $source, $dummy1, $is_subnet) = @args;
  $pathway{$pid}{name}     = $name;
  $pathway{$pid}{organism} = $organism;
  $pathway{$pid}{source}   = $source;
  $pathway{$pid}{dummy1}   = $dummy1;
  if ($is_subnet eq "Subnet") {
    $pathway{$pid}{issubnet} = 1;
  } else {
    $pathway{$pid}{issubnet} = 0;
  }
  $pathway_id = $pid;
}

######################################################################
sub DoCurators {
  my ($curators) = @_;

  for my $c (split(";", $curators, 1000)) {
    $c =~ s/^\s*'?\s*//;
    $c =~ s/\s*'?\s*$//;
    $c =~ s/\s+/ /g;
    if ($c ne "") {
      push @{ $pathway2curator{$pathway_id} }, $c;
    }
  }
}

######################################################################
sub SplitArgs {
  my ($s) = @_;

  my @args;
  while ($s) {
    if ($s =~ /^'([^']*)',/) {
      push @args, $1;
      $s =~ s/^'[^']*',//;
    } elsif ($s =~ /^'([^']*)'$/) {
      push @args, $1;
      $s =~ s/^'[^']*'//;
    } elsif ($s =~ /^(-?\d*),/) {
      push @args, $1;
      $s =~ s/^-?\d*,//;
    } elsif ($s =~ /^(-?\d*)$/) {
      push @args, $1;
      $s =~ s/^-?\d*$//;
    } else {
      print STDERR "can't parse args: $s\n";
      exit;
    }
  }
  for (@args) {
    s/^\s+//;
    s/\s+$//;
    s/\s+\&/\&/g;
    s/\&\s+/\&/g;
    s/\s*;\s*/;/g;
  }
  return \@args;
}

######################################################################
sub DoShownamed {
  my ($args) = @_;

  my ($op, @args) = @{ SplitArgs($args) };
  $op =~ s/\d*Gif$//;

  my ($xcoord, $ycoord) = ($args[4], $args[5]);
  if ($ycoord < 0) {
    print STDERR "## ycoord $ycoord < 0: " . join(",", @args) . "\n";
    return;
  }
  my $interaction_id = int($ycoord/100) + 1;
  $tmp_interactions[$interaction_id]{$xcoord} =
      join("\t", $op, FilterArgs(@args));

# all operators are followed by 13 arguments
# (although there appear to be some "references" that are followed by 14 args).
# arg 5: identifies the reaction (div 100)
# arg 4: appears to identify the left-to-right ordering of reactants

# Associates: args 1-4, 10-12 are null
#           : args 5-6  are integers
#           : arg  7    is 0
#           : arg  8    is 0
#           : arg  9    is -2
#           : arg 13    is 1

# Modification: arg 1   is a molecule name
#             : arg 2   is an external id: CA, CAS, LL, UP, simple UniProt id, in a variety of forms
#             : arg 3   is an activity-state label
#             : arg 4   is a ;-separated list of location,function labels
#             : arg 5-8 are integers (either arg 5 == arg 7, or arg 7 == 0)
#                                    (either arg 6 < arg 8, or arg 8 == 0)
#                                    (when arg 7 == 0, arg 8 also == 0)
#                                    (no obvious relation between arg 5 and arg 6)
#             : arg 9   is -1
#             : arg 10  is ;-separated font stuff
#             : arg 11  is ;-separated: first ;-expr is a &-phrase
#             : arg 12  is ;-separated: null; process term?; complex|compound|protein
#             : arg 13  is 1

# Complex: 
#        : arg  9   is 2
#        : arg 10   is font stuff


  if ($op eq "Modification"      || $op eq "Transcription"   || $op eq "Translocation") {
  } elsif ($op eq "NoModification" || $op eq "NoTranscription" || $op eq "NoTranslocation") {
  } elsif ($op eq "Associates") {
  } elsif ($op eq "MacroProcess" || $op eq "NoMacroProcess") {
  } elsif ($op eq "Pathway" || $op eq "NoPathway" || $op eq "Subnet" || $op eq "NoSubnet") {
## Subnet and Pathway
## Apparently, these have to have an agent or input molecule. So the
## name and alt id refer to this molecule. Perhaps the identity of the
## subnet or pathway would be in the macroprocess subtype field of
## MISC_INFO
  } elsif ($op eq "Complex") {
  } elsif ($op eq "Protein" || $op eq "Compound" || $op eq "RNA") {
  } elsif (defined $evidence_ops{$op}) {
  } elsif ($op eq "Family") {
  } elsif ($op eq "Cleave") {
  } elsif ($op eq "Reference") {
  } elsif ($op eq "Notes") {
  } else {
    print STDERR "unrecognized op $op\n";
  }
}

######################################################################
sub FilterArgs {
  my (@args) = @_;

  my @filtered_args = (
    $args[0],  # molname
    $args[1],  # external ids
    $args[2],  # activity
    $args[3],  # location & function labels
#    $args[4],   unused xcoord
#    $args[5],   unused ycoord
#    $args[6],   unused coord
#    $args[7],   unused coord
    $args[8],  # role cue
#    $args[9],   unused font stuff
    $args[10], # ptms or complex components
    $args[11]  # if an interaction-type node, then molecule  & other info
#    $args[12]  # unused constant 1 ??
  );

  return join("\t", @filtered_args);
}

######################################################################
1;
######################################################################
