package PWLabel;
require Exporter;
@ISA = qw(Exporter);
#@EXPORT = qw(
#);

use strict;

## lv:  label value id (integer)  -> label name \t label value name
##             where label value name is in standard (mixed case) spelling
## rlv: label name \t lowercase(label value name) -> label value id (integer)
## isa: label value id X label value id -> 1 | undef
      ## global isa relation (reflexive and transitive)
## risa: inverse of isa
## labelsort: label value id -> label id (integer) [i.e. the kind of label]
## labelname: label id (integer) -> label name (string)
## parent: label value id (integer) -> immediate parent label value id

######################################################################
sub new {
  my ($self, $db, $schema) = @_;
  my $x = {};
  if (defined $db) {
    GetLabelValues($x, $db, $schema);
    GetLabelValueParents($x, $db, $schema);
    GetLabelExtIds($x, $db, $schema);
    CreateIsA($x);
    ## don't need to CreateLabelSort
  }
  return bless $x;
}

######################################################################
sub AddBasicMolType {
  my ($self, $abbrev, $lvid) = @_;
  $self->{basicmoltype}{$abbrev}= $lvid;
}

######################################################################
sub AddSort {
  my ($self, $label_name, $label_id) = @_;
  $self->{labelname}{$label_id} = $label_name;
}

######################################################################
sub AllSorts {
  my ($self) = @_;

  return [ keys %{ $self->{labelname} } ];
}

######################################################################
sub AllSortPairs {
  my ($self) = @_;

  return $self->{labelname};
}

######################################################################
sub AllNamePairs {
  my ($self) = @_;

## need to return the standard (mixed-case) spellings
##  return $self->{rlv};
  my %temp;
  for my $lvid (values %{ $self->{rlv} }) {
    $temp{$self->{lv}{$lvid} } = $lvid;    
  }
  return \%temp;
}

######################################################################
sub AllUsedIds {
  my ($self) = @_;

  return [values %{ $self->{rlv} } ];
}

######################################################################
sub AllValuesInFamily {
  my ($self, $sort_name, $parent_label_name) = @_;
  
  my $lvid = $self->StringToLabelValue($sort_name, $parent_label_name);  
  return [ keys %{ $self->{risa}{$lvid} } ];
}

######################################################################
sub AddParent {
  my ($self, $child_id, $parent_id) = @_;

  $self->{parent}{$child_id}           = $parent_id;
  $self->{isa}{$child_id}{$parent_id}  = 1;
  $self->{isa}{$child_id}{$child_id}   = 1;
  $self->{isa}{$parent_id}{$parent_id} = 1;
}

######################################################################
sub AddNamePair {
  my ($self, $label_value_id, $label_value_name, $label_name) = @_;

  my $lvn = lc($label_value_name);
  $self->{lv}{$label_value_id} = "$label_name\t$label_value_name";
  $self->{rlv}{"$label_name\t$lvn"} = $label_value_id;
  $self->{isa}{$label_value_id}{$label_value_id} = 1;
}

######################################################################
sub ParentOf {
  my ($self, $child) = @_;

  return $self->{parent}{$child};
}

######################################################################
sub AncestorsOf {
  my ($self, $child) = @_;

  return [ keys %{ $self->{isa}{$child} } ];

  my %ancs;
  for (my $x = $child; ; $x = $self->ParentOf($x)) {
    $ancs{$x} = 1;
    if ($x eq $self->ParentOf($x)) {
      last;
    }
  }
  return [ keys %ancs ];
}

######################################################################
sub IsA {
  my ($self, $child, $parent) = @_;
  if (defined $self->{isa}{$child}{$parent}) {
    return 1;
  } else {
    return 0;
  }
}

######################################################################
sub FindLowestAncestor {
  my ($self, $child, $parent_set) = @_;

  ##
  ## child may or may not be a member of parent_set
  ##

  my ($p, $q, $i, $j, $ci, $cj, @candidates, %failures);

  for $p (@{ $parent_set }) {
    if ($p == $child) {
      next;
    }
    if ($self->IsA($child, $p)) {
      push @candidates, $p;
    } 
  }
  if (@candidates == 0) {
    return $p;
  } elsif (@candidates == 1) {
    return $candidates[0];
  } else {
    ##
    ## elminate any candidate that is a parent of another candidate
    ##
    for ($i = 0; $i < @candidates - 1; $i++) {
      for ($j = $i + 1; $j < @candidates; $j++) {
        $ci = $candidates[$i];
        $cj = $candidates[$j];
        if ($ci == $cj) {
          next;
        }
        if ($self->IsA($ci, $cj)) {
          $failures{$cj} = 1;
        } elsif ($self->IsA($cj, $ci)) {
          $failures{$ci} = 1;
        }
      }
    }
    for $q (@candidates) {
      if (not defined $failures{$q}) {
        return $q;
      }
    }
    ##
    ## should never get here !!!
    ##
    return $p;
  }
}

######################################################################
sub LabelSort {
  my ($self, $lvid) = @_;
  if (defined $self->{labelsort}{$lvid}) {
    return $self->{labelsort}{$lvid};
  } else {
    print STDERR "LabelSort: unrecognized label value id $lvid\n";
    return "";
  }
}

######################################################################
sub LabelString {
  my ($self, $sortid) = @_;
  if (defined $self->{labelname}{$sortid}) {
    return $self->{labelname}{$sortid};
  } else {
#    print STDERR "LabelString: unrecognized label value id $sortid\n";
    return "";
  }
}

######################################################################
sub StringToLabelValue {
  my ($self, $label_name, $label_value_name) = @_;

  my $lvn = lc($label_value_name);
  if (defined $self->{rlv}{"$label_name\t$lvn"}) {
    return $self->{rlv}{"$label_name\t$lvn"};
  } elsif ($label_value_name =~ /^GO:\d+$/ &&
      ($label_name eq "function" || $label_name eq "process-type" ||
      $label_name eq "location")) {
    if (defined $self->{extids}{$label_value_name}) {
      return $self->{extids}{$label_value_name}
    } else {
      return 0;
    }
  } else {
    return 0;
  }
}

######################################################################
sub LabelValueToString {
  my ($self, $lvid) = @_;
  if (defined $self->{lv}{$lvid}) {
    my ($label_name, $label_value_name) = split "\t", $self->{lv}{$lvid};
    return ($label_name, $label_value_name);
  } else {
#    print STDERR "request for unrecognized label value id $lvid\n";
#    return ("", "");
  }
}

######################################################################
sub BasicMolTypeCode {
  my ($self, $basic_type_str) = @_;
  if (defined $self->{basicmoltype}{$basic_type_str}) {
    return $self->{basicmoltype}{$basic_type_str};
  } else {
    print STDERR "unrecognized basic mol type string $basic_type_str\n";
    return 0;
  }
}

######################################################################
sub BasicMolType {
  my ($self, $label_value_id) = @_;

  my ($name);
  for $name (keys %{ $self->{basicmoltype} }) {
    if ($label_value_id == $self->{basicmoltype}{$name}) {
      return $name;
    }
  }
  return "??";
}

######################################################################
sub CreateIsA {
  my ($self) = @_;

  my ($c, $p, $g);

  my $done = 0;
  while (not $done) {
    $done = 1;
    for $c (keys %{ $self->{isa} }) {
      for $p (keys %{ $self->{isa}{$c} }) {
        for $g (keys %{ $self->{isa}{$p} }) {
          if (! defined $self->{isa}{$c}{$g}) {
            $done = 0;
            $self->{isa}{$c}{$g} = 1;
          }
        }
      }
    }
  }
  for my $c (keys %{ $self->{isa} }) {
    for my $p (keys %{ $self->{isa}{$c} }) {
      $self->{risa}{$p}{$c} = 1;
    }
  }
}

######################################################################
sub CreateLabelSort {
  my ($self) = @_;

  my %name2id;
  my ($label_name, $label_value_name, $pair);

  for my $label_id (keys %{ $self->{labelname} }) {
    $name2id{$self->{labelname}{$label_id}} = $label_id;
  }
  for my $label_value_id (keys %{ $self->{lv} }) {
    $pair = $self->{lv}{$label_value_id};
    ($label_name, $label_value_name) = split("\t", $pair);
    $self->{labelsort}{$label_value_id} = $name2id{$label_name};
  }
}


######################################################################
sub GetLabelValues {
  my ($self, $db, $SCHEMA) = @_;

  my ($sql, $stm);
  my ($label_value_id, $label_id, $label_name, $label_value_name);

  $sql = qq!
select
  lv.label_value_id, l.label_id, l.label_name, lv.label_value_name
from
  $SCHEMA.pw_label_value lv,
  $SCHEMA.pw_label l
where
  lv.label_id = l.label_id
!;

  $stm = $db->prepare($sql);
  if(not $stm) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "prepare call failed\n";
    die;
  }
  if(!$stm->execute()) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "execute call failed\n";
    die;
  }
  while (($label_value_id, $label_id, $label_name, $label_value_name) =
      $stm->fetchrow_array()) {
    $self->{labelsort}{$label_value_id} = $label_id;
    $self->{lv}{$label_value_id} = "$label_name\t$label_value_name";
##    $self->{rlv}{"$label_name\t$label_value_name"} = $label_value_id;
## temporary fix while juggling old/new label values that clash
##
    my $lvn = lc($label_value_name);
    if (defined $self->{rlv}{"$label_name\t$lvn"}) {
      my $id = $self->{rlv}{"$label_name\t$lvn"};
      if ($id < $label_value_id) {
        $self->{rlv}{"$label_name\t$lvn" . "[old]"} = $id;
        $self->{rlv}{"$label_name\t$lvn"} = $label_value_id;
      } elsif ($label_value_id < $id)  {
        $self->{rlv}{"$label_name\t$lvn" . "[old]"} = $label_value_id;
        $self->{rlv}{"$label_name\t$lvn"} = $id;
      }
    } else {
      $self->{rlv}{"$label_name\t$lvn"} = $label_value_id;
    }   
    $self->{labelname}{$label_id} = $label_name;
    if ($label_name eq "molecule-type") {
      if ($label_value_name eq "protein") {
        $self->{basicmoltype}{"PR"} = $label_value_id;
      } elsif ($label_value_name eq "compound") {
        $self->{basicmoltype}{"CM"} = $label_value_id;
      } elsif ($label_value_name eq "complex") {
        $self->{basicmoltype}{"CX"} = $label_value_id;
      } elsif ($label_value_name eq "rna") {
        $self->{basicmoltype}{"RN"} = $label_value_id;
      } elsif ($label_value_name eq "molecule-type") {
        $self->{basicmoltype}{"MO"} = $label_value_id;
      }
    }
  }
}

######################################################################
sub GetLabelExtIds {
  my ($self, $db, $SCHEMA) = @_;

  my ($sql, $stm);
  my ($lvid, $ext_id, $id_type);

  $sql = qq!
select
  p.label_value_id, p.ext_id, p.id_type
from
  $SCHEMA.pw_ext_lv_id p
where
  p.id_type = 'GO'
!;

  $stm = $db->prepare($sql);
  if(not $stm) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "prepare call failed\n";
    die;
  }
  if(!$stm->execute()) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "execute call failed\n";
    die;
  }
  while (($lvid, $ext_id, $id_type) =
      $stm->fetchrow_array()) {
    $self->{extids}{$ext_id} = $lvid;
    $self->{go}{$lvid}{$ext_id} = 1;
  }
}

######################################################################
sub AddGOTerm {
  my ($self, $lvid, $go) = @_;

  $self->{go}{$lvid}{$go} = 1;
  $self->{extids}{$go} = $lvid;
}

######################################################################
sub GOTermsFor {
  my ($self, $lvid) = @_;

  if (defined $self->{go}{$lvid}) {
    return [ keys %{ $self->{go}{$lvid} } ];
  } else {
    return [];
  }
}

######################################################################
sub GetLabelValueParents {
  my ($self, $db, $SCHEMA) = @_;

  my ($sql, $stm);
  my ($child_value_id, $parent_value_id);
  my ($c, $p, $g, $done);

  $sql = qq!
select
  p.child_value_id, p.parent_value_id
from
  $SCHEMA.pw_label_value_parent p
!;

  $stm = $db->prepare($sql);
  if(not $stm) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "prepare call failed\n";
    die;
  }
  if(!$stm->execute()) {
    print STDERR "$sql\n";
    print STDERR "$DBI::errstr\n";
    print STDERR "execute call failed\n";
    die;
  }
  while (($child_value_id, $parent_value_id) =
      $stm->fetchrow_array()) {
    $self->{isa}{$child_value_id}{$parent_value_id} = 1;
    $self->{isa}{$child_value_id}{$child_value_id} = 1;
    $self->{isa}{$parent_value_id}{$parent_value_id} = 1;
    $self->{parent}{$child_value_id} = $parent_value_id;
  }

}

######################################################################
1;


