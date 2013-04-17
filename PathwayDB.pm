

# Copyright SRA International
#
# Distributed under the OSI-approved BSD 3-Clause License.
# See http://ncip.github.com/pathway-interaction-database/LICENSE.txt for details.


package PathwayDB;
require Exporter;
@ISA = qw(Exporter);
#@EXPORT = qw(
#);

use strict;
use Pathway;
use DBI;

use constant ORACLE_LIST_LIMIT => 500;
use constant MAX_LONG_LEN       => 16384;

my $INCOMING;
my $OUTGOING;
my $pathway_type;
my $subnet_type;

######################################################################
sub new {
  my ($self, $db, $schema, $source_id_list, $evidence_code_list, $pw, $lv) = @_;
  my $x = {};
  $x->{db} = $db;
  $x->{pw} = $pw;
  $x->{lv} = $lv;
  $x->{schema} = $schema;
  $db->{LongReadLen} = MAX_LONG_LEN;
  for my $src_id (split(",", $source_id_list)) {
    $x->{source_id}{$src_id} = 1;
  }

  # if ($evidence_code_list eq '') {
  #  $evidence_code_list = 'NIL,IC,IFC,IMP,IPI,RGE,IAE,IDA,IGI,IOS,RCA,TAS';
  #}

  for my $code (split(",", $evidence_code_list)) {
    $x->{evidence_code}{$code} = 1;
  }
  InitializeSource($x);
  InitializePartOf($x);
  InitializeFamilyMember($x);
  $pathway_type = $lv->StringToLabelValue("process-type", "pathway");
  $subnet_type = $lv->StringToLabelValue("process-type", "subnet");
  return bless $x;
}

######################################################################
sub InitializeSource {
  my ($self) = @_;

  my $db     = $self->{db};
  my $pw     = $self->{pw};
  my $lv     = $self->{lv};
  my $SCHEMA = $self->{schema};

  my ($sql, $stm);
  my ($source_id, $source_name);
  my ($atom_id, $pathway_id, $evidence_code);
  my ($source_spec, $evidence_spec);
  if (defined $self->{source_id}) {
    $source_spec = 1;
  }
  if (defined $self->{evidence_code}) {
    $evidence_spec = 1;
  }

  $sql = qq!
select
  s.source_id,
  s.source_name
from
  $SCHEMA.pw_source s
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
  while (($source_id, $source_name) = $stm->fetchrow_array()) {
    $pw->AddSource($source_id, $source_name);
  }


  $sql = qq!
select
  a.atom_id,
  a.atom_source_id
from
  $SCHEMA.pw_atom a
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
  while (($atom_id, $source_id) = $stm->fetchrow_array()) {
#    $pw->AddAtomSource($atom_id, $source_id);
    if ($source_spec && (defined $self->{source_id}{$source_id})) {
      $self->{atomsinsource}{$atom_id} = 1;    
    }
  }

  $sql = qq!
select
  v.atom_id,
  v.evidence_code
from
  $SCHEMA.pw_evidence v
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
  while (($atom_id, $evidence_code) = $stm->fetchrow_array()) {
#    $pw->AddAtomEvidence($atom_id, $evidence_code);
    if ($evidence_spec && (defined $self->{evidence_code}{$evidence_code})) {
      $self->{atomsinevidence}{$atom_id} = 1;    
    }
  }

if ($evidence_spec && (defined $self->{evidence_code}{NIL})) {  
  $sql = qq!
select 
  atom_id
from
  $SCHEMA.pw_atom
where
  atom_id not in 
    (select atom_id from $SCHEMA.pw_evidence)
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

  while (($atom_id) = $stm->fetchrow_array()) {
        $self->{atomsinevidence}{$atom_id} = 1;
  }

}
  
$sql = qq!
select
  p.pathway_id,
  p.pathway_source_id
from
  $SCHEMA.pw_pathway p
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
  while (($pathway_id, $source_id) = $stm->fetchrow_array()) {
    $pw->AddPathwaySource($pathway_id, $source_id);
  }


}

######################################################################
sub InitializePartOf {
  my ($self) = @_;

  my $db     = $self->{db};
  my $pw     = $self->{pw};
  my $lv     = $self->{lv};
  my $SCHEMA = $self->{schema};

  my ($sql, $stm);
  my ($part_mol_id, $whole_mol_id, $part_mol_type, $whole_mol_type,
      $start_site, $end_site);
  my ($mol_id, $name_type, $mol_name);
  my ($ext_mol_id, $id_type);

  ##
  ## Part-of and basic mol type
  ##

  $sql = qq!
select
  p.part_mol_id,
  p.whole_mol_id,
  mp.basic_mol_type,
  mw.basic_mol_type,
  p.start_site,
  p.end_site
from
  $SCHEMA.pw_mol_part p,
  $SCHEMA.pw_mol mp,
  $SCHEMA.pw_mol mw
where
      p.part_mol_id = mp.mol_id
  and p.whole_mol_id = mw.mol_id
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
  while (($part_mol_id, $whole_mol_id, $part_mol_type, $whole_mol_type,
      $start_site, $end_site)
      = $stm->fetchrow_array()) {
    $pw->AddMol($part_mol_id, $lv->BasicMolTypeCode($part_mol_type));
    $pw->AddMol($whole_mol_id, $lv->BasicMolTypeCode($whole_mol_type));
    $pw->AddMolPart($part_mol_id, $whole_mol_id, $start_site, $end_site);
  }

  ##
  ## mol name
  ##

  $sql = qq!
select
  n.mol_id,
  n.name_type,
  n.mol_name
from
  $SCHEMA.pw_mol_name n,
  $SCHEMA.pw_mol_part p
where
  n.mol_id = p.part_mol_id or n.mol_id = p.whole_mol_id
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
  while (($mol_id, $name_type, $mol_name)
      = $stm->fetchrow_array()) {
    $pw->AddMolName($mol_id, $name_type, $mol_name);
  }

  ##
  ## mol ext id
  ##

  $sql = qq!
select
  e.mol_id,
  e.id_type,
  e.ext_mol_id
from
  $SCHEMA.pw_ext_mol_id e,
  $SCHEMA.pw_mol_part p
where
  e.mol_id = p.part_mol_id or e.mol_id = p.whole_mol_id
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
  while (($mol_id, $id_type, $ext_mol_id)
      = $stm->fetchrow_array()) {
    $pw->AddMolExId($mol_id, $id_type, $ext_mol_id);
  }

}

######################################################################
sub InitializeFamilyMember {
  my ($self) = @_;

  my $db     = $self->{db};
  my $pw     = $self->{pw};
  my $lv     = $self->{lv};
  my $SCHEMA = $self->{schema};

  my ($sql, $stm);
  my ($member_mol_id, $family_mol_id, $member_mol_type, $family_mol_type);
  my ($mol_id, $name_type, $mol_name);
  my ($ext_mol_id, $id_type);

  ##
  ## Family-member and basic mol type
  ##

  $sql = qq!
select
  f.member_mol_id,
  f.family_mol_id,
  mm.basic_mol_type,
  mf.basic_mol_type
from
  $SCHEMA.pw_family_member f,
  $SCHEMA.pw_mol mm,
  $SCHEMA.pw_mol mf
where
      f.member_mol_id = mm.mol_id
  and f.family_mol_id = mf.mol_id
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
  while (($member_mol_id, $family_mol_id, $member_mol_type, $family_mol_type)
      = $stm->fetchrow_array()) {
    $pw->AddMol($member_mol_id, $lv->BasicMolTypeCode($member_mol_type));
    $pw->AddMol($family_mol_id, $lv->BasicMolTypeCode($family_mol_type));
    $pw->AddFamilyMember($member_mol_id, $family_mol_id);
  }

  ##
  ## Family-member modification
  ##

  my ($ptm_expression_id, $uniprot_acc, $mod_type, $residue, $pos);

  $sql = qq!
select
  m.family_mol_id,
  m.member_mol_id,
  m.ptm_expression_id,
  p.uniprot_acc,
  p.ptm_type,
  p.ptm_residue,
  p.ptm_pos
from
  $SCHEMA.pw_member_modification m,
  $SCHEMA.pw_ptm_expression p
where
      m.ptm_expression_id = p.ptm_expression_id
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
  while (($family_mol_id, $member_mol_id, $ptm_expression_id, 
      $uniprot_acc, $mod_type, $residue, $pos) = $stm->fetchrow_array()) {
    my ($label, $value)  = $lv->LabelValueToString($mod_type);
    $pw->AddMemberPTM($family_mol_id, $member_mol_id, $uniprot_acc,
        $pos, $residue, $mod_type, $value);
  }

  ##
  ## Family-member label
  ##

  my ($lvid);

  $sql = qq!
select
  m.family_mol_id,
  m.member_mol_id,
  m.label_value_id
from
  $SCHEMA.pw_member_label m
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
  while (($family_mol_id, $member_mol_id, $lvid) = $stm->fetchrow_array()) {
    $pw->AddMemberLabel($family_mol_id, $member_mol_id, $lvid);
  }

  ##
  ## mol name
  ##

  $sql = qq!
select
  n.mol_id,
  n.name_type,
  n.mol_name
from
  $SCHEMA.pw_mol_name n,
  $SCHEMA.pw_family_member f
where
  n.mol_id = f.member_mol_id or n.mol_id = f.family_mol_id
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

  while (($mol_id, $name_type, $mol_name)
      = $stm->fetchrow_array()) {
    $pw->AddMolName($mol_id, $name_type, $mol_name);
  }

  ##
  ## mol ext id
  ##

  $sql = qq!
select
  e.mol_id,
  e.id_type,
  e.ext_mol_id
from
  $SCHEMA.pw_ext_mol_id e,
  $SCHEMA.pw_family_member f
where
  e.mol_id = f.member_mol_id or e.mol_id = f.family_mol_id
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
  while (($mol_id, $id_type, $ext_mol_id)
      = $stm->fetchrow_array()) {
    $pw->AddMolExId($mol_id, $id_type, $ext_mol_id);
  }

}


######################################################################
######################################################################
##
##    spec by degree
##
######################################################################
######################################################################
sub AtomsOfDegree {
  my ($self, $direction, $degree, $prune_mol_hash,
      $prune_atom_hash) = @_;

  my $db     = $self->{db};
  my $pw     = $self->{pw};
  my $lv     = $self->{lv};
  my $SCHEMA = $self->{schema};

  if ($degree < 1) {
    return;
  }

  ## !!! Initial fill of $pw has been pruned for mols,
  ## but duplicate atoms have not been detected and pruned

  $INCOMING = $lv->StringToLabelValue("edge-type", "incoming-edge");
  $OUTGOING = $lv->StringToLabelValue("edge-type", "outgoing-edge");

  my ($atoms_in, $mols_in, %mols_out, $prune_mol_list);
  my ($layer, $atom_id);
  my (%consumed, %produced, $producedCount);

  if ($prune_mol_hash) {
    $prune_mol_list = join(",", keys %{ $prune_mol_hash });
  }

  $atoms_in = {};
  my $atomidlist;
  for $atom_id (@{ $pw->Atoms() }) {
    $$atoms_in{$atom_id} = 1;
    $atomidlist = $atomidlist . $atom_id;
  }

  FindInsAndOuts($self, $pw->Atoms(), \%produced, \%consumed);

  $producedCount = scalar(keys %produced);
  
  for ($layer = 1; $layer <= $degree; $layer++) {
#print STDERR "layer = $layer\n";
    if (scalar(keys %consumed) > 0 &&
        ($direction eq "back" || $direction eq "both")) {
      $mols_in = join(",", keys %consumed);
      undef %consumed;
      OneDegreeDirection($self, "back", $prune_mol_list,
          $prune_mol_hash, $prune_atom_hash,
          $mols_in, \%consumed, $atoms_in, $layer);
#print STDERR "direction = back, mols_in = $mols_in, mols_out = " .
#join(",", keys %consumed) . "\n";
    }

    if (scalar(keys %produced) > 0 &&
        ($direction eq "forward" || $direction eq "both")) {
      
      $mols_in = join(",", keys %produced);
      undef %produced;
      my $retval = OneDegreeDirection($self, "forward", $prune_mol_list,
          $prune_mol_hash, $prune_atom_hash,
          $mols_in, \%produced, $atoms_in, $layer);

      # return $retval;
#print STDERR "direction = forward, mols_in = $mols_in, mols_out = " .
#join(",", keys %produced) . "\n";
    }
  }
}

######################################################################
sub OneDegreeDirection {
  my ($self, $direction, $prune_mol_list, $prune_mol_hash,
      $prune_atom_hash,
      $mols_in, $mols_out, $atoms_in, $layer) = @_;

  my $db     = $self->{db};
  my $pw     = $self->{pw};
  my $lv     = $self->{lv};
  my $SCHEMA = $self->{schema};

  my ($atom_id, %new_atoms, %pro, %con);
  my $atomidlist;

  $atomidlist = NewAtomsOfDegree($self, $direction, $prune_mol_list,
      $prune_atom_hash, $mols_in,
      $atoms_in, \%new_atoms, $pw, $lv);
  # return "atomidlist:  $atomidlist";
  if (scalar(keys %new_atoms) == 0) {
    return 'none';
  }
  FillAtoms($self, join(",", keys %new_atoms));
  if ($prune_mol_hash && scalar(keys %{ $prune_mol_hash }) > 0) {
    $pw->PruneMol($prune_mol_hash);
  }
  FindInsAndOuts($self, [keys %new_atoms], \%pro, \%con);
  if ($direction eq "back" || $direction eq "both") {
    for (keys %con) {
      $$mols_out{$_} = 1;
    }
  }
  if ($direction eq "forward" || $direction eq "both") {
    for (keys %pro) {
      $$mols_out{$_} = 1;
    }
  }
  for $atom_id (keys %new_atoms) {
    $pw->SetLayer($atom_id, $layer);
    $$atoms_in{$atom_id} = 1;
  }
  
  return 'some';
}

######################################################################
sub FindInsAndOuts {
  my ($self, $atom_set, $produced, $consumed) = @_;

  my $db     = $self->{db};
  my $pw     = $self->{pw};
  my $lv     = $self->{lv};
  my $SCHEMA = $self->{schema};

  my ($atom_id, $edge, $label);

  for $atom_id (@{ $atom_set }) {
    for $edge (@{ $pw->Edges($atom_id) }) {
      for $label (@{ $pw->EdgeLabel($atom_id, $edge) }) {
        if ($lv->IsA($label, $INCOMING)) {
          $$consumed{$pw->EdgeMol($atom_id, $edge)} = 1;
        }
        if ($lv->IsA($label, $OUTGOING)) {
          $$produced{$pw->EdgeMol($atom_id, $edge)} = 1;
        }
      }
    }
  }
}

######################################################################
sub DirectionsMatch {
  my ($self, $spec, $a_lab, $b_lab) = @_;

  my $db     = $self->{db};
  my $pw     = $self->{pw};
  my $lv     = $self->{lv};
  my $SCHEMA = $self->{schema};

  ## point of view is from "a": "back" means from b to a
  ## (hence b should have an outgoing edge and a an incoming)

  if ($spec eq "both") {
    return 1;
  }

  my ($a_lv, $b_lv);
  my ($a_in, $a_out, $b_in, $b_out) = (0, 0, 0, 0);

  for $a_lv (@{ $a_lab }) {
    if ($lv->IsA($a_lv, $INCOMING)) {
      $a_in = 1;
    } elsif ($lv->IsA($a_lv, $OUTGOING)) {
      $a_out = 1;
    }
  }
  for $b_lv (@{ $b_lab }) {
    if ($lv->IsA($b_lv, $INCOMING)) {
      $b_in = 1;
    } elsif ($lv->IsA($b_lv, $OUTGOING)) {
      $b_out = 1;
    }
  }

  if ($spec eq "back") {
    if ($a_in && $b_out) {
      return 1;
    } else {
      return 0;
    }
  } elsif ($spec eq "forward") {
    if ($a_out && $b_in) {
      return 1;
    } else {
      return 0;
    }
  }

}

######################################################################
sub NewAtomsOfDegree {
  my ($self, $direction, $prune_mol_list, $prune_atom_hash,
      $mols_in, $atoms_in,
      $atoms_out) = @_;

  my $db     = $self->{db};
  my $pw     = $self->{pw};
  my $lv     = $self->{lv};
  my $SCHEMA = $self->{schema};
  my $pw2    = $self->{pw};

  ## atoms_in are already added to $pw

  my ($sql, $stm, $stm2);
  my ($atom_id, $edge_seq_id, $mol_id, $mol_lv, $edge_lv, $atomidlist, $ptm_id);
  my (%mol_inst, %edge_inst, %ptm_inst);
  my ($i, @mol_list, $list, $prune_clause);
  my ($ptm_expression_id, $uniprot_acc, $mod_type, $residue, $pos);

  if ($prune_mol_list) {
    $prune_clause = "and e.mol_id not in ($prune_mol_list)";
1  } else {
    $prune_clause = "";
  }

  ##
  ## First, fetch all atoms that share a mol_id with the seed list
  ##

  @mol_list = split ",", $mols_in;

  for($i = 0; $i < @mol_list; $i += ORACLE_LIST_LIMIT) {
 
    if(($i + ORACLE_LIST_LIMIT - 1) < @mol_list) {
      $list = join(",", @mol_list[$i..$i+ORACLE_LIST_LIMIT-1]);
    }
    else {
      $list = join(",", @mol_list[$i..@mol_list-1]);
    }
    $sql = qq!
select
  e.atom_id, e.edge_seq_id, e.mol_id, m.label_value_id, el.label_value_id
from
  $SCHEMA.pw_edge e,
  $SCHEMA.pw_mol_label m,
  $SCHEMA.pw_edge_label el
where
      e.atom_id = m.atom_id (+)
  and e.edge_seq_id = m.edge_seq_id (+)
  and e.atom_id = el.atom_id
  and e.edge_seq_id = el.edge_seq_id
  and e.mol_id in ($list)
  $prune_clause
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
    while (($atom_id, $edge_seq_id, $mol_id, $mol_lv, $edge_lv) =
        $stm->fetchrow_array()) {
      if ($pw->HasAtom($atom_id)) {
        next;
      }
#      if ((defined $self->{atomsinsource}) &&
      if ((defined $self->{source_id}) &&
          !(defined $self->{atomsinsource}{$atom_id})) {
        next;
      }
#      if ((defined $self->{atomsinevidence}) &&
      if ((defined $self->{evidence_code}) &&
          !(defined $self->{atomsinevidence}{$atom_id})) {
        next;
      }
      if ($prune_atom_hash && defined $$prune_atom_hash{$atom_id}) {
        next;
      }
      $mol_inst{$atom_id}{$edge_seq_id}{$mol_id}{$mol_lv} = 1;
      $edge_inst{$atom_id}{$edge_seq_id}{$edge_lv} = $mol_id;

      } 
    
  }

  my ($candidate_atom, $candidate_edge, $candidate_inst);
  my ($inst, $hit);

  for $atom_id (keys %edge_inst) {
    $hit = 0;
    for $edge_seq_id (keys %{ $mol_inst{$atom_id} }) {
      for $mol_id (keys %{ $mol_inst{$atom_id}{$edge_seq_id} }) {

        my $ptms = $ptm_inst{$atom_id}{$edge_seq_id}; 
         $inst = Pathway::MolInstToString($mol_id,
            [ keys %{ $mol_inst{$atom_id}{$edge_seq_id}{$mol_id} } ]);
 
        # $inst = $pw->NormalMolInstString($atom_id, $edge_seq_id);

        $sql = qq!
        select
          pep.ptm_type,
          pep.ptm_residue,
          pep.ptm_pos
        from
          $SCHEMA.pw_edge_modification pem,
          $SCHEMA.pw_ptm_expression pep
        where
          atom_id = $atom_id
          and edge_id = $edge_seq_id
          and pem.ptm_expression_id = pep.ptm_expression_id
        !;

        $stm2 = $db->prepare($sql);
        if(not $stm2) {
          print STDERR "$sql\n";
          print STDERR "$DBI::errstr\n";
          print STDERR "prepare call failed\n";
          die;
        }
        if(!$stm2->execute()) {
          print STDERR "$sql\n";
          print STDERR "$DBI::errstr\n";
          print STDERR "execute call failed\n";
          die;
        }
        my $ptm_string;
        my $first = 1;
        while (($mod_type, $residue, $pos) =
          $stm2->fetchrow_array()) {
          if ($first == 1) {
            $first = 0;
            $ptm_string = "$pos$residue$mod_type";
          } else {
            $ptm_string = $ptm_string . "," . "$pos$residue$mod_type"; 
          }
        }
  
        $inst = $inst . $ptm_string;
 
        for $candidate_atom (@{ $pw->MolUse($mol_id) }) {
          if (! defined $$atoms_in{$candidate_atom}) {
            next;
          }
          for $candidate_edge (@{ $pw->Edges($candidate_atom) }) {
              
            $candidate_inst = $pw->NormalMolInstString(
              $candidate_atom, $candidate_edge);
            # $atomidlist = $atomidlist .  "PTM:  $ptm_string Atom: $atom_id Edge: $edge_seq_id Cand: $candidate_atom $candidate_edge Inst: $inst Cand: $candidate_inst<BR>";
            if ($inst eq $candidate_inst) {
              if (DirectionsMatch($self, $direction,
                  $pw->EdgeLabel($candidate_atom, $candidate_edge),
                  [ keys %{ $edge_inst{$atom_id}{$edge_seq_id} } ])) {
                $hit = 1;
                # $atomidlist = $atomidlist . "HIT!";
              } else {
              }
            }
          }
        }
      }
    }
    if ($hit) {
      $$atoms_out{$atom_id} = 1;
    }
  }
  return $atomidlist;
}

######################################################################
######################################################################
##
##    spec by molecules
##
######################################################################
######################################################################


######################################################################
sub RecursiveGetContainingComplexes {
  my ($self, $seed_mol_list, $cumulative_mol_list, $added_mol_list) = @_;

  my $db     = $self->{db};
  my $pw     = $self->{pw};
  my $lv     = $self->{lv};
  my $SCHEMA = $self->{schema};

  my ($sql, $stm, $list, $complex_mol_id, $i);

  for($i = 0; $i < @{ $seed_mol_list }; $i += ORACLE_LIST_LIMIT) {
 
    if(($i + ORACLE_LIST_LIMIT - 1) < @{ $seed_mol_list }) {
      $list = join(",", @{ $seed_mol_list }[$i..$i+ORACLE_LIST_LIMIT-1]);
    }
    else {
      $list = join(",", @{ $seed_mol_list }[$i..@{ $seed_mol_list } - 1]);
    }

  $sql = qq!
select
  c.complex_mol_id
from
  $SCHEMA.pw_complex_component c
where
  c.component_mol_id in ($list)
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
    while (($complex_mol_id) = $stm->fetchrow_array()) {
      if (! defined $$cumulative_mol_list{$complex_mol_id}) {
        $$cumulative_mol_list{$complex_mol_id} = 1;
        $$added_mol_list{$complex_mol_id} = 1;
      }
    }
  }
}

######################################################################
sub GetContainingComplexes {
  my ($self, $mol_list_in, $cumulative_mol_list) = @_;

  my $db     = $self->{db};
  my $pw     = $self->{pw};
  my $lv     = $self->{lv};
  my $SCHEMA = $self->{schema};

  my $seed_mol_list = $mol_list_in;
  my %added_mol_list;
  for my $mol_id (@{ $seed_mol_list }) {
    $$cumulative_mol_list{$mol_id} = 1;
  }
  while (1) {
    RecursiveGetContainingComplexes($self, $seed_mol_list, $cumulative_mol_list,
        \%added_mol_list);
    if (scalar(keys %added_mol_list) == 0) {
      last;
    }
    $seed_mol_list = [ keys %added_mol_list ];
    undef %added_mol_list;
  }
}

######################################################################
sub AtomsOfMols {
  my ($self, $mol_list, $include_complex_uses) = @_;

  my $db     = $self->{db};
  my $pw     = $self->{pw};
  my $lv     = $self->{lv};
  my $SCHEMA = $self->{schema};

  my ($sql, $stm);
  my ($atom_id, %atoms, $i, $list);
  my @mol_list = split ",", $mol_list;

  if ($include_complex_uses) {
    my %cumulative_mol_list;
    GetContainingComplexes($self, \@mol_list, \%cumulative_mol_list);
    @mol_list = keys %cumulative_mol_list;
  }

  for($i = 0; $i < @mol_list; $i += ORACLE_LIST_LIMIT) {
 
    if(($i + ORACLE_LIST_LIMIT - 1) < @mol_list) {
      $list = join(",", @mol_list[$i..$i+ORACLE_LIST_LIMIT-1]);
    }
    else {
      $list = join(",", @mol_list[$i..@mol_list-1]);
    }

  $sql = qq!
select
  e.atom_id
from
  $SCHEMA.pw_edge e
where
  e.mol_id in ($list)
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
    while (($atom_id) = $stm->fetchrow_array()) {
      $atoms{$atom_id} = 1;
    }
  }

  FillAtoms($self, join(",", keys %atoms));

}

######################################################################
######################################################################
##
##    MacroProcess
##
######################################################################
######################################################################

######################################################################
sub FillMacroProcesses {
  my ($self, $lvid_list) = @_;

  my $db     = $self->{db};
  my $pw     = $self->{pw};
  my $lv     = $self->{lv};
  my $SCHEMA = $self->{schema};

  my ($sql, $stm);
  my ($lvid, $atom_id, %atom_ids);

  if ($lvid_list eq "") {
    return;
  }

  $sql = qq!
select
  atom_id
from
  $SCHEMA.pw_atom_condition
where
  label_value_id in ($lvid_list)
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
  while (($atom_id) = $stm->fetchrow_array()) {

    if ((defined $self->{source_id}) &&
        !(defined $self->{atomsinsource}{$atom_id})) {
      next;
    }
    if ((defined $self->{evidence_code}) &&
        !(defined $self->{atomsinevidence}{$atom_id})) {
      next;
    }

    $atom_ids{$atom_id} = 1;
  }

  $sql = qq!
select
  atom_id
from
  $SCHEMA.pw_atom_label
where
  label_value_id in ($lvid_list)
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
  while (($atom_id) = $stm->fetchrow_array()) {

    if ((defined $self->{source_id}) &&
        !(defined $self->{atomsinsource}{$atom_id})) {
      next;
    }
    if ((defined $self->{evidence_code}) &&
        !(defined $self->{atomsinevidence}{$atom_id})) {
      next;
    }

    $atom_ids{$atom_id} = 1;

  }

  $self->FillAtoms(join(",", keys %atom_ids));

}
######################################################################
sub oldFillMacroProcesses {
  my ($self, $name_list) = @_;

  my $db     = $self->{db};
  my $pw     = $self->{pw};
  my $lv     = $self->{lv};
  my $SCHEMA = $self->{schema};

  my ($sql, $stm);

  my ($id, %lvids, %atom_ids, $lvid_list, $atom_id);

  my $mps = $lv->AllValuesInFamily("process-type", "macroprocess");
  for my $mp (@{ $mps }) {
    my ($dummy, $name) = $lv->LabelValueToString($mp);
    for my $probe (@{ $name_list }) {
      if ($name =~ /$probe/i) {
        $lvids{$mp} =1;
      }
    }
  }

#  my $macro_process_id =
#    $lv->StringToLabelValue("process-type", "macroprocess");
#  for (@{ $name_list }) {
#    $id = $lv->StringToLabelValue("process-type", $_);
#    if ($lv->IsA($id, $macro_process_id)) {
#      $lvids{$id} = 1;
#    }
#  }

  $lvid_list = join(",", keys %lvids);

  if ($lvid_list eq "") {
    return;
  }

  $sql = qq!
select
  atom_id
from
  $SCHEMA.pw_atom_condition
where
  label_value_id in ($lvid_list)
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
  while (($atom_id) = $stm->fetchrow_array()) {

    if ((defined $self->{source_id}) &&
        !(defined $self->{atomsinsource}{$atom_id})) {
      next;
    }
    if ((defined $self->{evidence_code}) &&
        !(defined $self->{atomsinevidence}{$atom_id})) {
      next;
    }

    $atom_ids{$atom_id} = 1;
  }


  $sql = qq!
select
  atom_id
from
  $SCHEMA.pw_atom_label
where
  label_value_id in ($lvid_list)
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
  while (($atom_id) = $stm->fetchrow_array()) {

    if ((defined $self->{source_id}) &&
        !(defined $self->{atomsinsource}{$atom_id})) {
      next;
    }
    if ((defined $self->{evidence_code}) &&
        !(defined $self->{atomsinevidence}{$atom_id})) {
      next;
    }

    $atom_ids{$atom_id} = 1;

  }

  $self->FillAtoms(join(",", keys %atom_ids));

}

######################################################################
######################################################################
##
##    Atoms only
##
######################################################################
######################################################################


######################################################################
sub FillAtoms {
  my ($self, $atom_list) = @_;

  my $db     = $self->{db};
  my $pw     = $self->{pw};
  my $lv     = $self->{lv};
  my $SCHEMA = $self->{schema};

  my (%temp, $a);

#  if (defined $self->{atomsinsource}) {
  if (defined $self->{source_id}) {
    if ($atom_list eq "*") {
#      if (defined $self->{atomsinevidence}) {
      if (defined $self->{evidence_code}) {
        $atom_list = join(",", keys %{ $self->{atomsinsource} });
        for $a (split(",", $atom_list)) {
          if (defined $self->{atomsinevidence}{$a}) {
            $temp{$a} = 1;
          }
        }
        $atom_list = join(",", keys %temp);
      } else {
        $atom_list = join(",", keys %{ $self->{atomsinsource} });
      }
    } else {
      for $a (split(",", $atom_list)) {
        if (defined $self->{atomsinsource}{$a}) {
#          if (defined $self->{atomsinevidence}) {
          if (defined $self->{evidence_code}) {
            if (defined $self->{atomsinevidence}{$a}) {
              $temp{$a} = 1;
            }
          } else {
            $temp{$a} = 1;
          }
        }
      }
      $atom_list = join(",", keys %temp);
    }
#  } elsif (defined $self->{atomsinevidence}) {
  } elsif (defined $self->{evidence_code}) {
    if ($atom_list eq "*") {
      $atom_list = join(",", keys %{ $self->{atomsinevidence} });
    } else {
      for $a (split(",", $atom_list)) {
        if (defined $self->{atomsinevidence}{$a}) {
          $temp{$a} = 1;
        }
      }
      $atom_list = join(",", keys %temp);
    }
  }
  if ($atom_list eq "*") {
    GetEdgesForAtoms($self, "*");
    GetAtomsAndLabelsForAtoms($self, "*");
    GetMoleculesForAtoms($self, "*");
    GetEvidenceForAtoms($self, "*");
    GetReferencesForAtoms($self, "*");
    GetNotesForAtoms($self, "*");
  } else {
    my @atom_list = split(",", $atom_list);
    my ($i, $list);

    for($i = 0; $i < @atom_list; $i += ORACLE_LIST_LIMIT) {

      if(($i + ORACLE_LIST_LIMIT - 1) < @atom_list) {
        $list = join(",", @atom_list[$i..$i+ORACLE_LIST_LIMIT-1]);
      }
      else {
        $list = join(",", @atom_list[$i..@atom_list-1]);
      }

      GetEdgesForAtoms($self, $list);
      GetAtomsAndLabelsForAtoms($self, $list);
      GetMoleculesForAtoms($self, $list);
      GetPathwaysForAtoms($self, $list);
      GetEvidenceForAtoms($self, $list);
      GetReferencesForAtoms($self, $list);
      GetNotesForAtoms($self, $list);
    }
  }
}

######################################################################
sub GetAtomBasics {
  my ($self, $atom_id) = @_;

  my $db     = $self->{db};
  my $pw     = $self->{pw};
  my $lv     = $self->{lv};
  my $SCHEMA = $self->{schema};
  my ($organism, $atom_source_id, $source_name);
  my ($sql, $stm);

  $sql = qq!
select
  a.organism,
  a.atom_source_id,
  s.source_name
from
  $SCHEMA.pw_atom a,
  $SCHEMA.pw_source s
where
      a.atom_id = $atom_id
  and a.atom_source_id = s.source_id
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
  while (($organism, $atom_source_id, $source_name) =
      $stm->fetchrow_array()) {
    last;
  }
  $stm->finish();
  return ($organism, $atom_source_id, $source_name);
}

######################################################################
sub GetAbstractions {
  my ($self, $atoms) = @_;

  my $db     = $self->{db};
  my $pw = $self->{pw};
  my $lv = $self->{lv};
  my $SCHEMA = $self->{schema};

  if (keys %{ $atoms } < 1) {
    return;
  }
  my $abs = join(",", keys %{ $atoms });
  my ($sql, $stm);
  my ($atom_id, $pathway_id, $pathway_name, $ext_pathway_id);

  $sql = qq!
select distinct
    a.atom_id, p.pathway_id, p.pathway_name, p.ext_pathway_id
from
  $SCHEMA.pw_abstraction a,
  $SCHEMA.pw_pathway p
where
      a.atom_id in ($abs)
  and (p.pathway_id = a.pathway_id or
       p.ext_pathway_id = a.ext_pathway_id
      )
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
  while (($atom_id, $pathway_id, $pathway_name, $ext_pathway_id) =
      $stm->fetchrow_array()) {
    $pw->AddAbstraction($atom_id, $pathway_id, $pathway_name,
        $ext_pathway_id);
  }
}

######################################################################
sub GetAtomsAndLabelsForAtoms {
  my ($self, $atom_list) = @_;

  my $db     = $self->{db};
  my $pw     = $self->{pw};
  my $lv     = $self->{lv};
  my $SCHEMA = $self->{schema};

  my ($sql, $stm);
  my ($atom_id, $label_value_id, $atom_source_id, $negative_flag);
  my ($name, $value);

  my ($atom_clause);
  if ($atom_list eq "*") {
    $atom_clause = "";
  } else {
    $atom_clause = "and al.atom_id in ($atom_list)";
  }

  $sql = qq!
select
  al.atom_id, al.label_value_id, a.atom_source_id
from
  $SCHEMA.pw_atom_label al,
  $SCHEMA.pw_atom a
where
      a.atom_id = al.atom_id
$atom_clause
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
  my %abstractions;
  while (($atom_id, $label_value_id, $atom_source_id) =
      $stm->fetchrow_array()) {
    $pw->AddAtomLabel($atom_id, $label_value_id);
##    ($name, $value) = split /\t/, $lvhash{$label_value_id};
    ($name, $value) = $lv->LabelValueToString($label_value_id);
    if ($name eq "process-type") {
      $pw->AddAtom($atom_id, $label_value_id);
      if ($label_value_id eq $pathway_type ||
          $label_value_id eq $subnet_type) {
        $abstractions{$atom_id} = 1;
      }
      $pw->AddAtomSource($atom_id, $atom_source_id);
    }
  }
  if (keys %abstractions) {
    $self->GetAbstractions(\%abstractions);
  }


  if ($atom_list eq "*") {
    $atom_clause = "";
  } else {
    $atom_clause = "where al.atom_id in ($atom_list)";
  }

  ##
  ## Now conditions
  ##

  $sql = qq!
select
  al.atom_id, al.label_value_id, al.negative_flag
from
  $SCHEMA.pw_atom_condition al
$atom_clause
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
  while (($atom_id, $label_value_id, $negative_flag) =
      $stm->fetchrow_array()) {
    if ($negative_flag == 1) {
      $pw->AddAtomNegativeCondition($atom_id, $label_value_id);
    } else {
      $pw->AddAtomCondition($atom_id, $label_value_id);
    }
  }

}

######################################################################
sub GetEdgesForAtoms {
  my ($self, $atom_list) = @_;

  my $db     = $self->{db};
  my $pw     = $self->{pw};
  my $lv     = $self->{lv};
  my $SCHEMA = $self->{schema};

  my ($sql, $stm);
  my ($atom_id, $edge_seq_id, $mol_id, $edge_type);
  my ($label_value_id);
  my ($ptm_expression_id, $uniprot_acc, $mod_type, $residue, $pos);

  my ($atom_clause);
  if ($atom_list eq "*") {
    $atom_clause = "";
  } else {
    $atom_clause = "and e.atom_id in ($atom_list)";
  }

##!!! Assumes every edge is labeled with edge-type

  $sql = qq!
select
  e.atom_id, e.edge_seq_id, e.mol_id, el.label_value_id
from
  $SCHEMA.pw_edge e,
  $SCHEMA.pw_edge_label el,
  $SCHEMA.pw_label_value lv,
  $SCHEMA.pw_label l
where
      el.atom_id = e.atom_id
  $atom_clause
  and el.edge_seq_id = e.edge_seq_id
  and l.label_name = 'edge-type'
  and l.label_id = lv.label_id
  and lv.label_value_id = el.label_value_id
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
  while (($atom_id, $edge_seq_id, $mol_id, $edge_type) =
      $stm->fetchrow_array()) {
    if ($atom_id && $edge_seq_id && $mol_id && $edge_type) {
      $pw->AddEdge($atom_id, $edge_seq_id, $edge_type, $mol_id);
    }
  }

  ##
  ## Now edge labels
  ##
  if ($atom_list eq "*") {
  $sql = qq!
select
  el.atom_id, el.edge_seq_id, el.label_value_id
from
  $SCHEMA.pw_edge_label el
!;
  } else {
  $sql = qq!
select
  el.atom_id, el.edge_seq_id, el.label_value_id
from
  $SCHEMA.pw_edge_label el
  where el.atom_id in ($atom_list)
!;
  }

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
  while (($atom_id, $edge_seq_id, $label_value_id) =
      $stm->fetchrow_array()) {
    $pw->AddEdgeLabel($atom_id, $edge_seq_id, $label_value_id);
  }

  ##
  ## Now mol labels (i.e., labels on molecule instances)
  ##
  if ($atom_list eq "*") {
  $sql = qq!
select
  ml.atom_id, ml.edge_seq_id, ml.label_value_id
from
  $SCHEMA.pw_mol_label ml
!;
  } else {
  $sql = qq!
select
  ml.atom_id, ml.edge_seq_id, ml.label_value_id
from
  $SCHEMA.pw_mol_label ml
where
  ml.atom_id in ($atom_list)
!;
  }

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
  while (($atom_id, $edge_seq_id, $label_value_id) =
      $stm->fetchrow_array()) {
    if ($atom_id && $edge_seq_id && $label_value_id) {
      $pw->AddMolLabel($atom_id, $edge_seq_id, $label_value_id);
    }
  }

  ##
  ## Now ptms
  ## (for now, just a conjunction of PTM terms)
  ##

  if ($atom_list eq "*") {
    $sql = qq!
select
  m.atom_id,
  m.edge_id,
  p.ptm_expression_id,
  p.uniprot_acc,
  p.ptm_type,
  p.ptm_residue,
  p.ptm_pos
from
  $SCHEMA.pw_edge_modification m,
  $SCHEMA.pw_ptm_expression p
where
      m.ptm_expression_id = p.ptm_expression_id
    !;
  } else {
    $sql = qq!
select
  m.atom_id,
  m.edge_id,
  p.ptm_expression_id,
  p.uniprot_acc,
  p.ptm_type,
  p.ptm_residue,
  p.ptm_pos
from
  $SCHEMA.pw_edge_modification m,
  $SCHEMA.pw_ptm_expression p
where
      m.ptm_expression_id = p.ptm_expression_id
  and m.atom_id in ( $atom_list )
    !;
  }

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

  while (($atom_id, $edge_seq_id, $ptm_expression_id, $uniprot_acc,
      $mod_type, $residue,  $pos) = $stm->fetchrow_array()) {
    my ($label, $value)  = $lv->LabelValueToString($mod_type);
    $pw->AddEdgePTM($atom_id, $edge_seq_id, $uniprot_acc, $pos, $residue,
        $mod_type, $value);
  }

}

######################################################################
sub GetMoleculesForAtoms {
  my ($self, $atom_list) = @_;

  my $db     = $self->{db};
  my $pw     = $self->{pw};
  my $lv     = $self->{lv};
  my $SCHEMA = $self->{schema};

  my ($sql, $stm);
  my ($mol_id, $basic_mol_type, $ext_mol_id, $id_type, $mol_name, $name_type);
  my (%complexes);

  my ($atom_clause);
  if ($atom_list eq "*") {
    $atom_clause = "";
  } else {
    $atom_clause = "and e.atom_id in ($atom_list)";
  }

  ##
  ## basic mol type and external ids
  ##

  $sql = qq!
select unique
  m.mol_id, m.basic_mol_type, em.ext_mol_id, em.id_type
from
  $SCHEMA.pw_edge e,
  $SCHEMA.pw_ext_mol_id em,
  $SCHEMA.pw_mol m
where
      em.mol_id (+) = e.mol_id
  and m.mol_id = e.mol_id
  $atom_clause
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
  while (($mol_id, $basic_mol_type, $ext_mol_id, $id_type) =
      $stm->fetchrow_array()) {
    $pw->AddMol($mol_id, $lv->BasicMolTypeCode($basic_mol_type));
    if ($id_type ne "") {
      $pw->AddMolExId($mol_id, $id_type, $ext_mol_id);
    }
    ##
    ## complexes
    ##
    if ($basic_mol_type eq "CX") {
      $complexes{$mol_id} = 1;
    }
  }

  ##
  ## mol names
  ##

  $sql = qq!
select unique
  m.mol_id, m.mol_name, m.name_type
from
  $SCHEMA.pw_edge e,
  $SCHEMA.pw_mol_name m
where
      m.mol_id = e.mol_id
  $atom_clause
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
  while (($mol_id, $mol_name, $name_type) =
      $stm->fetchrow_array()) {
    $pw->AddMolName($mol_id, $name_type, $mol_name);
  }

  ##
  ## complexes
  ##

  GetComplexes($self, \%complexes);

}

######################################################################
sub GetPathwaysForAtoms {
  my ($self, $atom_list) = @_;

  my $db     = $self->{db};
  my $pw     = $self->{pw};
  my $lv     = $self->{lv};
  my $SCHEMA = $self->{schema};

  my ($sql, $stm);
  my ($atom_id, $pathway_id);

  if ($atom_list eq "*") {
    $sql = qq!
select
  atom_id, pathway_id
from
  $SCHEMA.pw_pathway_atom
    !;
  } else {
    $sql = qq!
select
  atom_id, pathway_id
from
  $SCHEMA.pw_pathway_atom
where
  atom_id in ($atom_list)
    !;
  }
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
  while (($atom_id, $pathway_id) = $stm->fetchrow_array()) {
#    if ((defined $self->{atomsinsource}) &&
    if ((defined $self->{source_id}) &&
        !(defined $self->{atomsinsource}{$atom_id})) {
      next;
    }
#    if ((defined $self->{atomsinevidence}) &&
    if ((defined $self->{evidence_code}) &&
        !(defined $self->{atomsinevidence}{$atom_id})) {
      next;
    }
    $pw->AddAtomPathway($atom_id, $pathway_id);
  }
}

######################################################################
sub GetEvidenceForAtoms {
  my ($self, $atom_list) = @_;

  my $db     = $self->{db};
  my $pw     = $self->{pw};
  my $lv     = $self->{lv};
  my $SCHEMA = $self->{schema};

  my ($sql, $stm);
  my ($atom_id, $evidence_code);

  if ($atom_list eq "*") {
    $sql = qq!
 select
   atom_id, evidence_code
 from
  $SCHEMA.pw_evidence
    !;
  } else {
    $sql = qq!
 select
  atom_id, evidence_code
 from
  $SCHEMA.pw_evidence
 where
  atom_id in ($atom_list)
    !;
  }
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
  while (($atom_id, $evidence_code) = $stm->fetchrow_array()) {
    $pw->AddAtomEvidence($atom_id, $evidence_code);
  }
}

######################################################################
sub GetReferencesForAtoms {
  my ($self, $atom_list) = @_;

  my $db     = $self->{db};
  my $pw     = $self->{pw};
  my $lv     = $self->{lv};
  my $SCHEMA = $self->{schema};

  my ($sql, $stm);
  my ($atom_id, $pmid);

  if ($atom_list eq "*") {
    $sql = qq!
 select
   atom_id, pmid
 from
  $SCHEMA.pw_references
    !;
  } else {
    $sql = qq!
 select
  atom_id, pmid
 from
  $SCHEMA.pw_references
 where
  atom_id in ($atom_list)
    !;
  }

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
  while (($atom_id, $pmid) = $stm->fetchrow_array()) {
    $pw->AddAtomReferences($atom_id, $pmid);
  }
}

######################################################################
sub GetNotesForAtoms {
  my ($self, $atom_list) = @_;

  my $db     = $self->{db};
  my $pw     = $self->{pw};
  my $lv     = $self->{lv};
  my $SCHEMA = $self->{schema};

  my ($sql, $stm);
  my ($atom_id, $note);

  if ($atom_list eq "*") {
    $sql = qq!
 select
   atom_id, note
 from
  $SCHEMA.pw_notes
    !;
  } else {
    $sql = qq!
 select
  atom_id, note
 from
  $SCHEMA.pw_notes
 where
  atom_id in ($atom_list)
    !;
  }

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
  while (($atom_id, $note) = $stm->fetchrow_array()) {
    $pw->AddAtomNotes($atom_id, $note);
  }
}

######################################################################
######################################################################
##
##    Whole pathway
##
######################################################################
######################################################################


######################################################################

sub FillPathway {
  my ($self, $pid) = @_;

  my $db     = $self->{db};
  my $pw     = $self->{pw};
  my $lv     = $self->{lv};
  my $SCHEMA = $self->{schema};

  GetPathwayInfo($self, $pid);
  GetEdgesForPathway($self, $pid);
  GetAtomsAndLabelsForPathway($self, $pid);
  GetMoleculesForPathway($self, $pid);
  GetCuratorsForPathway($self, $pid);

}

######################################################################
sub GetPathwayInfo {
  my ($self, $pid) = @_;

  my $db     = $self->{db};
  my $pw     = $self->{pw};
  my $lv     = $self->{lv};
  my $SCHEMA = $self->{schema};

  my ($sql, $stm);
  my ($pathway_name, $organism, $ext_pathway_id, $pathway_source_id,
      $subnet);


  $sql = qq!
select
  p.pathway_name, p.organism, p.ext_pathway_id,
  p.pathway_source_id, p.subnet
from
  $SCHEMA.pw_pathway p
where
      p.pathway_id = $pid
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
  $stm->bind_columns(\$pathway_name, \$organism, \$ext_pathway_id,
      \$pathway_source_id, \$subnet);
  while($stm->fetch) {
  }

  $pw->SetPathwayName($pid, $pathway_name);
  $pw->SetPathwayExId($pid, $ext_pathway_id);
  $pw->SetPathwayOrg($pid, $organism);
  $pw->SetPathwaySrcId($pid, $pathway_source_id);
  $pw->SetIsSubnet($pid, $subnet);

}

######################################################################
sub GetEdgesForPathway {
  my ($self, $pid) = @_;

  my $db     = $self->{db};
  my $pw     = $self->{pw};
  my $lv     = $self->{lv};
  my $SCHEMA = $self->{schema};

  my ($sql, $stm);
  my ($atom_id, $edge_seq_id, $mol_id, $edge_type);
  my ($label_value_id);
  my ($ptm_expression_id, $uniprot_acc, $mod_type, $residue, $pos);

##!!! Assumes every edge is labeled with edge-type

  $sql = qq!
select
  e.atom_id, e.edge_seq_id, e.mol_id, el.label_value_id
from
  $SCHEMA.pw_edge e,
  $SCHEMA.pw_pathway_atom pa,
  $SCHEMA.pw_edge_label el,
  $SCHEMA.pw_label_value lv,
  $SCHEMA.pw_label l
where
      pa.pathway_id = $pid
  and el.atom_id = pa.atom_id
  and pa.atom_id = e.atom_id
  and el.edge_seq_id = e.edge_seq_id
  and l.label_name = 'edge-type'
  and l.label_id = lv.label_id
  and lv.label_value_id = el.label_value_id
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
  while (($atom_id, $edge_seq_id, $mol_id, $edge_type) =
      $stm->fetchrow_array()) {

    if ((defined $self->{source_id}) &&
        !(defined $self->{atomsinsource}{$atom_id})) {
      next;
    }
    if ((defined $self->{evidence_code}) &&
        !(defined $self->{atomsinevidence}{$atom_id})) {
      next;
    }

    if ($atom_id && $edge_seq_id && $mol_id && $edge_type) {
      $pw->AddEdge($atom_id, $edge_seq_id, $edge_type, $mol_id);
    }
  }

  ##
  ## Now edge labels
  ##

  $sql = qq!
select
  el.atom_id, el.edge_seq_id, el.label_value_id
from
  $SCHEMA.pw_pathway_atom pa,
  $SCHEMA.pw_edge_label el
where
      pa.pathway_id = $pid
  and el.atom_id = pa.atom_id
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
  while (($atom_id, $edge_seq_id, $label_value_id) =
      $stm->fetchrow_array()) {

    if ((defined $self->{source_id}) &&
        !(defined $self->{atomsinsource}{$atom_id})) {
      next;
    }
    if ((defined $self->{evidence_code}) &&
        !(defined $self->{atomsinevidence}{$atom_id})) {
      next;
    }

    $pw->AddEdgeLabel($atom_id, $edge_seq_id, $label_value_id);
  }

  ##
  ## Now mol labels (i.e., labels on molecule instances)
  ##

  $sql = qq!
select
  ml.atom_id, ml.edge_seq_id, ml.label_value_id
from
  $SCHEMA.pw_pathway_atom pa,
  $SCHEMA.pw_mol_label ml
where
      pa.pathway_id = $pid
  and ml.atom_id = pa.atom_id
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
  while (($atom_id, $edge_seq_id, $label_value_id) =
      $stm->fetchrow_array()) {

    if ((defined $self->{source_id}) &&
        !(defined $self->{atomsinsource}{$atom_id})) {
      next;
    }
    if ((defined $self->{evidence_code}) &&
        !(defined $self->{atomsinevidence}{$atom_id})) {
      next;
    }

    if ($atom_id && $edge_seq_id && $label_value_id) {
      $pw->AddMolLabel($atom_id, $edge_seq_id, $label_value_id);
    }
  }

  ##
  ## Now ptms
  ## (for now, just a conjunction of PTM terms)
  ##

  $sql = qq!
select
  m.atom_id,
  m.edge_id,
  p.ptm_expression_id,
  p.uniprot_acc,
  p.ptm_type,
  p.ptm_residue,
  p.ptm_pos
from
  $SCHEMA.pw_edge_modification m,
  $SCHEMA.pw_ptm_expression p,
  $SCHEMA.pw_pathway_atom pa
where
      m.ptm_expression_id = p.ptm_expression_id
  and pa.atom_id = m.atom_id
  and pa.pathway_id = $pid
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
#print STDERR "$sql\n";
  while (($atom_id, $edge_seq_id, $ptm_expression_id,
      $uniprot_acc, $mod_type, $residue,
      $pos) = $stm->fetchrow_array()) {
#print STDERR "in fetch\n";
    my ($label, $value)  = $lv->LabelValueToString($mod_type);

    if ((defined $self->{source_id}) &&
        !(defined $self->{atomsinsource}{$atom_id})) {
      next;
    }
    if ((defined $self->{evidence_code}) &&
        !(defined $self->{atomsinevidence}{$atom_id})) {
      next;
    }

    $pw->AddEdgePTM($atom_id, $edge_seq_id, $uniprot_acc, $pos, $residue,
        $mod_type, $value);
  }

}

######################################################################
sub GetPathwayID {
  my ($self, $ext_pname) = @_;

  my $db     = $self->{db};
  my $pw     = $self->{pw};
  my $lv     = $self->{lv};
  my $SCHEMA = $self->{schema};

  my ($sql, $stm) = @_;
  my ($pid, $pname, $org, $external_pwid, $pw_srcid, $is_subnet);

  $sql = qq!
select
  p.pathway_id, p.pathway_name, p.organism, p.ext_pathway_id,
  p.pathway_source_id, p.subnet
from
  $SCHEMA.pw_pathway p
where
      p.ext_pathway_id = '$ext_pname'
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
  $stm->bind_columns(\$pid, \$pname, \$org, \$external_pwid, \$pw_srcid,
    \$is_subnet);
  while($stm->fetch) {
  }

  $pw->SetPathwayName($pid, $pname);
  $pw->SetPathwayExId($pid, $external_pwid);
  $pw->SetPathwayOrg($pid, $org);
  $pw->SetPathwaySrcId($pid, $pw_srcid);
  $pw->SetIsSubnet($pid, $is_subnet);

}

######################################################################
sub GetMoleculesForSource {
  my ($self, $source_id, $molmap) = @_;

  my $db     = $self->{db};
  my $pw     = $self->{pw};
  my $lv     = $self->{lv};
  my $SCHEMA = $self->{schema};

  my ($sql, $stm);
  my ($mol_id, $basic_mol_type, $ext_mol_id, $id_type, $mol_name, $name_type);
  my (%complexes);

  ##
  ## basic mol type and external ids
  ##

  $sql = qq!
select unique
  m.mol_id, m.basic_mol_type, em.ext_mol_id, em.id_type
from
  $SCHEMA.pw_ext_mol_id em,
  $SCHEMA.pw_mol m,
  $SCHEMA.pw_mol_source ms
where
      em.mol_id (+) = ms.mol_id
  and m.mol_id = ms.mol_id
  and ms.mol_source_id in ($source_id)
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
  while (($mol_id, $basic_mol_type, $ext_mol_id, $id_type) =
      $stm->fetchrow_array()) {
    $pw->AddMol($mol_id, $lv->BasicMolTypeCode($basic_mol_type));
    if ($id_type ne "") {
      $pw->AddMolExId($mol_id, $id_type, $ext_mol_id);
    }

    ##
    ## complexes
    ##
    if ($basic_mol_type eq "CX") {
      $complexes{$mol_id} = 1;
    }
  }

  ##
  ## mol names
  ##

  $sql = qq!
select unique
  m.mol_id, m.mol_name, m.name_type
from
  $SCHEMA.pw_mol_name m,
  $SCHEMA.pw_mol_source ms
where
      m.mol_id = ms.mol_id
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
  while (($mol_id, $mol_name, $name_type) =
      $stm->fetchrow_array()) {
    $pw->AddMolName($mol_id, $name_type, $mol_name);
  }

  ##
  ## complexes
  ##

  GetComplexes($self, \%complexes);

  ##
  ## Get original mol-ids (per mol source) if available
  ##

  my ($original_mol_id, $src_id);

  $sql = qq!
select unique
  m.original_mol_id, m.mol_id, m.mol_source_id
from
  $SCHEMA.pw_mol_map m
where
      m.mol_source_id in ($source_id) 
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
  while (($original_mol_id, $mol_id, $src_id) =
      $stm->fetchrow_array()) {
      $$molmap{$mol_id}{$original_mol_id} = $src_id;
  }
  
}

######################################################################
sub GetAllMolecules {
  my ($self) = @_;

  my $db     = $self->{db};
  my $pw     = $self->{pw};
  my $lv     = $self->{lv};
  my $SCHEMA = $self->{schema};

  my ($sql, $stm);
  my ($mol_id, $basic_mol_type, $ext_mol_id, $id_type, $mol_name, $name_type);
  my (%complexes);

  ##
  ## basic mol type and external ids
  ##

  $sql = qq!
select unique
  m.mol_id, m.basic_mol_type, em.ext_mol_id, em.id_type
from
  $SCHEMA.pw_ext_mol_id em,
  $SCHEMA.pw_mol m
where
      em.mol_id (+) = m.mol_id
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
  while (($mol_id, $basic_mol_type, $ext_mol_id, $id_type) =
      $stm->fetchrow_array()) {
    $pw->AddMol($mol_id, $lv->BasicMolTypeCode($basic_mol_type));
    if ($id_type ne "") {
      $pw->AddMolExId($mol_id, $id_type, $ext_mol_id);
    }

    ##
    ## complexes
    ##
    if ($basic_mol_type eq "CX") {
      $complexes{$mol_id} = 1;
    }
  }

  ##
  ## mol names
  ##

  $sql = qq!
select unique
  m.mol_id, m.mol_name, m.name_type
from
  $SCHEMA.pw_mol_name m
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
  while (($mol_id, $mol_name, $name_type) =
      $stm->fetchrow_array()) {
    $pw->AddMolName($mol_id, $name_type, $mol_name);
  }

  ##
  ## complexes
  ##

  GetComplexes($self, \%complexes);
  
}

######################################################################
sub GetMoleculesForPathway {
  my ($self, $pid) = @_;

  my $db     = $self->{db};
  my $pw     = $self->{pw};
  my $lv     = $self->{lv};
  my $SCHEMA = $self->{schema};

  my ($sql, $stm);
  my ($mol_id, $basic_mol_type, $ext_mol_id, $id_type, $mol_name, $name_type);
  my (%complexes);

  ##
  ## basic mol type and external ids
  ##

  $sql = qq!
select unique
  m.mol_id, m.basic_mol_type, em.ext_mol_id, em.id_type
from
  $SCHEMA.pw_edge e,
  $SCHEMA.pw_pathway_atom pa,
  $SCHEMA.pw_ext_mol_id em,
  $SCHEMA.pw_mol m
where
      pa.pathway_id = $pid
  and e.atom_id = pa.atom_id
  and em.mol_id (+) = e.mol_id
  and m.mol_id = e.mol_id
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
  while (($mol_id, $basic_mol_type, $ext_mol_id, $id_type) =
      $stm->fetchrow_array()) {
    $pw->AddMol($mol_id, $lv->BasicMolTypeCode($basic_mol_type));
    if ($id_type ne "") {
      $pw->AddMolExId($mol_id, $id_type, $ext_mol_id);
    }
    ##
    ## complexes
    ##
    if ($basic_mol_type eq "CX") {
      $complexes{$mol_id} = 1;
    }
  }

  ##
  ## mol names
  ##

  $sql = qq!
select unique
  m.mol_id, m.mol_name, m.name_type
from
  $SCHEMA.pw_edge e,
  $SCHEMA.pw_pathway_atom pa,
  $SCHEMA.pw_mol_name m
where
      pa.pathway_id = $pid
  and e.atom_id = pa.atom_id
  and m.mol_id = e.mol_id
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
  while (($mol_id, $mol_name, $name_type) =
      $stm->fetchrow_array()) {
    $pw->AddMolName($mol_id, $name_type, $mol_name);
  }

  ##
  ## complexes
  ##

  GetComplexes($self, \%complexes);
  
}

######################################################################
sub GetComplexes {
  my ($self, $complexes) = @_;

  my $db     = $self->{db};
  my $pw     = $self->{pw};
  my $lv     = $self->{lv};
  my $SCHEMA = $self->{schema};

  my (%temp, %new_complexes);

  for (keys %{ $complexes }) {
    $temp{$_} = 1;
  }

  while (1) {
    undef %new_complexes;
    GetComplexComponents($self, \%temp, \%new_complexes);
    if (scalar(keys %new_complexes) == 0) {
      last;
    }
    undef %temp;
    for (keys %new_complexes) {
      if (defined $$complexes{$_}) {
        ##
        ## actually, not necessarily an error: consider a pathway
        ## containing complex A = "B:C" and complex C = "D:E". But
        ## let's not visit it multiple times.
        ##
      } else {
        $$complexes{$_} = 1;
        $temp{$_} = 1;
      }
    }
  }

}


######################################################################
sub GetComplexComponents {
  my ($self, $complexes, $new_complexes) = @_;

  my $db     = $self->{db};
  my $pw     = $self->{pw};
  my $lv     = $self->{lv};
  my $SCHEMA = $self->{schema};

  ##
  ## on entry, the mol_id and mol_labels for the complex itself has
  ## already been added via $pw->AddMol
  ##

  my ($sql, $stm);
  my ($i, $list);
  my @mol_list = keys %{ $complexes };
  my ($cx_molid, $cp_seqid, $cp_molid, $cp_lvid);
  my ($mol_id, $basic_mol_type, $ext_mol_id, $id_type);
  my ($mol_name, $name_type);
  my (%sub_molids);
  my ($ptm_expression_id, $uniprot_acc, $pos, $residue, $mod_type, $value);

  for($i = 0; $i < @mol_list; $i += ORACLE_LIST_LIMIT) {
 
    if(($i + ORACLE_LIST_LIMIT - 1) < @mol_list) {
      $list = join(",", @mol_list[$i..$i+ORACLE_LIST_LIMIT-1]);
    }
    else {
      $list = join(",", @mol_list[$i..@mol_list-1]);
    }

  ##
  ## First get component ids and their labelings
  ##

    $sql = qq!
select unique
  cc.complex_mol_id,
  cc.component_seq_id,
  cc.component_mol_id,
  cl.label_value_id
from
  $SCHEMA.pw_complex_component cc, 
  $SCHEMA.pw_component_labeling cl
where
      cc.complex_mol_id = cl.complex_mol_id (+)
  and cc.component_seq_id = cl.component_seq_id (+)
  and cc.complex_mol_id in ($list)
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

    while (($cx_molid, $cp_seqid, $cp_molid, $cp_lvid) =
        $stm->fetchrow_array()) {
      $pw->AddComponent($cx_molid, $cp_seqid, $cp_molid);
      if ($cp_lvid ne "") {
        $pw->AddComponentLabel($cx_molid, $cp_seqid, $cp_lvid);
      }
      $sub_molids{$cp_molid} = 1;
    }

  ##
  ## Now ptms
  ## (for now, just a conjunction of PTM terms)
  ##

    $sql = qq!
select
  c.complex_mol_id,
  c.component_id,
  p.ptm_expression_id,
  p.uniprot_acc,
  p.ptm_type,
  p.ptm_residue,
  p.ptm_pos
from
  $SCHEMA.pw_component_modification c,
  $SCHEMA.pw_ptm_expression p
where
      c.ptm_expression_id = p.ptm_expression_id
  and c.complex_mol_id in ( $list )
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
#print STDERR "$sql\n";
    while (($cx_molid, $cp_seqid, $ptm_expression_id,
        $uniprot_acc, $mod_type, $residue,
        $pos) = $stm->fetchrow_array()) {
#print STDERR "in fetch\n";
      my ($label, $value)  = $lv->LabelValueToString($mod_type);
      $pw->AddComponentPTM($cx_molid, $cp_seqid, $uniprot_acc, $pos, $residue,
          $mod_type, $value);
    }
  }

  ##
  ## Now get basic mol type and external ids for component molecules
  ## (some of which may, in turn, complexes)
  ##

  @mol_list = keys %sub_molids;

  for($i = 0; $i < @mol_list; $i += ORACLE_LIST_LIMIT) {
 
    if(($i + ORACLE_LIST_LIMIT - 1) < @mol_list) {
      $list = join(",", @mol_list[$i..$i+ORACLE_LIST_LIMIT-1]);
    }
    else {
      $list = join(",", @mol_list[$i..@mol_list-1]);
    }

    $sql = qq!
select unique
  m.mol_id, m.basic_mol_type, em.ext_mol_id, em.id_type
from
  $SCHEMA.pw_ext_mol_id em,
  $SCHEMA.pw_mol m
where
      em.mol_id (+) = m.mol_id
  and m.mol_id in ($list)
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

    while (($mol_id, $basic_mol_type, $ext_mol_id, $id_type) =
        $stm->fetchrow_array()) {
      $pw->AddMol($mol_id, $lv->BasicMolTypeCode($basic_mol_type));
      if ($id_type ne "") {
        $pw->AddMolExId($mol_id, $id_type, $ext_mol_id);
      }
      ##
      ## complexes
      ##
      if ($basic_mol_type eq "CX") {
        $$new_complexes{$mol_id} = 1;
      }
    }
  }

  ##
  ## Finally, get the names of those component molecules
  ##

  for($i = 0; $i < @mol_list; $i += ORACLE_LIST_LIMIT) {
 
    if(($i + ORACLE_LIST_LIMIT - 1) < @mol_list) {
      $list = join(",", @mol_list[$i..$i+ORACLE_LIST_LIMIT-1]);
    }
    else {
      $list = join(",", @mol_list[$i..@mol_list-1]);
    }

    $sql = qq!
select unique
  m.mol_id, m.mol_name, m.name_type
from
  $SCHEMA.pw_mol_name m
where
  m.mol_id in ($list)
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
    while (($mol_id, $mol_name, $name_type) =
        $stm->fetchrow_array()) {
      $pw->AddMolName($mol_id, $name_type, $mol_name);
    }
  }

}


######################################################################
sub GetAtomsAndLabelsForPathway {
  my ($self, $pid) = @_;

  my $db     = $self->{db};
  my $pw     = $self->{pw};
  my $lv     = $self->{lv};
  my $SCHEMA = $self->{schema};

  my ($sql, $stm);
  my ($atom_id, $label_value_id, $atom_source_id, $negative_flag);
  my ($name, $value, $empty_evidence);

  my $atom_list;
  if (defined $self->{evidence_code}) {
  } else {
    $empty_evidence = 1;
  }
  my ($atom_clause);

  my @a_list = split ",", $atom_list;
  my $list;

  $sql = qq!
select
  al.atom_id, al.label_value_id, a.atom_source_id
from
  $SCHEMA.pw_pathway_atom pa,
  $SCHEMA.pw_atom_label al,
  $SCHEMA.pw_atom a
where
      pa.pathway_id = $pid
  and al.atom_id = pa.atom_id
  and al.atom_id = a.atom_id
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
  my %abstractions;
  while (($atom_id, $label_value_id, $atom_source_id) =
      $stm->fetchrow_array()) {
    if (($empty_evidence == 1) || ($self->{atomsinevidence}{$atom_id})) { 
    $pw->AddAtomLabel($atom_id, $label_value_id);
    $pw->AddAtomPathway($atom_id, $pid);
##    ($name, $value) = split /\t/, $lvhash{$label_value_id};
    ($name, $value) = $lv->LabelValueToString($label_value_id);
    if ($name eq "process-type") {
      $pw->AddAtom($atom_id, $label_value_id);
      if ($label_value_id eq $pathway_type ||
          $label_value_id eq $subnet_type) {
        $abstractions{$atom_id} = 1;
      }
      $pw->AddAtomSource($atom_id, $atom_source_id);
    }
    }
  }
  if (keys %abstractions) {
    $self->GetAbstractions(\%abstractions);
  }

  ##
  ## Now conditions
  ##

  $sql = qq!
select
  ac.atom_id, ac.label_value_id, ac.negative_flag
from
  $SCHEMA.pw_pathway_atom pa,
  $SCHEMA.pw_atom_condition ac
where
      pa.pathway_id = $pid
  and ac.atom_id = pa.atom_id
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
  while (($atom_id, $label_value_id, $negative_flag) =
      $stm->fetchrow_array()) {
    if ($negative_flag == 1) {
      $pw->AddAtomNegativeCondition($atom_id, $label_value_id);
    } else {
      $pw->AddAtomCondition($atom_id, $label_value_id);
    }
  }

  ##
  ## Now evidence
  ##

  my ($evidence_code);

  $sql = qq!
 select
  e.atom_id, e.evidence_code
 from
  $SCHEMA.pw_pathway_atom pa,
  $SCHEMA.pw_evidence e
 where
      pa.pathway_id = $pid
  and pa.atom_id = e.atom_id
    !;
#HERE
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
  while (($atom_id, $evidence_code) = $stm->fetchrow_array()) {
    if (($empty_evidence == 1) ||($self->{atomsinevidence}{$atom_id})) {
      $pw->AddAtomEvidence($atom_id, $evidence_code);
    }
  }

  ##
  ## Now references
  ##

  my ($pmid);

  $sql = qq!
 select
  r.atom_id, r.pmid
 from
  $SCHEMA.pw_pathway_atom pa,
  $SCHEMA.pw_references r
 where
      pa.pathway_id = $pid
  and pa.atom_id = r.atom_id
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
  while (($atom_id, $pmid) = $stm->fetchrow_array()) {
    if (($empty_evidence == 1) || ($self->{atomsinevidence}{$atom_id})) {
      $pw->AddAtomReferences($atom_id, $pmid);
    }
  }

  ##
  ## Now notes
  ##

  my ($note);

  $sql = qq!
 select
   n.atom_id, n.note
 from
  $SCHEMA.pw_pathway_atom pa,
  $SCHEMA.pw_notes n
 where
      pa.pathway_id = $pid
  and pa.atom_id = n.atom_id
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
  while (($atom_id, $note) = $stm->fetchrow_array()) {
    if (($empty_evidence == 1) || ($self->{atomsinevidence}{$atom_id})) {
      $pw->AddAtomNotes($atom_id, $note);
    }
  }
}
######################################################################
sub GetCuratorsForPathway {
  my ($self, $pid) = @_;

  my $db     = $self->{db};
  my $pw     = $self->{pw};
  my $lv     = $self->{lv};
  my $SCHEMA = $self->{schema};

  my ($sql, $stm);
  my ($curator, $role);

  $sql = qq!
 select unique
  curator, role
 from
  $SCHEMA.pw_curators
 where
  pathway_id = $pid
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
  while (($curator, $role) = $stm->fetchrow_array()) {
    if ($role eq 'C') {
      $pw->AddCurator($pid, $curator);
    } else {
      $pw->AddReviewer($pid, $curator);
    }
  }
}

######################################################################
1;
######################################################################
