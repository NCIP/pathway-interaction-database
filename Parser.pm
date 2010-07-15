package Parser;
require Exporter;
@ISA = qw(Exporter);
#@EXPORT = qw(
#  CreateIsA
#  CreateLabelSort
#);

use strict;
use PWLabel;
use Pathway;

use constant EOF => "%&%&%&%&";

my %opclass = (
    "complex-def"      => "MOLECULE_DEF",      ## redundant, covered by "molecule-def"
    "molecule-def"     => "MOLECULE_DEF",
    "pathway"          => "PATHWAY",
    "pathway-info"     => "PATHWAY_INFO",
    "protein"          => "MOLECULE",
    "rna"              => "MOLECULE",
    "compound"         => "MOLECULE",
    "complex"          => "MOLECULE",
    "reaction"         => "ATOM",
    "modification"     => "ATOM",
    "binding"          => "ATOM",
    "transcription"    => "ATOM",
    "translocation"    => "ATOM",
#    "cell-process"     => "ATOM",
#    "ubiquitination"   => "ATOM",
#    "replicative-senescence"  => "ATOM",
#    "dna-fragmentation" => "ATOM",
#    "dna-damage"       => "ATOM",
#    "dna-repair"       => "ATOM",
#    "cell-survival"    => "ATOM",
#    "cell-proliferation" => "ATOM",
#    "apoptosis"        => "ATOM",
    "macroprocess"     => "ATOM",
    "cell-process"     => "ATOM",   ## for now, synonym for macroprocess
    "agent"            => "EDGE",
    "input"            => "EDGE",
    "output"           => "EDGE",
    "inhibitor"        => "EDGE",
    "component"        => "COMPONENT",
    "location"         => "LABEL_VALUE",
    "activity-state"   => "LABEL_VALUE",
    "edge-type"        => "LABEL_VALUE",
    "reversible"       => "LABEL_VALUE",
    "process-condition"  => "LABEL_VALUE",
    "process-type"     => "LABEL_VALUE",   ## to specify kind of macroprocess
    "function"         => "LABEL_VALUE",   ## label on inputs/agents
    "EC"               => "NONUNIQUE_VALUE",
    "CA"               => "NONUNIQUE_VALUE",
    "LL"               => "NONUNIQUE_VALUE",
    "GO"               => "NONUNIQUE_VALUE",
    "KG"               => "NONUNIQUE_VALUE",
    "AS"               => "NONUNIQUE_VALUE",
    "atom-id"          => "UNIQUE_VALUE",
    "pathway-id"       => "UNIQUE_VALUE",
    "component-seq-id" => "UNIQUE_VALUE",
    "edge-seq-id"      => "UNIQUE_VALUE",
    "molecule-id"      => "UNIQUE_VALUE",
    "name"             => "UNIQUE_VALUE",
    "organism"         => "UNIQUE_VALUE",
    "source-id"        => "UNIQUE_VALUE",
    "source-pathway-id" => "UNIQUE_VALUE"
);

my @EXT_MOL_ID = (
    "EC",
    "LL",
    "KG",
    "CA",
    "GO"
);

my @MOL_NAME = (
    "AS"
);

my @LABEL = (
    "function",
    "edge-type",
    "reversible",
    "location",
    "activity-state",
    "process-type",
    "process-condition"
);

my %label_name_of_op = (
    "function"          => "function",
    "edge-type"         => "edge-type",
    "reversible"        => "reversible",
    "location"          => "location",
    "activity-state"    => "activity-state",
    "process-type"      => "process-type",
    "process-condition" => "process-type"
);

######################################################################
sub new {
  my ($self, $pw, $lv, $atomidbase, $pathwayidbase) = @_;

  my $x = {};
  $x->{pw} = $pw;
  $x->{lv} = $lv;
  $x->{ROOT} = 0;
  $x->{tokeni} = -1;
  $x->{atomidbase} = $atomidbase;
  $x->{pathwayidbase} = $pathwayidbase;
  SetProcessTypes($lv);
  return bless $x;
}

######################################################################
sub AtomIdMap {
  my ($self) = @_;

  return $self->{atomidmap};
}

######################################################################
sub PathwayIdMap {
  my ($self) = @_;

  return $self->{pathwayidmap};
}

######################################################################
sub SetProcessTypes {
  my ($lv) = @_;

  my ($pair, $label_name, $label_value_name);
##  my $cell_process_lvid =
##      $lv->StringToLabelValue("process-type", "cell-process");
  my $macroprocess_lvid =
      $lv->StringToLabelValue("process-type", "macroprocess");
  for $pair (keys %{ $lv->AllNamePairs }) {
    ($label_name, $label_value_name) = split(/\t/, $pair);
    if ($label_name eq "process-type" &&
        $lv->IsA($lv->StringToLabelValue($label_name, $label_value_name),
        $macroprocess_lvid) ) {
##        $cell_process_lvid) ) {
      if (! defined $opclass{$label_value_name}) {
        $opclass{$label_value_name} = "ATOM";
      }
    }    
  }
}

######################################################################
sub ListErrors {
  my ($self) = @_;

  return $self->{errors}
}

######################################################################
sub OpClass {
  my ($self, $op) = @_;

  if (defined $opclass{$op}) {
    return $opclass{$op};
  } elsif ($op =~ /^\$\d+$/) {
    return "LITERAL";
  } elsif ($op eq "(" || $op eq ")") {
   return $op;
  } else {
    push @{ $self->{errors} }, "#! illegal op: $op, on line " . $self->CurLine() . "\n";
  }
}

######################################################################
sub NextToken {
  my ($self) = @_;

  $self->{tokeni}++;
  if (! defined $self->{tokens}) {
    push @{ $self->{errors} }, "#! empty token list\n";
    return EOF;
  }
  if ($self->{tokeni} >= @{ $self->{tokens} }) {
    return EOF;
  } else {
    return $self->{tokens}[$self->{tokeni}];
  }
}

######################################################################
sub CurToken {
  my ($self) = @_;

  if (! defined $self->{tokens}) {
    push @{ $self->{errors} }, "#! empty token list\n";
    return EOF;
  }
  if ($self->{tokeni} >= @{ $self->{tokens} }) {
    return EOF;
  } else {
    return $self->{tokens}[$self->{tokeni}];
  }
}

######################################################################
sub CurLine {
  my ($self) = @_;

  if (! defined $self->{tokens}) {
    push @{ $self->{errors} }, "#! empty token list\n";
    return "end of file";
  }
  if ($self->{tokeni} >= @{ $self->{tokens} }) {
    return "end of file";
  } else {
    return $self->{lineno}[$self->{tokeni}];
  }
}

######################################################################
sub TokIsStringLit {
  my ($tok) = @_;
  return $tok =~ /^\$/;
}

######################################################################
sub FetchStringLit {
  my ($self, $stridx) = @_;
  $stridx =~ s/^\$//;
  return $self->{str_array}[$stridx];
}

######################################################################
sub PushToken {
  my ($self, $token) = @_;

  my $level = $self->{level};
  my $op = $self->OpClass($token);
  if ($op eq "") {
    push @{ $self->{errors} }, "#! unrecognized opclass $op at " . $self->CurLine . "\n"
  }
  $self->{opstack}[$level] = $op;
  $self->{tokstack}[$level] = $token;
  $self->{stacklineno}[$level] = $self->{lineno}[$self->{tokeni}];
  if ($token eq "pathway") {
    $self->{ROOT} = $level;
  }
}

######################################################################
sub PushLevel {
  my ($self, $token) = @_;
  $self->{level}++;
  $self->PushToken($token);
}

######################################################################
sub CheckLVID {
  my ($self, $lvid, $label, $value) = @_;
  if ($lvid == 0) {
    push @{ $self->{errors} }, "#! can't find label value id for $value of type $label in parse unit " .
        "ending at line " . $self->CurLine() . "\n";
  }
}

######################################################################
sub CheckStack {
  my ($self, @acceptable_ops) = @_;
  for my $op (@acceptable_ops) {
    if ($self->{opstack}[$self->{level}-1] eq $op) {
      return 1;
    }
  }
  push @{ $self->{errors} }, "#! ill-formed parse stack: " .
    "parent = " . $self->{opstack}[$self->{level}-1] . ", " .
    "child = " . $self->{opstack}[$self->{level}] . " " .
    "for parse unit ending at line " . $self->CurLine() . "\n";
  return 0;
}

######################################################################
sub PopLevel {
  my ($self) = @_;

  my $pw = $self->{pw};
  my $lv = $self->{lv};

  my $level  = $self->{level};
  my $local  = \% { $self->{data}{$level} };
  my $parent = \% { $self->{data}{$level-1} };
  my $root   = \% { $self->{data}{$self->{ROOT}} };

  if ($level <= 0) {
    push @{ $self->{errors} }, "#! Attempt to pop empty stack\n";
    return;
  }
  my ($molid, $op, $atomid, $edgeid, $pathwayid, $lvid);
  my ($edgetype);
  my ($moltype, $value);
  my ($compseqid, $compmolid);

  my $curop    = $self->{opstack}[$level];
  my $token    = $self->{tokstack}[$level];
  my $parentop = $self->{opstack}[$level-1];

  my $org = $self->{data}{$self->{ROOT}}{"organism"};

  if ($curop eq "LABEL_VALUE") {

    $$parent{label}{$token} = $$local{value};

  } elsif ($curop eq "UNIQUE_VALUE") {

    if ($token eq "molecule-id") {
      $$local{value} = $self->MappedMolId($$local{value});
    } elsif ($token eq "atom-id") {
      my $original_atom_id = $$local{value};
      $$local{value} = $self->MappedAtomId($$local{value});
      $self->{atomidmap}{$original_atom_id} = $$local{value};
    } elsif ($token eq "pathway-id") {
      my $original_pathway_id = $$local{value};
      $$local{value} = $self->MappedPathwayId($$local{value});
      $self->{pathwayidmap}{$original_pathway_id} = $$local{value};
    }

    $$parent{unique_value}{$token} = $$local{value};

  } elsif ($curop eq "NONUNIQUE_VALUE") {

    push @{ $$parent{nonunique_value}{$token} }, $$local{value};

  } elsif ($curop eq "EDGE") {

    if (! $self->CheckStack("ATOM")) {
      return 0;
    }

    $edgeid = $$local{unique_value}{"edge-seq-id"};
    $molid  = $$local{unique_value}{"molecule-id"};
    push @{ $$parent{edge} }, "$edgeid\t$molid";
    $$parent{edge_type}{$edgeid} = $token;
    for $op (keys %{ $$local{label} }) {
      $value = $$local{label}{$op};
      $$parent{edge_label}{$edgeid}{$op} = $value;
    }
    for $op (keys %{ $$local{mol_label} }) {
      $value = $$local{mol_label}{$op};
      $$parent{mol_label}{$edgeid}{$op} = $value;
    }

  } elsif ($curop eq "ATOM") {

    if (! $self->CheckStack("PATHWAY") ) {
      return 0;
    }

    $atomid = $$local{unique_value}{"atom-id"};
    for (@{ $$local{edge} }) {
      my ($edgeid, $molid) = split /\t/;
      $edgetype = $$local{edge_type}{$edgeid}; 
      $lvid = $lv->StringToLabelValue("edge-type", $edgetype);
      $self->CheckLVID($lvid, $op, $edgetype);
      $pw->AddEdge($atomid, $edgeid, $lvid, $molid);
      $pw->AddEdgeLabel($atomid, $edgeid, $lvid);
    }
    for $edgeid (keys %{ $$local{mol_label} }) {
      for $op (keys %{ $$local{mol_label}{$edgeid} }) {
        $value = $$local{mol_label}{$edgeid}{$op};
        $lvid = $lv->StringToLabelValue($label_name_of_op{$op}, $value);
        $self->CheckLVID($lvid, $op, $value);
        $pw->AddMolLabel($atomid, $edgeid, $lvid);
      }
    }
    for $edgeid (keys %{ $$local{edge_label} }) {
      for $op (keys %{ $$local{edge_label}{$edgeid} }) {
        $value = $$local{edge_label}{$edgeid}{$op};
        $lvid = $lv->StringToLabelValue($label_name_of_op{$op}, $value);
        $self->CheckLVID($lvid, $op, $value);
	$pw->AddEdgeLabel($atomid, $edgeid, $lvid);
      }
    }
    for $op (@LABEL) {
      if (defined $$local{label}{$op}) {
        $value = $$local{label}{$op};
        $lvid = $lv->StringToLabelValue($label_name_of_op{$op}, $value);
        $self->CheckLVID($lvid, $op, $value);
        if ($op eq "process-type") {
          if ($token ne "cell-process" && $token ne "macroprocess") {
            push @{ $self->{errors} }, "#! mis-use of \"process-type\" label " .
            "in parse unit ending at line " . $self->CurLine() . "\n";
          } else {
            $$local{label}{$op} = $lvid;
          }
        } elsif ($op eq "process-condition") {
          $pw->AddAtomCondition($atomid, $lvid);
        } else {
          $pw->AddAtomLabel($atomid, $lvid);
        }
      }
    }
    push @{ $$parent{atom} }, $atomid;


    ##
    ## Don't forget the process-type
    ##
    if (defined $$local{label}{"process-type"}) {
      $lvid = $$local{label}{"process-type"};
    } else {
      $lvid = $lv->StringToLabelValue("process-type", $token);
      $self->CheckLVID($lvid, $op, $token);
    }
    $pw->AddAtom($atomid, $lvid);
  } elsif ($curop eq "COMPONENT") {

    if (! $self->CheckStack("MOLECULE") ) {
      return 0;
    }

    $compseqid = $$local{unique_value}{"component-seq-id"};
    $compmolid = $$local{unique_value}{"molecule-id"};

    $$parent{component}{$compseqid} = $compmolid;

    for $op (keys %{ $$local{mol_label} }) {
      $$parent{submol_label}{$compseqid}{$op} = $$local{mol_label}{$op}
    }
    for $op (keys %{ $$local{label} }) {
      $$parent{comp_label}{$compseqid}{$op} = $$local{label}{$op}
    }

  } elsif ($curop eq "MOLECULE") {

    if (! $self->CheckStack("COMPONENT", "EDGE", "MOLECULE_DEF") ) {
      return 0;
    }

    $molid   = $$local{unique_value}{"molecule-id"};
    $moltype = $lv->StringToLabelValue("molecule-type", $token);

    $pw->AddMol($molid, $moltype, $org);

    $$parent{unique_value}{"molecule-id"} = $molid;

    for $compseqid (keys %{ $$local{component} }) {
      $compmolid = $$local{component}{$compseqid};
      $pw->AddComponent($molid, $compseqid, $compmolid);
    }

    for $compseqid (keys %{ $$local{comp_label} }) {
      for $op (keys %{ $$local{comp_label}{$compseqid} }) {
        $value = $$local{comp_label}{$compseqid}{$op};
        $lvid = $lv->StringToLabelValue($op, $value);
        $self->CheckLVID($lvid, $op, $value);
        $pw->AddComponentLabel($molid, $compseqid, $lvid);
      }
    }

    for $compseqid (keys %{ $$local{submol_label} }) {
      for $op (keys %{ $$local{submol_label}{$compseqid} }) {
        $value = $$local{submol_label}{$compseqid}{$op};
        $lvid = $lv->StringToLabelValue($op, $value);
        $self->CheckLVID($lvid, $op, $value);
        $pw->AddComponentLabel($molid, $compseqid, $lvid);
      }
    }

    for $op (@LABEL) {
      if (defined $$local{label}{$op}) {
        $$parent{mol_label}{$op} = $$local{label}{$op};
      }
    }

    for $op (@EXT_MOL_ID) {
      if (defined $$local{nonunique_value}{$op}) {
        for (@{ $$local{nonunique_value}{$op} }) {
          $pw->AddMolExId($molid, $op, $_);
        }
      }
    }

    for my $op (@MOL_NAME) {
      if (defined $$local{nonunique_value}{$op}) {
        for (@{ $$local{nonunique_value}{$op} }) {
          $pw->AddMolName($molid, $op, $_);
        }
      }
    }

  } elsif ($curop eq "MOLECULE_DEF") {


  } elsif ($curop eq "PATHWAY") {

    $pathwayid = $$local{unique_value}{"pathway-id"};

    ##
    ## If there was a PATHWAY_COMPONENT, this was already done
    ##
    $pw->AddPathwayId($pathwayid);
    $pw->SetPathwayName($pathwayid, $$local{unique_value}{"name"});
    $pw->SetPathwayExId($pathwayid,
        $$local{unique_value}{"source-pathway-id"});
    $pw->SetPathwayOrg($pathwayid, $$local{unique_value}{"organism"});
##    $pw->SetPathwaySrcId($pathwayid, "?");

    for (@{ $$local{atom} }) {
      $pw->AddAtomPathway($_, $pathwayid);
    }

  } elsif ($curop eq "PATHWAY_INFO") {

    if (! $self->CheckStack("PATHWAY") ) {
      return 0;
    }

    $pathwayid = $$local{unique_value}{"pathway-id"};
    $$parent{unique_value}{"pathway-id"} = $pathwayid;
    $$parent{unique_value}{"name"} = $$local{unique_value}{"name"};
    $$parent{unique_value}{"source-pathway-id"} =
        $$local{unique_value}{"source-pathway-id"};
    $$parent{unique_value}{"organism"} = $$local{unique_value}{"organism"};

    $pw->AddPathwayId($pathwayid);
    $pw->SetPathwayName($pathwayid, $$local{unique_value}{"name"});
    $pw->SetPathwayExId($pathwayid,
        $$local{unique_value}{"source-pathway-id"});
    $pw->SetPathwayOrg($pathwayid, $$local{unique_value}{"organism"});
##    $pw->SetPathwaySrcId($pathwayid, "?");

  }

  undef %{ $local };  
  $self->{level}--;
  return 1;
}


######################################################################
sub DoStr {
  my ($self, $line) = @_;

  my ($a, $b, $c);
  while ($line =~ /(.*)(\")([^"]*)(")(.*)/) {
    ($a,$b,$c) = ($1, $3, $5);
    $b =~ s/^\s+//;
    $b =~ s/\s+$//;
    if (($b ne "") && (not defined $self->{str}{$b})) {
      $self->{next_str}++;
      $self->{str}{$b} = $self->{next_str};
      $self->{str_array}[$self->{next_str}] = $b;
    }
    $line = $a . " " . "\$" . $self->{str}{$b} . " " . $c;
  }
  return $line;
}

######################################################################
sub Parse {
  my ($self) = @_;

  my $token = $self->NextToken();
  my $opclass;
  while ($token ne EOF) {
    $opclass = $self->OpClass($token);
    if ($opclass eq "LITERAL") {
      $self->{data}{$self->{level}}{value} = $self->FetchStringLit($token);
    } elsif ($opclass eq "(") {
      $token = $self->NextToken();
      $opclass = $self->OpClass($token);
      $self->PushLevel($token);
    } elsif ($opclass eq ")") {
      if (! $self->PopLevel()) {
        return 0;
      }
    }
    $token = $self->NextToken();
  }
  if ($self->{level} != 0) {
    push @{ $self->{errors} }, "#! Mismatched parentheses; ending level = " .
        $self->{level} . "\n"
  }
  if (defined $self->{errors} && @{ $self->{errors} } > 0) {
    return 0;
  } else {
    return 1;
  }
}

######################################################################
sub CheckTokenStream {
  my ($self) = @_;

  my ($i, $prev, $succ, $tok, $err);
  my ($strprev, $strsucc, $strtok);

  my $lineno = $self->{lineno};
  if (! defined $self->{tokens}) {
    push @{ $self->{errors} }, "#! empty token list\n";
    return 0;
  }
  for ($i = 1; $i < @{ $self->{tokens} } - 1; $i++) {
    $prev = $self->{tokens}[$i-1];
    $tok =  $self->{tokens}[$i];
    $succ = $self->{tokens}[$i+1];
    if (TokIsStringLit($tok)) {
      $strtok = $self->FetchStringLit($tok);
      ##
      ## following condition to allow names with spaces to act as
      ## ATOM operators
      ##
      if ($prev eq "(" && $self->OpClass($strtok) eq "ATOM") {
        $tok                = $strtok;
        $self->{tokens}[$i] = $strtok;
      }
    } else {
      $strtok = $tok;
    }
    if (TokIsStringLit($prev)) {
      $strprev = $self->FetchStringLit($prev);
    } else {
      $strprev = $prev;
    }
    if (TokIsStringLit($succ)) {
      $strsucc = $self->FetchStringLit($succ);
    } else {
      $strsucc = $succ;
    }

    ##
    ##  ' ) X ) ', where X is not '(' or ')'
    ##
    if ($prev eq ")" && $succ eq ")" && $tok ne ")") {
      push @{ $self->{errors} }, "#! syntax error ')' at line $$lineno[$i]\n";
      $err++;
      next;
    }
    ##
    ##  ' X Y Z ', where X, Y, Z are not '(' or ')'
    ##
    if ($prev ne "(" && $prev ne ")" &&
        $succ ne "(" && $succ ne ")" &&
        $tok ne "(" && $tok ne ")" ) {
      push @{ $self->{errors} }, "#! syntax error '$strsucc' at line $$lineno[$i+1]\n";
      $err++;
      next;
    }
    ##
    ##  ' ( X ) '
    ##
    if ($prev eq "(" && $succ eq ")") {
      push @{ $self->{errors} }, "#! syntax error '$strsucc' at line $$lineno[$i+1]\n";
      $err++;
      next;
    }
    
    ##
    ##  ' ) X ( ', where X is not ')'
    ##
    if ($prev eq ")" && $succ eq "(" && $tok ne ")") {
      push @{ $self->{errors} }, "#! syntax error '$strtok' at line $$lineno[$i]\n";
      $err++;
      next;
    }
    
    ##
    ##  ' ( X ', where X is a string literal
    ##
    if ($prev eq "(" && TokIsStringLit($tok)) {
      push @{ $self->{errors} }, "#! syntax error '$strtok' at line $$lineno[$i]\n";
      $err++;
      next;
    }
    
    ##
    ##  ' X ) ', where X is not a string literal and not ')'
    ##
    if ($succ eq ")" && !TokIsStringLit($tok) && $tok ne ")") {
      push @{ $self->{errors} }, "#! syntax error '$strtok' at line $$lineno[$i]\n";
      $err++;
      next;
    }
    
  }
  return $err;
}

######################################################################
sub ReadFile {
  my ($self, $fh) = @_;

#  open(INF, $f) or die "Cannot open input file $f";

  my ($line, $lineno, $token);
  while (<$fh>) {
    $lineno++;
    chop;
    $line = $_;
    $line =~ s/;;.*//;
    $line = $self->DoStr($line);
    $line =~ s/\(/ ( /g;
    $line =~ s/\)/ ) /g;
    $line =~ s/\s+/ /g;
    $line =~ s/^\s+//;
    $line =~ s/\s+$//;
    if ($line =~ /^\s+$/) {
      next;
    }
    for $token (split(/ /, $line)) {
      push @{ $self->{tokens} }, $token;
      push @{ $self->{lineno} }, $lineno;
    }
  }
#  close INF;
  close $fh;
  if ($self->CheckTokenStream()) {
    push @{ $self->{errors} },  "#! quitting because of syntax errors\n";
    return 0;
  }
  return 1;
}

######################################################################
sub ReadData {
  my ($self, $data) = @_;

  my ($line, $lineno, $token);
  for (split("\n", $data)) {
    $lineno++;
    $line = $_;
    $line =~ s/;;.*//;
    $line = $self->DoStr($line);
    $line =~ s/\(/ ( /g;
    $line =~ s/\)/ ) /g;
    $line =~ s/\s+/ /g;
    $line =~ s/^\s+//;
    $line =~ s/\s+$//;
    if ($line =~ /^\s+$/) {
      next;
    }
    for $token (split(/ /, $line)) {
      push @{ $self->{tokens} }, $token;
      push @{ $self->{lineno} }, $lineno;
    }
  }
  close INF;
  if ($self->CheckTokenStream()) {
    push @{ $self->{errors} },  "#! quitting because of syntax errors\n";
    return 0;
  }
  return 1;
}

######################################################################
sub MappedMolId {
  my ($self, $molid) = @_;

  if (defined $self->{molmap}) {
    if (defined $self->{molmap}{$molid}) {
      return $self->{molmap}{$molid};
    } else {
      push @{ $self->{errors} },  "#! no map for molecule-id $molid\n";
      return 0;
    }
  } else {
    return $molid;
  }
}

######################################################################
sub MappedAtomId {
  my ($self, $atomid) = @_;

  if (defined $self->{atomidbase}) {
    my $x = $self->{atomidbase};
    $self->{atomidbase}++;
    return $x;
#    return $atomid + $self->{atomidbase};
  } else {
    return $atomid;
  }
}

######################################################################
sub MappedPathwayId {
  my ($self, $pathwayid) = @_;

##  if (defined $self->{pathwayidbase}) {
##    return $pathwayid + $self->{pathwayidbase};
## now it's a simple set pathway id, if non zero/empty string
  if ($self->{pathwayidbase}) {
    return $self->{pathwayidbase};
  } else {
    return $pathwayid;
  }
}

######################################################################
sub ReadMolMap {
  my ($self, $f) = @_;

  if ($f) {
    my ($local_mol_id, $global_mol_id);
    open(MAPF, $f) or die "Cannot open molecule map file $f";
    $self->{molmap} = {};
    while (<MAPF>) {
      chop;
      ($local_mol_id, $global_mol_id) = split /\t/;
      $self->{molmap}{$local_mol_id} = $global_mol_id;
    }
    close MAPF;
  }
}

######################################################################
1;
######################################################################
