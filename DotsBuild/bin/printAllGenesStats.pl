#!/usr/local/bin/perl

##Note that am presenting as the number of "genes" (for now represent as  clusters
## of RNAs)  Simplest to create a tmp table to hold the mapping btwn assemblies
## and genes to facilitate queries...

## Brian Brunk 11/28/00

use strict;
use Getopt::Long;
use GUS::ObjRelP::DbiDatabase;
use GUS::Common::GusConfig;

$| = 1;

my ($gusConfigFile,$verbose,$go_cvs_version,$createTmpTableOnly,$useExistingTmpTable);
my $dropTmpTable = 0;		#3drop the tmp table until learn how to keep it..

&GetOptions('gusConfigFile=s' => \$gusConfigFile,
	    'dropTmpTable!' => \$dropTmpTable,
            'verbose!' => \$verbose, 
	    'go_cvs_version=s' => \$go_cvs_version,
	    "createTmpTableOnly!" => \$createTmpTableOnly,
            'useExistingTmpTable' => \$useExistingTmpTable,
           );

#put all global vars here...
my($humAssSeqs,$musAssSeqs,$totalHumAssemblies,$totalMusAssemblies);
my $humNonSingletons = 0;
my $musNonSingletons = 0;
my $humDoTSGenes = 0;
my $musDoTSGenes = 0;
my $humGenScanPreds = 0;
my $musGenScanPreds = 0;
my $humGrailPreds = 0;
my $musGrailPreds = 0;
my $humGenesWithNR = 0;
my $musGenesWithNR = 0;
my $humAssWithNR = 0;
my $musAssWithNR = 0;
my $humDoTSWithNR = 0;
my $musDoTSWithNR = 0;
my $humGrailWithNR = 0;
my $musGrailWithNR = 0;
my $humGenScanWithNR = 0;
my $musGenScanWithNR = 0;
my $humGenesWithCellRoles = 0;
my $musGenesWithCellRoles = 0;
my $humDoTSWithCellRoles = 0;
my $musDoTSWithCellRoles = 0;
my $musAssWithCellRoles = 0;
my $humAssWithCellRoles = 0;
my $humGrailWithCellRoles = 0;
my $musGrailWithCellRoles = 0;
my $humGenScanWithCellRoles = 0;
my $musGenScanWithCellRoles = 0;
my $humGenesWithGOFunction = 0;
my $musGenesWithGOFunction = 0;
my $humDoTSWithGOFunction = 0;
my $musDoTSWithGOFunction = 0;
my $musAssWithGOFunction = 0;
my $humAssWithGOFunction = 0;
my $humGrailWithGOFunction = 0;
my $musGrailWithGOFunction = 0;
my $humGenScanWithGOFunction = 0;
my $musGenScanWithGOFunction = 0;
my %humDoTSRoles;
my %musDoTSRoles;
my $totalDoTShumRoles = 0;
my $totalDoTSmusRoles = 0;
my %humGrailRoles;
my %musGrailRoles;
my %humGenScanRoles;
my %musGenScanRoles;
my $totalGrailhumRoles = 0;
my $totalGenScanhumRoles = 0;
my $totalGrailmusRoles = 0;
my $totalGenScanmusRoles = 0;
my %humDoTSGOFun;
my %musDoTSGOFun;
my $totalDoTShumGOFun = 0;
my $totalDoTSmusGOFun = 0;
my %humGrailGOFun;
my %musGrailGOFun;
my %humGenScanGOFun;
my %musGenScanGOFun;
my $totalGrailhumGOFun = 0;
my $totalGenScanhumGOFun = 0;
my $totalGrailmusGOFun = 0;
my $totalGenScanmusGOFun = 0;
my $humConsAssemblies = 0;
my $musConsAssemblies = 0;
my $humConsDoTS = 0;
my $musConsDoTS = 0;
my $humDoTSConsGenes = 0;
my $musDoTSConsGenes = 0;
my %cellRoleIdMap;
my %goFunctionIdMap;

print STDERR "Establishing dbi login\n" if $verbose;

my $gusconfig = GUS::Common::GusConfig->new($gusConfigFile);

$dbh = DBI->connect($gusconfig->getDbiDsn(),
		    $gusconfig->getReadOnlyDatabaseLogin(),
		    $gusconfig->getReadOnlyDatabasePassword());


my $dbh = $db->getQueryHandle();
my $stmt;

# ##first create the tmp table that will be used to determine "genes"
print STDERR "creating tmp table\n";
##first check to see if it already exists....will drop it at end...
$stmt = $dbh->prepare("select table_name from all_tables where table_name = 'PRINTSTATSTMP'");
$stmt->execute();
my $tmpTableExists = 0;
while (my($tn) = $stmt->fetchrow_array()) {
  $tmpTableExists = 1;
}

my $rows = 0;

if (! $useExistingTmpTable) {
  ##should drop the table first if it exists...
  print STDERR "\nDropping PrintStatsTmp: ",$dbh->do("drop table PrintStatsTmp"),"\n" if $tmpTableExists;
  print STDERR "\nCreating PrintStatsTmp\n";

  my $sql = 
"create table PrintStatsTmp as
select r.gene_id,a.na_sequence_id,a.number_of_contained_sequences,
       1 as total_seqs,a.taxon_id,ts.aa_sequence_id
from Dots.RNA r, Dots.RNASequence rs, Dots.NAFeature f,
     Dots.Assembly a, Dots.TranslatedAAFeature tf,
     Dots.TranslatedAASequence ts
where r.rna_id = rs.rna_id
and rs.na_feature_id = f.na_feature_id
and a.na_sequence_id = f.na_sequence_id
and a.taxon_id in (8,14)
and tf.na_feature_id = f.na_feature_id
and ts.aa_sequence_id = tf.aa_sequence_id ";

  $rows = $dbh->do($sql);
  $dbh->commit();
  print STDERR "tmp table created..entered $rows rows\n";

  $sql =
"update PrintStatsTmp set total_seqs = 2
where gene_id in (
     select gene_id
     from PrintStatsTmp
     group by gene_id
     having sum(number_of_contained_sequences) > 1
)";

  print STDERR "  Upating total_seqs: ".$dbh->do($sql)." rows\n";

  $dbh->commit();

  $sql = 
"update printstatstmp set total_seqs = 2
where na_sequence_id in (
      select na_sequence_id
      from dots.assembly
      where contains_mrna = 1 )
and total_seqs = 1";

  ##want to updat total_seqs to two where sequence is singleton but mRNA
  print STDERR "  Upating total_seqs: ".$dbh->do($sql);
  $dbh->commit();
  print STDERR "  Updating number_of_contained_sequences: ".$dbh->do("update PrintStatsTmp set number_of_contained_sequences = 2 where number_of_contained_sequences > 1")." rows\n";
  $dbh->commit();
  ##create indexes to facillitate queries
  $dbh->do("create index tran_id_ind on printstatstmp (gene_id)");
  $dbh->do("create index na_id_ind on printstatstmp (na_sequence_id)");
  $dbh->do("create index aa_id_ind on printstatstmp (aa_sequence_id)");
  $dbh->commit();
  
} else {
  print STDERR "using existing PrintStatsTmp table and data\n";
}

if ( $createTmpTableOnly) {
  print STDERR "Creating TmpTable only....exiting\n";
  exit;
} 

##first get the number of Assembled Sequences
$stmt = $dbh->prepare("select a.taxon_id,sum(a.number_of_contained_sequences),count(*) from Dots.Assembly a group by a.taxon_id");
print STDERR "Retrieving number of AssemblySequences\n" ;
$stmt->execute();

while (my($taxon_id,$num,$count) = $stmt->fetchrow_array()) {
  if ($taxon_id == 8) {
    $humAssSeqs = $num;
    $totalHumAssemblies = $count;
  } elsif ($taxon_id == 14) {
    $musAssSeqs = $num;
    $totalMusAssemblies = $count;
  } else {
    print STDERR "ERROR: unknown taxon_id '$taxon_id' has $num AssemblySequences\n";
  }
}
$stmt = $dbh->prepareAndExecute("select taxon_id,count(*) from Dots.Assembly where number_of_contained_sequences > 1 group by taxon_id");
while (my($taxon_id,$nonSing) = $stmt->fetchrow_array()) {
  if ($taxon_id == 8) {
    $humNonSingletons = $nonSing;
  } elsif ($taxon_id == 14) {
    $musNonSingletons = $nonSing;
  } else {
    print STDERR "ERROR: unknown taxon_id '$taxon_id' has $nonSing nonSingletons\n";
  }
}
print STDERR "  Human: $humAssSeqs, $totalHumAssemblies, $humNonSingletons\n  Mouse: $musAssSeqs, $totalMusAssemblies, $musNonSingletons\n\n" ;

##now the number of DoTS genes....
$stmt = $dbh->prepare("select taxon_id,count(distinct gene_id) from PrintStatsTmp where total_seqs > 1 group by taxon_id");
print STDERR "Retrieving number of DoTS genes\n" ;
$stmt->execute();

while (my($taxon_id,$num) = $stmt->fetchrow_array()) {
  if ($taxon_id == 8) {
    $humDoTSGenes = $num;
  } elsif ($taxon_id == 14) {
    $musDoTSGenes = $num;
  } else {
    print STDERR "ERROR: unknown taxon_id '$taxon_id' has $num AssemblySequences\n";
  }
}
print STDERR "  Human: $humDoTSGenes\n  Mouse: $musDoTSGenes\n\n" ;

##now want to get the number of consistent DoTS alignments on genomic sequence...only for human...
$stmt = $dbh->prepare("select taxon_id,count(distinct gene_id)
 from PrintStatsTmp t, Dots.ConsistentAlignment c
 where t.total_seqs > 1
 and c.transcript_na_sequence_id = t.na_sequence_id
 and c.is_consistent = 1
 group by taxon_id");
print STDERR "Retrieving number of Consistent DoTS genes\n";
$stmt->execute();

while (my($taxon_id,$num) = $stmt->fetchrow_array()) {
  if ($taxon_id == 8) {
    $humDoTSConsGenes = $num;
  } elsif ($taxon_id == 14) {
    $musDoTSConsGenes = $num;
  } else {
    #		print STDERR "ERROR: unknown taxon_id '$taxon_id' has $num AssemblySequences\n";
  }
}
print STDERR "  Human: $humDoTSConsGenes\n  Mouse: $musDoTSConsGenes\n\n" ;

##now the number of consistent DoTS and Assemblies...
$stmt = $dbh->prepare("select taxon_id,number_of_contained_sequences,count(*)
 from PrintStatsTmp t, Dots.ConsistentAlignment c
 where c.transcript_na_sequence_id = t.na_sequence_id
 and c.is_consistent = 1
 group by taxon_id,number_of_contained_sequences");
print STDERR "Retrieving number of Consistent DoTS sequences\n" ;
$stmt->execute();

while (my($taxon_id,$num_seqs,$count) = $stmt->fetchrow_array()) {
  if ($taxon_id == 8) {
    if ($num_seqs == 2) {
      $humConsAssemblies = $count;
    }
    $humConsDoTS += $count;	##total number...
  } elsif ($taxon_id == 14) {
    if ($num_seqs == 2) {
      $musConsAssemblies = $count;
    }
    $musConsDoTS += $count;	##total number...
  } else {
    print STDERR "ERROR: unknown taxon_id '$taxon_id' has $count consistent\n";
  }
}
print STDERR "  Human: $humConsDoTS, $humConsAssemblies \n  Mouse: $musConsDoTS, $musConsAssemblies\n\n" ;

## Now determine the "known" Genes ...
print STDERR "Determining number of DoTS genes with NR neigbors:\n" ;
$stmt = $dbh->prepare("select pst.taxon_id,count(distinct pst.gene_id) from PrintStatsTmp pst, Dots.Similarity s where s.query_table_id = 56 and s.query_id = pst.na_sequence_id and s.subject_table_id = 83 and pst.total_seqs > 1 group by pst.taxon_id");
$stmt->execute();

while (my($taxon_id,$totGenes) = $stmt->fetchrow_array()) {
  if ($taxon_id == 8) {
    $humGenesWithNR = $totGenes;
  } elsif ($taxon_id == 14) {
    $musGenesWithNR = $totGenes;
  } else {
    print STDERR "ERROR: unknown taxon_id '$taxon_id' has $totGenes Assembly neighbors\n";
  }
}
print STDERR "  Human: $humGenesWithNR\n  Mouse: $musGenesWithNR\n\n" ;

##now assemblies ...
print STDERR "Determining number of DoTS assemblies with NR neigbors:\n" ;
$stmt = $dbh->prepare("select a.taxon_id,a.number_of_contained_sequences,count(distinct a.na_sequence_id) from PrintStatsTmp a, Dots.Similarity s where s.query_table_id = 56 and s.query_id = a.na_sequence_id and s.subject_table_id = 83 group by a.taxon_id,a.number_of_contained_sequences");
$stmt->execute();

while (my($taxon_id,$num_seqs,$num) = $stmt->fetchrow_array()) {
  if ($taxon_id == 8) {
    $humAssWithNR = $num unless $num_seqs == 1;
    $humDoTSWithNR += $num;
  } elsif ($taxon_id == 14) {
    $musAssWithNR = $num unless $num_seqs == 1;
    $musDoTSWithNR += $num;
  } else {
    print STDERR "ERROR: unknown taxon_id '$taxon_id' has $num Assembly neighbors\n";
  }
}
print STDERR "  Human: Assem: $humAssWithNR, DoTS: $humDoTSWithNR\n  Mouse: Assem: $musAssWithNR, DoTS: $musDoTSWithNR\n\n" ;


print STDERR "Determining number of DoTS genes with GoFunctions:\n" ;
$stmt = $dbh->prepare("select pst.taxon_id,count(distinct pst.gene_id) from PrintStatsTmp pst, Dots.AASequenceGOFunction scr, Dots.GOFunction f where pst.total_seqs > 1 and scr.aa_sequence_id = pst.aa_sequence_id and f.go_function_id = scr.go_function_id and f.go_cvs_version = '$go_cvs_version' group by pst.taxon_id");
$stmt->execute();

while (my($taxon_id,$totGenes) = $stmt->fetchrow_array()) {
  if ($taxon_id == 8) {
    $humGenesWithGOFunction = $totGenes;
  } elsif ($taxon_id == 14) {
    $musGenesWithGOFunction = $totGenes;
  } else {
    print STDERR "ERROR: unknown taxon_id '$taxon_id' has $totGenes GOFunction\n";
  }
}
print STDERR "  Human: $humGenesWithGOFunction Mouse: $musGenesWithGOFunction\n\n" ;

#DoTS assemblies and singletons
print STDERR "Determining number of DoTS Assemblies with GoFunctions:\n" ;
$stmt = $dbh->prepare("select pst.taxon_id,pst.number_of_contained_sequences,count(distinct pst.na_sequence_id) from PrintStatsTmp pst, Dots.AASequenceGOFunction scr, Dots.GOFunction f where scr.aa_sequence_id = pst.aa_sequence_id and f.go_function_id = scr.go_function_id and f.go_cvs_version = '$go_cvs_version' group by pst.taxon_id,pst.number_of_contained_sequences");

$stmt->execute();

while (my($taxon_id,$num_seqs,$num) = $stmt->fetchrow_array()) {
  if ($taxon_id == 8) {
    $humAssWithGOFunction = $num unless $num_seqs == 1;
    $humDoTSWithGOFunction += $num;
  } elsif ($taxon_id == 14) {
    $musAssWithGOFunction = $num unless $num_seqs == 1;
    $musDoTSWithGOFunction += $num;
  } else {
    print STDERR "ERROR: unknown taxon_id '$taxon_id' has $num GOFunction\n";
  }
}
print STDERR "  Human: Assem: $humAssWithGOFunction, DoTS: $humDoTSWithGOFunction, Mouse: Assem: $musAssWithGOFunction, DoTS: $musDoTSWithGOFunction\n\n" ;

# ##now lets generate the GOFunction dist...
##first DoTS Genes..
print STDERR "Determining breakdown of GOFunctions DoTS Genes:\n";
$stmt = $dbh->prepare("select pst.taxon_id,f.name,count(distinct pst.gene_id)
      from PrintStatsTmp pst, Dots.AASequenceGOFunction scr, Dots.GOFunction f
      where pst.total_seqs > 1
      and scr.aa_sequence_id = pst.aa_sequence_id
      and f.go_function_id = scr.go_function_id
      and f.minimum_level = 1
  and f.go_cvs_version = '$go_cvs_version'
      group by pst.taxon_id,f.name");
$stmt->execute();

while (my($taxon_id,$name,$count) = $stmt->fetchrow_array()) {
  next if $name eq 'obsolete' || $name =~ /unknown/;
  if ($taxon_id == 8) {
    $humDoTSGOFun{$name} = $count;
    $totalDoTShumGOFun += $count;
  } elsif ($taxon_id == 14) {
    $musDoTSGOFun{$name} = $count;
    $totalDoTSmusGOFun += $count;
  } else {
    print STDERR "unknown taxon_id $taxon_id for $name, $count\n";
  }
}
foreach my $a (keys%humDoTSGOFun) {
  print STDERR "  $a: human=$humDoTSGOFun{$a}, mus=$musDoTSGOFun{$a}\n";
}

##get the current date...
#Wed Nov 29 10:00:01 EST 2000
my @date = split(/\s+/,`date`);
my $date = "$date[1] $date[2], $date[-1]";

##now need to get the hashes mapping cell roles and gofunctions to ids....
# $stmt = $dbh->prepare("select name,cell_role_id from CellRole where source = 'TIGR'");
# $stmt->execute();
# while (my($name,$id) = $stmt->fetchrow_array()){
# 	$cellRoleIdMap{$name} = $id;
# }

my $goQuery = "select name,go_function_id from Dots.GOFunction where minimum_level = 1 and go_cvs_version = '$go_cvs_version'";

print STDERR "\n$goQuery\n\n" if $verbose;

$stmt = $dbh->prepare($goQuery);
$stmt->execute();
while (my($name,$id) = $stmt->fetchrow_array()) {
  $goFunctionIdMap{$name} = $id;
}

##print out the statistics...
&printStats();

if ($dropTmpTable) {
  $dbh->do("drop table PrintStatsTmp"); 
  $dbh->commit();
}

$stmt->finish();
$dbh->disconnect();

sub printStats {

  print <<endOfStats;
<HTML>
<TABLE BGCOLOR="#ffffff" WIDTH="100%">
<TR><TD>
<FONT SIZE="+1"><I><B>1. Total input sequences</B></I></FONT><IMG SRC="/images/dblue1.gif" WIDTH="150" HEIGHT="2"><BR><BR>
</TD></TR>
<TR><TD>
<TABLE CELLPADDING="1" CELLSPACING="0" BORDER="0" BGCOLOR="#000000" ALIGN="left">

<TR><TD>
<TABLE CELLPADDING="3" CELLSPACING="1" BORDER="0" BGCOLOR="#ffffff" WIDTH="100%">
<TR>
  <TH ROWSPAN="2">&nbsp;</TH>
  <TH BGCOLOR="#cc0000">
     <FONT FACE="helvetica,sans-serif" SIZE="+1" COLOR="#ffffff">DOTS</FONT></TH>
</TR>
<TR>
  <TH><FONT FACE="helvetica,sans-serif" SIZE="+1">(ESTs and mRNAs)</FONT></TH>
</TR>
<TR><TD ALIGN="center">Human</TD><TD ALIGN="center">$humAssSeqs</TD></TR>
<TR><TD ALIGN="center">Mouse</TD><TD ALIGN="center">$musAssSeqs</TD></TR>
</TABLE>

</TD></TR></TABLE>
</TD></TR>
<TR><TD>

<BR CLEAR="both"><BR><BR>

<FONT SIZE="+1"><I><B>2. Predicted and "known" genes</B></I></FONT><IMG SRC="/images/dblue1.gif" WIDTH="200" HEIGHT="2"><BR><BR>
</TD></TR>
<TR><TD>
<TABLE CELLPADDING="1" CELLSPACING="0" BORDER="0" BGCOLOR="#000000" ALIGN="left">

<TR><TD>
<TABLE CELLPADDING="3" CELLSPACING="1" BORDER="0" BGCOLOR="#ffffff">
<TR>
  <TH ROWSPAN="2">&nbsp;</TH>
  <TH BGCOLOR="#cc0000" COLSPAN="3">
    <FONT FACE="helvetica,sans-serif" SIZE="+1" COLOR="#ffffff">DOTS</FONT></TH>
</TR>
<TR>
  <TH><FONT FACE="helvetica,sans-serif" SIZE="+1">Total</FONT></TH>
  <TH><FONT FACE="helvetica,sans-serif" SIZE="+1">Assemblies</FONT></TH>
  <TH><FONT FACE="helvetica,sans-serif" SIZE="+1">"Genes"</FONT></TH>

<TR>
<TD ALIGN="left">Human</TD><TD ALIGN="center">$totalHumAssemblies</TD><TD ALIGN="center">$humNonSingletons</TD><TD ALIGN="center">$humDoTSGenes</TD>
</TR>

<TR>
<TD ALIGN="left">&nbsp;&nbsp;  Consistent alignment</TD><TD ALIGN="center">$humConsDoTS</TD><TD ALIGN="center">$humConsAssemblies</TD><TD ALIGN="center">$humDoTSConsGenes</TD>
</TR>

<TR>
<TD ALIGN="left">&nbsp;&nbsp;  "Known"</TD><TD ALIGN="center">$humDoTSWithNR</TD><TD ALIGN="center">$humAssWithNR</TD><TD ALIGN="center">$humGenesWithNR</TD>
</TR>

<TR>
<TD ALIGN="left">&nbsp;&nbsp;  GO function assigned</TD><TD ALIGN="center">$humDoTSWithGOFunction</TD><TD ALIGN="center">$humAssWithGOFunction</TD><TD ALIGN="center">$humGenesWithGOFunction</TD>
</TR>

<TR>
<TD ALIGN="left">Mouse</TD><TD ALIGN="center">$totalMusAssemblies</TD><TD ALIGN="center">$musNonSingletons</TD><TD ALIGN="center">$musDoTSGenes</TD>
</TR>
endOfStats

  if ($musConsDoTS) {		##have consistent DoTS mus alignments..
    print <<endMusCons
<TR>
<TD ALIGN="left">&nbsp;&nbsp;  Consistent alignment</TD><TD ALIGN="center">$musConsDoTS</TD><TD ALIGN="center">$musConsAssemblies</TD><TD ALIGN="center">$musDoTSConsGenes</TD>
</TR>
endMusCons
  }

  print <<endOfStats;

<TR>
<TD ALIGN="left">&nbsp;&nbsp;  "Known"</TD><TD ALIGN="center">$musDoTSWithNR</TD><TD ALIGN="center">$musAssWithNR</TD><TD ALIGN="center">$musGenesWithNR</TD>
</TR>

<TR>
<TD ALIGN="left">&nbsp;&nbsp;  GO function assigned</TD><TD ALIGN="center">$musDoTSWithGOFunction</TD><TD ALIGN="center">$musAssWithGOFunction</TD><TD ALIGN="center">$musGenesWithGOFunction</TD>
</TR>

</TD></TR></TABLE>
</TABLE>

</TD></TR>
<TR><TD>
<BR CLEAR="both"><BR>

<DIV ALIGN="left">
<B>DOTS Total:</B> all EST assemblies, including singletons<BR>
<B>DOTS Assemblies:</B> EST assemblies with >1 input sequence (i.e., excluding singletons)<BR>
<B>DOTS Genes:</B> Clusters (by similarity) of DoTS assemblies with >1 input sequence (i.e., excluding singletons)<BR>
<B>Consistent alignment:</B> DoTS sequences that align (SIM4) to UCSC Golden Path sequences in a manner consistent with being a processed transcript<BR>
<B>"Known" genes:</B> homology to a GenBank NR protein with BLAST p-value < 1E-5<BR>
<BR><BR>
</TD></TR>
<TR><TD>

endOfStats


  ##now make a table for the GOFunctions...
  print <<gofuns;
<FONT SIZE="+1"><I><B>3. GO function distributions</B></I></FONT><IMG SRC="/images/dblue1.gif" WIDTH="350" HEIGHT="2"><BR><BR>
</TD></TR>
<TR><TD>
<TABLE CELLPADDING="1" CELLSPACING="0" BORDER="0" BGCOLOR="#000000" ALIGN="left">

<TR><TD>
<TABLE CELLPADDING="3" CELLSPACING="1" BORDER="0" BGCOLOR="#ffffff">
<TR>
  <TH ROWSPAN="3">&nbsp;</TH>
  <TH BGCOLOR="#cc0000" COLSPAN="4">
    <FONT FACE="helvetica,sans-serif" SIZE="+1" COLOR="#ffffff">DOTS "Genes"</FONT></TH>
</TR>
<TR>
  <TH COLSPAN="2"><FONT FACE="helvetica,sans-serif" SIZE="+1">Human</FONT></TH>
  <TH COLSPAN="2"><FONT FACE="helvetica,sans-serif" SIZE="+1">Mouse</FONT></TH>
</TR>
<TR>
<TD ALIGN="center">Total</TD><TD ALIGN="center">Percentage</TD>
<TD ALIGN="center">Total</TD><TD ALIGN="center">Percentage</TD>
</TR>
gofuns

  foreach my $a (sort{$humDoTSGOFun{$b} <=> $humDoTSGOFun{$a}}keys%humDoTSGOFun) {
    &printDistRow($a,$humDoTSGOFun{$a},$totalDoTShumGOFun,$musDoTSGOFun{$a},
		  $totalDoTSmusGOFun,
		  {
		   'query' => 'goFunction','pm0' => $goFunctionIdMap{$a},
		   'method' => 'dots','organism' => 'both'});
  }

  print <<bottomOfStats;
</TD></TR></TABLE>
</TABLE>
</TD></TR>
<TR><TD>
<BR CLEAR="both">
<BR>
<FONT SIZE="+1" FACE="helvetica,sans-serif"><I><B>Statistics last updated $date</B></I></FONT>
</TD></TR>
</TABLE>

</HTML>

bottomOfStats

}

##note:  $link is hashref and must have the keys (query,method,organism,pm0 (value is id to be queried for))
sub printDistRow {
  my($name,$left,$lTotal,$right,$rTotal,$link) = @_;
  my $linkName = $name;
  if ($link && scalar(keys%$link) == 4) {
    $linkName = '<A href="http://www.allgenes.org/gc/servlet?page=query&rowsPerPage=50&run=submit';
    foreach my $key (keys %$link) {
      if ($key eq 'method') {
	$linkName .= '&pm1=' . ($link->{$key} =~ /dots/i ? '%3D+4' : '%3D+5');
      } elsif ($key =~ /org/) {
	$linkName .= '&pm2=' . ($link->{$key} =~ /both/i ? '%3E+7' : ($link->{$key} =~ /^h/i ? '%3D+8' : ' %3D+14'));
      } else {
	$linkName .= "\&$key=$link->{$key}";
      }
    }
    $linkName .= "\">$name</A>"
  }
  my ($dl,$ll) = &getSizes($left,$lTotal);
  my $ldkblue = "<IMG SRC=\"http://www.allgenes.org/images/blue1.gif\" HEIGHT=\"12\" WIDTH=\"$dl\">";
  my $lltblue = "<IMG SRC=\"http://www.allgenes.org/images/lblue1.gif\" HEIGHT=\"12\" WIDTH=\"$ll\">";
  my ($dr,$lr) = &getSizes($right,$rTotal);
  my $rdkblue = "<IMG SRC=\"http://www.allgenes.org/images/blue1.gif\" HEIGHT=\"12\" WIDTH=\"$dr\">";
  my $rltblue = "<IMG SRC=\"http://www.allgenes.org/images/lblue1.gif\" HEIGHT=\"12\" WIDTH=\"$lr\">";
  print "<TR><TD ALIGN=\"left\">$linkName</TD><TD ALIGN=\"center\">$left</TD><TD>$ldkblue$lltblue\&nbsp;", &getPercent($left,$lTotal),"\%</TD>\n";
  print "<TD ALIGN=\"center\">$right</TD><TD>$rdkblue$rltblue\&nbsp;", &getPercent($right,$rTotal),"\%</TD></TR>\n";
}

sub getPercent {
  my($a,$t) = @_;
  $t = $t ? $t : 1;
  my $p = $a/$t*100;
  $p =~ s/^(\d+\.\d).*/$1/;
  return $p;
}

sub getSizes {
  my($a,$t) = @_;
  $t = $t ? $t : 1;
  my $dark = int($a/$t*200) == 0 ? 1 : int($a/$t*200);
  return ($dark,200 - $dark);
}

1;

