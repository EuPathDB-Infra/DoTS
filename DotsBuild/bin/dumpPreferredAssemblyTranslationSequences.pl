#!@perl@

## dumps sequences from sequence table 
##note the sequence must be returned as thelast item

## Brian Brunk 01/05/2000

use strict;
use lib "$ENV{GUS_HOME}/lib/perl";
use Getopt::Long;
use GUS::ObjRelP::DbiDatabase;
use CBIL::Bio::SequenceUtils;
use GUS::Common::GusConfig;

my ($gusConfigFile,$debug,$verbose,$outFile,$minLength,$taxon);
&GetOptions("verbose!"=> \$verbose,
            "outputFile=s" => \$outFile,
            "minLength=i" => \$minLength,
	    "taxon=i" => \$taxon,
            "debug!" => \$debug,
            "gusConfigFile=s" => \$gusConfigFile);

if(!$taxon || !$outFile){
	die "usage: dumpSequencesFromFile.pl --outputFile <outfile> --taxon --verbose --debug --minLength <minimum length to output [1]> --gusConfigFile [\$GUS_CONFIG_FILE]\n";
}

##set the defaults
$minLength = $minLength ? $minLength : 1;

print STDERR "Establishing dbi login\n" if $verbose;
my $gusconfig = GUS::Common::GusConfig->new($gusConfigFile);

my $db = GUS::ObjRelP::DbiDatabase->new($gusconfig->getDbiDsn(),
					$gusconfig->getReadOnlyDatabaseLogin(),
					$gusconfig->getReadOnlyDatabasePassword,
					$verbose,0,1,
					$gusconfig->getCoreSchemaName,
					$gusconfig->getOracleDefaultRollbackSegment());

my $dbh = $db->getQueryHandle();

##want to be able to restart it....
my %done;
if(-e $outFile){
	open(F,"$outFile");
	while(<F>){
		if(/^\>(\S+)/){
			$done{$1} = 1;
		}
	}
	close F;
	print STDERR "Ignoring ".scalar(keys%done)." entries already dumped\n" if $verbose;
}

my %RNA;

my $sql = "select ri.rna_id,a.na_sequence_id from dots.assembly a,dots.rnafeature rf,dots.rnainstance ri where a.taxon_id = $taxon and a.na_sequence_id = rf.na_sequence_id and rf.na_feature_id = ri.na_feature_id";

print STDERR "$sql\n";

my $stmt = $dbh->prepareAndExecute($sql) || die "Cannot execute $sql\n";

while(my ($rna,$dt) = $stmt->fetchrow_array()){
  $RNA{$dt}=$rna unless $done{"DT.".$dt}==1;
}

$stmt->finish();

my $totRNA = scalar (keys %RNA);

print STDERR "total number of RNA $totRNA\n";

open(OUT,">>$outFile");
my $count = 0;
my @row;

foreach my $dt (keys %RNA) {
  my $idSQL = "select 'DT.$dt,length of predicted protein sequence = '||aas.length,aas.sequence from dots.protein p, dots.translatedaafeature taf, dots.proteininstance pi,dots.aasequenceimp aas where p.rna_id = $RNA{$dt} and p.protein_id = pi.protein_id and pi.is_reference = 1 and pi.aa_feature_id = taf.aa_feature_id and taf.aa_sequence_id = aas.aa_sequence_id";

  my $idStmt = $dbh->prepareAndExecute($idSQL);

  @row = $idStmt->fetchrow_array();

  $count++;
  print STDERR "Getting id for $count\n" if $count % 1000 == 0;
}

&printSequence();


sub printSequence{
  #	print STDERR "$gene_id,$na_id,$description,$number,$taxon,$assembly_id,$length,$seq_ver\n";
  my $sequence = pop(@row);
  if(length($sequence) < $minLength){
    print STDERR "ERROR: $row[0] too short: ",length($sequence),"\n";
    return;
  }
	my $defline = "\>".join(' ',@row)."\n";
	print OUT $defline . CBIL::Bio::SequenceUtils::breakSequence($sequence,60);
}

