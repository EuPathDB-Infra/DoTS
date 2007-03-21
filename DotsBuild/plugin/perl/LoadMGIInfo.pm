package DoTS::DotsBuild::Plugin::LoadMGIInfo;

@ISA = qw( GUS::PluginMgr::Plugin );
use GUS::PluginMgr::Plugin;

use strict;
use GUS::Model::SRes::DbRef;


$| = 1;

sub new {
  my ($class) = @_;
  my $self = {};
  bless($self,$class);

  my $purpose = <<PURPOSE;
Add information from downloaded files to MGI ids in DbRef
PURPOSE

  my $purposeBrief = <<PURPOSE_BRIEF;
Add MGI information from files
PURPOSE_BRIEF

  my $notes = <<NOTES;
NOTES

  my $tablesAffected = <<TABLES_AFFECTED;
DbRef
TABLES_AFFECTED

  my $tablesDependedOn = <<TABLES_DEPENDED_ON;
DbRef
TABLES_DEPENDED_ON

  my $howToRestart = <<RESTART;
RESTART

  my $failureCases = <<FAIL_CASES;
FAIL_CASES

  my $documentation ={ purpose => $purpose,
                       purposeBrief => $purposeBrief,
                       notes => $notes,
                       tablesAffected => $tablesAffected,
                       tablesDependedOn => $tablesDependedOn,
                       howToRestart => $howToRestart,
                       failureCases     => $failureCases };

  my $argsDeclaration =
    [
     integerArg({name => 'test_number',
                 descr => 'number of iterations for testing',
                 constraintFunc => undef,
                 reqd => 0,
                 isList => 0
		}),
     stringArg({name => 'infoFile',
                descr => 'file containing additional info from MGI',
                constraintFunc => undef,
                reqd => 1,
                isList => 0
               }),
     stringArg({name => 'geneFile',
                descr => 'file from MGI containing Entrez Gene ids',
                constraintFunc => undef,
                reqd => 1,
                isList => 0
               }),
     integerArg({name => 'external_db_release_id',
                 descr => 'db rel id for MGI',
                 constraintFunc => undef,
                 reqd => 1,
                 isList => 0
                })
	 ];

    $self->initialize({requiredDbVersion => 3.5,
		       cvsRevision => '$Revision$',
		       name => ref($self),
		       argsDeclaration => $argsDeclaration,
		       documentation => $documentation});
    return $self;
}


sub run {
  my ($self)  = shift;

  my $testnum = $self->getArg('testnumber') if $self->getArg('testnumber');

  if ($testnum) {
    $self->log ("Testing on $testnum insertions into temp table");
  }

  my $infoFile = $self->getArg('infoFile');

  my $geneFile = $self->getArg('geneFile');

  my $external_db_release_id = $self->getArg('external_db_release_id');

  my $dataHash = $self->makeDataHash($infoFile, $geneFile, $testnum);

  my $stmt = $self->updateDbRef($dataHash, $external_db_release_id);

  return "$stmt";
}

sub makeDataHash {

    my ($self,$infoFile, $geneFile, $testnum) = @_;

    my %dataHash;

    my %entryHash;

    my $num = 0;

    open (FILE, $infoFile) || die "Can't open the info file\n";

    while (<FILE>) {
      chomp;
      my $line = $_;

      my @arr = split(/\t/,$line);

      my $id = $arr[0];
      my $chrom = $arr[4];
      my $cm = $arr[3];
      my $symbol = $arr[1];
      $cm = undef if (!($cm =~ /^[\d\.]+$/));
      my $descr = $arr[2];
      $descr =~ s/^(\s+)//;
      $descr =~ s/(\s+)$//;

      $entryHash{$id}++;

      if ($entryHash{$id} > 1) {
	$self->log ("Duplicate entries for $id");
      }

      $dataHash{$id}= [$chrom,$cm,$symbol,$descr];
      $num++;

      if ($testnum && $num >= $testnum) {
	last;
      }
    }

    close(FILE);

    open (GENE, $geneFile) || die "Can't open the Entrez Gene file\n";
    my $geneNum = 0;
    while (<GENE>) {
      chomp;
      my $line = $_;
      my @arr = split(/\t/,$line);
      next if ($arr[2] ne 'O' || !$arr[8]);
      my $id = $arr[0];
      my $gene = $arr[8];
      $geneNum++;
      $dataHash{$id}[4]=$gene;
    }
    close(GENE);

    $self->log ("$num MGI entries will be processed\nThere are $geneNum corresponding gene ids");
    return \%dataHash;
}

sub updateDbRef($dataHash) {
    my ($self,$dataHash,$external_db_release_id) = @_;

    my $num = 0;

    foreach my $id (keys %$dataHash) {
	my $chromosome = $dataHash->{$id}->[0];
	my $centimorgans = $dataHash->{$id}->[1];
	my $symbol = $dataHash->{$id}->[2];
	my $remark = $dataHash->{$id}->[3];
	my $secondary_identifier =  $dataHash->{$id}->[4];
	my $newDbRef = GUS::Model::SRes::DbRef->new({'primary_identifier'=>$id,'external_database_release_id'=>$external_db_release_id});
	$newDbRef->retrieveFromDB;

	if ($chromosome ne $newDbRef->get('chromosome')) {
	    $newDbRef->setChromosome($chromosome);
	}
	if ($centimorgans ne $newDbRef->get('centimorgans')) {
	    $newDbRef->setCentimorgans($centimorgans);
	}
	if ($symbol ne $newDbRef->get('gene_symbol')) {
	    $newDbRef->setGeneSymbol($symbol);
	}
	if ($remark ne $newDbRef->get('remark')) {
	    $newDbRef->setRemark($remark);
	}
	if ($secondary_identifier ne $newDbRef->get('secondary_identifier')) {
            $newDbRef->setSecondaryIdentifier($secondary_identifier);
	  }
	$num += $newDbRef->submit();
	$newDbRef->undefPointerCache();
    }

    my $stmt = "$num DbRef rows processed";
    $self->log ("$stmt");

    return $stmt;
}



1;
__END__

