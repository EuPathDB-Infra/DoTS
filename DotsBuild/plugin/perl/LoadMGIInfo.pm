package DoTS::DotsBuild::Plugin::LoadMGIInfo;

@ISA = qw(GUS::PluginMgr::Plugin);

use strict;
use GUS::Model::SRes::DbRef;

my $debug = 0;

$| = 1;

sub new {
    my ($class) = @_;
    my $self = {};
    bless($self,$class);

    my $usage = 'Plug_in to load additional information from MGI files to entries in DbRef';

    my $easycsp =
	[{o => 'testnumber',
	  t => 'int',
	  h => 'number of iterations for testing',
         },
	 {o => 'infoFile',
	  t => 'string',
	  h => 'file downloaded from MGI containing additional information for rows in DbRef',
         },
	 {o => 'geneFile',
          t => 'string',
          h => 'file downloaded from MGI containing Entrez Gene ids for rows in DbRef',
         },
	 {o => 'external_db_release_id',
	  t => 'string',
	  h => 'file downloaded from MGI containing additional information for rows in DbRef',
         }
	 ];

    $self->initialize({requiredDbVersion => 3.5,
		       cvsRevision => '$Revision$',  # cvs fills this in!
		       name => ref($self),
		       revisionNotes => 'make consistent with GUS 3.0',
		       easyCspOptions => $easycsp,
		       usage => $usage
		       });
    return $self;
}


sub run {

  my $M  = shift;
  my $ctx = shift;
  my $testnum;
  print STDERR $ctx->{cla}->{'commit'}?"***COMMIT ON***\n":"**COMMIT TURNED OFF**\n";
  if ($ctx->{'cla'}->{'testnumber'}) {
    print STDERR "Testing on $ctx->{'cla'}->{'testnumber'} insertions 
                into temp table\n" 
      if $ctx->{'cla'}->{'testnumber'};
    $testnum = $ctx->{'cla'}->{'testnumber'};
  }

  if (!$ctx->{'cla'}->{'infoFile'} || !$ctx->{'cla'}->{'external_db_release_id' || !$ctx->{'cla'}->{'geneFile'}}) {

    die "--infoFile --geneFile --external_db_release_id must be supplied\n";

  }
  my $infoFile = $ctx->{'cla'}->{'infoFile'};

  my $geneFile = $ctx->{'cla'}->{'geneFile'};

  my $external_db_release_id = $ctx->{'cla'}->{'external_db_release_id'};

  my $dataHash = &makeDataHash($infoFile, $geneFile, $testnum);

  &updateDbRef($dataHash, $external_db_release_id);    
}

sub makeDataHash {

    my ($infoFile, $geneFile, $testnum) = @_;

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
	print STDERR ("Duplicate entries for $id\n");
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

    print STDERR ("$num MGI entries will be processed\nThere are $geneNum corresponding gene ids\n");
    return \%dataHash;
}

sub updateDbRef($dataHash) {
    my ($dataHash,$external_db_release_id) = @_;

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
    print STDERR ("$num DbRef rows processed\n");
}



1;
__END__
=pod
=head1 Description
B<LoadMGIInfo> - a plug_in that udates rows in DbRef with information from MGI for C<ga> (GUS application) package.

=head1 Purpose
B<LoadMGIInfo> plug_in that updates informtion in DbRef rows.
=cut
