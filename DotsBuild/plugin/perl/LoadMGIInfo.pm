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
	 {o => 'inputfile',
	  t => 'string',
	  h => 'file downloaded from MGI containing additional information for rows in DbRef',
         }
	 {o => 'external_db_release_id',
	  t => 'string',
	  h => 'file downloaded from MGI containing additional information for rows in DbRef',
         }
	 ];
    
    $self->initialize({requiredDbVersion => {},
		       cvsRevision => '$Revision$',  # cvs fills this in!
		     cvsTag => '$Name$', # cvs fills this in!
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

    if (!$ctx->{'cla'}->{'inputfile'}) {

	die "--inputfile must be supplied\n";

    }
    my $inputfile = $ctx->{'cla'}->{'inputfile'};

    my $external_db_release_id = $ctx->{'cla'}->{'external_db_release_id'};
    
    my $dataHash = &makeDataHash($inputfile, $testnum);
    
    &updateDbRef($dataHash, $external_db_release_id);    
}

sub makeDataHash {
    
    my ($inputfile, $testnum) = @_;
    
    my %dataHash;
    
    my %entryHash;

    my $num = 0;
    
    open (FILE, $inputfile) || die "Can't open the input file\n"; 
    while (<FILE>) {
	if ((/-/) && (/^[\s-]+$/)) {
	    $_ = <FILE>;
	    last;
	}
    }
    
    while (<FILE>) {
	last if (/^\s*$/);
	
	if (/^ (MGI:\d+)\s+(\d+|UN|X|Y|XY|MT)\s+([\d\.]+|syntenic|N\/A|NULL)\s+(\S.*\S|\S)\s+([.\s]+)$/) {
	    my $id = $1;
	    my $chrom = $2;
	    my $cm = $3;
	    my $symbol = $4;
	    $cm = undef if (!($cm =~ /^[\d\.]+$/));
	    
	    my $descr = $5;
	    $descr =~ s/^(\s+)//;
	    $descr =~ s/(\s+)$//;
	    
	    $entryHash{$id}++;
	    
	    if ($entryHash{$id} > 1) {
		print STDERR ("Duplicate entries for $id\n");
	    }
	    
	    $dataHash->{$id}= [$chrom,$cm,$symbol,$descr];

	    $num++;

	    if ($num >= $testnum) {
		
		last;
	    }		  
	    
	} else {
	    print STDERR ("Unable to parse '$_' \n");
	}
    }
    close(FILE);
    
    print STDERR ("$num MGI entries will be processed\n");
    
    return \%dataHash;
}

sub updateDbRef($dataHash) {
    
    my ($dataHash,$external_db_release_id) = @_;
    
    
    my $num = 0;
    
    foreach my $primary_identifier (%$dataHash) {
	
	my $chromosome = $dataHash->{$id}->[0];
	my $centimorgans = $dataHash->{$id}->[1];
	my $secondary_identifier = $dataHash->{$id}->[2];
	my $remark = $dataHash->{$id}->[3];
	
	my $newDbRef = GUS::Model::SRes::DbRef->new({'primary_identifier'=>$primary_identifier,'external_db_release_id'=>$external_db_release_id});
	
	$newDbRef->retrieveFromDB;
	
	if ($chromosome ne $newDbRef->get('chromosome')) {
	    
	    $newDbRef->setChromosome($chromosome);
	}
	if ($centimorgans ne $newDbRef->get('centimorgans')) {
	    
	    $newDbRef->setCentimorgans($centimorgans);
	}
	if ($secondary_identifier ne $newDbRef->get('secondary_identifier')) {
	    
	    $newDbRef->setSecondaryIdentifier($secondary_identifier);
	}
	if ($remark ne $newDbRef->get('remark')) {
	    
	    $newDbRef->setRemark($remark);
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
