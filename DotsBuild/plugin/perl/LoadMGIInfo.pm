########################################################################
##LoadMGIInfo.pm 
##
##This is a ga plug_in to add information to DbRef from the MRK_List2.sql.rpt 
##and MRK_Synonym.sql.rpt files downloaded from 
##ftp://www.informatics.jax.org/pub/informatics/reports/
##
##MGI external_db_id = 4893
##
##Created Nov. 21, 2002
##
##Deborah Pinney 
##
##algorithm_id=9195      
##algorithm_imp_id=11027      
########################################################################
package LoadMGIInfo;

use Objects::GUSdev::DbRef;


my $Cfg;
my $ctx;
my $count=0;
my $debug = 0;
$| = 1;

sub new {
    my $Class = shift;
    $Cfg = shift;
    return bless {}, $Class;
}

sub Usage {
    my $M   = shift;
    return 'Plug_in to populate the NRDBEntry table';
}

sub EasyCspOptions {
    my $M   = shift;
    {
	
	testnumber  => {
	    o => 'testnumber=i',
	    h => 'number of entries for testing',
	},
	
	inputfile   => {
	    o => 'inputfile=s',
	    h => 'file downloaded from MGI containing additional information for rows in DbRef'},
    }
}

sub Run {

    my $M  = shift;
    $ctx = shift;
    my $testnum;
    print STDERR $ctx->{'commit'}?"***COMMIT ON***\n":"**COMMIT TURNED OFF**\n";
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
    
    my $dataHash = &makeDataHash($inputfile, $testnum);
    
    &updateDbRef($dataHash);    
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
    
    my $dataHash = @_;
    
    my $external_db_id = 4893;
    
    my $num = 0;
    
    foreach my $primary_identifier (%$dataHash) {
	
	my $chromosome = $dataHash->{$id}->[0];
	my $centimorgans = $dataHash->{$id}->[1];
	my $secondary_identifier = $dataHash->{$id}->[2];
	my $remark = $dataHash->{$id}->[3];
	
	my $newDbRef = DbRef->new({'primary_identifier'=>$primary_identifier,'external_db_id'=>$external_db_id});
	
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
B<LoadMGIInfo> plug_in that updates inforamtion in DbRef rows.
=cut
