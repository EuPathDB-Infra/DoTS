package DoTS::DotsBuild::Plugin::LoadEPCRInfo;

@ISA = qw(GUS::PluginMgr::Plugin);

use strict;
use GUS::Model::DoTS::RHMarker;

sub new {
	my ($class) = @_;
	
	my $self = {};
	bless($self, $class);

	my $usage = 'Plug_in to load e-pcr primer from rh.sts file to dots.RHMarker';

	my $easycsp =
		[{o => 'testnumber',
	  	  t => 'int',
		  h => 'number of iterations for testing',
         	 },					
		 { o => 'inputfile',
		   i  => 'string',
		   h => 'rh.sts file containing e-pcr primer',
		 },
		 {o => 'outputfile',
		  t => 'string',
		  h => 'not matched e-pcr information',
         	 } ];
	$self->initialize({requiredDbVersion => {},
			       cvsRevision => '$Revision$',  # cvs fills this in!
			       cvsTag => '$Name$', # cvs fills this in!
 		               name => ref($self),
		               revisionNotes => 'make consistent with GUS 3.0',
		               easyCspOptions => $easycsp,
		               usage => $usage});
	
	return $self;
		
}

sub run {
	my $self = shift;
	my $testnum;
	$self -> log($self->getArgs()->{'commit'} ? "***COMMIT ON***\n":"**COMMIT TURNED OFF**\n"; 
	if ($self->getArgs()->{'testnumber'}) {
		$testnum = $self->getArgs()->{'testnumber'};
		print STDERR "Testing on $testnum on the RHMarker table";
	}
	my $inputfile = $self->getArgs()->{'inputfile'} || die "must supply input file \n";
	my $outputfile = $self->getArgs()->{'outputfile'} || die "must supply output file \n";
	$outputfile = "> " . $outputfile;
	my $epcrHash = $self->getEPCR();

	$self->updateDbRHMarker($epcrHash, $outputfile);
}

sub getEPCR(){
	my ($inputfile) = @_;

	my %epcrHash;
	
	my $num = 0;

	open (FILE, $inputfile) || die "Can't open the input file \n";
	
	while(<FILE>){
		chomp;
		my $line = $_;
		
		my @arr = split(/\t/,	$line);
		
		my $id = $arr[0];
		my $forward = $arr[1];
		my $backward = $arr[2];

		$epcrHash{$id} = [$forward, $backward];		
		$num++;
		if ($testnum $$ $num >= $testnum){
			last;
		}
	}

	close(FILE);
	
	print STDERR("$num entries will be made as hash");

	return \%epcrHash;
}

sub updateDbRHMarker($epcrHash, $outputfile){
	my ($epcrHash, $outputfile) = @_;

	my $num = 0;
	
	open(OUT, $outputfile) || die "Can't open the output file \n";

	foreach my $id (keys %$epcrHash) {
		my $forward = $epcrHash->{$id}->[0];
		my $backward = $epcrHash->{$id}->[1];

		my $newDbRHMarker = GUS::Model::DoTS::RHMarker->new({'rh_marker_id'=>$id, 'taxon_id'=>'8'});

		$newDbRHMarker->retrieveFromDB;

		if ($id eq $newDbRHMarker->get('rh_marker_id')){
			$newDbRHMarker->setForwardPrimer($forward);					
			$newDbRHMarker->setBackwardPrimer($backward);			
		}else{
			print(OUT $id, "\t", $forward, "\t", $backward, "\n");
		}
		
		$num += $newDbRHMarker->submit();
		$newDbRHMarker->undefPointerCache();					
	}

	print STDERR("$num DBRHMarker rows processed \n");
	close(OUT);
}

1;

=pod
=head1 Description
B<LoadEPCRInfo> - a plug_in that updates forward_primer and backward_primer in DbRHMarker with information from rh.sts file

=head1 Purpose
B<LoadEPCRInfor> plug_in that updates information in DbRHMarker rows.
=cut  