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
	 {o => 'inputfile',
	  t => 'string',
	  h => 'file for ePCR primer',
         },
	 {o => 'outputfile',
	  t => 'string',
	  h => 'file to write unmatched information',
         }
	 ];

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
	$self -> log($self->getArgs()->{'commit'} ? "***COMMIT ON***\n" : "**COMMIT TURNED OFF**\n"); 	
	
	my $epcrHash = $self->getEPCR();

	$self->updateDbRHMarker($epcrHash);
}

sub getEPCR(){

	my ($self) = @_;

	my $inputfile = $self->getArgs()->{'inputfile'} || die "must supply input file \n";
	print STDERR "The af inputfile is $inputfile\n";

	my $testnum;
	if ($self->getArgs()->{'testnumber'}) {
		$testnum = $self->getArgs()->{'testnumber'};
		print STDERR "Testing on $testnum on the DoTS.RHMarker table\n";
	}	

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
		my $len = $arr[3];
   		#$len = substr($len, 0, length($len)-1);
		$len =~ s/\s//;

		$epcrHash{$id} = [$forward, $backward, $len];		
		$num++;
		if ($testnum && $num >= $testnum){
			last;
		}
	}

	close(FILE);
	
	print STDERR "$num entries will be made as hash\n";

	return \%epcrHash;
}

sub updateDbRHMarker {
	my ($self, $epcrHash) = @_;

	my $num = 0;
	my $matchno = 0;
	my $len_upd_no = 0;

	my $outputfile = $self->getArgs()->{'outputfile'} || die "must supply output file \n";
	$outputfile = "> " . $outputfile;
	
	open(OUT, $outputfile) || die "Can't open the output file \n";

	foreach my $id (keys %$epcrHash) {
		my $forward = $epcrHash->{$id}->[0];
		my $backward = $epcrHash->{$id}->[1];
		my $len = $epcrHash->{$id}->[2];
		my $newDbRHMarker = GUS::Model::DoTS::RHMarker->new({'rh_marker_id'=>$id});

		$newDbRHMarker->retrieveFromDB;
		# my $temp;
		# $temp = $newDbRHMarker->get('rh_marker_id');
		# print STDERR " id = $id and temp rh_id = $temp\n";
		if ($id eq $newDbRHMarker->get('rh_marker_id') && length($forward) > 1){
			$newDbRHMarker->setForwardPrimer($forward);					
			$newDbRHMarker->setReversePrimer($backward);			
			if ($len ne $newDbRHMarker->get('product_length')){ 
				$newDbRHMarker->setProductLength($len);	
				$len_upd_no ++;
			}
			$matchno ++;
		}else{
			print(OUT $id, "\t", $forward, "\t", $backward, "\t", $len, "\n");
		}
		
		$num += $newDbRHMarker->submit();
		$newDbRHMarker->undefPointerCache();					
	}

	print STDERR "$num DBRHMarker rows processed \n";
	print STDERR "Matched rows are $matchno\n";
	print STDERR "Number of rows of updated product_length is $len_upd_no\n";
	close(OUT);
}

1;

=pod
=head1 Description
B<LoadEPCRInfo> - a plug_in that updates forward_primer and backward_primer in DbRHMarker with information from rh.sts file

=head1 Purpose
B<LoadEPCRInfor> plug_in that updates information in DbRHMarker rows.
=cut  