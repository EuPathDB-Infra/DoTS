package DoTS::DotsBuild::Plugin::LoadEPCR;


@ISA = qw(GUS::PluginMgr::Plugin);
use strict;
use GUS::Model::DoTS::EPCR;
use FileHandle;
use DirHandle;

sub new {
  my ($class) = @_;

  my $self = {};
  bless($self,$class);

  my $usage = 'Parses epcr output files for EPCR table in GUS. No notion of updates currentl';

  my $easycsp =
    [
     {o => 'idcol',
      h => 'Name of column from NASequenceImp to search for external ids',
      t => 'string', 
      d => 'string1'
     },
     {o => 'file',
      t => 'string',
      h => 'File name to load data from',
     },
     {o => 'dir',
      t => 'string',
      h => 'Working directory',
     },
     {o => 'log',
      t => 'string',
      h => 'Log file location (full path).',
     },
     {o => 'maptableid',
      t => 'string',
      h => 'Id for mapping table, ??? for plasmo, 2782 for genemap ',
     },
     {o => 'seqsubclassview',
      t => 'string',
      h =>  'Subclass view: Assembly for DoTS, ExternalNASequence for others',
      d => 'Assembly',
     },
     {o => 'testnumber',
      t => 'int',
      h => 'number of iterations for testing',
     },
     {o => 'start',
      t => 'string',
      h => 'Line number to start entering from',
      d => '1',
     },
    ];

  $self->initialize({requiredDbVersion => {},
		     cvsRevision => '$Revision$', # cvs fills this in!
		     cvsTag => '$Name$', # cvs fills this in!
		     name => ref($self),
		     revisionNotes => 'make consistent with GUS 3.0',
		     easyCspOptions => $easycsp,
		     usage => $usage
		    });

  return $self;
}

my $debug = 0;

sub run {
	my $self   = shift;
	my $args = $self->getArgs;
	my $algoInvo = $self->getAlgInvocation;

  print $args->{'commit'} ? "***COMMIT ON***\n" : "***COMMIT TURNED OFF***\n";
  print "Testing on $args->{'testnumber'}\n" if $args->{'testnumber'};

	############################################################
	# Put loop here...remember to undefPointerCache()!
	############################################################
	my ($idcol,$fh,$dir,$file,$count,$log,$maptableid, $seqsubclassview,$is_reversed);
	$idcol = $args->{ idcol };
	$dir = $args->{ dir };
	$file = $args->{ file };
	$fh = new FileHandle("$file", "<") || die "File $dir$file not open";
	$log = new FileHandle("$args->{ log }", ">") || die "Log not open";
        $maptableid = $args->{maptableid};
	$seqsubclassview = $args->{seqsubclassview};	
	while (my $l = $fh->getline()) {
		$count++;
		
		## Only interested sequence id, location, position, and map id 
		my ($seq_id, $loc, $order, $map_id) = split(/\s+/,$l); 
		$order =~ s/\(|\)//g;
                if ($order eq "+")
		{
		$is_reversed = 0;
		}
		else {
		$is_reversed = 1;
		}
		## Start seq_id; parese file until it is reached
		next if (defined $args->{ start } && $args->{ start } > $count);
		## Split location into start and stop
		my @locs = split /\.\./, $loc;
                

		## Build entry
		my $sgm = GUS::Model::DoTS::EPCR->new( {'start_pos' => $locs[0], 'stop_pos' => $locs[1]});
			$sgm->setIsReversed($is_reversed);
			$sgm->setMapId($map_id);
		
		## Grab na_sequence_id
		my $naid;
			$seq_id =~ s/^DT\.//;
			$naid = $seq_id;
 print STDERR "NA_ID: $naid \n" ;



		if (defined $naid){
			$sgm->setNaSequenceId($naid);
		}else {
			$log->print("No na_seq_id: $l\n");
			$algoInvo->undefPointerCache();
			next;
		}
                $sgm->setMapTableId($maptableid);
		$sgm->setSeqSubclassView($seqsubclassview);

		## Submit entry
			$log->print("SUBMITTED: $seq_id\t$map_id\t$count\n");
			$sgm->submit();

		$algoInvo->undefPointerCache();
		last if (defined $args->{testnumber} && $count >= $args->{testnumber});
		print STDERR "ENTERED: $count \n" if ($count % 1000 == 0);
	} #end while()

	print STDERR "ENTERED: $count\n";
	return "Entered $count.";
}

1;

__END__

=pod
=head1 Description
B<LoadEPCR> - Parses epcr output files for EPCR table in GUS. No notion of updates currently

=head1 Purpose
Parses epcr output files for EPCR table in GUS. No notion of updates currently

=cut
