###################################
# Program: LoadEPCR                # 
# Author: Shannon McWeeney         # 
###################################

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
     {o => 'maptable',
      t => 'string',
      h => 'Name of mapping source table',
      d => 'GeneMap',
     },
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
      h => 'Id for mapping table, 471 for plasmo, 366 for genemap ',
      d => '366',
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

sub Run {
	my $M   = shift;
	my $ctx = shift;

  print $ctx->{'commit'} ? "***COMMIT ON***\n" : "***COMMIT TURNED OFF***\n";
  print "Testing on $ctx->{'testnumber'}\n" if $ctx->{'testnumber'};

	############################################################
	# Put loop here...remember to undefPointerCache()!
	############################################################
	my ($idcol,$fh,$dir,$file,$count,$log,$maptable, $maptableid, $seqsubclassview,$is_reversed);
	$idcol = $ctx->{ cla }->{ idcol };
	$dir = $ctx->{ cla }->{ dir };
	$file = $ctx->{ cla }->{ file };
	$fh = new FileHandle("$file", "<") || die "File $dir$file not open";
	$log = new FileHandle("$ctx->{ cla }->{ log }", ">") || die "Log not open";
	$maptable = $ctx->{cla}->{maptable};
        $maptableid = $ctx->{cla}->{maptableid};
	$seqsubclassview = $ctx->{cla}->{seqsubclassview};	
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
		next if (defined $ctx->{ cla }->{ start } && $ctx->{ cla }->{ start } > $count);
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
			$ctx->{ self_inv }->undefPointerCache();
			next;
		}
                $sgm->setMapTableId($maptableid);
		$sgm->setSeqSubclassView($seqsubclassview);

		## Submit entry
			$log->print("SUBMITTED: $seq_id\t$map_id\t$count\n");
			$sgm->submit();

		$ctx->{ self_inv }->undefPointerCache();
		last if (defined $ctx->{cla}->{testnumber} && $count >= $ctx->{cla}->{testnumber});
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
