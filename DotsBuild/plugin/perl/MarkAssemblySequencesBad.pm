############################################################
## Change Package name....
############################################################
package MarkAssemblySequencesBad;

use strict;

############################################################
# Add any specific objects (GUSdev::) here
############################################################
use Objects::GUSdev::AssemblySequence;

sub new {
  my $Class = shift;
  return bless {}, $Class
}

sub Usage {
  my $M   = shift;
  return 'Marks AssemblySequences bad from Chimera files generated from Clustering';
}

############################################################
# put the options in this method....
############################################################
sub EasyCspOptions {
  my $M   = shift;
  {

    #		test_opt1 => {
    #									o => 'opt1=s',
    #									h => 'option 1 for test application',
    #									d => 4,
    #									l => 1,	ld => ':',
    #									e => [ qw( 1 2 3 4 ) ],
    #								 },

 testnumber        => {
                       o => 'testnumber=i',
                       h => 'number of iterations for testing',
                      },
                        
 filename          => {
                       o => 'filename=s',
                       h => 'name of chimera file',
                      },
                        
 processed_category  => {
                         o => 'processed_category=s',
                         h => 'processed_category',
                         e => [ qw ( chimera repeat vector low_quality ) ],
                        },
 regex_id            => {
                         o => 'regex_id',
                         t => 'string',
                         h => 'regular expression for pulling out assembly_sequence_id',
                        },
                      }
}

sub Run {
  my $M   = shift;
  my $ctx = shift;

  my $accession;
  my $deepChimera = 0;
  my @bad;

  print $ctx->{'commit'} ? "COMMIT ON\n" : "COMMIT TURNED OFF\n";
  print "Testing on $ctx->{'testnumber'}\n" if $ctx->{'testnumber'};

  open(F,"$ctx->{'filename'}") || die "file $ctx->{'filename'} not found\n";
  if ($ctx->{'processed_category'} eq 'chimera') {
    while (<F>) {
      if (/^\>(\S+)/) {
				##record here...
        push(@bad,$accession) if $deepChimera;
        $deepChimera = 0;
        $accession = $1;
      } elsif (/Chimera\s\d+:\s1,.*numberLeft=(\d+),\snumRight=(\d+)/) {
        $deepChimera = 1 if ($1 >= 2 && $2 >= 2);
      } elsif (/^Cluster_\d+:.*Chimeras\s\((.*)\)/) {
        push(@bad,split(', ',  $1));
      }
    }
  } else {
    while (<F>) {
      if ($ctx->{cla}->{regex_id}){
        if(/$ctx->{cla}->{regex_id}/){
          push(@bad,$1);
        }
      }elsif (/\>(\S+)/) {
        push(@bad,$1);
      }
    }
  }

  close F;
  print STDERR "marking ".scalar(@bad)." AssemblySequences as $ctx->{'processed_category'}\n";

  my $count = 0;
  foreach my $ass_seq_id (@bad) {

    last if ($ctx->{'testnumber'} && $count >= $ctx->{'testnumber'});

    my $ass = AssemblySequence->new({'assembly_sequence_id' => $ass_seq_id});
    if ($ass->retrieveFromDB()) {
      $ass->setHaveProcessed(1) unless $ass->getHaveProcessed() == 1;
      $ass->setProcessedCategory($ctx->{'processed_category'}) unless $ass->getProcessedCategory() eq $ctx->{cla}->{processed_category};

      ##if low_quality want to set quality_start and quality_end so is < 50 bp...so  doesn't  turn up in an assembly..
      if($ctx->{cla}->{processed_category} eq 'low_quality'){  
        $ass->setQualityStart(1) unless $ass->getQualityStart == 1;
        $ass->setQualityEnd(1) unless $ass->getQualityEnd == 1;
      }
      $ass->resetAssemblySequence();  ##sets sequence_start and end appropriately..etc...
      $count++ if $ass->hasChangedAttributes();
      $ass->submit();
    } else {
      print "ERROR: unable to retrieve $ass_seq_id from AssemblySequence table\n";
    }
    print "$count: have_processed = 1, processed_category = $ctx->{cla}->{processed_category} for $ass_seq_id\n" if $count % 100 == 0;
		
    ##following must be in loop to allow garbage collection...
    $ctx->{'self_inv'}->undefPointerCache();
  }

  my $ret = "Marked $count AssemblySequences as $ctx->{'processed_category'}";
  print "\n$ret\n";
  return $ret; 
}

1;

__END__

