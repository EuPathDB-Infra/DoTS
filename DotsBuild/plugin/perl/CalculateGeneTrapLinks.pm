package DoTS::DotsBuild::Plugin::CalculateGeneTrapLinks;

@ISA = qw(GUS::PluginMgr::Plugin);

use strict;

use CBIL::Bio::BLAST2::BLAST2;
use GUS::Model::DoTS::GeneTrapAssembly;

my $VERSION = '$Revision$'; $VERSION =~ s/Revision://; $VERSION =~ s/\$//g; $VERSION =~ s/ //g;

sub new {
    my ($class) = @_;
    
    my $self = {};
    bless($self,$class);
    
    my $usage = 'Calculate the correspondence between gene trap tag sequences and DoTS assemblies';
    
    my $easycsp =
	[{o => 'external_db_release',
	  t => 'int',
	  h => 'external_db_release_id of the gene trap tag sequences.',
         },
	 {o => 'blast_dir',
	  t => 'string',
	  h => 'Directory containing BLASTN output files from searching tag sequences against DoTS assemblies.',
	 },
	 {o => 'min_pct_id',
	  t => 'float',
	  h => 'Minimum percent identity for a match/HSP to be loaded.',
	  d => 90.0,
         },
	 {o => 'min_pct_len',
	  t => 'float',
	  h => "Minimum percent of the tag's length for a match/HSP to be loaded.",
	  d => 50.0,
	 },
         {o => 'logfile',
	  t => 'string',
          h => 'log used for output and restart',
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
    my $M = shift;
    my $ctx = shift;
    my $dbh = $ctx->{'self_inv'}->getQueryHandle();
    
    my $blastDir = $ctx->{cla}->{blast_dir};
    my $extDbRel = $ctx->{cla}->{external_db_release};
    my $minPctId = $ctx->{cla}->{min_pct_id};
    my $minPctLen = $ctx->{cla}->{min_pct_len};
    my $logfile = $ctx->{cla}->{logfile};

    # Query for and cache the na_sequence_id, name, and length of each tag
    #

    my $restartHash = &getRestartHash($logfile) if (-s $logfile);
    my $tagSeqs = &getTagSequences($dbh, $extDbRel);
    my $nTags = scalar(@$tagSeqs);
    my $nLinks = 0;
    my $failcnt = 0;
    # Examine BLAST results for each tag
    #
    foreach my $t (@$tagSeqs) {
	my $naSeqId = $t->{na_sequence_id};
	my $srcId = $t->{source_id};
	my $len = $t->{length};

	my $bfile = $blastDir . "/$naSeqId-musdots.blastn";

	if (!-e $bfile) {
	    print STDERR "ERROR - could not find $bfile\n";
	}
        eval {
	  $nLinks += &processBLASTResults($t, $bfile, $minPctId, $minPctLen,$restartHash);
	  };
	&handleFailure($t, $failcnt, $ctx->{cla}, $@) if ($@);
    }

    my $summary = "Generated $nLinks links for $nTags gene trap tag sequence(s).";
    print STDERR $summary, "\n";
    return $summary;
}
sub getRestartHash {
    my ($logfile) = @_;
    my %restartHash;

    open (LOG, $logfile) || die "Can't open $logfile\n";

    while (<LOG>) {
      chomp;
      if ($_ =~ /for\s(\S+)$/) {
	$restartHash{$1}=1;
      }
    }

    close (LOG);
    return \%restartHash;
}

sub getTagSequences {
    my($dbh, $extDbRel) = @_;
    my $seqs = [];

    my $q = ("select na_sequence_id, source_id, length " .
	     "from dots.ExternalNASequence " .
	     "where external_database_release_id = $extDbRel");

    my $sth = $dbh->prepare($q);

    $sth->execute();
    while (my $h = $sth->fetchrow_hashref('NAME_lc')) {
	my %h2 = %$h;
	push(@$seqs, \%h2);
    }
    $sth->finish();

    return $seqs;
}

sub processBLASTResults {
    my($tag, $bfile, $minPct, $minLenPct,$restartHash) = @_;
    my $tagId = $tag->{na_sequence_id};
    my $tagSrcId = $tag->{source_id};
    my $restart = "0 (done in previous run - skipping)\n";
    return $restart if (${$restartHash->{$tagSrcId}} == 1);
    my $tagLen = $tag->{length};
    my $minMatchLen = ($minLenPct / 100.0) * $tagLen;

    my $b = CBIL::Bio::BLAST2::BLAST2::parseBLAST2output("cat $bfile |");
    my $ns = $b->getNumSbjcts();
    my $numMatches = 0;
    my $numMeetingCriteria = 0;

    # Loop over subjects
    #
    for (my $s = 0;$s < $ns;++$s) {
	my $sbj = $b->getSbjct($s);

	my $nh = $sbj->getNumHSPs();

	# Loop over HSPs; worry about multiple HSPs per subject later
	#
	for (my $h = 0;$h < $nh;++$h) {
	    ++$numMatches;
	    my $bestHit = (($s == 0) && ($h == 0)) ? 1 : 0;
	    my $hsp = $sbj->getHSP($h);
	    my $ids = $hsp->{'identities'};
	    my $len = $hsp->{'length'};
	    my $pct = int(($ids/$len) * 100.0);

	    my $meetsCriteria = (($pct >= $minPct) && ($len >= $minMatchLen));

	    if ($meetsCriteria) {
	      ++$numMeetingCriteria;
	      my($dotsId) = ($sbj->{'description'} =~ /(\d+)\s\[Mus musculus\]/);

	      #print "Match: $tagSrcId ($tagId) against $dotsId ";
	      #print " $pct% over $len/$tagLen bp bestHit=$bestHit\n";

	      my $obj = GUS::Model::DoTS::GeneTrapAssembly->new({
		    'tag_na_sequence_id' => $tagId,
		    'assembly_na_sequence_id' => $dotsId,
		    'is_best_match' => $bestHit,
		    'match_start' => $hsp->{s_start},
		    'match_end' => $hsp->{s_end},
		    'is_reversed' => $hsp->isReversed() ? 1 : 0,
		    'percent_identity' => $pct
		});

		$obj->submit();
		$obj->undefPointerCache();
	    }
	}
    }

    print "Loaded ${numMeetingCriteria}/${numMatches} match(es) for $tagSrcId\n";
    print STDERR "WARNING - no BLAST hits for $tagSrcId\n" if ($numMatches == 0);
    print STDERR "WARNING - no qualifying hits found for $tagSrcId\n" if ($numMeetingCriteria == 0);

    return $numMeetingCriteria;
}

sub handleFailure {
  my ($t, $failCnt, $cla, $errMsg) = @_;
  my $naSeqId = $t->{na_sequence_id};
  my $srcId = $t->{source_id};

  my $failTol = 100;

  die "More than 100 entries failed.  Aborting."
    if ($failCnt++ > $failTol);

  print "Failed: $naSeqId, $srcId\n$errMsg\n";
}


# ----------------------------------------------------------
# Generate Algorithm/AlgorithmImplementation
# ----------------------------------------------------------

if ($0 !~ /ga$/i) {

    my $usg = 'Calculate the correspondence between gene trap tag sequences and DoTS assemblies.';
    my $name = $0; $name =~ s/\.pm$//; $name =~ s/^.+\///;
    my $md5 = `/usr/bin/md5sum $0`;
    chomp $md5;
    $md5 =~ s/^(\S+).+/$1/;

    print <<XML;
<Algorithm xml_id="1001">
  <name>$name</name>
  <description>$usg</description>
</Algorithm>

<AlgorithmImplementation xml_id="1002" parent="1001">
  <version>$VERSION</version>
  <executable>$0</executable>
</AlgorithmImplementation>
XML
}

1;

