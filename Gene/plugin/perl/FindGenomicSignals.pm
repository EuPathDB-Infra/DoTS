package DoTS::Gene::Plugin::FindGenomicSignals;
@ISA = qw( GUS::PluginMgr::Plugin);

use strict 'vars';

use CBIL::Util::Disp;
use CBIL::Util::PropertySet;

use GUS::PluginMgr::Plugin;

use GUS::Model::DoTS::BLATAlignment;

sub new {
    my ($class) = @_;
    my $self = {};
    bless($self,$class);

    my $purposeBrief = 'look for splice signal, polyA signal, and polyA tracks around genomic alignments of DTs.';
  
    my $purpose = <<PURPOSE; 
Post process of BlatAlignments of DTs:
1) splice signal: GT-AG (in constituent DoTS assemblies)
2) polyA signal type (AATAAA or ATTAAA and the reverse complements)
3) polyA track
PURPOSE

    my $tablesDependedOn = [['DoTS::BlatAlignment','Genomic alignment of DTs'],['DoTS::VirtualSequence', 'Genomic sequence']];
    my $tablesAffected = [];

    my $howToRestart = <<RESTART; 
Cannot be restarted, yet. 
RESTART

    my $failureCases = <<FAILURE_CASES;
DBMS run out of temp space.
FAILURE_CASES

    my $notes = <<NOTES;
=pod

=head2 F<General Description>

Plugin scan genomic alignments of DTs, examine genomic sequences around exon/intron junctions and 3' ends to in search for splice signals, polyA signals, and polyA tracks. 

=head1 AUTHORS

Y Thomas Gan

=head1 COPYRIGHT

Copyright, Trustees of University of Pennsylvania 2004. 

=cut

NOTES
    
    my $documentation = {purpose=>$purpose, purposeBrief=>$purposeBrief,tablesAffected=>$tablesAffected,tablesDependedOn=>$tablesDependedOn,howToRestart=>$howToRestart,failureCases=>$failureCases,notes=>$notes};
    my $argsDeclaration  =
    [
     stringArg({name => 'skip',
		descr => 'comma separated ids of chromosomes for which to skip process (for restart)',
		constraintFunc=> undef,
		reqd  => 0,
		isList => 1, 
	    }),
    
     stringArg({name => 'temp_login',
		descr => 'login for temp table usage',
		constraintFunc=> undef,
		reqd  => 1,
		isList => 0, 
	    }),

    integerArg({name => 'taxon_id',
		descr => 'taxon id',
		constraintFunc=> undef,
		reqd  => 1,
		isList => 0 
		}),

    integerArg({name => 'genome_db_rls_id',
		descr => 'genome external database release id',
		constraintFunc=> undef,
		reqd  => 1,
		isList => 0 
		})
    ];

    $self->initialize({requiredDbVersion => {},
		     cvsRevision => '$Revision$',
		     cvsTag => '$Name$',
		     name => ref($self),
		     revisionNotes => '',
		     argsDeclaration => $argsDeclaration,
		     documentation => $documentation
		    });
    return $self;
}

my $SIG_OPTS = { polya_signal_search_range => 30,
                 polya_track_range_cdna => 10,
                 polya_track_range_genomic => 20,
		 polya_track_min_length => 17,
		 polya_track_min_percent => 88 };
my $debug = 0;

sub run {
  my $self = shift;
    
  $self->logAlgInvocationId();
  $self->logCommit();
  $self->logArgs();

  my $dbh = $self->getQueryHandle();
  my $tempLogin = $self->getArgs->{'temp_login'};
  my $genomeId = $self->getArgs->{'genome_db_rls_id'};
  my @skip = split(/,/, $self->getArgs->{'skip'});
  my %skip_chroms; foreach (@skip) { $skip_chroms{$_} = 1; }
  my $taxonId = $self->getArgs->{'taxon_id'};
  my $t_taxon_id = $taxonId;
  my $t_table_id = 245;
  my $q_taxon_id = $taxonId;
  my $q_table_id = 56;

  # create the result table if it is not there yet, clean old result if there
  my @cols = &cleanOrCreateTempTable($dbh, $tempLogin, $genomeId, scalar(@skip));

  # get all chromosome ids
  my ($chrom_ids, $chrom_names) = &getChromosomeIdNames($dbh, $genomeId);

  # all alignments
  my ($tot, $c) = (0, 0);

  my @done_chrs;
  foreach my $cid (@$chrom_ids) {
      print "# done chromosome ids: ", join(',', @done_chrs), "\n" if scalar(@done_chrs) > 0;
      if ($skip_chroms{$cid}) {
	  print "# chromosome id $cid:  skip\n";
	  push @done_chrs, $cid;
	  next;
      } else {
	  print "# chromosome id $cid:  processing\n";
      }
      my $bids = $self->getBlatIds($dbh, $cid, $t_taxon_id, $t_table_id, $q_taxon_id, $q_table_id, $tempLogin);

      # HACK: ORACLE run out of memory periodically, entailing restarts.
      # redo the dbh for every chrom to see what happens
      $dbh = $self->getDb->makeNewHandle(); $dbh->{'LongReadLen'} = 4000000;

      my ($sql, $sth);
      foreach my $bid (@$bids) {
	  $sql = " SELECT target_na_sequence_id, is_reversed, query_na_sequence_id,"
	      . " target_start, target_end, tstarts, blocksizes"
	      . " FROM DoTS.BlatAlignment WHERE blat_alignment_id = $bid";
	  $sth = $dbh->prepareAndExecute($sql);
	  my $r = $sth->fetchrow_hashref('NAME_lc');
	  die "unexpected: no blat alignment for id: $bid\n" unless $r;
	  my $ba_info = getBlatAlignmentInfo($dbh, $r);
	  $ba_info->{blat_alignment_id} = $bid;
	  my $signals = getSignals($dbh, $ba_info, \@cols);
	  unless (++$tot % 1000) {
	      print "# processed $tot blat alignment entries\n";
	  }

	  next unless $signals;

	  $signals->{$cols[0]} = $genomeId;
	  $signals->{$cols[1]} = $cid;
	  &saveSignals($dbh, $signals, \@cols, $tempLogin);

	  unless (++$c % 200) {
	      print "# found $c informative blat alignment signals entries\n";
	      $dbh->commit();
	  }
      }
      $sth->finish if defined $sth;
      push @done_chrs, $cid;
      print "sleeping for 1 minutes\n" unless $chrom_names->{$cid} =~ /random/i;
      sleep 60 unless $chrom_names->{$cid} =~ /random/i;
      $dbh->commit();
  }
  print "# ALL DONE: processed $c blat alignment signal entries\n\n";
}

#----------------
#
# sub-routines
#
#----------------
sub saveSignals {
  my ($dbh, $signals, $cols, $tmpLogin) = @_;

  my $tmp = uc($tmpLogin);

  my @vals;
  foreach (@$cols) { push @vals, $signals->{$_}; }
  my $sql = "INSERT INTO ${tmp}.BlatAlignmentSignals(" . join(',', @$cols). ")"
          . "VALUES (" . join(',', @vals) . ")";
  $dbh->sqlexec($sql);
}

sub getBlatIds {
    my ($self, $dbh, $cid, $t_taxon_id, $t_table_id, $q_taxon_id, $q_table_id, $tmpLogin) = @_;

    my $sql =<<BA_SQL;
SELECT  blat_alignment_id
FROM DoTS.BlatAlignment
WHERE target_na_sequence_id = $cid
AND target_taxon_id = $t_taxon_id
AND target_table_id = $t_table_id
AND query_taxon_id = $q_taxon_id
AND query_table_id = $q_table_id
AND blat_alignment_id not in
(select blat_alignment_id from ${tmpLogin}.BlatAlignmentSignals
 where target_na_sequence_id = $cid)
BA_SQL

    $self->log("running $sql ...");
    my $sth = $dbh->prepareAndExecute($sql);
    my @bids;
    while (my ($bid) = $sth->fetchrow_array) { push @bids, $bid; }
    $sth->finish;
    \@bids;
}

sub cleanOrCreateTempTable {
  my ($dbh, $tmpLogin, $extDbId, $isRestart) = @_;

  my $tmp = uc($tmpLogin);

  my @cols = ('target_external_db_release_id',
              'target_na_sequence_id',
              'blat_alignment_id',
              'splice_signal_score',
              'splice_signal_score_opposite',
              'polya_signal_score',
              'polya_signal_score_opposite',
              'has_polya_track',
              'has_polya_track_opposite');

  my $sql = "select count(*) from all_tables "
          . "where table_name = 'BLATALIGNMENTSIGNALS' and owner = '${tmp}'";
  my $sth = $dbh->prepareAndExecute($sql);
  my ($c) = $sth->fetchrow_array;
  if ($c) {
      if ($isRestart) {
          print "# looks like it is a restart, do not clean table\n";
      } else {
          print "# cleaning up ${tmp}.BlatAlignmentSignals table...\n";
          $dbh->sqlexec("delete ${tmp}.BlatAlignmentSignals where target_external_db_release_id = $extDbId");
	  $dbh->commit();
      }
  } else {
    $sql = <<BAS_SQL;
CREATE TABLE ${tmp}.BlatAlignmentSignals(
    $cols[0] NUMBER(5),
    $cols[1] NUMBER(12),
    $cols[2] NUMBER(10),
    $cols[3] NUMBER(4),
    $cols[4] NUMBER(4),
    $cols[5] NUMBER(2),
    $cols[6] NUMBER(2),
    $cols[7] NUMBER(1),
    $cols[8] NUMBER(1),
    PRIMARY KEY ($cols[2])
    /* FOREIGN KEY ($cols[2]) REFERENCES DoTS.BlatAlignment */
)
BAS_SQL

    print "# creating ${tmp}.BlatAlignmentSignals table...\n";
    $dbh->sqlexec($sql);
    $dbh->sqlexec("GRANT SELECT ON ${tmp}.BlatAlignmentSignals to PUBLIC");
  }

  @cols;
}

sub getSignals {
  my ($dbh, $ba_info, $cols) = @_;

  my $pas_range = $SIG_OPTS->{polya_signal_search_range};
  my $pat_len = $SIG_OPTS->{polya_track_min_length};
  my $pct_a = $SIG_OPTS->{polya_track_min_percent};

  my ($cn_dbr,$cn_chr,$cn_bid, $cn_spl,$cn_splo, $cn_pas,$cn_paso, $cn_pat,$cn_pato) = @$cols;

  # all the pieces of info needed:
  my $rev = $ba_info->{is_reversed};
  my $strand = $rev ? '-' : '+';
  my $exon_seq = $ba_info->{exon_seq};
  my $gflank_l = $ba_info->{genomic_flank_seq_l};
  my $gflank_r = $ba_info->{genomic_flank_seq_r};
  my $intron_flanks = $ba_info->{intron_flank_seqs};

  my ($splsig, $pas_type, $has_pat) = (0, 0, 0);

  # splice signals
  my ($spl, $splo) = &spliceSignalScores($strand, $intron_flanks);

  # polyA signal
  my ($pas, $paso) = &polyaSignalScores($strand, $exon_seq, $pas_range);

  # polyA track
  my ($pat, $pato) = &hasPolyaTracks($strand, $gflank_l, $gflank_r, $pat_len, $pct_a);

  unless ($spl || $splo || $pas || $paso || $pat || $pato) { return undef; } 

  my $bid = $ba_info->{$cn_bid};
  return { $cn_bid   => $bid,
           $cn_spl   => $spl,
           $cn_splo  => $splo,
           $cn_pas   => $pas,
           $cn_paso  => $paso,
           $cn_pat   => $pat,
           $cn_pato  => $pato, };
}

sub getBlatAlignmentInfo {
  my ($dbh, $ba_info) = @_;

  # intron flanking sequences [-1,+3]..[-3,+1]
  my ($exon_bps, $intron_bps) = (1, 3); # look a little further (1 bp both sides at both ends)
  my $intron_flanks = &getIntronFlankSeqs($dbh, $ba_info, $exon_bps, $intron_bps);
  $ba_info->{intron_flank_seqs} = $intron_flanks;

  # genomic flanking sequences
  my ($gflank_l, $gflank_r) = &getGenomicFlankSeqs($dbh, $ba_info);
  $ba_info->{genomic_flank_seq_l} = $gflank_l;
  $ba_info->{genomic_flank_seq_r} = $gflank_r;

  # get exon sequence
  my $eseq = &getExonSequence($dbh, $ba_info);
  $ba_info->{exon_seq} = $eseq;


  $ba_info;
}

sub getExonSequence {
    my ($dbh, $ba_info) = @_;

    my $chr_id = $ba_info->{target_na_sequence_id};
    my $tss = $ba_info->{tstarts};
    my $bss = $ba_info->{blocksizes};

    my @tss = split(/,/, $tss);
    my @bss = split(/,/, $bss);
    my $bc = scalar(@bss);
    die "unexpected blat blocks info $tss and $bss" unless scalar(@tss) == $bc;

    my $seq;
    for (my $i=0; $i<$bc; $i++) {
	my ($s, $len) = ($tss[$i]+1, $bss[$i]);
	my $exon_seq = &getGenomicSeq($dbh, $chr_id, $s, $len);
	$seq .= $exon_seq;
    }

    $seq;
}

sub hasPolyaTracks {
    my ($strand, $gflank_l, $gflank_r, $pat_len, $pct_a) = @_;

    my ($pat, $pato) = (0, 0);


    my $ra = &isPoly('a', $gflank_r, $pat_len, $pct_a);
    my $lt = &isPoly('t', $gflank_l, $pat_len, $pct_a);
    if (($strand eq '+' && $ra) || ($strand eq '-' && $lt)) {
        $pat = 1;
    } elsif (($strand eq '+' && $lt) || ($strand eq '-' && $ra)) {
        $pato = 1;
    }

    if($debug) {
        print "DEBGU: strand=$strand, pat_len=$pat_len, pct_a=$pct_a\n";
        print "DEBUG: flankL=$gflank_l\n";
        print "DEBUG: flankR=$gflank_r\n";
        print "DEBUG: ra=$ra, lt=$lt, pat=$pat, pato=$pato\n";
        die "DEBUG: found polyA track!\n" if $pat or $pato;
    }

    ($pat, $pato);
}

sub isPoly {
    my ($what, $seq, $window_size, $pct) = @_;

    my $answer = 0;
    my $len = length $seq;
    for (my $i=0; $i<$len-$window_size; $i++)
    {
        my $sub_str = substr($seq, $i, $window_size);
        my $c = ($sub_str =~ s/$what//gi);
        my $o = ($sub_str =~ s/a|t|c|g//gi);
        next unless ($c+$o); 
        if (100 * $c/($c+$o) >= $pct) {
            $answer = 1;
            last;
        }
    }
    $answer;
}

sub polyaSignalScores {
    my ($strand, $seq, $pas_range) = @_;

    my ($pas, $paso) = (0, 0);

    my $len = length $seq;
    my $end_seq1 = substr($seq, 0, $pas_range);
    my $end_seq2 = substr($seq, $len - $pas_range - 1, $pas_range);

    if (($strand eq '+' && $end_seq2 =~ /aataaa/i) ||
	($strand eq '-' && $end_seq1 =~ /tttatt/i))
    {
	$pas = 2;
    } elsif (($strand eq '+' && $end_seq2 =~ /attaaa/i) ||
	($strand eq '-' && $end_seq1 =~ /tttaat/i))
    {
	$pas = 1;
    } elsif (($strand eq '+' && $end_seq1 =~ /aataaa/i) ||
	($strand eq '-' && $end_seq2 =~ /tttatt/i))
    {
	$paso = 2;
    } elsif (($strand eq '+' && $end_seq1 =~ /attaaa/i) ||
	($strand eq '-' && $end_seq2 =~ /tttaat/i))
    {
	$paso = 1;
    }

    ($pas, $paso);
}


sub spliceSignalScores {
  my ($strand, $intron_flanks) = @_;

  my $s = 0;
  if ($strand eq '+') {
    foreach (@$intron_flanks) { $s++ if /gt.+ag/i; }
  } else {
    foreach (@$intron_flanks) { $s++ if /ct.+ac/i; }
  }

  my $so = 0;
  if ($strand eq '-') {
    foreach (@$intron_flanks) { $so++ if /gt.+ag/i; }
  } else {
    foreach (@$intron_flanks) { $so++ if /ct.+ac/i; }
  }

  ($s, $so);
}

sub getIntronFlankSeqs {
    my ($dbh, $ba_info, $exon_bps, $intron_bps) = @_;

    my $chr_id = $ba_info->{target_na_sequence_id};
    my $tss = $ba_info->{tstarts};
    my $bss = $ba_info->{blocksizes};
    my $len = $exon_bps + $intron_bps;


    my @tss = split(/,/, $tss);
    my @bss = split(/,/, $bss);
    my $bc = scalar(@bss);
    die "unexpected blat blocks info $tss and $bss" unless scalar(@tss) == $bc;

    my $intron_flanks = [];
    for (my $i=0; $i<$bc-1; $i++) {
      my $s1 = $tss[$i] + $bss[$i] + 1 - $exon_bps;
      my $s2 = $tss[$i+1] + 1 - $intron_bps;
      next if $s1+$len >= $s2;
      my $fs1 = &getGenomicSeq($dbh, $chr_id, $s1, $len);
      my $fs2 = &getGenomicSeq($dbh, $chr_id, $s2, $len);
      push @$intron_flanks, "$fs1$fs2";
    }

    $intron_flanks;
}

sub getGenomicFlankSeqs {
    my ($dbh, $ba_info, $bps_i, $bps_o) = @_;

    my $chr_id = $ba_info->{target_na_sequence_id};
    my $cs = $ba_info->{target_start};
    my $ce = $ba_info->{target_end};

    my $bps_i = $SIG_OPTS->{polya_track_range_cdna};
    my $bps_o = $SIG_OPTS->{polya_track_range_genomic};
    my $gflank_l = &getGenomicSeq($dbh, $chr_id, $cs - $bps_o, $bps_o + $bps_i);
    my $gflank_r = &getGenomicSeq($dbh, $chr_id, $ce - $bps_i - 1, $bps_i + $bps_o);

    ($gflank_l, $gflank_r);
}

sub getGenomicSeq {
    my ($dbh, $chr_id, $s, $len) = @_;

    my $sql = "select substr(sequence, $s, $len) "
	    . "from DoTS.VirtualSequence "
	    . "where na_sequence_id = $chr_id";
    my $sth = $dbh->prepareAndExecute($sql);
    my $seq;
    if ($seq = $sth->fetchrow_array) { ; }
    die "could not get $len bp genomic seq for chr $chr_id at $s" unless defined $seq;
    $sth->finish;

    $seq;
}

sub getChromosomeIdNames {
  my ($dbh, $ext_db_rel_id, $skip_chroms) = @_;

  my @chrom_ids; 
  my %chrom_nams;
  my $sql = "select na_sequence_id, chromosome from DoTS.VirtualSequence "
    . "where external_database_release_id = $ext_db_rel_id order by chromosome_order_num asc";
  my $sth = $dbh->prepareAndExecute($sql);
  while (my ($sid, $nam) = $sth->fetchrow_array) {
    if($skip_chroms->{$sid}) {
      print "# chromosome id $sid: skip\n";
    } else {
      push @chrom_ids, $sid;
      $chrom_nams{$sid} = $nam;
    }
  }

  (\@chrom_ids, \%chrom_nams);
}

1;
