package DoTS::Gene::Plugin::CreateGenomeDotsGene;
@ISA = qw( GUS::PluginMgr::Plugin);

use strict 'vars';

use CBIL::Util::Disp;
use CBIL::Util::PropertySet;

use GUS::PluginMgr::Plugin;

use DoTS::Gene::GenomeAlignmentSelector;
use DoTS::Gene::GenomeAlignmentMerger;

sub new {
    my ($class) = @_;
    my $self = {};
    bless($self,$class);

    my $purposeBrief = 'create genome-based Dots Genes from genomic alignments of DTs';
  
    my $purpose = $purposeBrief;

    my $tablesDependedOn = [['DoTS::BlatAlignment','genomic alignments of DTs']];

    my $tablesAffected = [];

    my $howToRestart = <<RESTART; 
Cannot be restarted, yet. 
RESTART

    my $failureCases;

    my $notes = <<NOTES;
=pod

=head2 F<General Description>

Create genome-based Dots Genes from genomic  alignments of DTs. 
Merger alignments 
1) with overlap of exon bases on the same strand
2) within a distance, eg 500kb, and linked by 5'-3' EST pair(s)
3) very close to each other, say <75bp apart

=head1 AUTHORS

Y Thomas Gan

=head1 COPYRIGHT

Copyright, Trustees of University of Pennsylvania 2004. 

=cut

NOTES
    
    my $documentation = {purpose=>$purpose, purposeBrief=>$purposeBrief,tablesAffected=>$tablesAffected,tablesDependedOn=>$tablesDependedOn,howToRestart=>$howToRestart,failureCases=>$failureCases,notes=>$notes};
    my $argsDeclaration  =
    [
     integerArg({name => 'taxon_id',
		 descr => 'taxon id',
		 constraintFunc=> undef,
		 reqd  => 1,
		 isList => 0 
		 }),
 
    integerArg({name => 'genome_db_rls_id',
		descr => 'genome external dabtabase release id',
		constraintFunc=> undef,
		reqd  => 1,
		isList => 0 
		}),
 
     stringArg({name => 'temp_login',
		descr => 'login for temp table usage',
		constraintFunc=> undef,
		reqd  => 1,
		isList => 0 
		}),

     stringArg({name => 'est_pair_cache',
		descr => 'table that caches info about est clone pairs',
		constraintFunc=> undef,
		reqd  => 0,
		default => 'EstClonePair',
		isList => 0 
		}),

     stringArg({name => 'blat_signals_cache',
		descr => 'table that caches info about genomic signals of blat alignments',
		constraintFunc=> undef,
		reqd  => 0,
		default => 'BlatAlignmentSignals',
		isList => 0
		}),

     integerArg({name => 'initial_gdg_id',
		 descr => 'the initial id number to use for genome dots gene (if restart, get max in db + 1)',
		 constraintFunc=> undef,
		 reqd => 0,
		 default => 1,
		 isList => 0 
		 }),

     stringArg({name => 'skip_chrs',
		descr => 'comma separated chromosomes for which analysis should be skipped, for whole genome analysis in restart mode',
		constraintFunc=> undef,
		reqd => 0,
		isList => 1, 
	    }),

     stringArg({name => 'chr',
		descr => 'chromosome on which to do analysis',
		constraintFunc=> undef,
		reqd => 0,
		isList => 0, 
	    }),

     integerArg({name => 'start',
		descr => 'start of chromosome region in which to do analysis',
		constraintFunc=> undef,
		reqd => 0,
		isList => 0 
		}),

     integerArg({name => 'end',
		descr => 'end of chromosome region in which to do analysis',
		constraintFunc=> undef,
		reqd => 0,
		isList => 0 
		}),

     integerArg({name => 'splice_minimum',
		descr => 'minimum length of maximum intron required for a DT to be used in gene construction',
		constraintFunc=> undef,
		reqd => 0,
		default => 47,
		isList => 0 
		}),

     booleanArg({name => 'contains_mrna',
		descr => 'flag to indicate whether a DT has to contain mRNA to be used in gene construction',
		reqd => 0,
		default => 0,
		}),

     booleanArg({name => 'exclude_singleton',
		descr => 'flag to indicate whether to exclude singletons from analysis',
		reqd => 0,
		default => 0,
		}),

     integerArg({name => 'overlap_merge',
		descr => 'criteria to trigger merge by exon base overlap (eg 10)',
		constraintFunc=> undef,
		reqd => 0,
		default => 1,
		isList => 0, 
	    }),

     integerArg({name => 'clonelink_merge',
		descr => 'distance criteria to trigger merge by EST clone pair linkage',
		constraintFunc=> undef,
		reqd => 0,
		default => 500000,
		isList => 0, 
	    }),

     integerArg({name => 'proximity_merge',
		descr => 'distance criteria to trigger merge by alignment proximity',
		constraintFunc=> undef,
		reqd => 0,
		default => 75,
		isList => 0, 
	    }),

     booleanArg({name => 'test',
		descr => 'run plugin in a test region',
		reqd => 0,
		default => 0,
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

sub run {
    my $self = shift;
    
    $self->logAlgInvocationId();
    $self->logCommit();
    $self->logArgs();

    my $dbh = $self->getQueryHandle();
    my $tempLogin = $self->getArg('temp_login');

    $self->log("Clean out or create temp tables to hold genome dots gene analysis result...");
    my @tmpMeta = &createTempTables($dbh, $tempLogin) unless $self->getArg('skip_chrs');

    $self->log("Creating genome-based DoTS genes...");
    my ($coords, $skip_chrs) = $self->getRegionSelections();
    my @done_chrs = @$skip_chrs;

    my $qSel = $self->getBaseQuerySelector();
    my $tSel = $self->getBaseTargetSelector();
    my $optSel = $self->getOptionalQuerySelector();
    my $mrg_criteria = $self->getMergeCriteria();
    my $gdgId = $self->getArg('initial_gdg_id');

    foreach my $coord (@$coords) {
	$tSel->{target_na_sequence_id} = $coord->{chr_id};
	$tSel->{target_start} = $coord->{start} if $coord->{start};
	$tSel->{target_end} = $coord->{end} if $coord->{end};
	my $aln_sel = DoTS::Gene::GenomeAlignmentSelector->new($dbh, $qSel, $tSel, $optSel);
	my $srt_alns = $aln_sel->getSortedAlignments();

	my $gdg_mrg = DoTS::Gene::GenomeAlignmentMerger->new($srt_alns, $mrg_criteria, $dbh);
	my $genes = $gdg_mrg->getCompositeGenomeFeatures();
	$self->log("number of genes in this region: " . scalar(@$genes));
	foreach (@$genes) {
	    $self->saveGene($_, $gdgId++, \@tmpMeta);
	}
	$dbh->commit();
	push @done_chrs, $coord->{chr};
	$self->log("completed/skipped chromosomoes: " . join(', ', @done_chrs));
    }
}

##########################

sub getBaseQuerySelect {
    my $self = shift;
    my $taxonId = $self->getArg('taxon_id');
    return { query_taxon_id => $taxonId, query_table_id => 56 };
}

sub getBaseTargetSelect {
    my $self = shift;
    my $taxonId = $self->getArg('taxon_id');
    my $genomeId = $self->getArg('genome_db_rls_id');
    return { target_taxon_id => $taxonId, target_table_id => 245,
	     target_external_db_release_id => $genomeId };
}

sub getOptionalQuerySelect {
    my $self = shift;
    my $containsMrna = $self->getArg('contains_mrna');
    my $spliceMin = $self->getArg('splice_minimum');
    my $excSgl = $self->getArg('exclude_singleton');
    my $basc = $self->getArg('temp_login') . '.' . $self->getArg('blat_signals_cache');

    return { contains_mrna => $containsMrna, splice_minimum => $spliceMin,
	     exclude_singleton => $excSgl, blat_signals_cache => $basc };
}

sub getMergeCriteria {
    my $self = shift;
    my $oM = $self->getArg('overlap_merge');
    my $pM = $self->getArg('proximity_merge');
    my $cM = $self->getArg('clonelink_merge');
    my $epc = $self->getArg('temp_login') . '.' . $self->getArg('est_pair_cache');
    return { overlap_merge => $oM, proximity_merge => $pM,
	     clonelink_merge => $cM, est_pair_cache => $epc };
}

sub getRegionSelections {
    my ($self) = @_;

    my $dbh = $self->getQueryHandle();
    my $genomeId = $self->getArg('genome_db_rls_id');

    if ($self->getArg('test')) {
	my $chr = '1';
	my $chr_id = &getChromId($dbh, $genomeId, $chr);
	return [{chr=>$chr, chr_id=>$chr_id, start=>1e6, end=>11e6}];
    }

    my $chr = $self->getArg('chr'); $chr =~ s/^chr//i;
    my $start = $self->getArg('start');
    my $end = $self->getArg('end');

    die "start or end set while chr is not defined" if (defined($start) || $end) && !$chr;
    if ($chr) {
	my $chr_id = &getChromId($dbh, $genomeId, $chr);
	return ([{ chr_id=> $chr_id, chr => $chr, start => $start, end => $end }], []);
    } else {
	my $skipChrs = $self->getArg('skip_chrs'); $skipChrs =~ s/chr//gi;
	my @skipChrs = split(/,/, $skipChrs);
	my @chroms = &getChroms($dbh, $genomeId, \@$skipChrs);
	return (\@chroms, \@skipChrs);
    }
}

sub createTempTables {
  my ($dbh, $tmpLogin) = @_;

  my $gdg_tab = 'GenomeDotsGene';
  my @gdg_cols = ('genome_dots_gene_id', 'taxon_id', 'genome_external_db_release_id',
		  'chromosome', 'chromosome_start', 'chromosome_end', 'strand', 'gene_size',
		  'number_of_exons', 'min_exon', 'max_exon', 'min_intron', 'max_intron',
		  'exonstarts', 'exonends', 'deprecated', 'est_plot_score',
		  'number_of_splice_signals', 'has_human_mouse_orthology',
		  'polya_signal_type', 'has_polya_track', 'number_of_est_libraries',
		  'number_of_est_clones', 'number_of_est_p53pairs', 'confidence_score',
		  'number_of_rnas', 'contains_mrna', 'max_orf_length', 'max_orf_score', 'min_orf_pval', 'gene_id');
  my @gdg_types = ('NOT NULL NUMBER(10)', 'NOT NULL NUMER(10)', 'NUMER(10)',
		   'VARCHAR2(32)', 'NUMBER(12)', 'NUMBER(12)', 'CHAR(1)', 'NUMBER(8)',
		   'NUMBER(4)', 'NUMBER(8)', 'NUMBER(8)', 'NUMBER(8)', 'NUMBER(8)',
		   'CLOB', 'CLOB', 'NUMBER(1)', 'NUMBER(4)',
		   'NUMBER(4)', 'NUMBER(1)',
		   'NUMBER(1)', 'NUMBER(1)', 'NUMBER(4)',
		   'NUMBER(6)', 'NUMBER(4)', 'NUMBER(3)',
		   'NUMBER(4)', 'NUMBER(1)', 'NUMBER(8)', 'FLOAT(126)', 'FLOAT(126)', 'NUMBER(10)');
  my @gdg_csts = ("alter table ${tmpLogin}.$gdg_tab add constraint PK_" . uc($gdg_tab) . " primary key ($gdg_cols[0])",
		  "create index " . uc($gdg_tab) . "_IND01 on ${tmpLogin}.$gdg_tab($gdg_cols[1],$gdg_cols[2])"); 

  my $gdt_tab = 'GenomeDotsTranscript';
  my @gdt_cols = ('taxon_id', 'genome_external_db_release_id', $gdg_cols[0], 'blat_alignment_id', 'na_sequence_id');
  my @gdt_types = ('NOT NULL NUMER(10)', 'NUMER(10)',$gdg_types[0], 'NUMBER(10)', 'NUMBER(10)');
  my @gdt_csts = ("alter table ${tmpLogin}.$gdt_tab add constraint " . uc($gdg_tab)
		  . "_FK01 foreign key ($gdt_cols[3]) references ${tmpLogin}.$gdg_tab ($gdg_cols[0])",
		  "create index " . uc($gdt_tab) . "_IND01 on ${tmpLogin}.$gdt_tab($gdt_cols[1],$gdt_cols[2])"); 

#  my $g2s_tab = 'GdgSdgMap';
#  my @g2s_cols = ('taxon_id', 'genome_external_db_release_id', $gdg_cols[0], 'similarity_dots_gene_id', 'selected');
#  my $g2s_types = ('NOT NULL NUMER(10)', 'NUMER(10)', $gdg_types[0], 'NUMBER(10)', 'NUMBER(1)');
#  my @g2s_csts = ("alter table ${tmpLogin}.$g2s_tab add constraint " . uc($g2s_tab)
#		  . "_FK01 foreign key ($g2s_cols[2]) references ${tmpLogin}.$gdg_tab ($gdg_cols[0])",
#		  "create index " . uc($g2s_tab) . "_IND01 on ${tmpLogin}.$g2s_tab($g2s_cols[0],$g2s_cols[1])"); 

#  $dbh->do("drop table ${tmpLogin}.$g2s_tab");
  $dbh->do("drop table ${tmpLogin}.$gdt_tab");
  $dbh->do("drop table ${tmpLogin}.$gdg_tab");
  
  &createTable($dbh, $tmpLogin, $gdg_tab, \@gdg_cols, \@gdg_types, \@gdg_csts);
  &createTable($dbh, $tmpLogin, $gdt_tab, \@gdt_cols, \@gdt_types, \@gdt_csts);
 # &createTable($dbh, $tmpLogin, $g2s_tab, \@g2s_cols, \@g2s_types, \@g2s_csts);

  return ($gdg_tab, \@gdg_cols, \@gdg_types, $gdt_tab, \@gdt_cols, \@gdt_types);
}

sub createTable {
    my ($dbh, $schema, $tab, $cols, $types, $constraints) = @_;

    my $num_cols = scalar(@$cols);
    die "number of columns and number of column types mismatch" if scalar(@$types) != $num_cols;

    my $sql = "CREATE TABLE ${schema}.$tab(\n";
    for (my$i=0; $i<$num_cols; $i++) {
	$sql .= $cols->[$1] . ' ' . $types->[$i] . ($i < $num_cols - 1 ? ',' : '') . "\n";
    }
    print STDERR "# creating ${schema}.$tab table...\n";
    $dbh->sqlexec($sql);
    $dbh->sqlexec("GRANT SELECT ON ${schema}.$tab to PUBLIC");
    foreach (@$constraints) { $dbh->sqlexec($_); }
}

sub getChromId {
  my ($dbh, $ext_db_rel_id, $chr) = @_;

  my $sql = "select na_sequence_id from DoTS.VirtualSequence "
    . "where external_database_release_id = $ext_db_rel_id and chromosome = '$chr'";
  my $sth = $dbh->prepareAndExecute($sql);
  my ($chr_id) = $sth->fetchrow_array;

  die "no chr_id found for $chr and external database release $ext_db_rel_id\n" unless $chr_id;
  $chr_id;
}

sub getChroms {
  my ($dbh, $ext_db_rel_id, $skipChrs) = @_;

  my $sql = "select na_sequence_id, chromosome from DoTS.VirtualSequence "
    . "where external_database_release_id = $ext_db_rel_id "
    . "and chromosome not in (" . join(', ', map { "'" . $_ . "'" } @$skipChrs). ")";
  my $sth = $dbh->prepareAndExecute($sql);

  my @chroms;
  while (my ($chr_id, $chr) = $sth->fetchrow_array) {
      push @chroms, { chr_id => $chr_id, chr => $chr };
  }

  return @chroms;
}

sub saveGene {
    my ($self, $id, $gene, $tmpMeta) = @_;

    my ($gdg_tab, $gdg_cols, $gdg_types, $gdt_tab, $gdt_cols, $gdt_types) = @$tmpMeta;

    my $dbh = $self->getQueryHandle();
    my $taxonId = $self->getArg('taxon_id');
    my $genomeId = $self->getArg('genome_db_rls_id');
    my $tmpLogin = $self->getArg('temp_login');

    my $sql = "insert into ${tmpLogin}.$gdg_tab ($gdg_cols->[0], $gdg_cols->[1], "
	. "$gdg_cols->[2], $gdg_cols->[3], $gdg_cols->[4], $gdg_cols->[5], "
	. "$gdg_cols->[6], $gdg_cols->[7], $gdg_cols->[8], $gdg_cols->[9], "
	. "$gdg_cols->[10], $gdg_cols->[11], $gdg_cols->[12], $gdg_cols->[13], "
	. "$gdg_cols->[14]) values ($id, $taxonId, $genomeId, " . "'" . $gene->getChrom . "',"
	. $gene->getChromStart . "," . $gene->getChromEnd . ",'" . $gene->getStrand . "',"
	. $gene->getTotalSpanSize . "," . $gene->getNumberOfSpans . ","
	. $gene->getMinSpanSize . "," . $gene->getMaxSpanSize . ","
	. $gene->getMinInterspanSize. "," . $gene->getMaxInterspanSize. ","
	. "'" . $gene->getSpanStarts . "'," . $gene->getSpanEnds . "')";
    $dbh->sqlexec($sql);

    my @transcripts = $gene->getConstituents();
    foreach my $t (@transcripts) {
	my $bid = $t->getId();
        my $dots = $t->getAnnotationProperties()->{na_sequence_id};

	my $sql = "insert into ${tmpLogin}.$gdt_tab "
	    . "($gdg_cols->[0], $gdg_cols->[1], $gdg_cols->[2], $gdg_cols->[3], $gdg_cols->[4]) "
	    . "values ($taxonId, $genomeId, $id, $bid, $dots)";
	$dbh->sqlexec($sql);
    }
}

1;
