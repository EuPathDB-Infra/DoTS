#!@perlLocation@

#------------------------------------------------------------------
# blastAll.pl
#
# BLAST all gene trap insertion tag sites against musDoTS.
#
# Created: Tue Oct  2 12:53:07 EDT 2001
#
# Jonathan Crabtree
#modified 12/03/02 Deborah Pinney
#------------------------------------------------------------------

$| = 1;

use strict;
use GUS::Common::Sequence;
use Getopt::Long;

my ($BLASTN,$SEQFILE,$MUSDOTS,$TARGETDIR);
&GetOptions("blastn=s"=> \$BLASTN,
            "seqfile=s" => \$SEQFILE, 
            "musdots=s" => \$MUSDOTS,
            "targetdirlogin=s" => \$TARGETDIR);

if(!$BLASTN || !$SEQFILE || !$MUSDOTS || !$TARGETDIR){
    die "usage: blastAll.pl --blastn full path of blastn file to use --seqfile fasta file of all the gene tag sequences  --musdots fasta file of the final dots sequences --targetdirlogin directory to which the blastn output will be written\n";}


my $blastSeq = sub {
    my $s = shift;
    my $dl = $s->{'defline'}; $dl =~ s/^>//;
    
    my $fsa = &GUS::Common::Sequence::toFasta($s->{'sequence'}, $dl, 80);
    my $sfile = "/tmp/$$.fsa";

    open(FF, ">$sfile");
    print FF $fsa;
    close(FF);

    my $trapId = undef;

    if ($s->{'defline'} =~ /- (\d+)/) {
	$trapId = $1;
    } else {
	die "Unable to parse " . $s->{'defline'};
    }

    my $tfile = $TARGETDIR . "/$trapId-musdots.blastn";

    if (-e $tfile) {
	print STDERR "$tfile already exists - skipping this sequence\n";
    }

    my $cmd = "$BLASTN $MUSDOTS $sfile >$tfile";

    print "BLASTing $trapId versus $MUSDOTS\n";
#    print "$cmd\n";
    system($cmd);
    unlink $sfile;
};

my $nseqs = &GUS::Common::Sequence::parseFasta($SEQFILE, $blastSeq, 0, 0);
print "Processed $nseqs sequences\n";
