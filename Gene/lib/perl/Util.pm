#!/usr/bin/perl

# -------------------------------------------------------
# Util.pm
#
# Created: Tue May 21 12:07:23 EDT 2002
#
# Jonathan Crabtree
# -------------------------------------------------------

package Util;

use strict;
use DBI;

# -------------------------------------------------------
# Configuration
# -------------------------------------------------------

my $DISPLAY_URL = 'http://www.cbil.upenn.edu/cgi-bin/dev/mouse-chr5/displayRegion.pl';
my $UCSC_BROWSER_URL = 'http://genome.ucsc.edu/cgi-bin/hgTracks';

# -------------------------------------------------------
# Subroutines
# -------------------------------------------------------

$Util::GUS_DBH = undef;
$Util::GUSDEV_DBH = undef;

# Get a readonly database login.
#
sub getLogin {
    return &getGusdevLogin;
}
sub getGusLogin {
    return &_getLogin(0);
}
sub getGusdevLogin {
    return &_getLogin(1);
}
sub _getLogin {
    my ($dev) = @_;

    my $host = ($dev ? 'erebus' : 'nemesis');
    my $sid = ($dev ? 'gusdev' : 'gus');
    my $user = ($dev ? 'gusdevreadonly' : 'gusreadonly');
    my $pass = ($dev ? '' : '');

    if ($dev) {
	if (!defined($Util::GUSDEV_DBH)) {
	    $Util::GUSDEV_DBH = DBI->connect("dbi:Oracle:host=$host;sid=$sid", $user, $pass);
	    $Util::GUSDEV_DBH->{'LongReadLen'} = 4000000;
	}
	return $Util::GUSDEV_DBH;
    } else {
	if (!defined($Util::GUS_DBH)) {
	    $Util::GUS_DBH = DBI->connect("dbi:Oracle:host=$host;sid=$sid", $user, $pass);
	    $Util::GUS_DBH->{'LongReadLen'} = 4000000;
	}	
	return $Util::GUS_DBH;
    }
}


# Run query that returns at most one row and return the results
#
sub runOneRowQuery {
    my($sql) = @_;
    my $result = undef;
    my $dbh = &getLogin();
    my $sth = $dbh->prepare($sql);
    $sth->execute();

    if (my $h = $sth->fetchrow_hashref('NAME_lc')) {
	my %copy = %$h;
	$result = \%copy;
    }
    $sth->finish();
    return $result;
}

# Run a query and display the results
#
sub displayQuery {
    my($sql, $formatStr, $highlightFn) = @_;
    my $dbh = &getLogin();
    my $sth = $dbh->prepare($sql);
    my $n = 0;

    $sth->execute();

    while (my @a = $sth->fetchrow_array()) {

	# Make link to sequence browser
	#
	my $last = scalar(@a);
	my $naSeqId = $a[--$last];
	my $isConsist = $a[--$last];
	my $tEnd = $a[--$last];
	my $tStart = $a[--$last];
	my $url = $DISPLAY_URL . "?naSeqId=$naSeqId&start=$tStart&end=$tEnd&width=800";

	if (defined($highlightFn) && &$highlightFn(@a)) {
	    print "<B>";
	    printf($formatStr, @a);
	    print "  <A HREF=\"$url\">view</A>";
	    print "</B>\n";
	} else {
	    printf($formatStr, @a);
	    print "  <A HREF=\"$url\">view</A>";
	    print "\n";
	}
	++$n;
    }

    $sth->finish();

    print "\n$n rows returned\n";
}

# Get basic information about the target genomic sequence.
#
sub getSeqInfo {
    my($naSeqId) = @_;

    my $sql = ("select vs.na_sequence_id, vs.source_id, vs.taxon_id, vs.length, ed.name, t.scientific_name " .
	       "from VirtualSequence vs, ExternalDatabase ed, Taxon t " .
	       "where vs.na_sequence_id = $naSeqId " .
	       "and vs.external_db_id = ed.external_db_id " .
	       "and vs.taxon_id = t.taxon_id ");

    return &Util::runOneRowQuery($sql);
}

# Get na seq id for given chr.
#
sub getNaSeqId {
    my($ext_did, $chr) = @_;

    my $sql = ("select na_sequence_id from VirtualSequence " .
	       "where external_db_id = $ext_did and chromosome = '$chr'");
    return &Util::runOneRowQuery($sql)->{'na_sequence_id'};
}

# Get the URL to display a particular region in the UCSC browser.
#
sub getUcscBrowserUrl {
    my($db, $chrom, $start, $end, $pixwidth) = @_;
    return ($UCSC_BROWSER_URL . 
	    "?db=$db" .
	    "&position=${chrom}:${start}-${end}" .
	    "&pix=${pixwidth}");
}

1;

