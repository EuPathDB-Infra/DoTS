#!@perl@

use strict;
use lib "$ENV{GUS_HOME}/lib/perl";
use Getopt::Long;

my ($inFiles,$outFile);
&GetOptions("inFiles=s"=>\$inFiles,"outFile=s"=>\$outFile);

my @fileArr = `ls $inFiles`;
chomp @fileArr;

my $num = scalar @fileArr;

my $numCat =  0;

foreach my $file (@fileArr) {
  `cat $file >> $outFile`;
  $numCat++;
}

print "total number file = $num : number of files concatenated = $numCat\n";

  










