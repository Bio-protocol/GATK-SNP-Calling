#!/usr/bin/perl -w

use strict;
use warnings;

# input parameters:
# a sample name followed by a set of alignment log files

#11/16/2012 12:20:22 Fri
#/mnt/03/BSLE/alignment/PE1/B73-1a.sam	Total reads in the SAM output	701354
#/mnt/03/BSLE/alignment/PE1/B73-1a.sam	Reads could be mapped	690618
#/mnt/03/BSLE/alignment/PE1/B73-1a.sam	Passing criteria reads	530453
#/mnt/03/BSLE/alignment/PE1/B73-1a.sam	Unmapped reads	23694
#---- Parse criteria -----
#minimum number of identical bases	30
#mismatched value (include INDEL penalty) out of matched bases	5 100
#INDEL penalty used to calculate mismatched value	2
#maximum length of tail allowed at each side out of matched bases	3 75
#maximum length of the gap allowed	5000
#maximum mapping locations of each read	1

if ($#ARGV < 1) {
	print STDERR "Parameters of one sample name followed by log files are needed\n";
	exit;
}

my $total = 0;
my $mapped = 0;
my $confident_mapped = 0;
my $unmapped = 0;

foreach my $align_file (sort {$a cmp $b} @ARGV[1..$#ARGV]) {
	open (IN, $align_file) || die;
	while (<IN>) {
		chomp;
		if (/(\S+)\tTotal reads in the SAM output\t(\d+)/) { $total += $2; }
		if (/Reads could be mapped\t(\d+)/) { $mapped += $1; }
		if (/Passing criteria reads\t(\d+)/) { $confident_mapped += $1; }
		if (/Unmapped reads\t(\d+)/) { $unmapped += $1; }
	}
	close IN;
}

# output
#print "File\tTotal_Reads\tMapped\tMapped_perc\tConfidentMapped\tConfidentMapped_perc\tUnmapped\tUnmapped_perc\n";
#print "$ARGV[0]\t$total";
#printf ("\t%d\t%.1f", $mapped, $mapped/$total);
#printf ("\t%d\t%.1f", $confident_mapped, $confident_mapped/$total);
#printf ("\t%d\t%.1f", $unmapped, $unmapped/$total);
#print "\n";

# print alignment statistics:
my $dataset = $ARGV[0];
#$dataset =~ s/\.sam.*//g;
print "$dataset\tTotal reads in the SAM output\t$total\n";
print "$dataset\tReads could be mapped\t$mapped\n";
print "$dataset\tPassing criteria reads\t$confident_mapped\n";
print "$dataset\tUnmapped reads\t$unmapped\n";

# # print out the criteria:
# print STDERR "#---- Parse criteria -----\n";
# print STDERR "#minimum number of identical bases\t$identical\n";
# print STDERR "#mismatched value (include INDEL penalty) out of matched bases\t@mismatches\n";
# my $indel_penalty_value = INDEL_PENALTY;
# print STDERR "#INDEL penalty used to calculate mismatched value\t$indel_penalty_value\n";
# print STDERR "#maximum length of tail allowed at each side out of matched bases\t@maxtail\n";
# print STDERR "#maximum length of the gap allowed\t$maxgap\n";
# print STDERR "#minimum mapping score of each read\t$mappingscore\n";

