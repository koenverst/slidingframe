#!/usr/bin/perl

## This script performs an iterative analysis of BLASTs
## to find frames/windows that are conserved in a given set
## of genomes, but absent in another set.


use strict;
use warnings;
use diagnostics;
use File::Copy qw(copy);
use Getopt::Long;
use Bio::SeqIO;
use Bio::Seq;

##############################
## Variables
##############################

my $showhelp =			 0;		#Equals 1 when user indicates to display help text
my $verbose =			 0;		#Equals 1 when user indicates to run in verbose mode
my $blastversion =		"";		#Local version of blastn. Only needed to check if blastn is available

my $queryfilename =		"";		#Filename for the query sequence that will be blasted
my $outputfilename =		"";		#Filename for the output file for the user
my $tmpresultsfile =		"";		#File that will hold the frames until the next genome is analysed
my $countfilename =		"";		#Name of the file that will show the number of frames left after each iteration
my $refgenomefile =		"";		#Filename of the genome that will act as reference
my $outgroupdbname =		"";		#Name of the BLAST db with the outgroup sequences
my $ingroupfilename =		"";		#File that contains the filename of ingroup sequence that is processed

# Initialize critical values when interpreting BLAST-results
my $idOut =			 0;		#id-percentage of out-group sequences
my $idIn =			 0;		#id-percentage of in-group sequences
my $coverageOut =		 0;		#Query coverage of out-group sequences
my $coverageIn =		 0;		#Query coverage of in-group sequences

my $windowsize =		 0;		#Number of nt to be BLAST'ed
my $shiftsize =			 0;		#Number of nt to shift the frame after each iteration
my $blastcmd = 			"";		#bash command to run blast
my $blastresult =		"";		#Result that is returned by blastn
my $windowsequence =		"";
my $genomesanalysed =	 	 0;		#Number of genomes analysed/number of iterations
my $framesleft =		 0;		#Number of frames that are left after an iteration


##############################
## Arguments
##############################

## These arguments are extracted from the command line.

GetOptions(		"help"			        => \$showhelp,
				"verbose"		=> \$verbose,
				"framesize=i"	        => \$windowsize,
				"shift=i"		=> \$shiftsize,
				"reference=s"	        => \$refgenomefile,
				"output=s"		=> \$outputfilename,
				"outgroupdb=s"	        => \$outgroupdbname,
				"ingrouplist=s"	        => \$ingroupfilename,
				"idoutgroup=i"		=> \$idOut,
				"idingroup=i"		=> \$idIn,
				"coverageoutgroup=i"	=> \$coverageOut,
				"coverageingroup=i"	=> \$coverageIn) 
				or die ("Error in command line arguments. Use -help to obtain more information.\n");
				
if ( $showhelp ) { help(); }


##############################
## Pre-processing
##############################

# If script runs in verbose mode: print start time
if ( $verbose ) { warn scalar(localtime()) . " Script started\n"; }

	
# Check if blastn is available on machine
open( BLASTVER, "-|", "blastn -version");
$blastversion = <BLASTVER>;
close( BLASTVER );
if (!$blastversion) {
	die "ERROR: BLAST (blastn) is not available on this machine.\n";
}

# Use process id to create temporary filenames (pid is used in case multiple instances
# of the script run simultaneously)
$queryfilename = "query-". $$ . ".fas";
$tmpresultsfile = "tmpresults-" . $$ . ".list";


## Load reference genome, chop into pieces and put in tab-separated format
$countfilename = "framecount.list";
$framesleft = chop_to_pieces($refgenomefile, $outputfilename, $windowsize, $shiftsize);
open(COUNT, ">$countfilename");
print COUNT "iteration" . "\t" . "filename" . "\t" . "frames\n";
print COUNT $genomesanalysed . "\t" . $refgenomefile . "\t" . $framesleft . "\n";
close(COUNT);


##############################
## Main program
##############################

# When in verbose mode: display info message
if ( $verbose ) { warn scalar(localtime()) . " Analysing out-group\n"; }
	
# BLAST outgroup first, because it will significantly reduce the number of frames
# that will be left to be analysed.
$blastcmd = "blastn -db $outgroupdbname -query $queryfilename -outfmt \"6 pident qcovs\" -num_alignments 1";

open(TMP, ">$tmpresultsfile");
open(FRAMES, "<$outputfilename");


$framesleft = 0;

# For every frame in the file
while ( my $line = <FRAMES> )  {
	open(QUERY, ">$queryfilename");
	print QUERY "\>query\n";
	print QUERY return_sequence($line);
	close(QUERY);
	
	open(BLAST, "-|", $blastcmd);
	$blastresult = <BLAST>;
	close(BLAST);
	
	#Compare results to the cut-off values provided by the user
	#Checking if BLAST gives any result at all is duplicated in the sub
	if ( $blastresult ) {
		if ( compare_blast("out", $blastresult, $idOut, $coverageOut) ) {
			#If compare_blast() evalueates to TRUE then it means it meets the user's cut-offs
			print TMP $line;
			$framesleft++;
		}
	} else {
		print TMP $line;
		$framesleft++;
	}
}

close(TMP);
close(FRAMES);;

# Write number of frames to count-file
$genomesanalysed++;
open(COUNT, ">>$countfilename");
print COUNT $genomesanalysed . "\t" . $outgroupdbname . "\t" . $framesleft . "\n";
close(COUNT);

# Tmp results file becomes actual results file
copy $tmpresultsfile, $outputfilename;

# BLAST ingroup
# Using BLAST subject (FASTA file)

# When in verbose mode: display info message
if ( $verbose ) { 	warn scalar(localtime()) . " Analysing in-group\n"; }

open(INGROUP, "<", $ingroupfilename);
# For every file that is in this list, and thus should be analysed as in-group genome
while ( $ingroupfilename = <INGROUP> ) {

	chomp($ingroupfilename);

	# When in verbose mode: display info message
	if ( $verbose ) { 	warn scalar(localtime()) . " Analysing in-group: " . $ingroupfilename  . "\n"; }

	$blastcmd = "blastn -subject \"" .  $ingroupfilename . "\" -query \"" .  $queryfilename . "\" -outfmt \"6 pident qcovs\" -num_alignments 1";
	
	open(TMP, ">", $tmpresultsfile);
	open(FRAMES, "<", $outputfilename);

	$framesleft = 0;

	#For every frame in the file
	while ( my $line = <FRAMES> )  {
		open(QUERY, ">$queryfilename");
		print QUERY "\>query\n";
		print QUERY return_sequence($line);
		close(QUERY);
	
		open(BLAST, "-|", $blastcmd);
		$blastresult = <BLAST>;
		close(BLAST);
	
		#Compare results to the cut-off values provided by the user
		if ( $blastresult ) {
			if ( compare_blast("in", $blastresult, $idIn, $coverageIn) ) {
				#If compare_blast() evalueates to TRUE then it means it meets the user's cut-offs
				print TMP $line;
				$framesleft++;
			}
		}
	}

	close(TMP);
	close(FRAMES);
	
	# After every genome of the in-group: write number of frames to the count file
	$genomesanalysed++;
	open(COUNT, ">>$countfilename");
	print COUNT $genomesanalysed . "\t" . $ingroupfilename . "\t" . $framesleft . "\n";
	close(COUNT);
	
	# Tmp results file becomes actual results file
	copy $tmpresultsfile, $outputfilename;
}
close(INGROUP);

if ( $verbose ) { warn scalar(localtime()) . " Iterations finished\n"; }


##############################
## Cleanup
##############################

# Remove temporary files
system('rm ' . $queryfilename);
system('rm ' . $tmpresultsfile);


# If script runs in verbose mode: print end time
if ( $verbose ) { warn scalar(localtime()) . " Script finished\n"; }

exit;

##############################
## Subroutines
##############################

# This subroutine takes a fasta file and for each entry performs the following routine:
# Take the first n bases (=frame) and write to the output-file. Shift the frame k nucleotides to the right and repeat.
sub chop_to_pieces {

	my $genome =		   "";			#Holds the reference genome sequence
	my $entryseq =		   "";			#Holds the sequence of the current entry of the ref genome
	my $windowstart =	    0;			#Start position of current frame
	my $lastwindowstart =	0;
	my $entrylength =	    0;			#Length of the current fasta-entry in nt
	my $numberofframes =	0;			#Total number of frames
	
	my ($refgenomefile, $outputfilename, $windowsize, $windowshift) = @_;
	
	
	$genome = Bio::SeqIO->new(-file => $refgenomefile, -format => 'fasta');
	open(CHOPPED, ">$outputfilename");
	while ( $entryseq = $genome->next_seq ) {
		$windowstart = 1;
		$entrylength = $entryseq->length();
	
		#For each window in this fasta entry
		while ( $windowstart < $entrylength ) {
			if ( $windowstart + $windowsize < $entrylength ) {
				print CHOPPED $entryseq->display_id . "\t" . $windowstart . "\t" . $entryseq->subseq($windowstart, $windowstart + $windowsize) . "\n";
			} else {
				#If total entry length is smaller than windowsize -> use entire entry
				if ( $entrylength < $windowsize ) {
					print CHOPPED $entryseq->display_id . "\t" . $windowstart . "\t" . $entryseq->seq . "\n";
				} else {
				#If total entry length is larger than windowsize -> use last possible window
					$lastwindowstart = $entrylength - $windowsize +1;
					print CHOPPED $entryseq->display_id . "\t" . $lastwindowstart . "\t" . $entryseq->subseq($lastwindowstart, $entrylength) . "\n";
				}
			}
			$numberofframes++;
			$windowstart = $windowstart + $windowshift;
		} #End while-loop: no more windows in this fasta entry
	} #End while-loop: no more entries in this fasta file
	close(CHOPPED);
	
	return $numberofframes;
	
}


sub return_sequence {

	my ( $line ) = @_;
	my @frame = split("\t", $line);
	return $frame[2];

}


sub compare_blast {

	my ($inout, $blastresult, $identity, $querycoverage) = @_;
	
	my @results = split("\t", $blastresult);
	
	## For out-group
	if ($inout eq "out") {
		if (!$blastresult) {
			return 1;
		} else {
			if (($results[0] < $identity) | ($results[1] < $querycoverage) ) {
				return 1;
			} else {
				return 0;
			}
		}
	}
	
	## For in-group
	if ($inout eq "in") {
		if (($results[0] >= $identity) & ($results[1] >= $querycoverage) ) {
			return 1;
		} else {
			return 0;
		}
	}
}


sub help {
	print "\nNAME\n";
	print "\tIterative script for a sliding frame analysis of genome sequences.\n";

	print "\nSYNOPSIS\n";
	print "\tperl iterative-slidingwindow.pl -verbose -framesize [nt] -shift [nt]\n";
	print "\t-reference [FILE] -output [FILE] -outgroupdb [dbname] -ingrouplist [FILE]\n";
	print "\t-idoutgroup [%] -idingroup [%] -coverageoutgroup [%] -coverageingroup [%]\n";
	print "\t-help\n\n";

	print "DESCRIPTION\n";
	print "\t-verbose\n";
	print "\t\tDisplay information about the progress of the script.\n";

	print "\t-framesize [nt]\n";
	print "\t\tNumber of nucleotides that should be BLAST-ed.\n";

	print "\t-shift [nt]\n";
	print "\t\tNumber of nucleotides the frame should shift after a BLAST is performed.\n";

	print "\t-reference [FILE]\n";
	print "\t\tFilename of the reference genome in fasta format. This is the genome where\n";
	print "\t\tthe frame will slide over.\n";

	print "\t-output [FILE]\n";
	print "\t\tName of the outpout file. The file format needs to be described.\n";

	print "\t-outgroupdb [dbname]\n";
	print "\t\tName of the BLAST database that contains the outgroup genome sequences.\n";
	print "\t\tThese are requested in a BLAST db format because it will speed up the analysis\n";
	print "\t\tin larger data sets.\n";

	print "\t-ingrouplist [FILE]\n";
	print "\t\tFile that contains a list with filenames of the in-group genomes. The files that\n";
	print "\t\tare pointed towards in this list need to be in FASTA format. The script will\n";
	print "\t\tanalyse all the genomes in the list.\n";

	print "\t-idoutgroup [%]\n";
	print "\t\tMaximum percentage of identity that is allowed for a match in the out-group database\n";
	print "\t\tas a BLAST result.\n";
	print "\t\tFor the out-group, this script will allow frames that have an identity score lower than the\n";
	print "\t\tscore in -idoutgroup OR a coverage score lower than the percentage given in -coverageoutgroup.\n";

	print "\t-idingroup [%]\n";
	print "\t\tMinimum percentage of identity that is needed for a match in an in-group genome\n";
	print "\t\tas a BLAST result.\n";
	print "\t\tFor the in-group, this script will allow frames that have an identity score higher than the\n";
	print "\t\tscore in -idingroup AND a coverage score higher than the percentage given in -coverageingroup.\n";

	print "\t-coverageoutgroup [%]\n";
	print "\t\tMaximum alignment coverage that is allowed for a match in the out-group database as a\n";
	print "\t\tBLAST result.\n";
	print "\t\tFor the out-group, this script will allow frames that have an identity score lower than the\n";
	print "\t\tscore in -idoutgroup OR a coverage score lower than the percentage given in -coverageoutgroup.\n";

	print "\t-coverageingroup [%]\n";
	print "\t\tMinimum alignment coverage that is needed for a match in the in-group genomes\n";
	print "\t\tas a BLAST result.\n";
	print "\t\tFor the in-group this script will allow frames that have an identity score higher than the\n";
	print "\t\tscore in -idingroup AND a coverage score higher than the percenage given in -coverageingroup.\n";

	print "\t-help\n";
	print "\t\tPrint this help text.\n";

	print "\nDEPENDENCIES\n";
	print "\tBLASTn and BioPerl should be installed on your system. BLASTn should be added to \$PATH.\n";

	print "\nAUTHOR\n";
	print "\tThis script was written by Koen M. Verstappen at Utrecht University, the Netherlands (2016).\n";

	exit;
}

## Copyright 2017 Utrecht University

## Licensed under the Apache License, Version 2.0 (the "License");
## you may not use this file except in compliance with the License.
## You may obtain a copy of the License at

##     http://www.apache.org/licenses/LICENSE-2.0

## Unless required by applicable law or agreed to in writing, software
## distributed under the License is distributed on an "AS IS" BASIS,
## WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
## See the License for the specific language governing permissions and
## limitations under the License.
