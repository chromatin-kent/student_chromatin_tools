#!/usr/bin/perl
# Written: Nick Kent 11th May 2022 to process .bed and .int files from Rossi et al data
# Modified: Nick Kent 22nd Jan 2024 to be more logical to use
# Modified : Nick Kent 14th Nov 2024 - changed Posix date stamp to perl localtime() so that script will run on both
# linux and Windows machines.
#
# USAGE:- perl Rossi_interval_strand_caller_vxxx.plx
#
# Rossi et al. (2021) ChEx peak locations have no strandedness unlike yeast genes and other 
# features. If your Rossi .bed file contains locations that may be expected to bind 
# and regulate genes via a promoter-proximal sequence or within an ORF you might what to give those locations
# the same strandedness as the genes they are associated with. This script will help you do that.
#
# This script is a modifcation of interval_comparison.pl. It takes an interval file listing
# TSS minus some bp, introns, genic/ORF regions or whatever and asks whether a ChEx peak location in a
# Rossi et al (2021) .bed file occurs within any of those intervals. If the answer is "yes", then the
# peak is recoded with the strand identity (F or R) associated with the matching interval. 
#
# The script takes two input files in tab-delimited format (place these in the same directory as the script):
#
# (i) Rossi .bed file with six columns, the last being the useless "dot".
# (ii) An .int file containg a list of intervals - you'll have to ask Nick for any custom files. These 
# intervals have a strand designation and a gene ID. You must make sure the chr IDs match between the two 
# input files.
#
#
# The Script outputs a new .bed file in a time-stamped output folder:
# This new .bed file is almost identical to the original but now contains the gene name relating to
# the matching interval and a strand designation.
#
# The new stranded.bed file can them be used as input to RossiStrandedbed_vs_PartNsgr_to_CFD_vxxx.plx
# 
# The script also outputs at command line the total number of sites which found a matching interval. 
#
################################################################################

use strict;
use warnings;
use Cwd;



################################################################################
################################################################################
# MAIN PROGRAM
################################################################################
################################################################################
my $datestring = localtime();
$datestring =~ tr/ /_/;
$datestring =~ tr/:/-/;
print "\n\n*********************************************";
print "\n***** Add strand direction to Rossi .bed ******";
print "\n*********************************************\n";
print "Started:\t".$datestring."\n";

# define some variables

my $output_dir = $datestring."_stranded_bed_out";

unless (-e cwd."/$output_dir" or mkdir cwd."/$output_dir"){
		die "Unable to create $output_dir\n"
};

my $interval_indir_path =cwd;
my $peak_indir_path =cwd;
my $outdir_path =cwd."/$output_dir";
my $no = 0;
my $interval_outfile_name; 
my $infile_interval;
my $infile_peak;
my $outfile; 
my @line_interval;
my @line_peak;
my @files_interval;
my @files_peak;
my $interval_size;
my $peak_size;

################################################################################
# Read in the interval file values to five arrays
################################################################################

# store input file name in an array
opendir(DIR,$interval_indir_path) || die "Unable to access file at: $interval_indir_path $!\n";

@files_interval = readdir(DIR);

# process the input file within indir_path
foreach $infile_interval (@files_interval){    

    # ignore hidden files and only get those ending .int
    if (($infile_interval !~ /^\.+/) && ($infile_interval =~ /.*\.int/)){
        
        
print "Found, and processing, $infile_interval \n";

open(IN, "$interval_indir_path/$infile_interval")
            || die "Unable to open $infile_interval: $!";
        
        # define the arrays to store required values from infile
        my @interval_chr; # the interval chrn
		my @interval_id; # an informative interval gene name
        my @interval_start; # the start position of the interval
		my @interval_end; # the end position of the interval
		my @interval_strand; # Strand associated with interval
	
	
	
	# loop through infile to get values
        while(<IN>){
           
           
	    chomp;

            # split line by delimiter and store elements in an array
            @line_interval = split('\t',$_);

            # store the columns we want in five new arrays
            
	    push(@interval_id,$line_interval[1]);
	    push(@interval_strand,$line_interval[4]);
	    push(@interval_chr,$line_interval[0]);
	    push(@interval_start,$line_interval[2]);
	    push(@interval_end,$line_interval[3]);
	    
            
        }

	# close in file handle
        close(IN);
	closedir(DIR);
        
        # store size of bin array
        $interval_size = @interval_chr;


print "Contains: $interval_size intervals\n";


################################################################################
# Read in the peak file values to four arrays
################################################################################


opendir(DIR,$peak_indir_path) || die "Unable to access file at: $peak_indir_path $!\n";

@files_peak = readdir(DIR);

# process the input file within indir_path
foreach $infile_peak (@files_peak){    

    # ignore hidden files and only get those ending .sgr
    if (($infile_peak !~ /^\.+/) && ($infile_peak =~ /.*\.bed/)){
        
       

print "Found, and processing, $infile_peak \n";

open(IN, "$peak_indir_path/$infile_peak")
            || die "Unable to open $infile_peak: $!";
        
        # define three new arrays to store required values from the Rossi .bed infile
       
        my @peak_chr; # chr ID
        my @peak_start; # Start of ChEx peak
		my @peak_score; # Score for ChEx peak

        
        # loop through infile to get values
        while(<IN>){

            chomp;

            # split line by delimiter and store elements in an array
            @line_peak = split('\t',$_);

            # store the columns we want in the four new arrays
           
			push(@peak_chr,$line_peak[0]);
            push(@peak_start,$line_peak[1]);
            push(@peak_score,$line_peak[4]);
        }
        
        # close in file handle
        close(IN);
	closedir(DIR);

# store size of bin array
        $peak_size = @peak_chr;

print "Contains: $peak_size ChEx peak values\n";


######################################################################################
# The .bed output file
######################################################################################

# define outfile name and set ending to .bed
		$interval_outfile_name = substr($infile_peak,0,-4)."_nearestORF_stranded";
        $outfile = $interval_outfile_name;
        $outfile .= '.bed';

# try and open output file
        open(OUT,"> $outdir_path/$outfile")
             || die "Unable to open $outfile: $!";
        
print "Have just created $outfile\n";
        
# some counter variables
my $interval_count = 0; # interval site number
my $peak_count = 0; # Peak site number
my $hit_count = 0; # A counter for positive site/peak hits
my $match = 0; # running score of total peaks within this interval


# a very inefficient way of searching for matches - but it works
until ($interval_count == $interval_size){

	until ($peak_count == $peak_size){ 


	# this increments the peak counter if the chrn values do not match
	if ($interval_chr[$interval_count] ne $peak_chr[$peak_count]){

		$peak_count++;

        # this tests whether or not a peak bin occurs within an interval 

	}elsif ($peak_start[$peak_count]>=($interval_start[$interval_count]) && 
			$peak_start[$peak_count]<=($interval_end[$interval_count]) && 
			$interval_chr[$interval_count] eq $peak_chr[$peak_count]){

		# this writes out the 6 column tab-delimited .bed
		print(OUT 
		
		$peak_chr[$peak_count]."\t".
		$peak_start[$peak_count]."\t".
		($peak_start[$peak_count]+1)."\t".
		$interval_id[$interval_count]."\t".
		$peak_score[$peak_count]."\t".
		$interval_strand[$interval_count]."\n");
		
		$hit_count++;
		$peak_count++;

	}else{
		$peak_count ++;
	}
}
	
	if ($hit_count==0){

		$interval_count++;
		$peak_count = 0;

	}else{
		$match += $hit_count;
		$interval_count++;
		$hit_count = 0;
		$peak_count = 0;
	}


}
	
       
        # close out file handle
        close(OUT);

	
      

	print "Out of $peak_size peaks, $match occured within the specified intervals\n";

}}}}
       
