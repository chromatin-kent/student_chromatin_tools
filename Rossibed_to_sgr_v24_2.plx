#!/usr/bin/perl
# Written: Nick Kent, 28th Jan 2024
# Last updated:
# Adapted from Interval_speedy_wrangler.pl: Nick Kent 23rd Jan 2024
# Modified: Nick Kent 27th Nov 2024 - to use perl Localtime() rather than Posix

#
# This script converts Rossi et al., .bed files into .sgr format for input into a "...to_CFD" plotter.
# It takes an interval file templated by a SacCer3 whole genome .sgr file and uses this to convert
# Rossi et al .bed file peaks into the .sgr format via a wildly inefficient brute force interval comparison. 
#
# The script takes two types of input files in tab-delimited format:
#
# 1.The interval file can be any SacCer3 .sgr file, but given an .int file end.
# 
# 2.Any number of standard Rossi et al., SacCer3.bed files describing the factors you are trying to match to the "target"/
#
# The Script .sgr format file listing the matches where Rossi peaks fall within .sgr file bins. The script also
# outputs at command line the total number of peaks which found a matching interval (effectively calling bins
# with more than one peak as the number greater than the input value). 
#
#
# WARNING: this is relatively untested and was written on the fly in about 30 mins; it is appallingly inefficient.
################################################################################

use strict;
use warnings;
use Cwd;
use List::Util;


###############################################################################
# SET VARIABLES HERE AS REQUIRED
#
# $bin_window should match the bin width of the input .sgr format .int file
###############################################################################

my $bin_width = 10;

################################################################################
################################################################################
# MAIN PROGRAM
################################################################################
################################################################################

my $datestring = localtime();
$datestring =~ tr/ /_/;
$datestring =~ tr/:/-/;
print "\n\n*********************************************";
print "\n********** Bed to sgr? You fool! ***********";
print "\n*********************************************\n";
print "Started:\t".$datestring."\n";

#temp new stuff
# define some variables

my $output_dir = $datestring."_bed2sgr_out";

unless (-e cwd."/$output_dir" or mkdir cwd."/$output_dir"){
		die "Unable to create $output_dir\n"
};

my $interval_indir_path = cwd;
my $peak_indir_path = cwd;
my $outdir_path =cwd."/$output_dir";
my $infile_interval;
my $infile_peak;
my $outfile; 
my @line_interval;
my @line_peak;
my @files_interval;
my @files_peak;
my $interval_size;
my $peak_size;
my $descriptor;
my $line_count;



################################################################################
# Read in the interval file values to arrays
################################################################################

# store input file name in an array
opendir(DIR,$interval_indir_path) || die "Unable to access file at: $interval_indir_path $!\n";

@files_interval = readdir(DIR);

# process the input file within indir_path
foreach $infile_interval (@files_interval){    

    # ignore hidden files and only get those ending .int
    if (($infile_interval !~ /^\.+/) && ($infile_interval =~ /.*\.int/)){
     $descriptor = substr($infile_interval,0, -4);   
        
print "Found, and processing, $infile_interval \n";

open(IN, "$interval_indir_path/$infile_interval")
            || die "Unable to open $infile_interval: $!";
        
        # define the arrays to store required values from infile
    my @interval_chr; # the interval chrn
    my @interval_start; # the start position of the interval
	my @interval_end; # the end position of the interval
	my @interval_value; # default is zero
	
	
	
	# loop through infile to get values
        while(<IN>){
           
           
	    chomp;

            # split line by delimiter and store elements in an array
            @line_interval = split('\t',$_);

            # store the columns we want in four new arrays
            
	    push(@interval_chr,$line_interval[0]);
	    push(@interval_start,$line_interval[1]);
	    push(@interval_end,$line_interval[1]+$bin_width);
	    push(@interval_value,0);	    
            
        }

	# close in file handle
        close(IN);
	closedir(DIR);
        
        # store size of bin array
        $interval_size = @interval_chr;


print "Contains: $interval_size intervals\n";


opendir(DIR,$peak_indir_path) || die "Unable to access file at: $peak_indir_path $!\n";

@files_peak = readdir(DIR);

# process the input file within indir_path
foreach $infile_peak (@files_peak){    

    # ignore hidden files and only get those ending .bed
    if (($infile_peak !~ /^\.+/) && ($infile_peak =~ /.*\.bed/)){
        
       

print "Found, and processing, $infile_peak \n";

open(IN, "$peak_indir_path/$infile_peak")
            || die "Unable to open $infile_peak: $!";
        
        # define three new arrays to store required values from .sgr infile
       
        my @peak_chr; # the .bed file chrn
        my @peak_pos; # the ChIP exo peak position in bp
        
        # loop through infile to get values
        while(<IN>){

            chomp;

            # split line by delimiter and store elements in an array
            @line_peak = split('\t',$_);

            # store the columns we want in the four new arrays
           
	    push(@peak_chr,$line_peak[0]);
        push(@peak_pos,$line_peak[1]);
        
        }
        
        # close in file handle
        close(IN);


# store size of bin array
        $peak_size = @peak_chr;

print "Contains: $peak_size peak values\n";


######################################################################################
# The output file
######################################################################################

        
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

	}elsif ($peak_pos[$peak_count]>=($interval_start[$interval_count]) && 
			$peak_pos[$peak_count]<=($interval_end[$interval_count]) && 
			$interval_chr[$interval_count] eq $peak_chr[$peak_count]){

			
		# Increment the frequency value at this [$interval_count]
		$interval_value[$interval_count]++;
		
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
	

	print "Out of $peak_size peaks, $match occured within the specified intervals\n";     
     	#define outfile name and set ending to .sgr
        $outfile = $infile_peak;
        $outfile .= '.sgr';
		
# try and open output file
        open(OUT,"> $outdir_path/$outfile")
             || die "Unable to open $outfile: $!";
        
print "Have just created $outfile\n"; 
$line_count=0;

	while($line_count<$interval_size){

	  # print out .sgr format
		print(OUT 
		$interval_chr[$line_count]."\t".
		$interval_start[$line_count]."\t".
		$interval_value[$line_count]."\n");
		
		$line_count++;
	}
	
	#clear the interval values back to zero
	$_ =0 for @interval_value;
	
# close out file handle
        close(OUT);
 

}

}
 
		closedir(DIR); 
}
 
}
      
