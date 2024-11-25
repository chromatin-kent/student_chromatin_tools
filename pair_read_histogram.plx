#!/usr/bin/perl
# Written by Nick Kent, Aug 2011
# Fiddled with: Nick Kent, Oct 2016
# Fiddled with: Nick Kent, Nov 2018 - removed stupid ID flag logic
# Fiddled with: Nick Kent, Oct 2019 - use Cwd
# Fiddled with: Nick Kent, Apr 2021 - tidy commenting
#
# USAGE:- perl pair_read_histogram.plx
#
# This script takes takes a paired read SAM format file and produces a frequency 
# distribution of the ISIZE values in a user specificed range and resolution for 
# plotting as a histogram next to an EtBr-stained gel or whatever.
#
# You could apply this to a whole .sam output file from Bowtie, or to a single .txt
# output file from chrgrep.sh. You will have to specify the expected file-ending
# at line 63.
#
# NOTE: the script will process up to the first 10,000,000 reads (unless you change
# the value at line 88. Going higher than this will hog masses of memory and make your
# coffee go cold.
#
################################################################################
use strict;
use warnings;
use Math::Round;
use Cwd;
################################################################################
# SET THE VARIABLES BELOW AS REQUIRED
# $indir_path   - The directory containing the .SAM file(s) to be processed
# $outdir_path  - The directory to store the output file(s)
# $top		- The upper ISIZE boundary
# $bin_width    - The histogram bin width value in bp
################################################################################

my $inA_indir_path =cwd."/in";
my $outdir_path =cwd."/out";
my $top = 1000;
my $bin_width = 10;

################################################################################
# MAIN PROGRAM - you should not need to edit below this line
################################################################################

# define some variables

my $infile_A;
my @line_A;
my @files_A;

################################################################################
# Count SAM file IDs 
################################################################################

# store input file name in an array
opendir(DIR,$inA_indir_path) || die "Unable to access file at: $inA_indir_path $!\n";

@files_A = readdir(DIR);

# process the input file within indir_path
foreach $infile_A (@files_A){    

    # ignore hidden files and only get those with the correct ending
    if (($infile_A !~ /^\.+/) && ($infile_A =~ /.*\.sam/)){
        
        
print "Found, and processing, $infile_A \n";

open(IN, "$inA_indir_path/$infile_A")
            || die "Unable to open $infile_A: $!";
        
        # define arrays and variables to store required values from infile
       
		my $id_counter = 0; # the SAM ID counter
		my @isize = 0; # the array of ISIZE (end-to-end distance) values

	
	
	# loop through infile to get values
        while(<IN>){
           
           
	    chomp;

            # split line by delimiter and store elements in an array
            @line_A = split('\t',$_);

            # test for number of lines and store the values we want
	    if ($id_counter >= 10000000){
	    last;
	    }
	    elsif ($line_A[0] !~ "@" && $line_A[8] > 0){
	    $id_counter ++;
	    push (@isize, $line_A[8]);
	    }

           print "Current count is $id_counter\r"; 
        }

	# close in file handle
        close(IN);
        
       


print "\n\nProcessed: $id_counter aligned reads\n";



#################################################################################
# Tally the histogram and produce output file
#################################################################################

        # define outfile name from infile name
       my $outfile = substr($infile_A,0,-4)."_pair_histogram";
        $outfile .= '.txt';
        
		
		
		# Define the number of bins 
		my$bin_no = (int($top/$bin_width))+1;
		
		# Define the distribution frequency array
		my @dist_freq;
		my $i=0;
		
		# Fill the frequency distribution "bins" with 1's. This prevents div 0 errors.
		for ($i=0; $i<$bin_no; $i++){
			push (@dist_freq, 1);
			}
			
		# Reset the incrementor and define the bin hit variable
		$i=0;
		my $bin_hit = 0;
		
		# The tally counter 
		while ($i < $id_counter){
			$bin_hit = int($isize[$i]/$bin_width);
			$dist_freq[$bin_hit] ++;
			$i ++;
			}
			
		
        # try and open output file
        open(OUT,"> $outdir_path/$outfile")
             || die "Unable to open $outfile: $!";
        
        
            # print required data to output file
	    print(OUT "Input file ".$infile_A." contained ".$id_counter." aligned reads with ISIZE greater than 0\n\n");
	    print(OUT "Bin:"."\t"."Freq:"."\t"."log Freq:"."\n");

			for ($i=0; $i<$bin_no; $i++){
			
            print(OUT ($i*$bin_width)."\t".$dist_freq[$i]."\t".(log($dist_freq[$i])/log(10))."\n");
			}
            
        }
        # close out file handle
        close(OUT);
    
}
