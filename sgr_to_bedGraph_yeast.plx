#!/usr/bin/perl
# Written: Nick Kent, 17th July 2013
# Last updated: Nick Kent, 29th June 2016
# Modified: Nick Kent, 20th Dec 2018 - changed directory finder

# USAGE:- perl sgr_to_bedGraph.plx
#

# 
################################################################################

use strict;
use warnings;
use Cwd;

################################################################################
# SET THE VARIABLES BELOW AS REQUIRED
# 
# $bin_size - this is the input bin size in the .sgr - usually 10 (bp).
# $end_pos - where you want the end of your chr - usually 0.
#
################################################################################

my $bin_size = 10;
my $end_pos = 0;

################################################################################
################################################################################
# MAIN PROGRAM
################################################################################
################################################################################

# define some variables
my $sgr_indir_path =cwd."/in";
my $outdir_path =cwd."/out";
my $infile_sgr;
my $sgr_outfile;
my @line_sgr;
my @files_sgr;
my $sgr_size;


################################################################################
# Read in the .sgr file values to  enormous arrays
################################################################################


opendir(DIR,$sgr_indir_path) || die "Unable to access file at: $sgr_indir_path $!\n";

@files_sgr = readdir(DIR);

# process the input file within sgr_indir_path
foreach $infile_sgr (@files_sgr){    

    # ignore hidden files and only get those ending .sgr
    if (($infile_sgr !~ /^\.+/) && ($infile_sgr =~ /.*\.sgr/)){
        
       # define outfile name from infile name
        $sgr_outfile = substr($infile_sgr,0,-4);
        $sgr_outfile .= '.bedGraph';
       # $sgr_outfile .= '.bed';
        
print "Found, and processing, $infile_sgr \n";

open(IN, "$sgr_indir_path/$infile_sgr")
            || die "Unable to open $infile_sgr: $!";
        
        # define three new arrays to store the .sgr values from infile
        my @sgr_chr;
        my @sgr_bin;
        my @sgr_freq;
        
        # loop through infile to get values
        while(<IN>){

            chomp;

            # split line by delimiter and store elements in an array
            @line_sgr = split('\t',$_);

            # store the columns we want in the three new arrays
			
            push(@sgr_chr,$line_sgr[0]);
            push(@sgr_bin,$line_sgr[1]);
            push(@sgr_freq,($line_sgr[2]));
        }
        
        # close in file handle
        close(IN);
	

# store size of bin array
        $sgr_size = @sgr_freq;

print "Contains a whopping: $sgr_size bin values\n";


######################################################################################
# The output file
######################################################################################


# try and open the .sgr output file
        open(OUT,"> $outdir_path/$sgr_outfile")
             || die "Unable to open $sgr_outfile: $!";
        
print "Have just created $sgr_outfile\n";
print (OUT
      "track type=bedGraph name=".$sgr_outfile."\n");

# a counter variables
my $count = 0; # Counter for each line ID

until ($count == $sgr_size){ #until 1

$end_pos = $sgr_bin[$count]+$bin_size;

print(OUT 
		$sgr_chr[$count]."\t".
		$sgr_bin[$count]."\t".
		$end_pos."\t".
		$sgr_freq[$count]."\n");
		
		$count++;

} #until 1 closer


 # close .sgr out file handle
        close(OUT);


}}
