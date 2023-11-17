#!/usr/bin/perl
# Written: Nick Kent, 12th Nov 2012
# Modified: Nick Kent, 12th Nov 2012
# Modified: Nick Kent, 22nd Feb 2019
# Modified: Nick Kent, 6th April 2022
# Tidied: Nick Kent, 17th Nov 2023
# USAGE:- perl Rossi_bed_chr_changer_v22_xxx.plx
#
# This script takes chr1, chr2 format .bed files from the Rossi et al., 2021 GitHub
# in its home dir and converts to chrI, chrII format for the SacCer3 build of the yeast genome.
# Output will be in a new time-stamped directory.
# 
################################################################################

use strict;
use warnings;
use Cwd;
use POSIX qw(strftime);


################################################################################
# RESET THE VARIABLES BELOW IF REQUIRED
# $bed_indir_path   - The directory containing the .bed file
# $outdir_path  - The directory to store the flipped .bed output file
# $bed_outfile_name - User-defined name for the .bed output file
# %chr_pairs - a hash containing the chr conversion criteria - check it's right!!!
#
################################################################################
my $datestring =strftime "%F_%H-%M-%S", localtime;
print "\n\n*********************************************";
print "\n*********** Time to change your bed?!********";
print "\n*********************************************\n";
print "Started:\t".$datestring."\n";
my $output_dir = $datestring."_bed_change_out";

unless (-e cwd."/$output_dir" or mkdir cwd."/$output_dir"){
		die "Unable to create $output_dir\n"
};

my $bed_indir_path =cwd;
my $outdir_path =cwd."/$output_dir";
my %chr_pairs = (
'chr1'=>'chrI',
'chr2'=>'chrII',
'chr3'=>'chrIII',
'chr4'=>'chrIV',
'chr5'=>'chrV',
'chr6'=>'chrVI',
'chr7'=>'chrVII',
'chr8'=>'chrVIII',
'chr9'=>'chrIX',
'chr10'=>'chrX',
'chr11'=>'chrXI',
'chr12'=>'chrXII',
'chr13'=>'chrXIII',
'chr14'=>'chrXIV',
'chr15'=>'chrXV',
'chr16'=>'chrXVI',
'chr17'=>'chrM',
'chr18'=>'2micron');

################################################################################
################################################################################
# MAIN PROGRAM
################################################################################
################################################################################

# define some variables 


my $infile_bed;
my $bed_outfile;
my @line_bed;
my @files_bed;
my $bed_size;


################################################################################
# Read in the .bed file values to four enormous arrays
################################################################################


opendir(DIR,$bed_indir_path) || die "Unable to access file at: $bed_indir_path $!\n";

@files_bed = readdir(DIR);

# process the input file within bed_indir_path
foreach $infile_bed (@files_bed){    

    # ignore hidden files and only get those ending .bed
    if (($infile_bed !~ /^\.+/) && ($infile_bed =~ /.*\.bed/)){
        
       # define outfile name from infile name
        $bed_outfile = substr($infile_bed,0,-4)."_SacCer3";
        $bed_outfile .= '.bed';

print "Found, and processing, $infile_bed \n";

open(IN, "$bed_indir_path/$infile_bed")
            || die "Unable to open $infile_bed: $!";
        
        # define  new arrays to store the .bed values from infile
        my @bed_chr;
        my @bed_pos_start;
        my @bed_pos_end;
		my @bed_id;
        my @bed_score;
		my @bed_dot;
        
        # loop through infile to get values
        while(<IN>){

            chomp;

            # split line by delimiter and store elements in an array
            @line_bed = split('\t',$_);

            # store the columns we want in the  new arrays and rename the chr column
			# this is cumbersome - should have done something with subst $_ - in another life
			

            push(@bed_chr,$chr_pairs{$line_bed[0]});         
            push(@bed_pos_start,$line_bed[1]);
			push(@bed_pos_end,($line_bed[2]));
			push(@bed_id,($line_bed[3]));
			push(@bed_score,($line_bed[4]));
			push(@bed_dot,($line_bed[5]));
        }
        
        # close in file handle
        close(IN);
	

# store size of bin array
        $bed_size = @bed_pos_start;

print "Contains: $bed_size bin values\n";


######################################################################################
# The output files
######################################################################################


# try and open the .bed output file
        open(OUT,"> $outdir_path/$bed_outfile")
             || die "Unable to open $bed_outfile: $!";
        
print "Have just created $bed_outfile\n";

# a counter variables
my $count = 0; # Counter for each line ID

until ($count == $bed_size){ #until 1

print(OUT 
		$bed_chr[$count]."\t".
		$bed_pos_start[$count]."\t".
		$bed_pos_end[$count]."\t".
		$bed_id[$count]."\t".
		$bed_score[$count]."\t".
		$bed_dot[$count]."\n");
		
		$count++;

} #until 1 closer


 # close .bed out file handle
        close(OUT);


}}
