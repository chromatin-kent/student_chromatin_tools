#!/usr/bin/perl
# Written: Nick Kent, 10th Dec 2010
# Last updated: Nick Kent, 13 Aug 2015
# Last updated: Nick Kent, 19 Jun 2019 to use Cwd
# Last updated: Nick Kent, 13th Jan 2026 to work for FYP2026 data
# Last updated: Nick Kent, 19th Jan 2026 to work for FYP2026 .bed files
# USAGE:- perl Stranded_bed_vs_Anysgr_2_CFD_Compare_Stats_vxxx.plx
#
# This script takes one or several .site files containing lists of sites (these
# could be TSSs or TF sites or whatever you want) and compares them with a series
# of whole-genome, Partn .sgr files representing a WT/mutant comparison - referred
# to as "A.sgr" and "B.sgr". It will take one size class or a spectrum of Partn sizes.
#
# The input .bed file should have six columns: chrn; Site start pos;Site end pos; ID string;Score;strand.
#
# The script matches up A.sgr and B.sgr files on the basis of the Partn value present
# in the filename. It then performs a cumulative frequency calculation
# over each bin in the specified window surrounding each site feature. All the cumulative
# frequencies are then averaged over the entire window to generate  normalization
# factors for each paired dataset.
#
# The script then revisits the A and B .sgr data, scales each individual frequency value
# according to the respective nomalization factor (above) and pushes the
# individual values from each bin into an array. The arrays for the A and B .sgrs are
# then fed into the Statistics::DependantTTest, Statistics::Descriptive, Statistics::
# Test::WilcoxonRank and Statistics::RankCorrelation modules
# to output a range of values including p-values for paired T and Wilcoxon Mann-Whitney
# tests between the two scaled/normalized samples at each bin. 
# 
# In testing, the script is clearly working well to superimpose both datasets, and the
# T-test p-values seem to visually track obvious differences in the profiles. Note that
# a Student's T test is probably not valid for this type of data because the distribution
# of values is non-gaussian/normal and skewed. The non-parametric WMW is valid but weak.
# Remember a low p-value simply implies that the medians of the two datasets come from different
# distributions - the result says nothing specific about chomatin structure per se. Also remember
# that the WMW p-values will be nonsensical in low read number regions because the median
# values will moslty be zeroes. So although the idea of creating a "p-value landscape" seems
# wonderful it might, realistically, be un-obtainable with this data. Finally, remember that you
# should consider multiple testing correction. One approach here, is to use the Chunker.plx script
# to randomly subsample your .site lists to bring n values down. A second approach here is to apply
# the Bonferroni correction - e.g. if the number of sites n = 500, and you are interested
# in significance at the 0.05 level, then your calculated p-value should be <=0.05/500 <=0.0001 (<=alpha/n).
#
# The calculation of Spearman's tau and Kendall's rho values is highly compute intensive. Running
# this script with several thousand sites (e.g. TSS) might take several hours.
# 
################################################################################

use strict;
use warnings;
use Cwd;
use Statistics::DependantTTest;
use Statistics::Descriptive;
use Statistics::Distributions;
use Statistics::Test::WilcoxonRankSum;
use Statistics::RankCorrelation;


################################################################################
# SET THE VARIABLES BELOW AS REQUIRED
# $siteID_indir_path - The directory containing the site list file
#
# $sgr_A_indir_path   - The directory containing the full genome Partn.sgr file A
#
# $sgr_B_indir_path   - The directory containing the full genome Partn.sgr file B
# IMPORTANT: Both A and B .sgr files MUST be identical in bin numbers for each chr
#
# $outdir_path  - The directory to store the .sgr and .cfd output files
#
# $descriptor - A string which defines the outfile name.
#
# @particles - An array which contains a comma delimited list of the particle sizes you want
# the script to look for (you can put the full spectrum of size files in and the script will
# output for all f them allowing you to plot not only surface landscapes of the individual
# data sets, but also p-value landscapes.
#
# $bin_window - For 5bp binned data this is set to 120;the script will
# extract read frequency values from bins -600bp (120*5bp) to +600bp surrounding
# the site (potentially four nucleosomes-worth each side).
# 
# $bin_size - For 5bp bins this is set to 5. Alter the increment counters at Lines
# 487 to 489 to 3 to match 15bp bin output or to 1 in you want every bin in the CFD
#
# $scale_factor - Applied to .sgr B to correct for differences in read depth with 
# .sgr B
#
# $output_scale - Only relevant if using 5bp binned data.
################################################################################

my $siteID_indir_path =cwd."/Sites_in";
my $sgr_A_indir_path =cwd."/A_in";
my $sgr_B_indir_path =cwd."/B_in";
my $outdir_path =cwd."/Stats_out";
my $descriptor="WT_vs_isw1_CTDSer2";
my @particles = ("Part450"); #comma delimit
my $particle_no = @particles;
my $bin_window = 120;
my $bin_size = 10;
my $output_scale =1;

################################################################################
################################################################################
# MAIN PROGRAM
################################################################################
################################################################################

# define some variables

my $cwd = getcwd;
my $infile_sgr_A;
my $infile_sgr_B;
my $infile_siteID;
my $stats_outfile;
my @line_siteID; 
my @line_sgr_A;
my @line_sgr_B;
my @files_sgr_A;
my @files_sgr_B;
my @files_siteID;
my $sgr_A_size;
my $sgr_B_size;
my $siteID_size;
my $sgrA_file_count;
my $sgrB_file_count;
my $A_counter = 0;
my $B_counter = 0;
my $part_counter = 0;
my $F_siteID_size;
my $R_siteID_size;
my %bin_map; 
my $chr_count;

################################################################################
# Read the sgrA and B filenames into an array
################################################################################
opendir(DIR,$sgr_A_indir_path) || die "Unable to access file at: $sgr_A_indir_path $!\n";

@files_sgr_A = readdir(DIR);
$sgrA_file_count = @files_sgr_A;

print "There are $sgrA_file_count files in the SgrA directory. The .sgr files are called:\n";
for ($A_counter=0;$A_counter< $sgrA_file_count; $A_counter++){
  unless (($files_sgr_A[$A_counter] =~ /^\.+/) && ($files_sgr_A[$A_counter] !~ /.*\.sgr/)){
  print "$files_sgr_A[$A_counter]\n\n"
  }}
  closedir(DIR);

opendir(DIR,$sgr_B_indir_path) || die "Unable to access file at: $sgr_B_indir_path $!\n";

@files_sgr_B = readdir(DIR);
$sgrB_file_count = @files_sgr_B;

print "There are $sgrB_file_count files in the SgrB directory. The .sgr files are called:\n";
for ($B_counter=0;$B_counter< $sgrB_file_count; $B_counter++){
  unless (($files_sgr_B[$B_counter] =~ /^\.+/) && ($files_sgr_B[$B_counter] !~ /.*\.sgr/)){
  print "$files_sgr_B[$B_counter]\n\n"
  }}
  closedir(DIR);

	

################################################################################
# Get ded list and write to an array - from .bed format with four columns:
# chrn;siteID; site position; F/R 
################################################################################

# store input file name in an array
opendir(DIR,$siteID_indir_path) || die "Unable to access file at: $siteID_indir_path $!\n";

@files_siteID = readdir(DIR);

# process the input file within siteID_indir_path
foreach $infile_siteID (@files_siteID){    

    # ignore hidden files and only get those ending .bed
    if (($infile_siteID !~ /^\.+/) && ($infile_siteID =~ /.*\.bed/)){
        
        
print "Found, and processing, $infile_siteID \n\n";

open(IN, "$siteID_indir_path/$infile_siteID")
            || die "Unable to open $infile_siteID: $!";
        
        # define arrays to store site chromosome no., position and strand 
        my @site_chr;
		my @site_pos;
		my @site_strand;
		my @F_site_chr;
		my @F_site_pos;
		my @R_site_chr;
		my @R_site_pos;
	
	# loop through infile to get values
        while(<IN>){

	    chomp;

	    # split line by delimiter and store elements in an array
            @line_siteID = split('\t',$_);

            # store the required chrn, positions and strands in the three
			# arrays
			
			if($line_siteID[5] =~ "F") { # infile if 1
			
				push(@F_site_chr,$line_siteID[0]);
				push(@F_site_pos,$line_siteID[1]);
				push(@site_chr,$line_siteID[0]);
				push(@site_pos,$line_siteID[1]);
				push(@site_strand,$line_siteID[5]);
			
			}elsif ($line_siteID[5] =~ "R"){ # infile elsif 1
			
				push(@R_site_chr,$line_siteID[0]);
				push(@R_site_pos,$line_siteID[1]);
				push(@site_chr,$line_siteID[0]);
				push(@site_pos,$line_siteID[1]);
				push(@site_strand,$line_siteID[5]);
			}else{
			
				print "Failed to match strand at $line_siteID[0], $line_siteID[3]\n";
			
			} #infile if 1 closer

        }

	# close in file handle
        close(IN);
	closedir(DIR);
        
        # store sizes of the arrays
        $siteID_size = @site_pos;
        $F_siteID_size = @F_site_pos;
	$R_siteID_size = @R_site_pos;


print "Contains at total of: $siteID_size site IDs\n";
print "consisting of: $F_siteID_size forward strand site IDs; $R_siteID_size reverse
 strand site IDs\n\n";



################################################################################
# Read in the .sgr A file values to three arrays
################################################################################


opendir(DIR,$sgr_A_indir_path) || die "Unable to access file at: $sgr_A_indir_path $!\n";


# take each particle in the list and try to match it with an SgrA file
for ($part_counter=0; $part_counter<$particle_no; $part_counter++){ 
print "The current size range is $particles[$part_counter]\n\n";  

$A_counter = 0;
until ($files_sgr_A[$A_counter] =~ $particles[$part_counter]){
	  
	  $A_counter++;
	  if ($A_counter > 10){
		print "Crikey! It's all gone pear-shaped trying to match a Part size to SgrA\n";
		exit;}
	  }

print "Found, and processing, $files_sgr_A[$A_counter] \n";

open(IN, "$sgr_A_indir_path/$files_sgr_A[$A_counter]")
            || die "Unable to open $files_sgr_A[$A_counter]: $!";
        
        # define three new arrays to store the .sgr values from infile
        my @sgr_A_chr;
        my @sgr_A_bin;
	my @sgr_A_freq;
        
        # loop through infile to get values
        while(<IN>){

            chomp;

            # split line by delimiter and store elements in an array
            @line_sgr_A = split('\t',$_);

            # store the columns we want in the three new arrays
            push(@sgr_A_chr,$line_sgr_A[0]);
            push(@sgr_A_bin,$line_sgr_A[1]);
	    push(@sgr_A_freq,$line_sgr_A[2]);
        }
        
        # close in file handle
        close(IN);
		

# store size of bin array
        $sgr_A_size = @sgr_A_freq;

print "which contains a whopping: $sgr_A_size bin values\n\n";


#######################################################################################
# BUILD THE BIN MAP - using the Sgr A file
#######################################################################################

my $map_count = 0; # a counter variable

# Set bottom
$bin_map{$sgr_A_chr[$map_count]} = 0;

$map_count ++;

# scan through the @sgr_chr array and mark the bins where each new chromsomome starts

until ($map_count == $sgr_A_size){
  
      if ($sgr_A_chr[$map_count] ne $sgr_A_chr[$map_count-1]){
      
      $bin_map{$sgr_A_chr[$map_count]} = $map_count;
      $map_count ++;
      
      }
      else{
      
	  $map_count ++;
	  
	  }

}
# output the number of chromosome types found as the number of hash keys.
$chr_count = keys %bin_map;



#######################################################################################



#################################################################################
# Use sgrA values to determine a CFD nomalization value
#################################################################################

# define some arrays that will be reset during each iteration
my @A_F_cfd_bin;
my @A_F_cfd_freqsum;
my @A_R_cfd_bin;
my @A_R_cfd_freqsum; 


# The F strand calculations:

# some counter variables
my $site_count = 0; # Counter for each site ID
my $bin_count = 0; # Counter .sgr bin numbers
my $cfd_count = 0; # Counter for the cfd arrays
my $top_limit = 0; # A top limit for $bin_window
my @A_F_out_chr; # F .sgr output array for chr
my @A_F_out_bin; # F .sgr output array for bin pos
my @A_F_out_freq; # F .sgr output array for read freq
my $A_F_out_size = 0; # Size of F .sgr output arrays
my $i=0; # An iterator variable


until ($site_count == $F_siteID_size){ #until 1


# Use %bin_map to jump to correct region of sgr arrays

$bin_count = (int($F_site_pos[$site_count]/$bin_size) + $bin_map{$F_site_chr[$site_count]}) - 3;

# this looks mad, but it allows me to recycle all the code from the last version, and takes up any
# rounding slack which would come from different $bin_size values

# find an .sgr bin which contains the current site
	until ($F_site_chr[$site_count] eq $sgr_A_chr[$bin_count] &&
	       $F_site_pos[$site_count] >= $sgr_A_bin[$bin_count] &&
	       $F_site_pos[$site_count] <  $sgr_A_bin[$bin_count +1]){ #until 2

			$bin_count ++;

	} #until 2 closer

# now that we've found the match, let's write values to the output files

# set the bin_counter BACK $bin_window places and set the $top_limit

$bin_count -= $bin_window;
$top_limit = $bin_count + ($bin_window*2);

# Better test to see if match is close to ends of a chromosome. If so, the reported
# bins and read freqs will be chaemeric - we don't want this so we will ditch such matches

	if($F_site_chr[$site_count] ne $sgr_A_chr[$bin_count] ||
	   $F_site_chr[$site_count] ne $sgr_A_chr[$top_limit]){ #if 1

	print "Can't output forward strand values for $F_site_chr[$site_count] site: $F_site_pos[$site_count]\n";

	     } #if 1 closer
	else { #else 1


# Push the chrn, bin and freq values to the F .sgr arrays and add values to  F cfd freq array

	until ($bin_count == $top_limit+1){ #until 3

		push (@A_F_out_chr,$sgr_A_chr[$bin_count]);
		push (@A_F_out_bin,$sgr_A_bin[$bin_count]);
		push (@A_F_out_freq,$sgr_A_freq[$bin_count]);

		$A_F_cfd_freqsum[$cfd_count] += $sgr_A_freq[$bin_count];

		$bin_count ++;
		$cfd_count ++;


	} #until 3 closer
    } #else 1 closer

$cfd_count = 0;
$bin_count = 0;
$site_count ++;

} #until 1 closer

$A_F_out_size = @A_F_out_chr;


# The R strand calculations:

# reset the counter variables and define some more arrays
$site_count = 0; # Counter for each site ID
$cfd_count = 0; # Counter for the cfd arrays
$bin_count = 0;
my @A_R_out_chr; # R .sgr output array for chr
my @A_R_out_bin; # R .sgr output array for bin pos
my @A_R_out_freq; # R .sgr output array for read freq
my $A_R_out_size = 0; # Size of F .sgr output arrays

until ($site_count == $R_siteID_size){ #until 1


# Use %bin_map to jump to correct region of sgr arrays

$bin_count = (int($R_site_pos[$site_count]/$bin_size) + $bin_map{$R_site_chr[$site_count]}) - 3;


# find an .sgr bin which contains the current site
	until ($R_site_chr[$site_count] eq $sgr_A_chr[$bin_count] &&
	       $R_site_pos[$site_count] >= $sgr_A_bin[$bin_count] &&
	       $R_site_pos[$site_count] <  $sgr_A_bin[$bin_count +1]){ #until 2

			$bin_count ++;

	} #until 2 closer

# now that we've found the match, let's write values to the output files

# set the bin_counter BACK $bin_window places and set the $top_limit

$bin_count -= $bin_window;
$top_limit = $bin_count + ($bin_window*2);

# Better test to see if match is close to ends of a chromosome. If so, the reported
# bins and read freqs will be chaemeric - we don't want this so we will ditch such matches

	if($R_site_chr[$site_count] ne $sgr_A_chr[$bin_count] ||
	   $R_site_chr[$site_count] ne $sgr_A_chr[$top_limit]){ #if 1

	print "Can't output reverse strand values for $R_site_chr[$site_count] site: $R_site_pos[$site_count]\n";

	     } #if 1 closer
	else { #else 1


# Push the chrn, bin and freq values to the R .sgr arrays and add values to  R cfd freq array


	until ($bin_count == $top_limit+1){ #until 3

		push (@A_R_out_chr,$sgr_A_chr[$bin_count]);
		push (@A_R_out_bin,$sgr_A_bin[$bin_count]);
		push (@A_R_out_freq,$sgr_A_freq[$bin_count]);

		$A_R_cfd_freqsum[$cfd_count] += $sgr_A_freq[$bin_count];

		$bin_count ++;
		$cfd_count ++;


	} #until 3 closer
    } #else 1 closer

$cfd_count = 0;
$bin_count = 0;
$site_count ++;

} #until 1 closer

$A_R_out_size = @A_R_out_chr;


# At last, the CFD calculation

# Set counter variables and define new arrays
$bin_count = 0;
$cfd_count = 0;
my $A_cfd_sum = 0; # a sum of sums for normalizing the data
my $A_norm_factor = 0; # calced from $cfd_sum
my $A_R_cfd_count = $bin_window*2;
my @A_FandR_cfd; # array to hold summed F and R strand CFD values
my @A_R_cfd; # array to hold ordered R strand CFD values


$bin_count -= $bin_window;

	until ($bin_count == $bin_window+1){ #until 4

		# re-order reverse strand cfd freqsum values
		push (@A_R_cfd, $A_R_cfd_freqsum[$A_R_cfd_count]); 
		
		# calculate summed value for both F and R cfd freqsums
		push (@A_FandR_cfd, $A_F_cfd_freqsum[$cfd_count] + $A_R_cfd_freqsum[$A_R_cfd_count]);
	
		$bin_count ++;
		$cfd_count ++;
		$A_R_cfd_count --;

		} #until 4 closer

		
# Need to find average read values over bin_window to normalize data

$A_cfd_sum += $_ for @A_FandR_cfd;
$A_norm_factor = ($A_cfd_sum/(($bin_window*2)+1))/100000; # order of magnitude here is arbitrary

# reset counters once more
$bin_count = (0-$bin_window);
$cfd_count = 0;

########################################################################################################

print "A file normalization factor is: $A_norm_factor\n\n";

########################################################################################################


################################################################################
# Match the correct SgrB file and read values into arrays
################################################################################

opendir(DIR,$sgr_B_indir_path) || die "Unable to access file at: $sgr_B_indir_path $!\n";

$B_counter = 0;
until ($files_sgr_B[$B_counter] =~ $particles[$part_counter]){
	  
	  $B_counter++;
	  if ($B_counter > 10){
		print "Crikey! It's all gone pear-shaped trying to match a Part size to SgrB\n";
		exit;}
	  }
	  
print "Found, and processing, $files_sgr_B[$B_counter] \n";

open(IN, "$sgr_B_indir_path/$files_sgr_B[$B_counter]")
            || die "Unable to open $files_sgr_B[$B_counter]: $!";
        
        # define three new arrays to store the .sgr values from infile
        my @sgr_B_chr;
        my @sgr_B_bin;
		my @sgr_B_freq;
        
        # loop through infile to get values
        while(<IN>){

            chomp;

            # split line by delimiter and store elements in an array
            @line_sgr_B = split('\t',$_);

            # store the columns we want in the three new arrays
            push(@sgr_B_chr,$line_sgr_B[0]);
            push(@sgr_B_bin,$line_sgr_B[1]);
	    push(@sgr_B_freq,$line_sgr_B[2]);
        }
        
        # close in file handle
        close(IN);
		closedir(DIR);

# store size of bin array
        $sgr_B_size = @sgr_B_freq;

print "which contains a whopping: $sgr_B_size bin values\n\n";


#################################################################################
# Use sgrB values to determine a CFD nomalization value
#################################################################################

# define some arrays that will be reset during each iteration
my @B_F_cfd_bin;
my @B_F_cfd_freqsum;
my @B_R_cfd_bin;
my @B_R_cfd_freqsum; 


# The F strand calculations:

# some counter variables
$site_count = 0; # Counter for each site ID
$bin_count = 0; # Counter .sgr bin numbers
$cfd_count = 0; # Counter for the cfd arrays
$top_limit = 0; # A top limit for $bin_window
my @B_F_out_chr; # F .sgr output array for chr
my @B_F_out_bin; # F .sgr output array for bin pos
my @B_F_out_freq; # F .sgr output array for read freq
my $B_F_out_size = 0; # Size of F .sgr output arrays
$i=0; # An iterator variable


until ($site_count == $F_siteID_size){ #until 1

# Use %bin_map to jump to correct region of sgr arrays

$bin_count = (int($F_site_pos[$site_count]/$bin_size) + $bin_map{$F_site_chr[$site_count]}) - 3;

# find an .sgr bin which contains the current site
	until ($F_site_chr[$site_count] eq $sgr_B_chr[$bin_count] &&
	       $F_site_pos[$site_count] >= $sgr_B_bin[$bin_count] &&
	       $F_site_pos[$site_count] <  $sgr_B_bin[$bin_count +1]){ #until 2

			$bin_count ++;

	} #until 2 closer

# now that we've found the match, let's write values to the output files

# set the bin_counter BACK $bin_window places and set the $top_limit

$bin_count -= $bin_window;
$top_limit = $bin_count + ($bin_window*2);

# Better test to see if match is close to ends of a chromosome. If so, the reported
# bins and read freqs will be chaemeric - we don't want this so we will ditch such matches

	if($F_site_chr[$site_count] ne $sgr_B_chr[$bin_count] ||
	   $F_site_chr[$site_count] ne $sgr_B_chr[$top_limit]){ #if 1

	print "Can't output forward strand values for $F_site_chr[$site_count] site: $F_site_pos[$site_count]\n";

	     } #if 1 closer
	else { #else 1


# Push the chrn, bin and freq values to the F .sgr arrays and add values to  F cfd freq array

	until ($bin_count == $top_limit+1){ #until 3

		push (@B_F_out_chr,$sgr_B_chr[$bin_count]);
		push (@B_F_out_bin,$sgr_B_bin[$bin_count]);
		push (@B_F_out_freq,$sgr_B_freq[$bin_count]);

		$B_F_cfd_freqsum[$cfd_count] += $sgr_B_freq[$bin_count];

		$bin_count ++;
		$cfd_count ++;


	} #until 3 closer
    } #else 1 closer

$cfd_count = 0;
$bin_count = 0;
$site_count ++;

} #until 1 closer

$B_F_out_size = @B_F_out_chr;


# The R strand calculations:

# reset the counter variables and define some more arrays
$site_count = 0; # Counter for each site ID
$cfd_count = 0; # Counter for the cfd arrays
$bin_count = 0;
my @B_R_out_chr; # R .sgr output array for chr
my @B_R_out_bin; # R .sgr output array for bin pos
my @B_R_out_freq; # R .sgr output array for read freq
my $B_R_out_size = 0; # Size of F .sgr output arrays

until ($site_count == $R_siteID_size){ #until 1

# Use %bin_map to jump to correct region of sgr arrays

$bin_count = (int($R_site_pos[$site_count]/$bin_size) + $bin_map{$R_site_chr[$site_count]}) - 3;

# find an .sgr bin which contains the current site
	until ($R_site_chr[$site_count] eq $sgr_B_chr[$bin_count] &&
	       $R_site_pos[$site_count] >= $sgr_B_bin[$bin_count] &&
	       $R_site_pos[$site_count] <  $sgr_B_bin[$bin_count +1]){ #until 2

			$bin_count ++;

	} #until 2 closer

# now that we've found the match, let's write values to the output files

# set the bin_counter BACK $bin_window places and set the $top_limit

$bin_count -= $bin_window;
$top_limit = $bin_count + ($bin_window*2);

# Better test to see if match is close to ends of a chromosome. If so, the reported
# bins and read freqs will be chaemeric - we don't want this so we will ditch such matches

	if($R_site_chr[$site_count] ne $sgr_B_chr[$bin_count] ||
	   $R_site_chr[$site_count] ne $sgr_B_chr[$top_limit]){ #if 1

	print "Can't output reverse strand values for $R_site_chr[$site_count] site: $R_site_pos[$site_count]\n";

	     } #if 1 closer
	else { #else 1


# Push the chrn, bin and freq values to the R .sgr arrays and add values to  R cfd freq array


	until ($bin_count == $top_limit+1){ #until 3

		push (@B_R_out_chr,$sgr_B_chr[$bin_count]);
		push (@B_R_out_bin,$sgr_B_bin[$bin_count]);
		push (@B_R_out_freq,$sgr_B_freq[$bin_count]);

		$B_R_cfd_freqsum[$cfd_count] += $sgr_B_freq[$bin_count];

		$bin_count ++;
		$cfd_count ++;


	} #until 3 closer
    } #else 1 closer

$cfd_count = 0;
$bin_count = 0;
$site_count ++;

} #until 1 closer

$B_R_out_size = @B_R_out_chr;


# At last, the CFD calculation

# Set counter variables and define new arrays
$bin_count = 0;
$cfd_count = 0;
my $B_cfd_sum = 0; # a sum of sums for normalizing the data
my $B_norm_factor = 0; # calced from $cfd_sum
my $B_R_cfd_count = $bin_window*2;
my @B_FandR_cfd; # array to hold summed F and R strand CFD values
my @B_R_cfd; # array to hold ordered R strand CFD values


$bin_count -= $bin_window;

	until ($bin_count == $bin_window+1){ #until 4

		# re-order reverse strand cfd freqsum values
		push (@B_R_cfd, $B_R_cfd_freqsum[$B_R_cfd_count]); 
		
		# calculate summed value for both F and R cfd freqsums
		push (@B_FandR_cfd, $B_F_cfd_freqsum[$cfd_count] + $B_R_cfd_freqsum[$B_R_cfd_count]);
	
		$bin_count ++;
		$cfd_count ++;
		$B_R_cfd_count --;

		} #until 4 closer

		
# Need to find average read values over bin_window to normalize data

$B_cfd_sum += $_ for @B_FandR_cfd;
$B_norm_factor = ($B_cfd_sum/(($bin_window*2)+1))/100000; # remember the order of magnitude is arbitrary

# reset counters once more
$bin_count = (0-$bin_window);
$cfd_count = 0;

########################################################################################################
print "B file normalization factor is: $B_norm_factor\n\n";


#########################################################################################################

print "The sgr files contained values for: $chr_count chromosomes\n\n";

######################################################################################
# Generate hit-lists of site bin matches for .sgr_A and .sgr_B
######################################################################################

# some counter variables and arrays
$site_count = 0; # Counter for each site ID
$bin_count = 0; # Counter .sgr bin numbers
my $bottom_limit = 0; # A bottom limit for $bin_window
$top_limit = 0; # A top limit for $bin_window
my $A_hitlist_size = 0;
my $B_hitlist_size = 0;
my @A_hitlist; # Array for $bincounts matching sites for the A_.sgr
my @A_hitlist_strand; # Strand record for above
my @B_hitlist; # Array for $bincounts matching sites for the B_.sgr
my @B_hitlist_strand; # Strand record for above

###############################################################################
# The .sgr_A hitlist
###############################################################################

until ($site_count == $siteID_size){ #until 1

# Use %bin_map to jump to correct region of sgr arrays

$bin_count = (int($site_pos[$site_count]/$bin_size) + $bin_map{$site_chr[$site_count]}) - 3;

# find an .sgr_A bin which contains the current site
	until ($site_chr[$site_count] eq $sgr_A_chr[$bin_count] &&
	       $site_pos[$site_count] >= $sgr_A_bin[$bin_count] &&
	       $site_pos[$site_count] <  $sgr_A_bin[$bin_count +1]){ #until 2

		$bin_count ++;
			

	} #until 2 closer

# set the bottom_limit BACK $bin_window places and set the $top_limit FORWARD $bin_window places

$bottom_limit = $bin_count - $bin_window;
$top_limit = $bin_count + $bin_window;

# Test to see if match is close to ends of a chromosome. If so, the reported
# read freqs will be chaemeric - we don't want this so we will ditch such matches

	if($site_chr[$site_count] ne $sgr_A_chr[$bottom_limit] ||
	   $site_chr[$site_count] ne $sgr_A_chr[$top_limit]){ #if 1

	print "Just double-checking...No, still can't output values for $site_chr[$site_count] site: $site_pos[$site_count] for $infile_sgr_A \n";

	     } #if 1 closer
	else { #else 1

		push (@A_hitlist, $bin_count);
		push (@A_hitlist_strand, $site_strand[$site_count]);
	
		
		} #else 1 closer

$site_count ++;
$bin_count = 0;
} #until 1 closer


$A_hitlist_size = @A_hitlist;
$site_count = 0; 
$bin_count = 0;


print "Matched $A_hitlist_size sites to $files_sgr_A[$A_counter]\n\n";

############################################################################
# The .sgr_B hitlist
##############################################################################


until ($site_count == $siteID_size){ #until 1

# Use %bin_map to jump to correct region of sgr arrays

$bin_count = (int($site_pos[$site_count]/$bin_size) + $bin_map{$site_chr[$site_count]}) - 3;

# find an .sgr_B bin which contains the current site
	until ($site_chr[$site_count] eq $sgr_B_chr[$bin_count] &&
	       $site_pos[$site_count] >= $sgr_B_bin[$bin_count] &&
	       $site_pos[$site_count] <  $sgr_B_bin[$bin_count +1]){ #until 2

			$bin_count ++;

	} #until 2 closer

# set the bottom_limit BACK $bin_window places and set the $top_limit FORWARD $bin_window places

$bottom_limit = $bin_count - $bin_window;
$top_limit = $bin_count + $bin_window;
	
# Test to see if match is close to ends of a chromosome. If so, the reported
# read freqs will be chaemeric - we don't want this so we will ditch such matches

	if($site_chr[$site_count] ne $sgr_B_chr[$bottom_limit] ||
	   $site_chr[$site_count] ne $sgr_B_chr[$top_limit]){ #if 1

	print "Just double-checking...No, still can't output values for $site_chr[$site_count] site: $site_pos[$site_count] for $infile_sgr_B\n";

	     } #if 1 closer
	else { #else 1

		push (@B_hitlist, $bin_count);
		push (@B_hitlist_strand, $site_strand[$site_count]);
		
		} #else 1 closer

$site_count ++;
$bin_count = 0;
} #until 1 closer

$B_hitlist_size = @B_hitlist;
$site_count = 0; 
$bin_count = 0;

print "Matched $B_hitlist_size sites to $files_sgr_B[$B_counter]\n\n";



######################################################################################
# The output file
######################################################################################



# define outfile name and set correct endings


		
	$stats_outfile = $descriptor."_with_".substr($infile_siteID,0,-4)."_".$particles[$part_counter]."_SubStats";
        $stats_outfile .= '.txt';
		
		# try and open the stats output file
        open(OUT,"> $outdir_path/$stats_outfile")
             || die "Unable to open $stats_outfile: $!";
        
print "Just creating $stats_outfile\n\n";

# print outfile headers

print (OUT "Comparison of $files_sgr_A[$A_counter] and $files_sgr_B[$B_counter] using $infile_siteID\n");
print (OUT "Bin\tA mean\tA SD\tB mean\tB SD\tt-value\tt-testDF\tt-test_p-value\tWMWp_value\tSpearman's rho\tKendall's tau\n");
		
# Define strand-specific bin marker variables

my $range_counter = 0;
my $F_bin_marker = (0 - $bin_window); # e.g. -40
my $R_bin_marker = $bin_window; # e.g. +40		
		
# changed this to catch all $output_scale values		
until ($range_counter >= (($bin_window*2)+1)){ #until 1

my @A_data;
my @B_data;



	##########################################################################	
	# Produce the list of sgr_A frequency values	
	##########################################################################
	$site_count = 0;
	until ($site_count == $A_hitlist_size){ # until 

	# extract a bin count from the A hitlist
	$bin_count = $A_hitlist[$site_count];

		#test for strand and push frequency value from correct strand-specific position
		if($A_hitlist_strand[$site_count] =~"F"){ #if 
	
			push (@A_data, $sgr_A_freq[$bin_count + $F_bin_marker]/$A_norm_factor);

			} #if closer
		else{ #else
		
			push (@A_data, $sgr_A_freq[$bin_count + $R_bin_marker]/$A_norm_factor);
			
			} #else closer
				
	$site_count ++;
	} #until  closer
	###############################################################################	
	# Produce the list of sgr_B frequency values	
	###############################################################################
	$site_count = 0;
	
	until ($site_count == $B_hitlist_size){ # until 

	# extract a bin count from the B hitlist
	$bin_count = $B_hitlist[$site_count];

		#test for strand and push frequency value from correct strand-specific position
		if($B_hitlist_strand[$site_count] =~"F"){ #if 
	
			push (@B_data, $sgr_B_freq[$bin_count + $F_bin_marker]/$B_norm_factor);

			} #if closer
		else{ #else
		
			push (@B_data, $sgr_B_freq[$bin_count + $R_bin_marker]/$B_norm_factor);
			
			} #else closer
				
	$site_count ++;
	} #until  closer	
		
	###############################################################################
	# Do the stats
	###############################################################################
	print (OUT $F_bin_marker*$bin_size."\t");
	
	# descriptive stats for A_data
	my $stat = Statistics::Descriptive::Full->new();
	$stat -> add_data(@A_data);
	my $mean_A = $stat -> mean();
	print (OUT $mean_A."\t");
	my $sd_A = $stat -> standard_deviation();
	print (OUT $sd_A."\t");
	$stat -> clear();
	
	# descriptive stats for B_data
	$stat -> add_data(@B_data);
	my $mean_B = $stat -> mean();
	print (OUT $mean_B."\t");
	my $sd_B = $stat -> standard_deviation();
	print (OUT $sd_B."\t");
	$stat -> clear();
	
	# paired T-test
	my $t_test = Statistics::DependantTTest->new();
	$t_test -> load_data('before',@A_data);
	$t_test -> load_data('after',@B_data);
	my($t_value,$deg_freedom) = $t_test->perform_t_test('before','after');
	my($p_value) = Statistics::Distributions::tprob($deg_freedom,$t_value);
	print (OUT $t_value."\t".$deg_freedom."\t".$p_value."\t");
	
	# Wilcoxon Rank Sum test (note - this algorithm uses a normalization approximation
	# because the n value is high. I have checked that it outputs numerically identical
	# (to 7s.f.) by running in exact mode)
	
	my $wilcox_test = Statistics::Test::WilcoxonRankSum ->new();
	$wilcox_test -> load_data (\@A_data, \@B_data);
	my $prob = $wilcox_test -> probability();
	print (OUT $prob."\t");
	
	# Spearman's rho and Kendall's tau
	my $rank_corr = Statistics::RankCorrelation->new(\@A_data,\@B_data);
	my $Spear_rho = $rank_corr ->spearman;
	my $Kend_tau = $rank_corr ->kendall;
	print (OUT $Spear_rho."\t".$Kend_tau."\n");
	
	$F_bin_marker += $output_scale;
	$R_bin_marker -= $output_scale;
	$range_counter += $output_scale;
		
		
} #until 1 closer
 # close out file handle
        close(OUT);


}}}
