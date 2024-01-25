#!/usr/bin/perl
# Written: Nick Kent, 23rd Aug 2011
# Last updated: Nick Kent, 24th Apr 2012
# Last updated: Nick Kent, 8th Dec 2019; Use Cwd,
# Revised : Nick Kent, 24th Jan 2024 to handle unstranded Rossi et al., .bed files
# USAGE:- perl read the filename innit?
#
# FUNCTION:
# This script takes Rossi et al., .bed files containing a list of ChIP-exo-seq features (these
# and compares it with MNase -seq sequence read frequency distribution data in 
# whole-genome, Partn .sgr files. It then outputs: (i) CUMULATIVE FREQUENCY DISTRIBUTION
# values over a user-specified bin range centered on, and surrounding, the sites;
# (ii) The individual frequency distribution values for all bins (normalised) and
# formatted for input into Cluster 3.0 for k-means/HC analysis; (iii) F and R strand 
# specific .sgr files which show the raw traces for each of the loci which went up to 
# make the CFD.
#
# Sites close to chromosome ends, which would not yield the full range of data are ignored,
# as usual, but reported at the command line.
#
#
# INPUT AND OUTPUT (all tab-delimited):
# The input .bed files should have Rossi et al., format (6 columns ending in a dot).
#
# The input .sgr files should have three columns: chrn; bin pos; paired-read mid-point freq.
#
# The output CFD.txt file has an input file header and column headers and returns 
# 5 columns: Bin (relative to Site); F strand cumulative freq; R strand cumulative freq;
# summed F+R cumulative freq; normalised cumulative freq. The idea is to plot the first
# and last columns as a line graph to produce a TREND GRAPH for the data. Each bins F+R
# frequencies are normalised to the average F+R frequency for the entire bin window.
# 
# The output C3.txt file is a mix of the Site ID and .sgr data: it outputs a tab-
# delimited table of locally-normalized dyad frequency values for every bin position
# at and surrounding each input site. This text file can then be directly imported
# into Cluster 3.0 (http://bonsai.hgc.jp/~mdehoon/software/cluster/software.htm); The 
# raw data from this file can also be extracted on a column-by-column basis to perform
# a correlation analysis between two data sets (run this script on each .sgr dataset).
# (All this means that here are separate $C3_norm_factor and $CFD_norm_factor variables.
# The latter gives an average CFD value of 1.0, wheras the former scales
# the numbers upwards to help Cluster 3.0)
#
# Note: The script would handles F and R strand data separately. It has been kludged to
# ignore this but if you change it back note this:If you give it all F (or
# all R) strand sites it will work just fine, however, it will also throw a load of
# uninitialised variable warnings at the command line. If you find this upsetting,
# stick a # in front of use warnings (below)
# 
# For development see: Kent et al.,(2011) Chromatin particle spectrum analysis: a 
# method for comparative chromatin structure analysis using paired-end mode 
# next-generation DNA sequencing. NAR 39: e26.
################################################################################

use strict;
use warnings;
use Cwd;
use Math::Round;
use List::Util;
use POSIX qw(strftime);

################################################################################
# SET THE VARIABLES BELOW AS REQUIRED
#
# $bin_window - number of bins surrounding the site of interest. E.g. if you set
# this to 40 then you will get 40 bins either side of your site - 400bp if you
# were using 10bp binned data.
#
# $bin_size - binning interval of .sgr file in base pairs.
#
# $output_scale - controls how many bins are included in the output file. If set
# to 1 you will get every bin (use this). Set to 3 to output only every third bin
# in the series.You can use this feature to scale output files derived from input 
# .sgr data with different bin intervals.
################################################################################

my $bin_window = 40;
my $bin_size = 10;
my $output_scale = 1;

################################################################################
################################################################################
# MAIN PROGRAM
################################################################################
################################################################################
my $datestring =strftime "%F_%H-%M-%S", localtime;
print "\n\n*********************************************";
print "\n*********** Generate k-means input ************";
print "\n*********************************************\n";
print "Started:\t".$datestring."\n";



# define some variables
my $output_dir = $datestring."_Cluster3-kmeans_input_out";

unless (-e cwd."/$output_dir" or mkdir cwd."/$output_dir"){
		die "Unable to create $output_dir\n"
};
my $sgr_indir_path =cwd;
my $siteID_indir_path =cwd;
my $outdir_path =cwd."/$output_dir";


my $descriptor;
my $cwd = getcwd;
my $infile_sgr;
my $infile_siteID;
my $cfd_outfile;
my $c3_outfile;
my $F_sgr_outfile;
my $R_sgr_outfile;
my @line_siteID; 
my @line_sgr;
my @files_sgr;
my @files_siteID;
my $sgr_size;
my $F_siteID_size;
my $R_siteID_size;
my %bin_map;
my $chr_count;

################################################################################
# Get site list and write info to arrays - from Rossi .bed format
################################################################################

# store input file name in an array
opendir(DIR,$siteID_indir_path) || die "Unable to access file at: $siteID_indir_path $!\n";

@files_siteID = readdir(DIR);

# process the input file within siteID_indir_path
foreach $infile_siteID (@files_siteID){    

    # ignore hidden files and only get those ending .bed
    if (($infile_siteID !~ /^\.+/) && ($infile_siteID =~ /.*\.bed/)){

    $descriptor = substr($infile_siteID,0, -4);        
        
print "Found, and processing, $infile_siteID \n";

open(IN, "$siteID_indir_path/$infile_siteID")
            || die "Unable to open $infile_siteID: $!";
        
        # define strand-specific arrays to store site ID, chromosome no., and position
		# NOTE - Rossi data does not include strand specificity BUT have retained strand-specificity
		# function here so that you can modify this script in the future if .bed data does!		
		my @F_site_ID;
		my @F_site_chr;
		my @F_site_pos;

		my @R_site_ID;
		my @R_site_chr;
		my @R_site_pos;
	
	# loop through infile to get values
        while(<IN>){

	    chomp;

	    # split line by delimiter and store elements in an array
            @line_siteID = split('\t',$_);

			# ROSSI KLUDGE
            # store the required chrn, position in two pairs of strand-specific
			# arrays
			
			if($line_siteID[0] =~ "chr"){ # infile if 1
			#altered here to force numeric/string and to take Rossi .bed data as a F strand ID
				push(@F_site_ID,$line_siteID[3]);
				push(@F_site_chr,"".$line_siteID[0]);
				push(@F_site_pos,$line_siteID[1]);
				
			}elsif($line_siteID[3] =~ "ReverseStrandInfo"){
			#altered here to force numeric/string - with Rossi data this should fail to match anything
				push(@R_site_ID,$line_siteID[1]);
				push(@R_site_chr,"".$line_siteID[0]);
				push(@R_site_pos,$line_siteID[2]);
				
			}else{
			
				print "Failed to match strand at $line_siteID[0], $line_siteID[1], $line_siteID[2]\n";
			
			} #infile if 1 closer

        }

	# close in file handle
        close(IN);
	closedir(DIR);
        
        # store sizes of the arrays
        $F_siteID_size = @F_site_pos;
		$R_siteID_size = @R_site_pos;


print "Contains: $F_siteID_size forward strand site IDs; $R_siteID_size reverse
 strand site IDs\n";


################################################################################
# Read in the .sgr file values to three enormous arrays
################################################################################


opendir(DIR,$sgr_indir_path) || die "Unable to access file at: $sgr_indir_path $!\n";

@files_sgr = readdir(DIR);

# process the input file within sgr_indir_path
foreach $infile_sgr (@files_sgr){   

# define some arrays that will be reset during each iteration
my @F_cfd_bin;
my @F_cfd_freqsum;
my @R_cfd_bin;
my @R_cfd_freqsum; 

    # ignore hidden files and only get those ending .sgr
    if (($infile_sgr !~ /^\.+/) && ($infile_sgr =~ /.*\.sgr/)){
        
       

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
	    push(@sgr_freq,$line_sgr[2]);
        }
        
        # close in file handle
        close(IN);
	

# store size of bin array
        $sgr_size = @sgr_freq;

print "Contains a whopping: $sgr_size bin values\n";

#######################################################################################
# BUILD THE BIN MAP
#######################################################################################

my $map_count = 0; # a counter variable

# Set bottom
$bin_map{$sgr_chr[$map_count]} = 0;

$map_count ++;

# scan through the @sgr_chr array and mark the bins where each new chromsomome starts

until ($map_count == $sgr_size){
  
      if ($sgr_chr[$map_count] ne $sgr_chr[$map_count-1]){
      
      $bin_map{$sgr_chr[$map_count]} = $map_count;
      $map_count ++;
      
      }
      else{
      
	  $map_count ++;
	  
	  }

}
# output the number of chromosome types found as the number of hash keys.
$chr_count = keys %bin_map;

print "The sgr file contains values for: $chr_count chromosomes\n";


#######################################################################################
# FORWARD STRAND .sgr calculation:
#######################################################################################

# some counter variables
my $site_count = 0; # Counter for each site ID
my $bin_count = 0; # Counter .sgr bin numbers
my $cfd_count = 0; # Counter for the cfd arrays
my $top_limit = 0; # A top limit for $bin_window
my @F_out_chr; # F .sgr output array for chr
my @F_out_bin; # F .sgr output array for bin pos
my @F_out_freq; # F .sgr output array for read freq
my @F_out_ID; # F .sgr output array for relevant site ID string
my $F_out_size = 0; # Size of F .sgr output arrays
my $i=0; # An iterator variable

until ($site_count == $F_siteID_size){ #until 1


# Use %bin_map to jump to correct region of sgr arrays

$bin_count = (int($F_site_pos[$site_count]/$bin_size) + $bin_map{$F_site_chr[$site_count]}) - 3;

# this looks mad, but it allows me to recycle all the code from the last version, and takes up any
# rounding slack which would come from different $bin_size values




# find an .sgr bin which contains the current site
	until ($F_site_chr[$site_count] eq $sgr_chr[$bin_count] &&
	       $F_site_pos[$site_count] >= $sgr_bin[$bin_count] &&
	       $F_site_pos[$site_count] <  $sgr_bin[$bin_count +1]){ #until 2

			$bin_count ++;

	} #until 2 closer

# now that we've found the match, let's write values to the output files

# set the bin_counter BACK $bin_window places and set the $top_limit

$bin_count -= $bin_window;
$top_limit = $bin_count + ($bin_window*2);

# Better test to see if match is close to ends of a chromosome. If so, the reported
# bins and read freqs will be chaemeric - we don't want this so we will ditch such matches

	if($F_site_chr[$site_count] ne $sgr_chr[$bin_count] ||
	   $F_site_chr[$site_count] ne $sgr_chr[$top_limit]){ #if 1

	print "Can't output forward strand values for $F_site_chr[$site_count] site: $F_site_pos[$site_count]\n";

	     } #if 1 closer
	else { #else 1


# Push the chrn, bin and freq values to the F .sgr arrays, push an ID tag, and add values to  F cfd freq array

	until ($bin_count == $top_limit+1){ #until 3

		push (@F_out_chr,$sgr_chr[$bin_count]);
		push (@F_out_bin,$sgr_bin[$bin_count]);
		push (@F_out_freq,$sgr_freq[$bin_count]);
		push (@F_out_ID,$F_site_ID[$site_count]);

		$F_cfd_freqsum[$cfd_count] += $sgr_freq[$bin_count];

		$bin_count ++;
		$cfd_count ++;


	} #until 3 closer
    } #else 1 closer

$cfd_count = 0;
$bin_count = 0;
$site_count ++;

} #until 1 closer

$F_out_size = @F_out_chr;


		
####################################################################################		
# REVERSE STRAND .sgr calculation:
#####################################################################################

# reset the counter variables and define some more arrays
$site_count = 0; # Counter for each site ID
$cfd_count = 0; # Counter for the cfd arrays
$bin_count = 0;
my @R_out_chr; # R .sgr output array for chr
my @R_out_bin; # R .sgr output array for bin pos
my @R_out_freq; # R .sgr output array for read freq
my @R_out_ID; # R .sgr output array array for relevant site ID string
my $R_out_size = 0; # Size of F .sgr output arrays

until ($site_count == $R_siteID_size){ #until 1

# Use %bin_map to jump to correct region of sgr arrays

$bin_count = (int($R_site_pos[$site_count]/$bin_size) + $bin_map{$R_site_chr[$site_count]}) - 3;


# find an .sgr bin which contains the current site
	until ($R_site_chr[$site_count] eq $sgr_chr[$bin_count] &&
	       $R_site_pos[$site_count] >= $sgr_bin[$bin_count] &&
	       $R_site_pos[$site_count] <  $sgr_bin[$bin_count +1]){ #until 2

			$bin_count ++;

	} #until 2 closer

# now that we've found the match, let's write values to the output files

# set the bin_counter BACK $bin_window places and set the $top_limit

$bin_count -= $bin_window;
$top_limit = $bin_count + ($bin_window*2);

# Better test to see if match is close to ends of a chromosome. If so, the reported
# bins and read freqs will be chaemeric - we don't want this so we will ditch such matches

	if($R_site_chr[$site_count] ne $sgr_chr[$bin_count] ||
	   $R_site_chr[$site_count] ne $sgr_chr[$top_limit]){ #if 1

	print "Can't output reverse strand values for $R_site_chr[$site_count] site: $R_site_pos[$site_count]\n";

	     } #if 1 closer
	else { #else 1


# Push the chrn, bin and freq values to the R .sgr arrays and add values to  R cfd freq array


	until ($bin_count == $top_limit+1){ #until 3

		push (@R_out_chr,$sgr_chr[$bin_count]);
		push (@R_out_bin,$sgr_bin[$bin_count]);
		push (@R_out_freq,$sgr_freq[$bin_count]);
		push (@R_out_ID,$R_site_ID[$site_count]);

		$R_cfd_freqsum[$cfd_count] += $sgr_freq[$bin_count];

		$bin_count ++;
		$cfd_count ++;


	} #until 3 closer
    } #else 1 closer

$cfd_count = 0;
$bin_count = 0;
$site_count ++;

} #until 1 closer

$R_out_size = @R_out_chr;


###################################################################################################
# THE CFD and normalization calculation:
###################################################################################################



# Set counter variables and define new arrays
$bin_count = 0;
$cfd_count = 0;
my $cfd_sum = 0; # a sum of sums for normalizing the data
my $CFD_norm_factor = 0; # calced from $cfd_sum
my $C3_norm_factor = 0; # calced from $cfd_sum
my $R_cfd_count = $bin_window*2;
my @FandR_cfd; # array to hold summed F and R strand CFD values
my @R_cfd; # array to hold ordered R strand CFD values
my @normalized_cfd;

$bin_count -= $bin_window;

	until ($bin_count == $bin_window+1){ #until 4
	
		# ROSSI KLUDGE - very, very bad Perl!
		no warnings 'uninitialized';

		# re-order reverse strand cfd freqsum values
		push (@R_cfd, $R_cfd_freqsum[$R_cfd_count]); 
		
		# calculate summed value for both F and R cfd freqsums
		push (@FandR_cfd, $F_cfd_freqsum[$cfd_count] + $R_cfd_freqsum[$R_cfd_count]);
	
		$bin_count ++;
		$cfd_count ++;
		$R_cfd_count --;

		} #until 4 closer

		
# Need to find average read values over bin_window to normalize data

$cfd_sum += $_ for @FandR_cfd;
$C3_norm_factor = ($cfd_sum/(($bin_window*2)+1))/1000000;
$CFD_norm_factor = ($cfd_sum/(($bin_window*2)+1));

#############################################################################################
#############################################################################################
# The output files
############################################################################################
###########################################################################################

# define outfile names and set correct endings

        $F_sgr_outfile = substr($infile_sgr,0,-4)."_".$descriptor."_F";
        $F_sgr_outfile .= '.sgr';
		
	$R_sgr_outfile = substr($infile_sgr,0,-4)."_".$descriptor."_R";
        $R_sgr_outfile .= '.sgr';
		
	$cfd_outfile = substr($infile_sgr,0,-4)."_".$descriptor."_CFD";
        $cfd_outfile .= '.txt';

	$c3_outfile = substr($infile_sgr,0,-4)."_".$descriptor."_C3";
        $c3_outfile .= '.txt';

#############################################################################################
# The CFD file
##############################################################################################

# try and open the .cfd output file
        open(OUT,"> $outdir_path/$cfd_outfile")
             || die "Unable to open $cfd_outfile: $!";
        


# reset counters once more
$bin_count = (0-$bin_window);
$cfd_count = 0;

# print a header for the CFD.txt file so you can read it in Excel

print (OUT "Values from $cfd_outfile\n");
print (OUT "CFD sum: $cfd_sum\n");
print (OUT "Normalization Factor: $CFD_norm_factor\n");

# print column headers

print (OUT "Bin"."\t"."F Freq"."\t"."R Freq"."\t"."Comb Freq"."\t"."Norm Freq"."\n");

# print data values
		until ($bin_count >= $bin_window+1){ #until 5
		# ROSSI KLUDGE - very, very bad Perl!
		no warnings 'uninitialized';
		print(OUT 
		$bin_count*$bin_size."\t".
		$F_cfd_freqsum[$cfd_count]."\t".
		$R_cfd[$cfd_count]."\t".
		$FandR_cfd[$cfd_count]."\t".
		$FandR_cfd[$cfd_count]/$CFD_norm_factor."\n");
	
		$bin_count += $output_scale;
		$cfd_count += $output_scale;
		

		} #until 5 closer
		
 # close .cfd out file handle
        close(OUT);
print "Have just created $cfd_outfile\n";

#############################################################################################
# The .sgr files
#############################################################################################

# try and open the forward strand .sgr output file
        open(OUT,"> $outdir_path/$F_sgr_outfile")
             || die "Unable to open $F_sgr_outfile: $!";
        


# Output the F .sgr file
	for ($i=0; $i < $F_out_size; $i++){

	      print(OUT 
		$F_out_chr[$i]."\t".
		$F_out_bin[$i]."\t".
		$F_out_freq[$i]."\n");

} # for closer

 # close .sgr out file handle
        close(OUT);
print "Have just created $F_sgr_outfile\n";

# try and open the reverse strand .sgr output file
        open(OUT,"> $outdir_path/$R_sgr_outfile")
             || die "Unable to open $R_sgr_outfile: $!";
        


# Output the R .sgr file
# The reverse strand .sgr frequency values are shown as negative numbers to aid visualisation
	for ($i=0; $i < $R_out_size; $i++){

	      print(OUT 
		$R_out_chr[$i]."\t".
		$R_out_bin[$i]."\t".
		(0 - $R_out_freq[$i])."\n");

} # for closer


 # close .sgr out file handle
        close(OUT);
print "Have just created $R_sgr_outfile\n";

################################################################################################
# The Cluster 3.0 readable file
################################################################################################
	# try and open the .cfd output file
        open(OUT,"> $outdir_path/$c3_outfile")
             || die "Unable to open $c3_outfile: $!";

# write the Cluster 3.0 headers

print (OUT "ID_string"); # this may have to be YORF - and the IDs in YAL001C format

$bin_count = (0-$bin_window);

		until ($bin_count >= $bin_window+1){ 

		print(OUT 
		"\t".$bin_count*$bin_size);
		$bin_count += $output_scale;
		}

print (OUT "\n");


my $group_counter =0; # another counter for the groups of values
 
# write F strand IDs and values

	until ($group_counter >= $F_out_size){
	      
	      print (OUT $F_out_ID[$group_counter]); # print a row header

	      #write tab-delimited row of $bin_window/2 +1 values
	      for ($i=0; $i<($bin_window*2)+1; $i++){

		  print (OUT "\t".round($F_out_freq[$group_counter]/$C3_norm_factor));

		  $group_counter ++;
	      }

	      print (OUT "\n"); # finish row with a newline

	}


 
# write R strand IDs and values
$group_counter =0;

#reverse the required R arrays

@R_out_freq = reverse @R_out_freq;
@R_out_ID = reverse @R_out_ID; 



	until ($group_counter >= $R_out_size){
	      
	      print (OUT $R_out_ID[$group_counter]); # print a row header

	      #write tab-delimited row of $bin_window/2 +1 values
	      for ($i=0; $i<($bin_window*2)+1; $i++){

		  print (OUT "\t".round($R_out_freq[$group_counter]/$C3_norm_factor));

		  $group_counter ++;
	      }

	      print (OUT "\n"); # finish row with a newline

	}

 
 # close c3.txt out file handle
        close(OUT);
print "Have just created $c3_outfile\n";
   	    	


}}}}
