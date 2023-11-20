#!/usr/bin/perl
# Written by Nick Kent, Apr 2021
# Fiddled by Nick Kent 10th April 2021 to fix Read 2 bug

#
# USAGE:- perl SAM2PartN_SAM.plx
#
# This script will take a Bowtie 1 CPSA paired read alignment .sam format file
# and will select just the alignments of a particular chromatin particle size.
#
# User specifies a "PartN" value in bp - e.g. 150 for a nucleosome, 50 for a TF.
# User specifies a window +/- PartN as a fraction of 1 - 0.2 is eqivalent to +/- 20%
#
# The output is in .SAM format (complete with headers) and can then be wazzed into Danpos
# for peakmarking and dodgy stats. If you are going to use Danpos - this process *is* important;
# Danpos struggles valiently with CPSA data but makes cleaner nucleosome calls if you size-
# select first!
#
# NOTE: this script is the bastard love-child of ...SAM2PartN_sgr and SAM_parser.
################################################################################
use strict;
use warnings;
use Cwd;
################################################################################
# SET THE VARIABLES BELOW AS REQUIRED
# $indir_path   - The directory containing the .txt files to be processed
# $outdir_path  - The directory to store the .sgr output files
# $bin_width    - The histogram bin width value in bp
# $partn 	- The "chromatin particle size" e.g. 150 for a nucleosome
# $pwind	- The "particle window" e.g. 0.2 = +/- 20% of $partn
################################################################################

my $inA_indir_path =cwd."/in";
my $outdir_path =cwd."/out";
my $partn = 150;
my $pwind = 0.2;

################################################################################
################################################################################
# MAIN PROGRAM - you should not need to edit below this line
################################################################################
################################################################################
print "start:\t", `date`."\n";

# define some variables

my $infile_A;
my @line_A;
my @line_B;
my @files_A;
my @SAM_LN;
my @SAM_SN;
my @SAM_chr_id;
my @SAM_chr_length;
my @SAM_out; # this will be a whopper



# define lowwer and upper values of ISIZE as a percentage from $pwin

my $ISIZE_low = ($partn - ($pwind * $partn));
my $ISIZE_high =($partn + ($pwind * $partn));

################################################################################
# Check for and analyse SAM headers and add to output
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

	
	# loop through top of infile to get header values
        while(<IN>){
           
	    chomp;

            # split line by delimiter and store elements in an array
            @line_A = split('\t',$_);
            my $line_A_size = @line_A;

            # test for the correct headers and extract chr id and chr length
	    if ($line_A_size >3){
	    last;
				}
				
	    elsif ($line_A[0] eq '@HD'){
		print "Reassuringly, have found a SAM header, and found the following chromosome IDs and lengths: \n";
		push (@SAM_out, $_);
					}
					
	    elsif ($line_A[0] eq '@SQ'){
		@SAM_SN = split(':',$line_A[1]);
		@SAM_LN = split(':',$line_A[2]);
		push (@SAM_chr_id, $SAM_SN[1]);
		push (@SAM_chr_length, $SAM_LN[1]);
		push (@SAM_out, $_);
		print "Chromosome ID: $SAM_SN[1], Length: $SAM_LN[1] bp \n";
			}
			
		elsif ($line_A[0] eq '@PG'){
        push (@SAM_out, $_);
		}
	    else {
		print "Not sure this is a normal SAM file. Will stop so that you can re-check \n";
		exit;

	      }
        }

	# close in file handle
        close(IN);
	
	my $chr_list_size = @SAM_chr_id;
	
		
################################################################################
# Grab PartN data from SAM
################################################################################



	# define outfile name from infile name
        my $outfile = substr($infile_A,0,-4)."_Part".$partn;
        $outfile .= '.sam';

        open(IN, "$inA_indir_path/$infile_A")
            || die "Unable to open $infile_A: $!";
            
        my $hit = 0;
        my $probe = 0;

        
        # loop through infile to get values
        while(<IN>){
        
        chomp;
            
            # split line by delimiter and store elements in an array
            @line_B = split('\t',$_);
            # test for headers
             if($line_B[0] ne '@HD' && $line_B[0] ne '@SQ' && $line_B[0] ne '@PG'){
             
            # test for Read2 required from previous hit
             if ($line_B[8] < 0 && $hit == 1){
                  push (@SAM_out, $_);
                  $hit =0;
                  }
             
            
			# test for match within pwindow
                if ($line_B[8] > $ISIZE_low && $line_B[8] < $ISIZE_high && $hit == 0){
                   push (@SAM_out, $_);
				   $probe ++;
				   $hit = 1;
				   print "Currently working and have found $probe data points so far.\r";
				   
				   
					      }
				      }
           
        }

        # close in file handle
        close(IN);
  
        my $out_size = @SAM_out;
            



        # try and open output file
        open(OUT,"> $outdir_path/$outfile")
             || die "Unable to open $outfile: $!";
        
        
            # print required data in tab-delimited format to output file
			for (my$i=0; $i<$out_size; $i++){
			
            print(OUT $SAM_out[$i]."\n");
			}

            print "\nJust output: $outfile\n";

   # close out file handle
        close(OUT);
}
}
print "end:\t", `date`."\n";
