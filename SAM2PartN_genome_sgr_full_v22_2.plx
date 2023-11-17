#!/usr/bin/perl
#################################################################################
#WARNING: FOR BI3008 2022 Project (Yeast SacCer3) aligned reads ONLY
#################################################################################

use strict;
use warnings;
use Math::Round;
use Cwd;
use POSIX qw(strftime);

#################################################################################
# Written by Nick Kent, Dec 2017
# Fiddled with: Nick Kent, Dec 20th 2017 to add more processing information and time stamps.
# Fiddled with: Nick Kent, Dec 21st 2017 to load SAM data into memory.
# Fiddled with: Nick Kent, Jan 4th 2018 to process multiple PartN values
# Fiddled with: Nick Kent, Jan 5th 2018 to solve end-of-data run-off error
# Fiddled with: Nick Kent, Oct 23rd 2019 to use cwd
# Fiddled with: Nick Kent, Mar 29th 2022 to run under Gomphus environment for BI3008 and to
# create time-stamped output folders using POSIX strftime.
# Fiddled with: Nick Kent, April 1st 2022 to deposit £1M in your bank account
# Fiddled with: Nick Kent, April 2nd 2022 to convert long RefSeq chr ID to short chr form
#################################################################################
# USAGE:- perl SAM2PartN_genome_sgr_full_vxxx.plx
#
# This script takes a SORTED* Bowtie 1 CPSA paired read alignment .sam format file and will 
# generate full genome.sgr files (all chromosomes concatenated) containing 3MA-smoothed read 
# mid-point frequency values for nuclease-resistant chromatin particle reads.
#
# Place the input .sam file in the same folder as the script; the script will generate a 
# date specific out folder in which your .sgr files will appear.
#
# This script will just run on a 16Gb memory Laptop, but probably requires a high memory
# Linux machine for anything bigger than yeast. Typically you will need RAM eqivalent to the 0.75 * size 
# of your .sam file. There is also a limit to the number of reads that can be processed,
# and the number of genomic bins. These limits are both = 2^31 - the max size of a perl array. This
# script was designed for model eukaryote (yeast, Drosophila, Arabidopsis) and prokaryotic
# genomes, but may not work for mammalian size genomes.
#
# User specifies a list of  "PartN" value in bp - e.g. 150 for a nucleosome, 50 for a TF.
# User specifies a window +/- PartN as a fraction of 1 (i.e. 0.2 is eqivalent to +/- 20%)
#
#
# * To sort use: ./samtools sort Input.sam -o Input_sorted.sam
#
# NOTE: this script is a cut-n-shut of pair_read_histo.plx, histogram.plx and SiteWriterCFD.plx
#
# To do (never): MORE ERROR TRAPPING (e.g. non-sorted .sam; general fuck-wittery;); rename test variables.
################################################################################
# SET THE VARIABLES BELOW AS REQUIRED
# $bin_width    - The histogram bin width value in bp (10bp is good for most things)
# @partn 	- A bracketed, comma delimited list of "chromatin particle sizes" e.g. (50,150,300)
# $pwind	- The "particle window" e.g. 0.2 = +/- 20% of $partn
# %chr_pairs - Hash that specifically converts long NCBI RefSeq Id to chrI, chr II
################################################################################

my $bin_width = 10;
my @partn = (50,75,100,125,150,175,200,225,250,275,300,325,350,375,400,425,450);
my $pwind = 0.2;
my %chr_pairs = (
'gi|330443391|ref|NC_001133.9|'=>'chrI',
'gi|330443482|ref|NC_001134.8|'=>'chrII',
'gi|330443489|ref|NC_001135.5|'=>'chrIII',
'gi|330443520|ref|NC_001136.10|'=>'chrIV',
'gi|330443531|ref|NC_001137.3|'=>'chrV',
'gi|330443543|ref|NC_001138.5|'=>'chrVI',
'gi|330443578|ref|NC_001139.9|'=>'chrVII',
'gi|330443590|ref|NC_001140.6|'=>'chrVIII',
'gi|330443595|ref|NC_001141.2|'=>'chrIX',
'gi|330443638|ref|NC_001142.9|'=>'chrX',
'gi|330443667|ref|NC_001143.9|'=>'chrXI',
'gi|330443681|ref|NC_001144.5|'=>'chrXII',
'gi|330443688|ref|NC_001145.3|'=>'chrXIII',
'gi|330443715|ref|NC_001146.8|'=>'chrXIV',
'gi|330443743|ref|NC_001147.6|'=>'chrXV',
'gi|330443753|ref|NC_001148.4|'=>'chrXVI',
'gi|6226515|ref|NC_001224.1|'=>'chrM',
'gi|11466067|ref|NC_001398.1|'=>'2micron');


################################################################################
################################################################################
# MAIN PROGRAM - you should not need to edit below this line
################################################################################
################################################################################

my $datestring =strftime "%F_%H-%M-%S", localtime;
print "\n\n*********************************************";
print "\n*********** Welcome to SAM2PartN ************";
print "\n*********************************************\n";
print "Started:\t".$datestring."\n";

# set paths, outfolder name and define some variables
my $output_dir = $datestring."_sgr_out";

unless (-e cwd."/$output_dir" or mkdir cwd."/$output_dir"){
		die "Unable to create $output_dir\n"
};

my $inA_indir_path =cwd;
my $outdir_path =cwd."/$output_dir";
my $infile_A;
my @line_A;
my @line_B;
my @files_A;
my @SAM_LN;
my @SAM_SN;
my @SAM_chr_id;
my @SAM_chr_length;
my @test_ID;
my @test_pos;
my @test_ISIZE;
my $read_counter = 0;
my %SAM_map; 
my $SAM_data_count;
my $partn_size = @partn;
my $partsize;
my $mapsize;
my $arguement_string;



################################################################################
# Find Chr IDs and gather specific data from .sam
################################################################################

# store input file name in an array
opendir(DIR,$inA_indir_path) || die "Unable to access file at: $inA_indir_path $!\n";

@files_A = readdir(DIR);

# process the input file within indir_path
foreach $infile_A (@files_A){    

    # ignore hidden files and only get those with the correct ending
    if (($infile_A !~ /^\.+/) && ($infile_A =~ /.*\.sam/)){
    
    
# While we're at it, let's print some useful information
print "Frequency distributions will be binned in $bin_width bp intervals \n";
print "\nThe following nuclease-protection/chromatin-particle/PartN sizes have been specified:\n ";
    
    foreach $partsize (@partn){
	print "$partsize bp\n ";
    
			}
print "\nParticle size window will be +/-".($pwind*100)." percent.\n";
print "\nFound, and processing, $infile_A \n\n";


open(IN, "$inA_indir_path/$infile_A")
            || die "Unable to open $infile_A: $!";
        
       
	
	# loop through top of infile to get header values
        while(<IN>){
           
	    chomp;

            # split line by delimiter and store elements in an array
            @line_A = split('\t',$_);
            my $line_A_size = @line_A;
           

            # test for the correct headers and extract chr id and chr length
            # load three columns of data into huge arrays
	    
				
	    if ($line_A[0] eq '@HD'){
		print "Found the SAM header and the following chromosome IDs and lengths: \n";
					}
			
	    elsif ($line_A[0] eq '@SQ'){
		@SAM_SN = split(':',$line_A[1]);
		@SAM_LN = split(':',$line_A[2]);
		push (@SAM_chr_id, $chr_pairs{$SAM_SN[1]}); # changed here to hash chr ID
		push (@SAM_chr_length, $SAM_LN[1]);
		print "Chromosome ID: $chr_pairs{$SAM_SN[1]}, Length: $SAM_LN[1] bp \n"; # changed here to hash chr ID
			}
	    elsif ($line_A[0] eq '@PG'){
		print "End of the SAM header.\n";
        print "Starting to read SAM data into memory - this might take a while; please be patient.\n";
					}
			
	    elsif ($line_A_size >3 && $line_A[8]>0){
		push (@test_ID,$chr_pairs{$line_A[2]}); # changed here to hash chr ID
		push (@test_pos,$line_A[3]);
		push (@test_ISIZE,$line_A[8]);
		$read_counter ++;
				}
				
	  elsif ($line_A_size >3 && $line_A[8]<=0){
				}
				
	  else{
		print "\n Not sure this is a normal SAM file. Will stop so that you can re-check \n";
		exit;
	  }
        }

	# close in file handle
        close(IN);
	
	my $chr_list_size = @SAM_chr_id;
	my $data_list_size = @test_ID;
	
	
		
#######################################################################################
# BUILD AN ARRAY MAP
#######################################################################################

print "\n\nIndexing all the data according to chromosome ID\n";
my $map_count = 0; # a counter variable

# Set bottom
$SAM_map{$test_ID[$map_count]} = 0;

$map_count ++;

# scan through the @test_ID array and mark the places where each new chromsomome starts

until ($map_count == $data_list_size){
  
      if ($test_ID[$map_count] ne $test_ID[$map_count-1]){
      
      $SAM_map{$test_ID[$map_count]} = $map_count;
      $map_count ++;
      
      }
      else{
      
	  $map_count ++;
	  
	  }

}
# output the number of chromosome types found as the number of hash keys.
$SAM_data_count = keys %SAM_map;
$mapsize = $map_count;

print "The data contained values corresponding to: $SAM_data_count chromosome(s)\n\n";



################################################################################
# Cycle through each PartN size
################################################################################
foreach $partsize (@partn){

my $total_particles;
# define outfile name from infile name NOTE this changes original multiple file output to concatenation
my $outfile = substr($infile_A,0,-4)."_Part".$partsize."_".$bin_width;
$outfile .= '.sgr';

# try and open output file
        open(OUT,"> $outdir_path/$outfile")
             || die "Unable to open $outfile: $!";
        


# define lower and upper values of ISIZE as a percentage from $pwind
my $ISIZE_low = ($partsize - ($pwind * $partsize));
my $ISIZE_high =($partsize + ($pwind * $partsize));


################################################################################
# Plot histogram of PartN midpoint position for each chromosome
################################################################################

my $chr_counter = 0; # a counter variable
$map_count =0; #reset $map_count


until ($chr_counter == $chr_list_size){
  
	# set top bin of histogram
	my $top = $SAM_chr_length[$chr_counter];
	
        
        # define new array to store required dyad (paired read mid-point) position values
        my @dyad_pos;
        my $probe = 0;

        
        # use bin map to retrieve values
        
        $map_count = $SAM_map{$SAM_chr_id[$chr_counter]}; 
        
        while($SAM_chr_id[$chr_counter] eq $test_ID[$map_count]){
	    #test for end of data
            if ($map_count>=$mapsize -1){
            last;
            }      
	    #test for match within pwindow
            if ($test_ISIZE[$map_count] > $ISIZE_low && $test_ISIZE[$map_count] < $ISIZE_high){
                       
		push(@dyad_pos,($test_pos[$map_count] + ($test_ISIZE[$map_count] * 0.5)));
		$probe ++;
		$map_count ++;
		   
	      }
	  else{
		$map_count ++;
		}
				

        
   }
 

        # Tally counter to plot histogram
		
		my $dyadarray_size= @dyad_pos;
		
		# Define the number of bins for the relevant chromosome
		my $bin_no = (int($top/$bin_width))+1;
		
		# Define the distribution frequency array
		my @dist_freq;
		my $i=0;
		
		# Fill the frequency distribution "bins" with zeroes
		for ($i=0; $i<$bin_no; $i++){
			push (@dist_freq, 0);
			}
			
		# Reset the incrementor and define the bin hit variable
		$i=0;
		my $bin_hit = 0;
		
		# The tally counter 
		while ($i < $dyadarray_size){
			$bin_hit = int($dyad_pos[$i]/$bin_width);
			$dist_freq[$bin_hit] ++;
			$i ++;
			}
		
		# Calculate the 3 bin moving average
		my @moving_ave;
		my $ma = 0;
		my $count = 1;
		push (@moving_ave,0);
		
		while ($count<$bin_no-1){
			$ma = (($dist_freq[$count--] + $dist_freq[$count] + $dist_freq[$count++])/3);
			push (@moving_ave,$ma);
			$count ++;
			}
			push (@moving_ave,0);
        
            # print required data in tab-delimited format to output file
            # NK modifies to output chrn, bin and ma only
			for ($i=0; $i<$bin_no; $i++){
			
            print(OUT $SAM_chr_id[$chr_counter]."\t".($i*$bin_width)."\t".round($moving_ave[$i])."\n");
			}
            $chr_counter ++;
			$total_particles += $probe;
            #print "Output to: $outfile having found $probe data points for $partsize bp particles.\n";
        }
        
    print "Output to file having found $total_particles data points for $partsize bp particles.\n";     
         
 }

      
       
}

   # close out file handle
        close(OUT);
 
}

$datestring =strftime "%F_%H-%M-%S", localtime;
print "Finished:\t".$datestring." Wait for prompt - Perl is sluicing the chunk of RAM you just used.\n";