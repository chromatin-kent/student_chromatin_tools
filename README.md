# student_chromatin_tools
A suite of stand-alone Perl scripts for students to create and process Kent Lab MNase-seq/CPSA data. NOTE: The core suggested scripts for the 2024 FYPs (complete with suggested folder structure) can be downloaded in one go as Script_suite.zip

Rather than write huge app-like structures, Nick tends to write short scripts which perform single tasks. This modularity is designed to get you to think about the principle of the pipeline by inspecting the intermediate data as you go. You can, of course, write your own wrapper shell script to force the data from .fastq -> bedGraph if you want. You could also just use Danpos in which case you'll want to look at SAM2PartN_SAM.plx as a gateway to size-selecting the MNase size distribution data.

1. You should open each script with a code editor before use. This will allow you to read the comment header, which describes what the script does and what input and output to expect. You will also have the opportunity to set some variables; instead of taking values from @ARGV and using command line flags, these Perl scripts force you to set variables by opening the actual script and editing it. The idea here is to get students to view the code...and ultimately to learn to write it better than Nick does. The other reason for doing this is because Nick wants to get a functioning script ASAP without arsing around writing all the @ARGV handling and --help text you'd normally expect with SPades or whatever. Once you have edited the script, remember to Save.

2. You should place each script in its own working directory to keep your data and mind tidy. Check the explcit paths in the script. Some of them are set to take input files from one folder and deposit outputs in another (check the variables!); in this case you will need to create, and correctly name, these folders manually in your current working directory (cwd) before running the script for the first time. Others just need input files to be in the same directory as the script and will helpfully create a time-stamped output directory for you.

Kent Lab Software Disclaimer:
THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. 

IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, 
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO: PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE OF
DATA, LOSS OF PROFITS; BUSINESS INTERRUPTION; DEGREE FLUNKING; LAUNCH OF 
STRATEGIC NUCLEAR WEAPONS) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE, 
IGNORANCE, RANK INSANITY OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
