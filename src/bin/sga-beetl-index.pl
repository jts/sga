#! /usr/bin/perl
# Run Illumina's BEETL BWT builder on a set of sequence reads
# and convert the output to sga's index format.
use strict;
use Getopt::Long;
use Cwd;
use File::Temp;
use File::Basename;

# Program paths
my $SGA_BIN = "~/bin/sga-0.9.14";
my $BEETL_BIN = "/nfs/team118/js18/software/BEETL-tc/BCRext";

# Options
my $bHelp = 0;
my $tmpDir = ".";
GetOptions("help"     => \$bHelp,
           "tmp-dir=s" => \$tmpDir);

if($bHelp)
{
    usage();
    exit(0);
}

if(scalar(@ARGV) != 1)
{
    print STDERR "Error: A single input file must be provided\n\n";
    usage();
    exit(1);
}

# Parse the input filename and convert it to an absolute path
my $input = $ARGV[0];
my $abs_input = Cwd::abs_path($input);
print "Building index for: $abs_input\n";

# Make a temporary directory for BEETL to work in
my $beetl_dir = File::Temp::tempdir("beetl-working-XXXX", DIR => $tmpDir, CLEANUP => 1);
print "Working directory: $beetl_dir\n";

# Save the current directory and change to the temp
my $final_dir = Cwd::getcwd();
chdir($beetl_dir);

# Currently using the development version of BEETL
# Flatten the input file
my $flat_file = "flattened.txt";
my $bIsFastq = isFastq($abs_input);
my $ret = 0;
if($bIsFastq == 1)
{
    runCmd("cat $abs_input | awk 'NR % 2 == 0 && NR % 4 > 0' > $flat_file");
}
elsif($bIsFastq == 0)
{
    runCmd("cat $abs_input | awk 'NR % 2 == 0' > $flat_file");
}
else
{
    print STDERR "Error: Cannot determine if $abs_input is FASTQ or FASTA\n"; 

    # Set error status so no further commands are run
    $ret = 1;
}

# Run BEETL on the flattened input file
my $time_str = `date`;
print "Starting beetl at $time_str\n";
$ret = runCmd("$BEETL_BIN $flat_file > beetl.status") if($ret == 0);

# concatenate the beetl output files
my $beetl_bwt = "beetl.bwt";
$ret = runCmd("cat bwt-B* > $beetl_bwt") if($ret == 0);

# Run sga convert-beetl
$time_str = `date`;
print "Starting convert-beetl at $time_str\n";
$ret = runCmd("$SGA_BIN convert-beetl $beetl_bwt $abs_input") if($ret == 0);

# Copy the final files back to the original directory
my $base = File::Basename::basename($abs_input, (".fa", ".fastq"));
$ret = runCmd("mv $base* $final_dir") if($ret == 0);

# Change back to the original directory
chdir($final_dir);

sub usage
{
    print "Usage:\n";
    print "sga-beetl-index.pl FILE\n";
    print "Build a BWT of the input FILE using Illumina's BEETL program\n";
    print "Options:\n\n";
    print "--no-reverse          suppress construction of the reverse BWT\n";
    print "--tmp-dir=DIR         use DIR as the working directory for BEETL.\n";
    print "                      as BEETL is a disk-based algorithm, the performance\n";
    print "                      of the underlying storage system has a large impact on\n";
    print "                      the running time. For best results, use an SSD.\n";
    print "                      Alternatively, try local disk instead of shared\n";
    print "                      parallel storage like lustre\n";
}

sub isFastq
{
    my($file) = @_;
    open(F, $file);
    my $first = <F>;
    my $isFastq = -1;
    if($first =~ /^>/)
    {
        $isFastq = 0;
    }
    elsif($first =~ /^@/)
    {
        $isFastq = 1;
    }
    else
    {
        $isFastq = -1;
    }
    close(F);

    return $isFastq;
}

sub runCmd
{
    my($cmd) = @_;
    print $cmd . "\n";
    my $ret = system($cmd);
    if($ret != 0)
    {
        print STDERR "Error: Failed to run cmd: \"$cmd\"\n";
    }
    return $ret;
}

