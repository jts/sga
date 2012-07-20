#! /usr/bin/perl
# Run Illumina's BEETL BWT builder on a set of sequence reads
# and convert the output to sga's index format.
use strict;
use Getopt::Long;
use Cwd;
use File::Temp;
use File::Basename;

# Program paths
my $SGA_BIN = "/nfs/users/nfs_j/js18/work/devel/sga/src/build-lenny/SGA/sga";
my $BEETL_BIN = "/nfs/users/nfs_j/js18/work/devel/BEETL/Beetl";

# Options
my $bHelp = 0;
my $bNoConvert = 0;
my $tmpDir = ".";
GetOptions("help"     => \$bHelp,
           "tmp-dir=s" => \$tmpDir,
           "no-convert" => \$bNoConvert);

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

# Save the output name for the files
my $base = File::Basename::basename($abs_input, (".fa", ".fastq"));

# BCR expects FASTA input, convert if necessary
my $bcr_in_file = $abs_input;
my $bIsFastq = isFastq($abs_input);
my $ret = 0;
if($bIsFastq == 1)
{
    my $tmp_fasta = "bcr_input.fa";
    fastq2fasta($bcr_in_file, $tmp_fasta);
    $bcr_in_file = $tmp_fasta;
}

# Run BEETL on the flattened input file
my $time_str = `date`;
print "Starting beetl at $time_str\n";

# BCRext
$ret = runCmd("$BEETL_BIN ext -i $bcr_in_file -a > beetl.status") if($ret == 0);

# BCR
#$ret = runCmd("$BEETL_BIN bcr -i $bcr_in_file -o bcr.test > beetl.status") if($ret == 0);

# concatenate the beetl output files
my $beetl_bwt = $final_dir . "/" . $base . ".beetl.bwt";
$ret = runCmd("cat BCRext-B* > $beetl_bwt") if($ret == 0);

# Run sga convert-beetl
$time_str = `date`;
if(!$bNoConvert) {
    print "Starting convert-beetl at $time_str\n";
    $ret = runCmd("$SGA_BIN convert-beetl $beetl_bwt $abs_input") if($ret == 0);

    # Copy the final files back to the original directory
    $ret = runCmd("mv $base* $final_dir") if($ret == 0);
}

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

sub fastq2fasta
{
    my($in, $out) = @_;
    open(IN, $in) || die("Cannot read $in");
    open(OUT, ">$out") || die("Cannot write to $out");
    while(my $h = <IN>) {
        $h =~ s/^@/\>/;
        my $s = <IN>;
        print OUT $h;
        print OUT $s;
        $h = <IN>;
        $h = <IN>;
    }
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

