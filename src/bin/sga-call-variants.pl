#! /usr/bin/perl

use strict;
use Getopt::Long;
use File::Basename;
use File::stat;
use Time::localtime;

my $query_file = "";
my $control_file = "";
my $reference_file = "";
my $help = 0;
my $prefix = "sga-varcall";
my $threads = 0;
my $SGA_BIN = "~/work/devel/sga/src/build-lenny/SGA/sga";

# Variant calling parameters
# If these are -1, the defaults will be used
my $k = -1;
my $x = 5;
my $m = -1;

GetOptions("query-reads=s"         => \$query_file,
           "control-reads=s"       => \$control_file,
           "reference=s"           => \$reference_file,
           "threads=i"             => \$threads,
           "prefix=s"              => \$prefix,
           "kmer=i"                => \$k,
           "min-discovery-count=i" => \$x,
           "min-graph-count=i"     => \$m,
           "help"                  => \$help);

# Check for dependencies
checkDependency($SGA_BIN);

#
if($help) {
    printUsage();
    exit(0);
}

# A reference genome and set of query reads are required
if($reference_file eq "") {
    print STDERR "Error: a reference genome is required (--reference)\n";
    exit(1);
}

if(! -f $reference_file) {
    print STDERR "Error: the reference genome is not readable (file: $reference_file)\n";
    exit(1);
}

if($query_file eq "") {
    print STDERR "a set of input query reads is required (--query)\n";
    exit(1);
}

if(! -f $query_file) {
    print STDERR "Error: the query reads file is not readable (file: $query_file)\n";
    exit(1);
}

if($threads <= 0) {
    print STDERR "Error: you must specify the number of threads to use\n";
    exit(1);
}

# Open a log file to record the commands used
open(LOG, ">$prefix.log");

# Write the version number to the log
my $version_string = `$SGA_BIN --version | head -1`;
print LOG "#Calling variants with: $version_string";
print LOG "#Commands:\n";

#
# Index the reference genome
#

# Make a symlink to the reference genome then build the SGA index files
$reference_file = symlinkReference($reference_file);
indexReference($reference_file);

#
# Index the query reads
#
my $pp_query_file = preprocessReads($query_file);
indexReads($pp_query_file);

# Index the control reads, if they exist
my $pp_control_file = "";
if($control_file ne "") {
    $pp_control_file = preprocessReads($control_file);
    indexReads($pp_control_file);
}

#
# Run the caller
#
if($control_file eq "") {
    callVariantsSingle($pp_query_file, $reference_file);
} else {
    callVariantsCompare($pp_query_file, $pp_control_file, $reference_file);
}

# Cleanup
close(LOG);

#
# Utilities
#

# Call variants for a single file versus the reference
sub callVariantsSingle
{
    my($in_file, $reference) = @_;
    run("$SGA_BIN graph-diff -t $threads --debruijn -x $x --reference $reference --variant $in_file");
}

# Call variants by comparing two sets of reads
sub callVariantsCompare
{
    my($in_query_file, $in_control_file, $reference) = @_;
    run("$SGA_BIN graph-diff -t $threads -x $x --reference $reference --variant $in_query_file --base $in_control_file");
}

# Index a reference genome 
sub indexReference
{
    my($in_file) = @_;

    # Only index the reference if needed
    my ($base) = basename("$in_file", (".fa", ".fasta"));
    my $bwt_file = "$base.bwt";
    if(! -f $bwt_file || isOlderThan($bwt_file, $in_file)) {
        print "Starting to index reference genome\n";
        run("$SGA_BIN index $in_file > /dev/null");
        run("$SGA_BIN gen-ssa $in_file > /dev/null");
    }
}

# Index a set of reads
sub indexReads
{
    my($in_file) = @_;
    my($filename, $directories, $suffix) = fileparse(fileparse($in_file, ".gz"), (".fa", ".fastq", ".fq", ".fasta"));
    my $bwt_file = "$filename.bwt";
    if(! -f $bwt_file || isOlderThan($bwt_file, $in_file)) {
        print "Starting to index reads in file $in_file\n";
        run("$SGA_BIN index -a ropebwt --no-reverse --threads $threads $in_file");
    }
}

# Preprocess reads, compress the output. 
# Returns the name of the preprocessed reads.
sub preprocessReads
{
    my($in_file) = @_;
    my($filename, $directories, $suffix) = fileparse(fileparse($in_file, ".gz"), (".fa", ".fastq", ".fq", ".fasta"));
    if($suffix eq "") {
        print STDERR "Error: unrecognized file suffix ($suffix)\n";
        exit(1);
    }

    my $pp_out = "$filename.pp$suffix.gz";
    if( ! -f $pp_out || isOlderThan($pp_out, $in_file)) {
        run("$SGA_BIN preprocess -o $pp_out $in_file");
    }

    return $pp_out;
}

# Symlink the reference into the working directory
# unless it is already present. Returns the name of the symlinked file
sub symlinkReference
{
    my($reference_file) = @_;
    my($filename) = fileparse($reference_file);

    # Nothing to do if the reference is already in the CWD
    return $reference_file if($filename eq $reference_file);

    # Nothing to do if the link already exists
    return $filename if(-f $filename);

    # Make the link
    my $ret = symlink($reference_file, $filename);
    if($ret ne 1) {
        print STDERR "Unable to make symlink to $reference_file\n";
        exit(1);
    } else {
        return $filename;
    }
}

# Run a single command
sub run
{
    my($cmd) = @_;
    print LOG "$cmd\n";
    my $ret = system($cmd);
    if($ret != 0) {
        print STDERR "Error: command \"$cmd\" returned with non-zero exit code [$ret]\n";
        exit(1);
    }
}

# Returns true if file1 is older than file2
sub isOlderThan
{
    my($file1, $file2) = @_;
    my $mtime1 = stat($file1)->mtime;
    my $mtime2 = stat($file2)->mtime;

    if($mtime1 == "") {
        print STDERR "Error: could not stat $file1\n";
        exit(1);
    }

    if($mtime2 == "") {
        print STDERR "Error: could not stat $file2\n";
        exit(1);
    }

    return $mtime1 < $mtime2;
}

# Print out a help message
sub printUsage
{
    print "sga-call-variants.pl: run the sga alignment-free variant calling pipeline\n";
    print "\nExample of calling variants in a single high depth genome:\n";
    print "sga-call-variants.pl --reference human.grc37.fa --query-reads NA12878.fastq --threads 8\n";
    print "\nExample of calling mutations in a cancer from a sequenced tumor/normal pair\n";
    print "sga-call-variants.pl --reference human.grc37.fa --query-reads tumor.fastq --control-reads normal.fastq --threads 8\n";
    print "\n";
    print "Parameters:\n\n";
    print "        --help                    display this help and exit\n";
    print "        --prefix=STR              prefix the output files with STR [$prefix]\n";
    print "        --threads=N               use N computation threads [required]\n";
    print "        --reference=FILE          the reference genome is in FILE\n";
    print "        --query-reads=FILE        call variants using the reads in FILE\n";
    print "        --control-reads=FILE      when this file is given, do not call variants\n";
    print "                                  in the query reads when they are supported by the\n";
    print "                                  control reads\n";
    print "\n";
    print "Advanced Parameters:\n\n";
    print "        --kmer=K                  use K-mers when discovering variants\n";
    print "        --min-discovery-count=N   require N variant k-mer occurrencess to start assembly\n";
    print "        --min-graph-count=N       require N k-mer occurrences to use in graph\n";
}

# Check whether the programs are installed and executable
sub checkDependency
{
    while(my $program = shift) {
        my $ret = system("/bin/bash -c \"hash $program\"");
        if($ret != 0) {
            print STDERR "Could not find program $program. Please install it or update your PATH.\n";
            print STDERR "Return: $ret\n";
            exit(1);
        }
    }
}
