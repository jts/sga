#! /usr/bin/perl

use strict;
use Getopt::Long;

my $n = 5;
my $k = 99;
my $minLength = 200;
my $prefix = "";
my $numThreads = 1;
my $mind = -99; # minimum gap size to pass to abyss
my $mina = 30; # minimum alignment length to pass to abyss

# Filter the abyss distance est histogram to remove insert sizes
# with fewer than hist_min data points
my $hist_min = 3;

GetOptions("prefix=s" => \$prefix,
           "k=i"      => \$k,
           "n=i"      => \$n,
           "m=i"      => \$minLength,
           "mind=i"   => \$mind,
           "mina=i"   => \$mina,
           "t=i"      => \$numThreads);

checkDependency("abyss-fixmate");
checkDependency("DistanceEst");

my $bFail = 0;
if($prefix eq "")
{
    print("Error: A prefix must be specified\n");
    $bFail = 1;
}
my $bamFile = $ARGV[0];
if($bamFile eq "")
{
    print "Error: A bam file must be provided\n";
    $bFail = 1;
}

if($minLength <= 0)
{
    print "Error: minimum contig length (-m) must be greater than 0\n";
    $bFail = 1;
}

if($bFail)
{
    usage();
    exit(1);
}

my $cmd;

# fix-mate info
# use reformat.sh from BBmap bundle to convert from v1.4 CIGAR strings into v1.3, alternatively one could have used RTG-CORE (unlike RTG-TOOLS which is the stripped version lacking sammerge feature)
my $mypath = `which reformat.sh`;
my $cigar_convertor_cmd = "";
my $rtg_core_distdir = "";
if ($mypath ne '')
{
    $cigar_convertor_cmd = "reformat.sh sam=1.3 in=stdin.sam out=stdout.sam |";
}
else
{
    $mypath = `which rtg`;
    if ($mypath ne '')
    {
        $rtg_core_distdir = `dirname $mypath`;
        $cigar_convertor_cmd = "java -Djava.library.path=${rtg_core_distdir} -XX:ErrorFile=./hs_err_pid'$$'.log -Dreferences.dir=${rtg_core_distdir}/references -Dmodels.dir=${rtg_core_distdir}/models -Dtalkback=false -Xmx38g -jar ${rtg_core_distdir}/RTG.jar sammerge ${bamFile} --legacy-cigars --exclude-duplicates --exclude-unmapped --exclude-unplaced --output=- --no-gzip |";
    }
}
$cmd = "samtools view -h -F 0x800 -F 0x200 -F 0x400 -F 0x100 $bamFile | $cigar_convertor_cmd abyss-fixmate -h $prefix.tmp.hist | samtools view -h -Sb -o $prefix.diffcontigs.bam";
runCmd($cmd);

# DistanceEst is very slow when the learned insert size distribution is very wide
# To work around this, remove singleton data points from the histogram
runCmd("awk \'\$2 >= $hist_min\' $prefix.tmp.hist > $prefix.hist");

# sort 
$cmd = "samtools sort -\@$numThreads -o $prefix.diffcontigs.sorted.bam $prefix.diffcontigs.bam";
runCmd($cmd);

# distance est
my $mind_opt = "--mind $mind";
$cmd = "DistanceEst -s $minLength $mind_opt -n $n -k $k -j $numThreads -o $prefix.de $prefix.hist -l $mina $prefix.diffcontigs.sorted.bam";
runCmd($cmd);

sub usage
{
    print "sga-bam2de.pl - Make a distance estimate file from a bam file of reads aligned to contigs\n";
    print "sga-bam2de.pl -n N --prefix lib300 lib300.bam\n";
    print "Options:\n";
    print "                -n N             Minimum number of pairs required to consider two contigs linked\n";
    print "                -m LEN           Only find links between contigs with length at least LEN bp (default: 200)\n";
    print "                -t NUM           Use NUM threads for computing the distance estimates\n";
    print "                --prefix NAME    Use NAME as the prefix for the outfiles\n";
    print "                --mind D         Set the minimum distance estimate to test to be D. This should be a negative\n";
    print "                                 number if contigs are expected to overlap. Defaults to -99bp.\n";
    print "                --mina N         Set the minimum alignment length to be N.\n";
}

sub runCmd
{
    my($c) = @_;
    print "$c\n";
    system($c) and warn "Error: ", chomp($c), " failed with: $!\n" and die;
}

# Check whether the programs are found and executable
sub checkDependency
{
    while(my $program = shift) {
        my $ret = system("/bin/bash -c \"hash $program\"");
        if($ret != 0) {
            print STDERR "Could not find program $program: $!. Please install it or update your PATH.\n";
            print STDERR "Return: $ret\n";
            exit(1);
        }
    }
}
