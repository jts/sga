#! /usr/bin/perl
# Perform filtering of SGA variant calls

use strict;
use Getopt::Long;
use File::Basename;
use IPC::Cmd qw[can_run];

my $dbsnp_path = ""; # Filter variants against dbSNP at the given directory
my $sga_file = "";
my $dust_cutoff = 2.0;
my $strand_cutoff = 2.0;
my $hplen_cutoff = 7;
my $depth_cutoff = 0;
my $passed_only = 0;
my $min_af = 0;
my $tumor_bam = "";
my $normal_bam = "";
my $outname = "";
my $extra_dir = ""; # where all the dependencies live

GetOptions("dbsnp=s"       => \$dbsnp_path,
           "sga=s"         => \$sga_file,
           "extra-dir=s"   => \$extra_dir,
           "min-depth=i"   => \$depth_cutoff,
           "min-af=f"      => \$min_af,
           "passed-only"   => \$passed_only,
           "tumor-bam=s"   => \$tumor_bam,
           "normal-bam=s"  => \$normal_bam,
           "outname=s"     => \$outname);

die("The --sga option is mandatory") if($sga_file eq "");

# Set paths to dependencies
my $samtools_bin = "$extra_dir/samtools/samtools";
my $vcfsort_bin = "$extra_dir/vcflib/bin/vcfsort";
my $bgzip_bin = "$extra_dir/htslib/bgzip";
my $bcftools_bin = "$extra_dir/bcftools/bcftools";

# check dependencies
check_prerequisites($samtools_bin,
                    $vcfsort_bin,
                    $bgzip_bin,
                    $bcftools_bin);

my %non_dbsnp_sites;
my $using_dbsnp = 0;
if($dbsnp_path ne "") {
    $using_dbsnp = 1;
    load_nondbsnp_sites($sga_file, $dbsnp_path);
}

perform_filter($sga_file);

sub perform_filter
{
    my($in) = @_;
    my $base = basename($in, ".vcf");
    $outname = "$base.filters.vcf" if $outname eq "";

    open(OUT, ">$outname") || die("Cannot open $outname");
    open(IN, $in) || die("Cannot open $in");
    my $total_out = 0;
    my %seen_hash;

    while(<IN>) {
        # Print header lines
        if(/^#/) {
            print OUT;
            next;
        }

        chomp;
    
        # If this call is exactly the same as a previous call, ignore it
        my @f = split;
        my $key = sprintf("%s.%d.%s.%s", $f[0], $f[1], $f[3], $f[4]);
        if(defined($seen_hash{$key})) {
            next;
        } else {
            $seen_hash{$key} = 1;
        }

        my @filter_reason;
        my ($dust) = (/Dust=(\d+\.\d+)/);
        if($dust > $dust_cutoff) {
            push @filter_reason, "LowComplexity";
        }

        my ($af) = (/AF=(\d+\.\d+)/);
        if($af < $min_af) {
            push @filter_reason, "LowAlleleFreq";
        }

        my ($strand) = (/SB=(\d+\.\d+)/);
        if($strand > $strand_cutoff) {
            push @filter_reason, "StrandBias";
        }

        my ($hplen) = (/HPLen=(\d+)/);
        if($hplen > $hplen_cutoff) {
            push @filter_reason, "Homopolymer";
        }

        if($using_dbsnp && !defined($non_dbsnp_sites{$key})) { 
            push @filter_reason, "dbSNP";
        }

        if($tumor_bam ne "" && $normal_bam ne "") {
            my @bam_filters = filter_bam($_);
            push @filter_reason, @bam_filters;
        }

        if(scalar(@filter_reason) > 0) {
            $f[6] = join(";", @filter_reason);
        } else {
            $total_out++;
        }

        print OUT join("\t", @f) . "\n" if(!$passed_only || scalar(@filter_reason) == 0);
    }
    
    close(IN);
    close(OUT);

    print "$total_out calls passed all filters\n";
    return $outname;
}

sub filter_bam
{
    my($line) = $_;
    my @f = split(' ', $line);

    my $is_snv = length($f[3]) == 1 && length($f[4]) == 1;

    # Call mpileup
    my $cmd = sprintf("$samtools_bin mpileup -r %s:%d-%d %s %s", $f[0], $f[1], $f[1], $tumor_bam, $normal_bam);
    open(S, "$cmd|") || die("mpileup failed");
    my $mpl = <S>; # A single line is returned
    chomp $mpl;
    
    my @g = split(' ', $mpl);
    
    my $t_depth = $g[3];
    my $t_bases = $g[4];
    my $t_qual = $g[5];

    my $n_depth = $g[6];
    my $n_bases = $g[7];
    my $n_qual = $g[8];

    my @filters;
    push @filters, "LowNormalDepth" if($n_depth < 5);

    if($is_snv) {
        # Count the number of times the alt allele shows up in the normal reads
        my $alt_normal_count = 0;
        foreach my $b (split '', $n_bases) {
            ++$alt_normal_count if(uc($b) eq $f[4]);
        }
        print "$f[3] $f[4] $t_bases $n_bases $alt_normal_count\n";
        push @filters, "NormalEvidence" if $alt_normal_count >= 2;

        print "$t_bases $t_qual " . length($t_bases) . " " . length($t_qual) . "\n";

        my @base_array = split('', $t_bases);
        my @qual_array = split('', $t_qual);

        my $bi = 0;
        my $qi = 0;
        my @alt_quals;
        for(;$bi < scalar(@base_array); $bi++) {
            my $u_base = uc($base_array[$bi]);

            if($u_base eq $f[4]) {
                push(@alt_quals, ord($qual_array[$qi]) - 33);
            }

            $qi++ if($u_base eq 'A' || $u_base eq 'C' || $u_base eq 'G' || $u_base eq 'T');
        }
        
        push @filters, "NoAltEvidence" if(scalar(@alt_quals) == 0);

        my @sorted_quals = sort { $a <=> $b } @alt_quals;
        my $mi = scalar(@sorted_quals) / 2;
        my $median = $sorted_quals[$mi];
        if (scalar(@sorted_quals) % 2) {
                $median = $sorted_quals[ $mi ];
        } else {
                $median = ($sorted_quals[$mi-1] + $sorted_quals[$mi]) / 2;
        }
        print join(',', @sorted_quals) . "\n"; 
        print "Median: $median\n";
        push @filters, "LowQuality" if($median <= 15);
    }
    return @filters;
}

# Build a hash of calls that are NOT present in dbsnp
sub load_nondbsnp_sites
{
    my($in, $path) = @_;
    my $out = "$in.dbsnp_filtered.vcf";
    
    # sort, bgzip and tabix the file
    system("$vcfsort_bin $in > $in.tmp.sorted.vcf");
    system("$bgzip_bin -f $in.tmp.sorted.vcf");
    system("$bcftools_bin index -f $in.tmp.sorted.vcf.gz");
    
    # find the non-dbsnp sites
    open(SITES, "$bcftools_bin isec -C $in.tmp.sorted.vcf.gz $path |");
    while(<SITES>) {
        chomp;
        my @fields = split;
        my $site_key = "$fields[0].$fields[1].$fields[2].$fields[3]";
        $non_dbsnp_sites{$site_key} = 1;
    }
    close(SITES);

    # Cleanup tmp
    unlink("$in.tmp.sorted.vcf");
    unlink("$in.tmp.sorted.vcf.gz");
    unlink("$in.tmp.sorted.vcf.gz.tbi");
}

# check that each program can be run
sub check_prerequisites
{
    my (@programs) = @_;

    foreach my $program (@programs) {
        if(!can_run($program)) {
            print STDERR "Error: could not find program $program.\n";
            print STDERR "Please set the --extra option to the sga-extra directory.\n";
            exit(1);
        }
    }
}
