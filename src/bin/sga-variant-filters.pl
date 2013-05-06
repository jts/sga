#! /usr/bin/perl
# Perform filtering of SGA variant calls

use strict;
use Getopt::Long;
use File::Basename;

my $dbsnp_path = ""; # Filter variants against dbSNP at the given directory
my $sga_file = "";
my $dust_cutoff = 2.0;
my $strand_cutoff = 2.0;
my $hplen_cutoff = 7;

GetOptions("dbsnp=s" => \$dbsnp_path,
           "sga=s" => \$sga_file);

die("The --sga option is mandatory") if($sga_file eq "");
my %non_dbsnp_sites;
my $using_dbsnp = 0;
if($dbsnp_path ne "") {
    $using_dbsnp = 1;
    load_nondbsnp_sites($sga_file, $dbsnp_path);
}

filter_annotations($sga_file);

sub filter_annotations
{
    my($in) = @_;
    my $base = basename($in, ".vcf");
    my $outname = "$base.filters.vcf";
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

        if(scalar(@filter_reason) > 0) {
            $f[6] = join(";", @filter_reason);
        } else {
            $total_out++;
        }
        print OUT join("\t", @f) . "\n";
    }
    
    close(IN);
    close(OUT);

    print "$total_out calls passed all filters\n";
    return $outname;
}

# Build a hash of calls that are NOT present in dbsnp
sub load_nondbsnp_sites
{
    my($in, $path) = @_;
    my $out = "$in.dbsnp_filtered.vcf";
    
    # sort, bgzip and tabix the file
    system("cat $in | vcf-sort > $in.tmp.sorted.vcf");
    system("bgzip -f $in.tmp.sorted.vcf");
    system("vcf tabix -f -p vcf $in.tmp.sorted.vcf.gz");
    
    # find the non-dbsnp sites
    open(SITES, "vcf isec -C $in.tmp.sorted.vcf.gz $path |");
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
