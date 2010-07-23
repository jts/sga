package ParseDE;
use strict;

sub parseDERecord
{
    my ($line) = @_;
    chomp $line;
    my @fields = split(' ', $line);
    my $base = shift @fields;
    my @records;
    my $dir = 0;
    foreach my $linkRec (@fields)
    {
        if($linkRec eq ";")
        {
            ++$dir;
            next;
        }
        #print "B $base $linkRec\n";
        my ($tag, $dist, $n, $sd) = split(",", $linkRec);
        my $id = substr($tag, 0, -1);
        my $strand = substr($tag, -1, 1);

        my $rec;
        $rec->{ctg} = $id;
        $rec->{strand} = $strand;
        $rec->{dist} = $dist;
        $rec->{n} = $n;
        $rec->{sd} = $sd;
        $rec->{dir} = $dir;
        push @records, $rec;
    }
    return ($base, \@records);
}

1;
