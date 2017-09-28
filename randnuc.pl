#!/usr/bin/perl

# randnuc.pl
# generates a random sequence of given length based on a first-order Markov model
# see sample probability files for formatting

if($#ARGV!=2){
    die "usage: randnuc.pl sequence_length nucleotide_prob_file dinucleotide_prob_file\n";
}
my $count=$ARGV[0];

my $nfile=$ARGV[1];
open NFILE, "<$nfile" or die '$!';
my %nfreqs;
for my $i (1..4) { # A, C, G, T
    my $nline=<NFILE>;
    chomp $nline;
    my @nsplit=split(/\t/, $nline);
    $nfreqs{$nsplit[0]}=$nsplit[1];
}
close NFILE;

my $dnfile=$ARGV[2];
open DNFILE, "<$dnfile" or die '$!';
my %dnfreqs;
for my $i (1..16) { # AA, AC, AG, AT, CA, CC ...
    my $dnline=<DNFILE>;
    chomp $dnline;
    my @dnsplit=split(/\t/, $dnline);
    $dnfreqs{$dnsplit[0]}=$dnsplit[1];
}
close DNFILE;

sub firstNuc{ # choose first nucleotide randomly
    my $num=rand();
    if($num<$nfreqs{'A'}){return 'A';}
    elsif($num<$nfreqs{'A'}+$nfreqs{'C'}){return 'C';}
    elsif($num<$nfreqs{'A'}+$nfreqs{'C'}+$nfreqs{'G'}){return 'G';}
    else{return 'T';} # probability value of T is not used
}

sub nextNuc{ # choose next nucleotide randomly based on the last
    my $lastNuc=$_[0];
    my $num=rand();
    if($num<$dnfreqs{$lastNuc.'A'}){return 'A';}
    elsif($num<$dnfreqs{$lastNuc.'A'}+$dnfreqs{$lastNuc.'C'}){return 'C';}
    elsif($num<$dnfreqs{$lastNuc.'A'}+$dnfreqs{$lastNuc.'C'}+$dnfreqs{$lastNuc.'G'}){return 'G';}
    else{return 'T';} # probability values of AT, CT, GT, and TT are not used
}

my $letter=firstNuc();
for my $i (1..$count){
    print $letter;
    $letter=nextNuc($letter);
}
print "$letter\n";
