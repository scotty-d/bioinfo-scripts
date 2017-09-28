#!/usr/bin/perl

# vcf_substitute.pl
# changes positions in a genome assembly from reference to variants using a VCF file
# each chromosome in reference FASTA file should be on one line
# caution: reads entire chromosome into memory

if ($#ARGV!=1){
    die "usage: vcf_substitute.pl reference.fasta variants.vcf)";
}

my $ffile=$ARGV[0];
my $vfile=$ARGV[1];
open FASTA, "<$ffile" or die "$!";
  
my $seqname;
my @seq;
    
while(my $fline=<FASTA>){
    chomp $fline;
    if($fline=~/^>/){ # line is the sequence name
	$seqname=substr($fline, 1);
	next;
    } 
    else{ # line is sequence
	@seq=split(//,$fline);
	open VCF, "<$vfile" or die "$!";
	while(my $vline=<VCF>){
	    chomp $vline;
	    if($vline=~/^$seqname\s/){
		my @varr=split(/\t/,$vline);
		my $pos=$varr[1]-1;
		my @refsplit=split(//,$varr[3]);
	        my $reflen=scalar @refsplit;
		my @varsplit=split(//,$varr[4]);
	        my $varlen=scalar @varsplit;
		if($reflen==$varlen && $reflen==1){
	       	# simple substitution
		    $seq[$pos]=$varsplit[0];
		}
	        elsif($reflen<=$varlen){
       		# multiple substitution/insertion
		    foreach my $base (@refsplit) {
			$seq[$pos]=$varsplit[0];
			splice(@varsplit,0,1);
		    } continue { $pos++; }
		    if(scalar @varsplit>0) {
			# insertion; add whatever leftover bases to $pos
			my $str=join("",@varsplit);
			$pos--; # offsets last $pos++ to maintain position
			$seq[$pos] .= $str;
		    }
		 }
	         elsif($reflen>$varlen){
		    # deletion; populate variant with "-" to make it a substitution
		    while($varlen<$reflen){
			$varsplit[$varlen]="-";
			$varlen++;
		    }
		    foreach my $base (@refsplit) {
			$seq[$pos]=$varsplit[0];
			splice(@varsplit,0,1);
		    } continue { $pos++; }
		}
	    }   
	}
        close VCF;
	print ">$seqname\n";
	print join("",@seq),"\n";
    }
}
close FASTA;
