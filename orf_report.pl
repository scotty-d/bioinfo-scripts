#!/usr/bin/perl

# orf_report.pl
# prints the longest ORF in each sequence of a FASTA file (e.g. transcripts)
# includes nucleotide start and end positions
# uses the standard genetic code, start codon ATG/M
# each FASTA sequence must all be on one line

if($#ARGV!=0){ # only one file allowed at a time
    die "usage: orf_report.pl sequence_file.fasta\n";
}

my $file=$ARGV[0];
open FASTA, "<$file" or die "$!";

# standard genetic code
my %codon = ( TTT=>'F', TTC=>'F', TTA=>'L', TTG=>'L',
	      CTT=>'L', CTC=>'L', CTA=>'L', CTG=>'L',
	      ATT=>'I', ATC=>'I', ATA=>'I', ATG=>'M',
	      GTT=>'V', GTC=>'V', GTA=>'V', GTG=>'V',
	      TCT=>'S', TCC=>'S', TCA=>'S', TCG=>'S',
	      CCT=>'P', CCC=>'P', CCA=>'P', CCG=>'P',
	      ACT=>'T', ACC=>'T', ACA=>'T', ACG=>'T',
	      GCT=>'A', GCC=>'A', GCA=>'A', GCG=>'A',
	      TAT=>'Y', TAC=>'Y', TAA=>'*', TAG=>'*',
	      CAT=>'H', CAC=>'H', CAA=>'Q', CAG=>'Q',
	      AAT=>'N', AAC=>'N', AAA=>'K', AAG=>'K',
	      GAT=>'D', GAC=>'D', GAA=>'E', GAG=>'E',
	      TGT=>'C', TGC=>'C', TGA=>'*', TGG=>'W',
	      CGT=>'R', CGC=>'R', CGA=>'R', CGG=>'R',
	      AGT=>'S', AGC=>'S', AGA=>'R', AGG=>'R',
	      GGT=>'G', GGC=>'G', GGA=>'G', GGG=>'G' );

# convert from codon to aa
sub Translate{
    my $nuc=$_[0];
    my $aa;
    my $num=0;
    my $len=length($nuc);
    while ( $num < $len ) {
	my $sstr=uc substr($nuc,$num,3);
	$aa .= $codon{$sstr};
	$num+=3;
    }
    $aa;
}

my $seqname;
while(my $line=<FASTA>){
    chomp $line;
    if($line=~/^>/){ # name line
	# separate the name from the initial ">"
	$seqtemp=substr($line,1);
	@seqline=split(/  */, $seqtemp);
	$seqname=$seqline[0];
    }
    else{ # sequence line
	my @longest=(0,0,0);
	my @lstart=(0,0,0);
	my @lend=(0,0,0);
	my @lpep=("","","");
	print "$seqname\t";

	# translate each forward reading frame completely, then look for ORFs
	for(my $i=0; $i<=2; ++$i){
	    my $pepline=&Translate($line);
	    my @pepsplit=split(//,$pepline);

	    my $pep=""; # peptide sequence
	    my $run=0; # peptide length
	    my $start=$i; # ORF start position, based on reading frame

	    # scan through amino acids for start codon (M)
	    for my $num (0..$#pepsplit) {
		my $amino=$pepsplit[$num];

		if($run==0 && $amino eq 'M'){ # start codon
		    $run=1; # this is the beginning of an ORF
		    $start=($num+1)*3+$i;
                    # $num is the array index, which is the aa position minus 1.
		    # The nucleotide position is thus ($num+1)*3, offset by the reading frame $i.
		    $pep.=$amino; # i.e. M
		}

		elsif($run>0 && $amino eq '*'){ # stop codon (end of ORF)
		    if($run>$longest[$i]){
			$longest[$i]=$run;
			$lstart[$i]=$start;
			$lend[$i]=($num+1)*3+$i;
			$lpep[$i]=$pep;
		    }
		    $run=0; # end the run
		    $pep="";
		}

		elsif($run>0){ # other amino acid within ORF
		    $run++; # continue the run
		    $pep.=$amino;
		}
	    }
	    # continue with next reading frame
	    my $temp=substr($line,1);
	    $line=$temp;
	}

	# choose the longest ORF from the three reading frames
	# initialize to first frame
	my $longstr="$lstart[0]\t$lend[0]\t$lpep[0]\t";
	my $longestest=$longest[0]; # "longest of the longest"
	for my $j (1 .. 2) {
	    if($longest[$j]>$longestest){
		$longestest=$longest[$j];
		$longstr="$lstart[$j]\t$lend[$j]\t$lpep[$j]\t";
	    }
	}
	print "$longstr\n";
    }
}
