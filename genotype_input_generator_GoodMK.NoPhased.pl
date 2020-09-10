#!/usr/bin/perl;

### Program to change genotype code (0|0)(1|1) to alte/ref.
## Usage

#perl genotype_input_generator.pl vcf_file.vcf




#use strict;
use warnings;

$file = $ARGV[0];

$outF="matrix.$file.NoPhased.mtr";

open (OUT, "+>$outF") or die $!;

open(FF,$file);

$count1=0;

while(<FF>) {

	$li=$_;
	
	if($li =~ /^[#]/) { next;}
	
	@line =split("\t",$li);
	
	$totalId=$#line-8;
	chomp @line;	
	
	$chr=$line[0];
	$pos=$line[1];
	$ref=$line[3];
	$alt=$line[4];

	$arraout[0]="CP.$chr.$pos\t";

	##### Generating array 

	for($x=9; $x < $#line+1; ++$x) {
					$count2=$x-8;
					#print "VOY $x\n";	
					#print "REF $ref ALT $alt\n";
					#print "LL <$line[$x]>\n";
					if ($line[$x] eq "0|0") { $arraout[$count2]="$ref/$ref\t";}
					if ($line[$x] eq "0|1") { $arraout[$count2]="$ref/$alt\t";}
 					if ($line[$x] eq "1|0") { $arraout[$count2]="$ref/$alt\t";}
					if ($line[$x] eq "1|1") { $arraout[$count2]="$alt/$alt\t";}
					if($line[$x] eq "0") 	{$arraout[$count2]="$ref\t";}
					if($line[$x] eq "1")    {$arraout[$count2]="$alt\t";}
					#print "OO $arraout[$count2]\n";

					 #<STDIN
			}

	print OUT "@arraout\n"; 

}


close OUT;













