#!/usr/bin/perl;

### Program to convert VCF genotype codes (e.g., 0|0, 1|1) to REF/ALT allele format
### Usage:
# perl genotype_input_generator.pl vcf_file.vcf

# Turn on warnings to catch potential issues (but 'strict' is commented out)
# use strict;
use warnings;

# === Input file from command-line argument ===
$file = $ARGV[0];

# === Define output filename ===
$outF = "matrix.$file.NoPhased.mtr";

# === Open output file for writing (overwrite mode) ===
open (OUT, "+>$outF") or die $!;

# === Open input VCF file ===
open(FF, $file);

$count1 = 0;

# === Main loop over each line in the VCF ===
while (<FF>) {
    $li = $_;

    # Skip header lines starting with #
    if ($li =~ /^[#]/) { next; }

    # Split the VCF line into fields
    @line = split("\t", $li);
    chomp @line;

    # Total number of individuals (subtract metadata columns)
    $totalId = $#line - 8;

    # Extract variant information
    $chr = $line[0];
    $pos = $line[1];
    $ref = $line[3];
    $alt = $line[4];

    # Start building output row: e.g., CP.1.123456<TAB>
    $arraout[0] = "CP.$chr.$pos\t";

    # === Convert genotypes for each sample ===
    for ($x = 9; $x < $#line + 1; ++$x) {
        $count2 = $x - 8;

        # Translate genotype code to allele format
        if ($line[$x] eq "0|0") { $arraout[$count2] = "$ref/$ref\t"; }
        if ($line[$x] eq "0|1") { $arraout[$count2] = "$ref/$alt\t"; }
        if ($line[$x] eq "1|0") { $arraout[$count2] = "$ref/$alt\t"; }
        if ($line[$x] eq "1|1") { $arraout[$count2] = "$alt/$alt\t"; }

        # Haploid / unphased formats (e.g. mitochondrial)
        if ($line[$x] eq "0")   { $arraout[$count2] = "$ref\t"; }
        if ($line[$x] eq "1")   { $arraout[$count2] = "$alt\t"; }
    }

    # Print converted line to output
    print OUT "@arraout\n";
}

# === Close output file ===
close OUT;












