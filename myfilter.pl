#!/usr/bin/perl

#$sample= "sample_11";
#$input = $sample.".goodQual.snpSift.vcf";
#$output = $sample.".filtered.vcf";


$sample = $ARGV[0];
$input = $ARGV[1];
$output = $ARGV[2];

print $input."\n".$output."\n";
open (FIL,$input);
open (OUT,">$output");

while(<FIL>)
{
	$l = $_;
	@line_split = split("\t",$l);
	if($line_split[0] =~ /^#/) # print header line
	{
		print OUT $l;
		next;
	}
	if($line_split[0] =~ /^[0-9]/ || $line_split[0] =~ /^[XY]/) # filter chr1-chr22,X,Y
	{
		@info = split(";",$line_split[7]);	
		foreach (<@info>) # get MQ0,DP region
		{
			$info_split = $_;
			if($info_split =~ /^DP/) {$DP = $info_split;next;}
			if($info_split =~ /^MQ0/) {$MQ0 = $info_split;next;}
		}
		@DP_num = split("=",$DP);
		@MQ0_num = split("=",$MQ0);
		if( $MQ0_num[1] <= 4 || $MQ0_num[1] <= 0.1*$DP_num[1]) # filter MQ0 <= 4 or MQ0 <= 0.1*DP
		{
			@format_split = split(":",$line_split[-1]);
			if( $format_split[0] eq "0/1") # filter het only
			{
				@AD = split(",",$format_split[1]);
				if ( $AD[0]/($AD[1]+$AD[0]) <= 0.75) {print OUT $l;} # filter heterozygote Allele balance <= 0.75
			}
			else {print OUT $l}; # print hom
		} 
	}
}

close (FIL);
close (OUT);
