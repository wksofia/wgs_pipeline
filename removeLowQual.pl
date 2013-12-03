#!/usr/bin/perl
# data:		2012/09/26
# version:	1.0
# author:	Caibin,Lifengyu

$inFile = shift;
$outFile = shift or die "Usage: vcfFile output\n";


open IN, "<$inFile" or die "COF $inFile\n";
open OUT, ">$outFile" or die "CCF $outFile\n";

while ($line = <IN>) {
	if ($line=~/^#/){ #annotation
		print OUT $line;
		next;
	}
	
	 if (&checkQual($line)){
		print OUT $line if (&checkMut($line) ne 'ref');
	}

}


close OUT;
close IN;

sub checkQual{
$line = shift;
@l = split "\t" , $line;

die "Formate error:
$line
" if (@l<10);

return 0 if $l[6] ne 'PASS';
return 0 if !($l[7]=~/DP/);

@info = split ';',$l[7];
for $i(@info){
	if ($i=~/DP=/){#DP<5, bad quality
		($dp,$score)=split '=',$i;
		return 0 if $score < 5;
	}elsif ($i=~/MQ=/) {#MQ<20, bad quality
		($mq,$score)=split '=',$i;
		return 0 if $score < 20;
#	}elsif ($i=~/VQSLOD=/) {#VQSLOD<0, bad quality
#		($vq,$score)=split '=',$i;
#		return 0 if $score < 0;
	}else {
		next:
	}
}

@f = split ':',$l[8];
@s = split ':',$l[9];

$t = '';
for ($i=0;$i<@f;$i++){
	$t = $s[$i] if ($f[$i] eq 'AD');
}

@support = split ',',$t;

shift @support;
for (@support){#supported reads<3, bad quality
	return 0 if $_<3;
}

return 1;
}

sub checkMut{#if one site has INFO of DP, it should has GT. So wouldn't check it again here
$line = shift;
@l = split "\t",$line;

@f = split ':',$l[8];
@s = split ':',$l[9];

$t = '';
for ($i=0;$i<@f;$i++){
	$t = $s[$i] if ($f[$i] eq 'GT');
}

($a,$b) = split '/',$t;

return 'ref' if (($a+$b) == 0);
return 'homo' if ($a == $b);
return 'het';
}
