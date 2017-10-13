#!usr/bin/perl
#random_permutation_analysis
use warnings;

print "\nEnter number if iterations\n\n";
$iter = <STDIN>;
chomp $iter;
print "iter is $iter\n";

open INPUT2, 'c:\Users\tsouthal\Dropbox\Southall lab files\Projects\Cell specific chromatin accessibility profiling\Peak calling\ATAC-Seq_Peaks_rel5.bed' or die "\nCan't open file\n";   ##### change this for different file
@enhancers = <INPUT2>; 
$enhancersnum = @enhancers;

$iter_no = 0;

while($iter_no < $iter){



open (INPUT, 'C:\Users\tsouthal\Dropbox\Southall lab files\Projects\Cell specific chromatin accessibility profiling\Peak calling\GG_rep1_rel5_0.01FDR_50counts.gff') or die "\nCan't open file!\n";   ##### change this for different file
@file_array = <INPUT>;
$num = @file_array;



$ln1 = 0;

while($ln1 < $num){
@col = split(/\t/,$file_array[$ln1]);

$size = ($col[4] - $col[3]);

if($col[0] =~ m/^chr2L$/){push (@chr2L, "$size");}				
if($col[0] =~ m/^chr2R$/){push (@chr2R, "$size");}	
if($col[0] =~ m/^chr3L$/){push (@chr3L, "$size");}	
if($col[0] =~ m/^chr3R$/){push (@chr3R, "$size");}	
if($col[0] =~ m/^chr4$/){push (@chr4, "$size");}	
if($col[0] =~ m/^chrX$/){push (@chrX, "$size");}		

$ln1 = $ln1 + 1;
}

print "\nchr2L first line is $chr2L[0]\n";

################################
$ln1 = 0;
$num_2L = @chr2L;


while($ln1 < $num_2L){

$limit = (22968378 - $chr2L[$ln1]);

$random_start = int(rand($limit));
$end = $random_start + $chr2L[$ln1];

push (@peaks, "chr2L\t$random_start\t$end");
#print OUTPUT "chr2L\t.\trandom\t$random_start\t$end\t1\t.\t.\t.\n";
#print "chr2L\t.\trandom\t$random_start\t$end\t1\t.\t.\t.\n";

$ln1 = $ln1 + 1;
}

####################################

$ln1 = 0;

$num_2R = @chr2R;

while($ln1 < $num_2R){

$limit = (21142058 - $chr2R[$ln1]);

$random_start = int(rand($limit));
$end = $random_start + $chr2R[$ln1];

push (@peaks, "chr2R\t$random_start\t$end");
#print OUTPUT "chr2R\t.\trandom\t$random_start\t$end\t1\t.\t.\t.\n";

$ln1 = $ln1 + 1;
}

####################################
	
$ln1 = 0;

$num_3L = @chr3L;

while($ln1 < $num_3L){

$limit = (24530920 - $chr3L[$ln1]);


$random_start = int(rand($limit));
$end = $random_start + $chr3L[$ln1];

push (@peaks, "chr3L\t$random_start\t$end");
#print OUTPUT "chr3L\t.\trandom\t$random_start\t$end\t1\t.\t.\t.\n";

$ln1 = $ln1 + 1;
}

####################################
	
$ln1 = 0;

$num_3R = @chr3R;

while($ln1 < $num_3R){

$limit = (27894456 - $chr3R[$ln1]);

$random_start = int(rand($limit));
$end = $random_start + $chr3R[$ln1];

push (@peaks, "chr3R\t$random_start\t$end");
#print OUTPUT "chr3R\t.\trandom\t$random_start\t$end\t1\t.\t.\t.\n";

$ln1 = $ln1 + 1;
}

####################################	
	
$ln1 = 0;

$num_4 = @chr4;

while($ln1 < $num_4){

$limit = (1284582 - $chr4[$ln1]);

$random_start = int(rand($limit));
$end = $random_start + $chr4[$ln1];

push (@peaks, "chr4\t$random_start\t$end");
#print OUTPUT "chr4\t.\trandom\t$random_start\t$end\t1\t.\t.\t.\n";

$ln1 = $ln1 + 1;
}

####################################	
	
$ln1 = 0;

$num_X = @chrX;

while($ln1 < $num_X ){

$limit = (22421853 - $chrX[$ln1]);

$random_start = int(rand($limit));
$end = $random_start + $chrX[$ln1];

push (@peaks, "chrX\t$random_start\t$end");
#print OUTPUT "chrX\t.\trandom\t$random_start\t$end\t1\t.\t.\t.\n";

$ln1 = $ln1 + 1;
}

####################################

close INPUT;
@file_array = ();
@chr2L = ();			
@chr2R = ();	
@chr3L = ();	
@chr3R = ();	
@chr4 = ();	
@chrX = ();

$peaksnum = @peaks;

$state = 0;

$ln1 = 0;

while($ln1 < $enhancersnum){
	@col = split(/\t/,$enhancers[$ln1]);
	$chrom = $col[0]; $estart = $col[1]; $eend = $col[2]; chomp $eend; #print "\nChrom is $chrom, enhancer start is $estart and enhancer end is $eend";
	
$ln2 = 0;

while($ln2 < $peaksnum){
	@col2 = split(/\t/,$peaks[$ln2]);
	
	if(($chrom eq $col2[0]) && ((($estart >= $col2[1]) && ($estart <= $col2[2])) || (($eend >= $col2[1]) && ($eend <=$col2[2])) || (($estart <= $col2[1]) && ($eend >= $col2[2])))){
		
		$state = $state + 1; $ln2 = $peaksnum;
		
		#print "\nOverlap found! Peak chrom is $col2[0], peak start is $col2[3] and peak end is $col2[4]"; 
	}
		$ln2 = $ln2 + 1;

}
$ln1 = $ln1 + 1;
}

$percent = ($state / $enhancersnum) * 100;


open LOG, '>> d:\monte_carlo_for_ATAC.txt';   ####### change this for different output location
print LOG "$percent\n";
close LOG; 

print "\nSimulation number $iter_no done - percentage overlap was $percent";

@peaks = ();

$iter_no = $iter_no + 1;
}

exit;
	
	