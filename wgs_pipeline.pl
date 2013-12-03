# ============== info ==============
# code by ugoodlfy
# 1st version finish date: 2013/04/07
# ============== info ==============


use strict;
use Getopt::Long;



my($configFile)="wgs_pipeline.config";
GetOptions(
	'config=s' => \$configFile
);



my(%config)=getConfigHashTable($configFile);

our($sample)=$config{"Sample_Name"};
our($inputFolder)=$config{"Input_Folder"};
our($outputFolder)=$config{"Output_Folder"};
our($refSequence)=$config{"Reference_Sequence"};
our($refDbsnp)=$config{"Reference_Dbsnp"};
our($bwa)=$config{"Bwa"};

our($bwa_threads)=$config{"bwa_threads"};
our($bwa_trim)=$config{"bwa_trim"};

our($samtools)=$config{"Samtools"};
our($picardMerge)=$config{"Picard_MergeBam"};

our($Picard_MergeBam_validationStringency)=$config{"Picard_MergeBam_validationStringency"};
our($Picard_MergeBam_assumeSorted)=$config{"Picard_MergeBam_assumeSorted"};
our($Picard_MergeBam_useThreading)=$config{"Picard_MergeBam_useThreading"};

our($picardDuplicateRemove)=$config{"Picard_MarkDuplicates"};
our($picardDuplicateRemove_maxSequencesForDiskReadEndsMap)=$config{"Picard_MarkDuplicates_maxSequencesForDiskReadEndsMap"};
our($picardDuplicateRemove_removeDuplicates)=$config{"Picard_MarkDuplicates_removeDuplicates"};
our($picardDuplicateRemove_assumeSorted)=$config{"Picard_MarkDuplicates_assumeSorted"};

our($picardSummary)=$config{"Picard_CollectSummary"};
our($Queue)=$config{"Queue"};
our($DataProcessingPipeline)=$config{"Queue_DataProcessingPipeline"};
our($UnifiedGenotyper)=$config{"Queue_UnifiedGenotyper"};

our($gatk)=$config{"GATK"};
our($gatk_SnpRecalHapmap)=$config{"GATK_SnpRecalHapmap"};
our($gatk_SnpRecalOmni)=$config{"GATK_SnpRecalOmni"};
our($gatk_SnpRecalDbsnp)=$config{"GATK_SnpRecalDbsnp"};
our($gatk_IndelRecalMills)=$config{"GATK_IndelRecalMills"};

our($AnnovarDb)=$config{"Annovar_db"};

our($AnnovarConvertAnn)=$config{"Annovar_ConvertAnn"};
our($Annovar_Annotation)=$config{"Annovar_Annotation"};
our($snpSift)=$config{"SnpSift"};
our($vcf2bco)=$config{"Vcf2Bco"};
our($removeLowQual)=$config{"Filter_RemoveLowQual"};
our($myfilter)=$config{"Filter_myFilter"};
our($email)=$config{"Email"};

#print "$sample
#$inputFolder
#$outputFolder
#$refSequence
#$refDbsnp
#$bwa
#$bwa_threads
#$bwa_trim
#$samtools
#$picardMerge
#$Picard_MergeBam_validationStringency
#$Picard_MergeBam_assumeSorted
#$Picard_MergeBam_useThreading
#$picardDuplicateRemove
#$picardDuplicateRemove_maxSequencesForDiskReadEndsMap
#$picardDuplicateRemove_removeDuplicates
#$picardDuplicateRemove_assumeSorted
#$picardSummary
#$Queue
#$DataProcessingPipeline
#$UnifiedGenotyper
#$gatk
#$gatk_SnpRecalHapmap
#$gatk_SnpRecalOmni
#$gatk_SnpRecalDbsnp
#$gatk_IndelRecalMills
#$AnnovarDb
#$AnnovarConvertAnn
#$Annovar_Annotation
#$snpSift
#$vcf2bco
#$removeLowQual
#$myfilter
#$email
#";

print "All Start:";
print scalar localtime(),qq(\n);
print "\n---------------------------- All Start ----------------------------\n";

my(@waitingList)=("1");
my(@tmp);
our($projectFolder)="$outputFolder/$sample/"; # folder of the new sample


if(checkFolder($projectFolder)){mkdir($projectFolder,0755);makeStruture($projectFolder);}


chdir($projectFolder); #move to project folder
runbwaAln();
runbwaSampe();
runsamtoolsGenerateBam();
runpicardMerge();
runpicardDuplicateRemove();
@waitingList=(@waitingList,runpicardSummary());
runQueue();
@waitingList=(@waitingList,rungatkReadDepth());
rungatkGT();
my(@indedRecal)=rungatkIndelRecalModel();
rungatkSnpRecalModel();
rungatkApplySnpRecal();
checkJobFinish2(@indedRecal);
rungatkApplyIndelRecal();

@waitingList=(@waitingList,runAnnovar("$projectFolder/variant_calling/project.$sample.raw.vcf","$projectFolder/annotation/raw/"));
@waitingList=(@waitingList,runvcf2bco("$projectFolder/variant_calling/project.$sample.raw.vcf","$projectFolder/vcf2bco/raw/"));
@waitingList=(@waitingList,gatkEvaluation("$projectFolder/variant_calling/project.$sample.raw.vcf"));
runFilter("$projectFolder/variant_calling/$sample.recal.vcf","$projectFolder/filter/");
@waitingList=(@waitingList,runAnnovar("$projectFolder/filter/$sample.filtered.vcf","$projectFolder/annotation/filtered/"));
@waitingList=(@waitingList,runvcf2bco("$projectFolder/filter/$sample.filtered.vcf","$projectFolder/vcf2bco/filtered/"));
@waitingList=(@waitingList,gatkEvaluation("$projectFolder/filter/$sample.filtered.vcf"));

checkJobFinish2(@waitingList);


collectSummary();
#############change here, greeting file delete, and write "Your sample is finished!" 
system("mail -s $sample.report -a $sample.Summary $email < /wa/ugoodlfy/old/graduation_project/bak/greeting"); # tell me if finish
#system("echo Your sample analysis is finished| mail -s $sample.report -a $sample.Summary $email");
print "All finish:";
print scalar localtime(),qq(\n);
print "\n---------------------------- All Finish ----------------------------\n";
sub getConfigHashTable{
	my($configFile)=$_[0];
	my(%configHashTable);
	my(@tmp);
	my($x);
	open(CON,"<$configFile");
	while($x=<CON>){
		chomp($x);
		if($x !~ /^#/){
			@tmp=split("\t",$x);
			$configHashTable{$tmp[0]}=$tmp[1];
		}	
	}

	close(CON);	
	#foreach $x (keys(%configHashTable)) {print $x."\n";} # test !!!
	return(%configHashTable);
}
sub checkFolder{
	my($folder)=$_[0];
	if( -e $folder){
		print "$folder already exist! Rename the existed folder or rename the new sample!\n";
		exit 1; # the project folder already exist, exit with error
	}
	else{
		return 1;# pass check successfully
	}
	
}

sub makeStruture{
	my($base)=$_[0];

	mkdir("$base/bam",0755);
	mkdir("$base/bam/bwa",0755);
	mkdir("$base/bam/picard",0755);
	mkdir("$base/bam/samtools",0755);
	mkdir("$base/bam/Queue",0755);

	mkdir("$base/annotation",0755);
	mkdir("$base/annotation/raw",0755);
	mkdir("$base/annotation/filtered",0755);

	mkdir("$base/analysis",0755);
	mkdir("$base/analysis/readDepth",0755);
	mkdir("$base/analysis/rsVSnors",0755);

	mkdir("$base/vcf2bco",0755);
	mkdir("$base/vcf2bco/raw",0755);
	mkdir("$base/vcf2bco/filtered",0755);

	mkdir("$base/variant_calling",0755);
	mkdir("$base/filter",0755);
}

sub runbwaAln{
	#enter the specific path
#	print "\n========================================================\n";
#	print "runbwaAln start:";
#	print scalar localtime(),qq(\n);
	print "runbwaAln running...\n";
	chdir("$projectFolder/bam/bwa");
	
	# get fastq.gz file list in inputFolder
	opendir(TEMP,$inputFolder) || die "Cannot open $inputFolder";
	my(@files) = readdir TEMP;
	closedir(TEMP);
	my($x);
	my(@fqList);
	foreach $x (@files){
		push(@fqList,$x) if($x =~ /.*fastq.gz/); # input file in input folder
	}

	# qsub all fq files
	my($fq);
	my(@attr);
	my($command);
	unlink("bwaAln.jobID") if (-e "bwaAln.jobID"); # check if bwaAln.jobID exist , revise every time
	foreach $fq (@fqList){ 
		@attr=("aln","-t", $bwa_threads,"-q", $bwa_trim,$refSequence,$inputFolder.$fq);
		$command=join(" ","qsub -b y -wd \$PWD -e $fq.log -o $fq.sai",$bwa,@attr,">> bwaAln.jobID");
		my($recode)=system($command);
#		print "$recode\n$command\n";
	}


	# Check job status
	checkJobFinish("bwaAln.jobID");


	#go back to project path
	chdir($projectFolder);
#	print "runbwaAln end:";
#	print scalar localtime(),qq(\n);
}

sub extractJobIDList{
	my($jobListFile)=$_[0]; #input a .jobID file
	open(JOB,"<$jobListFile");# get all jobs ID and put them into @jobID
	my(@jobList)=<JOB>;
	close(JOB);
	my(@jobID);
	my(@f);
	my($x);
	my($job);
	my($info);
	
	#get job ID list
	foreach $x (@jobList){
		@f=split(" ",$x);
		push(@jobID,$f[2]);# the 3 field in the line is jobID
	}
	return @jobID;
}

sub checkJobFinish2{ # the function accept the a array as input 
       # Check job status
	my(@jobList)=@_;
	my($x);
	my($job);
	my($info);

	my($status)=1;#default value 1
	while($status){
		foreach $job (@jobList){
			$info=`qstat -j $job 2>&1`;
			if($info =~ "job_number"){
				$status=1;# still run
				next;
			}
			$status=0; # all finish
		}
		#print "status: $status\n";
		sleep(1);# check every 3 sec
	}
}

sub checkJobFinish{ # the function accept the a jobID file as input 
       # Check job status
	my($jobListFile)=$_[0]; #input a .jobID file
	open(JOB,"<$jobListFile");# get all jobs ID and put them into @jobID
	my(@jobList)=<JOB>;
	close(JOB);
	my(@jobID);
	my(@f);
	my($x);
	my($job);
	my($info);
	
	#get job ID list
	foreach $x (@jobList){
		@f=split(" ",$x);
		push(@jobID,$f[2]);# the 3 field in the line is jobID
	}

	my($status)=1;#default value 1
	while($status){
		foreach $job (@jobID){
			$info=`qstat -j $job 2>&1`;
			if($info =~ "job_number"){
				#print "still working\n";
				$status=1;# still run
				next;
			}
			$status=0; # all finish
		}
		sleep(1);# check every 3 sec
	}
}

sub runbwaSampe{
	#enter the specific path
#	print "\n========================================================\n";
#	print "runbwaSampe start:";
#	print scalar localtime(),qq(\n);
	print "runbwaSampe running...\n";
	my($workPath)="$projectFolder/bam/bwa/";
	chdir("$workPath");

	# get fastq.gz file list in inputFolder
	opendir(TEMP,$inputFolder) || die "Cannot open $inputFolder";
	my(@files) = readdir TEMP;
	closedir(TEMP);
	my($x);
	my(@fqList);
	foreach $x (@files){
		push(@fqList,$x) if($x =~ /.*fastq.gz/);
	}

	# For example:sample_15_ACAGTG_L004_R1_001.fastq.gz
	# Get prefix of every fq files like : sample_15_ACAGTG_L004_R1_001.fastq.gz
	my(%prefixHash);
	foreach $x (@fqList){
		my(@part)=split("_",$x);

		# to detect the "L001" in which part
		my($lanePartNum)=0;  
		my($lane);
		foreach $lane (@part){
			if($lane =~ /^L0/){last}# find out the lane part
			$lanePartNum+=1;
		}

		my($prefix)=join("_",@part[0..$lanePartNum]);
		$prefixHash{$prefix}="1"; # save in a hash table to avoid prefix duplicate
		#print $prefix,qq(\n);
	}
	
	#qsub all pair fq file
	my($lane);
	my($ID)=1;
	my($command);
	my(@prefixTable)=keys(%prefixHash);
	my(@attr);
	
	unlink("bwaSampe.jobID") if (-e "bwaSampe.jobID"); # check if bwaSampe.jobID exist , revise every time
	foreach $lane (@prefixTable){
		@attr=("sampe","-P","-r '\@RG\tID:$ID\tFuwai:BCHR\tPG:bwa\tPI:100\tPL:Illumina\tSM:$sample'",$refSequence,"$lane"."_R1_001.fastq.gz.sai","$lane"."_R2_001.fastq.gz.sai","$inputFolder/$lane"."_R1_001.fastq.gz","$inputFolder/$lane"."_R2_001.fastq.gz");
		$command=join(" ",$bwa,@attr,">$lane.sam");
		open(OUT,">sampeCommand.sh"); #write the command into a file , and then qsub it . Because qsub -b y can not handle a command contain '' 
		print OUT $command;
		close(OUT);
		$command=join(" ","qsub -wd \$PWD","sampeCommand.sh");
		system("$command >> bwaSampe.jobID");
		$ID+=1;
	}

	# Check job status
	checkJobFinish("bwaSampe.jobID");

	#go back to project path
	chdir($projectFolder);
#	print "runbwaSampe end:";
#	print scalar localtime(),qq(\n);
}

sub runsamtoolsGenerateBam{
	#enter the specific path
#	print "\n========================================================\n";
#	print "runsamtoolsConvertSamToBam start:";
#	print scalar localtime(),qq(\n);
	print "runsamtoolsConvertSamToBam running...\n";

	my($workPath)="$projectFolder/bam/samtools/";
	my($samPath)="$projectFolder/bam/bwa/";
	my($x);

	chdir($workPath);
	# get sam files list in bwa folder
	my(@samList);
	opendir(TEMP,$samPath) || die "Cannot open $samPath";
	my(@files) = readdir TEMP;
	closedir(TEMP);
	foreach $x (@files){
		push(@samList,$x) if($x =~ /sam$/);#end with sam
	}

	#qsub all sam to bam
	my(@attr);
	my($command);
	my($sam);
	unlink("samtoolsGB.jobID") if (-e "samtoolsGB.jobID"); # check if bwaSampe.jobID exist , revise every time
	foreach $sam (@samList){
		@attr=("view -bS","$samPath/$sam");
		$command=join(" ","qsub -b y -wd \$PWD -e $sam.log -o $sam.bam",$samtools,@attr,">> samtoolsGB.jobID");
		system($command);
#		print $command."\n";
	}

       # Check job status
	checkJobFinish("samtoolsGB.jobID");

	#qsub all sort job
	unlink("samtoolsSort.jobID") if (-e "samtoolsSort.jobID"); # check if bwaSampe.jobID exist , revise every time
	foreach $sam (@samList){
		@attr=("sort","$sam.bam","$sam.bam.sorted");
		$command=join(" ","qsub -b y -wd \$PWD",$samtools,@attr,">> samtoolsSort.jobID");
		system($command);
#		print $command."\n";
	}

       # Check job status
	checkJobFinish("samtoolsSort.jobID");

	#go back to project path
	chdir($projectFolder);
#	print "runsamtoolsConvertSamToBam end:";
#	print scalar localtime(),qq(\n);
}
sub runpicardMerge{
	#enter the specific path
#	print "\n========================================================\n";
#	print "runpicardMerge start:";
#	print scalar localtime(),qq(\n);
	print "runpicardMerge running...\n";

	my($workPath)="$projectFolder/bam/picard/";
	my($bamPath)="$projectFolder/bam/samtools/";
	#my($bamPath)="/wa/ugoodlfy/sequencing_201301/sample_11_newRef/bam/samtools/"; # test!!!
	my($x);

	chdir($workPath);
	# get bam files list in samtools folder
	my(@bamList);
	opendir(TEMP,$bamPath) || die "Cannot open $bamPath";
	my(@files) = readdir TEMP;
	closedir(TEMP);
	foreach $x (@files){
		push(@bamList,$x) if($x =~ /sorted.bam$/);#end with sorted bam files.Be cautious! samtools sort output end with sorted.bam , no sorted !!
		#push(@bamList,$x) if($x =~ /sorted.bam$/);#end with sorted bam files  test!!!
	}

	#qsub all sam to bam
	my(@attr);
	my(@inputList);
	my($command);

	# to adapt the input format of picard : INPUT=~.bam.sort
	foreach $x (@bamList){push(@inputList,"INPUT=$bamPath/$x")}

	#qsub merge bam
	unlink("picardMerge.jobID") if (-e "picardMerge.jobID"); # check if picardMerge.jobID exist , revise every time
	@attr=("VALIDATION_STRINGENCY=$Picard_MergeBam_validationStringency",@inputList,"OUTPUT=$sample.bam","ASSUME_SORTED=$Picard_MergeBam_assumeSorted","USE_THREADING=$Picard_MergeBam_useThreading");
	$command=join(" ","qsub -b y -wd \$PWD",$picardMerge,@attr,">> picardMerge.jobID");
#	print $command."\n";
	system("$command >> picardMerge.jobID");

	checkJobFinish("picardMerge.jobID");

	#go back to project path
	chdir($projectFolder);
#	print "runpicardMerge end:";
#	print scalar localtime(),qq(\n);
}

sub runpicardDuplicateRemove{
	#enter the specific path
#	print "\n========================================================\n";
#	print "runpicardDuplicateRemove start:";
#	print scalar localtime(),qq(\n);
	print "runpicardDuplicateRemove running...\n";

	my($workPath)="$projectFolder/bam/picard/";
	my($x);

	chdir($workPath);
	#qsub bam file to removed
	my(@attr);
	my($command);

	unlink("picardDuplicateRemove.jobID") if (-e "picardDuplicateRemove.jobID"); # check if picardDuplicateRemove.jobID exist , revise every time
	@attr=("VALIDATION_STRINGENCY=LENIENT","MAX_SEQUENCES_FOR_DISK_READ_ENDS_MAP=$picardDuplicateRemove_maxSequencesForDiskReadEndsMap","REMOVE_DUPLICATES=$picardDuplicateRemove_removeDuplicates","ASSUME_SORTED=$picardDuplicateRemove_assumeSorted","INPUT=$sample.bam","OUTPUT=$sample"."_duprmed.bam","METRICS_FILE=$sample"."_duprmedMetrics.info");
	$command=join(" ","qsub -b y -wd \$PWD",$picardDuplicateRemove,@attr,">> picardDuplicateRemove.jobID");
#	print $command."\n";
	system("$command >> picardDuplicateRemove.jobID");

	checkJobFinish("picardDuplicateRemove.jobID");

	#go back to project path
	chdir($projectFolder);
#	print "runpicardDuplicateRemove end:";
#	print scalar localtime(),qq(\n);
}

sub runpicardSummary{
	#enter the specific path
#	print "\n========================================================\n";
#	print "runpicardSummary start:";
#	print scalar localtime(),qq(\n);
	print "runpicardSummary running...\n";

	my($workPath)="$projectFolder/bam/picard/";
	my($x);

	chdir($workPath);
	#qsub bam file to removed
	my(@attr);
	my($command);

	@attr=("VALIDATION_STRINGENCY=LENIENT","REFERENCE_SEQUENCE=$refSequence","INPUT=$sample"."_duprmed.bam","OUTPUT=$sample"."_collect_AlignmentSummary");
	$command=join(" ","qsub -b y -wd \$PWD",$picardSummary,@attr,">> picardSummary.jobID");
#	print $command."\n";
	system("$command > picardSummary.jobID");
	my(@job)=extractJobIDList("picardSummary.jobID");
	#checkJobFinish("picardSummary.jobID");

	#go back to project path
	chdir($projectFolder);
#	print "runpicardSummary end:";
#	print scalar localtime(),qq(\n);
	return(@job);

}

sub runQueue{
	#enter the specific path
#	print "\n========================================================\n";
#	print "runQueue start:";
#	print scalar localtime(),qq(\n);
	print "runQueue running...\n";

	my($workPath)="$projectFolder/bam/Queue/";
	my($bamPath)="$projectFolder/bam/picard/";
	my($x);

	chdir($workPath);
	#qsub bam file to removed
	my(@attr);
	my($command);
##directory of DataProcessingPipeline.scala
#-S path of DataProcessingPipeline.scala
	@attr=("-S $DataProcessingPipeline","-i $bamPath/$sample"."_duprmed.bam","-R $refSequence","-D $refDbsnp","-jobRunner GridEngine","-run","-sg 10"); #-sg 10
	#@attr=("-S /wa/ugoodlfy/Queue_work/DataProcessingPipeline.scala","-i /data/software/bin/Queue-2.1-10/resources/exampleBAM.bam","-R /data/software/bin/Queue-2.1-10/resources/exampleFASTA.fasta","-D $refDbsnp","-jobRunner GridEngine","-run"); #test!!!!!!!!!!!!
	$command=join(" ",$Queue,@attr);
#	print $command."\n";
	system("$command > Queue.log"); # the command will not stop until it finish


	#go back to project path
	chdir($projectFolder);
#	print "runQueue end:";
#	print scalar localtime(),qq(\n);
}

sub rungatkReadDepth{
	#enter the specific path
#	print "\n========================================================\n";
#	print "rungatkReadDepth start:";
#	print scalar localtime(),qq(\n);
	print "rungatkReadDepth running...\n";

	my($workPath)="$projectFolder/analysis/readDepth/";
	my($bamPath)="$projectFolder/bam/Queue/";
	my($x);
	chdir($workPath);

	#qsub bam file to removed
	my(@attr);
	my($command);

	@attr=("-T DepthOfCoverage","-R $refSequence","-I $bamPath/project.$sample.clean.dedup.recal.bam","-o $sample.metrics","-ct 5 -ct 10 -ct 20");
	$command=join(" ","qsub -b y -wd \$PWD",$gatk,@attr);
#	print $command."\n";
	system("$command > DepthOfCoverage.jobID");

	my(@job)=extractJobIDList("DepthOfCoverage.jobID");

	#go back to project path
	chdir($projectFolder);
#	print "rungatkReadDepth end:";
#	print scalar localtime(),qq(\n);
	return(@job);

}

sub rungatkGT{
	#enter the specific path
#	print "\n========================================================\n";
#	print "rungatkGT start:";
#	print scalar localtime(),qq(\n);
	print "rungatkGT running...\n";

	my($workPath)="$projectFolder/variant_calling/";
	#my($workPath)="/wa/ugoodlfy/graduation_project/test"; # test!!!!!!!!!!!!!!!!!!!
	my($bamPath)="$projectFolder/bam/Queue/";
	my($x);

	chdir($workPath);

	#qsub bam file to removed
	my(@attr);
	my($command);

##-S path of UnifiedGenotyper.scala
	@attr=("-S $UnifiedGenotyper","-R $refSequence","-D $refDbsnp","-I $bamPath/project.$sample.clean.dedup.recal.bam","-jobRunner GridEngine","-run","-qsub");
#	@attr=("-T UnifiedGenotyper","-l INFO","-R $refSequence","-D $refDbsnp","-nt 2","-nct 2","-I $bamPath/project.$sample.clean.dedup.recal.bam","-glm BOTH","-o project.$sample.raw.vcf","--output_mode EMIT_VARIANTS_ONLY");
	#@attr=("-T UnifiedGenotyper","-l INFO","-R $refSequence","-D $refDbsnp","-I $bamPath/project.exampleBAM.clean.dedup.recal.bam","-glm BOTH","-o $sample.raw.vcf","--output_mode EMIT_VARIANTS_ONLY"); #test!!!!!!!!!!!!!!!!!!!!!!!!!!!
	$command=join(" ",$Queue,@attr);
#	$command=join(" ",$gatk,@attr);
#	print $command."\n";
	system("$command > ug_Queue.log"); # the command will not stop until it finish


	#go back to project path
	chdir($projectFolder);
#	print "rungatkGT end:";
#	print scalar localtime(),qq(\n);

}

sub rungatkSnpRecalModel{
	#enter the specific path
#	print "\n========================================================\n";
#	print "rungatkSnpRecalModel start:";
#	print scalar localtime(),qq(\n);
	print "rungatkSnpRecalModel running...\n";

	my($workPath)="$projectFolder/variant_calling/";
	my($x);
	chdir($workPath);

	#qsub bam file to removed
	my(@attr);
	my($command);

###-resource 
	@attr=("-T VariantRecalibrator","-R $refSequence","-input project.$sample.raw.vcf","-mode SNP","-resource:hapmap,known=false,training=true,truth=true,prior=15.0 $gatk_SnpRecalHapmap","-resource:omni,known=false,training=true,truth=false,prior=12.0 $gatk_SnpRecalOmni","-resource:dbsnp,known=true,training=false,truth=false,prior=6.0 $gatk_SnpRecalDbsnp","-an QD -an HaplotypeScore -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an DP","-recalFile $sample.SnpRecal.model","-tranchesFile $sample.SnpTranches.model","-rscriptFile $sample.SnpPlot.Rmodel");
	#@attr=("-T UnifiedGenotyper","-l INFO","-R $refSequence","-D $refDbsnp","-I $bamPath/project.exampleBAM.clean.dedup.recal.bam","-glm BOTH","-o $sample.raw.vcf","--output_mode EMIT_VARIANTS_ONLY"); #test!!!!!!!!!!!!!!!!!!!!!!!!!!!
	$command=join(" ","qsub -b y -wd \$PWD",$gatk,@attr);
#	print $command."\n";
	system("$command > gatkSnpRecalModel.jobID");
	checkJobFinish("gatkSnpRecalModel.jobID");


	#go back to project path
	chdir($projectFolder);
#	print "rungatkSnpRecalModel end:";
#	print scalar localtime(),qq(\n);

}

sub rungatkIndelRecalModel{
	#enter the specific path
#	print "\n========================================================\n";
#	print "rungatkIndelRecalModel start:";
#	print scalar localtime(),qq(\n);
	print "rungatkIndelRecalModel running...\n";

	my($workPath)="$projectFolder/variant_calling/";
	my($x);
	chdir($workPath);

	#qsub bam file to removed
	my(@attr);
	my($command);

########-resource
	@attr=("-T VariantRecalibrator","-R $refSequence","-input project.$sample.raw.vcf","-mode INDEL","--maxGaussians 4 -std 10.0 -percentBad 0.12","-resource:mills,known=true,training=true,truth=true,prior=12.0 $gatk_IndelRecalMills","-an QD -an FS -an HaplotypeScore -an ReadPosRankSum","-recalFile $sample.IndelRecal.model","-tranchesFile $sample.IndelTranches.model","-rscriptFile $sample.IndelPlot.Rmodel");
	$command=join(" ","qsub -b y -wd \$PWD",$gatk,@attr);
#	print $command."\n";
	system("$command > gatkIndelRecalModel.jobID");
	my(@jobList)=extractJobIDList("gatkIndelRecalModel.jobID");

	#go back to project path
	chdir($projectFolder);
#	print "rungatkIndelRecalModel end:";
#	print scalar localtime(),qq(\n);
	return(@jobList);

}

sub rungatkApplySnpRecal{
	#enter the specific path
#	print "\n========================================================\n";
#	print "rungatkApplySnpRecal start:";
#	print scalar localtime(),qq(\n);
	print "rungatkApplySnpRecal running...\n";

	my($workPath)="$projectFolder/variant_calling/";
	my($x);
	chdir($workPath);

	#qsub bam file to removed
	my(@attr);
	my($command);

	@attr=("-T ApplyRecalibration","-R $refSequence","-input project.$sample.raw.vcf","--ts_filter_level 99.0","-recalFile $sample.SnpRecal.model","-tranchesFile $sample.SnpTranches.model","-mode SNP","-o $sample.SnpRecal.vcf");
	$command=join(" ","qsub -b y -wd \$PWD",$gatk,@attr);
#	print $command."\n";
	system("$command > runApplySnpRecal.jobID"); 
	checkJobFinish("runApplySnpRecal.jobID"); 


	#go back to project path
	chdir($projectFolder);
#	print "rungatkApplySnpRecal end:";
#	print scalar localtime(),qq(\n);

}

sub rungatkApplyIndelRecal{
	#enter the specific path
#	print "\n========================================================\n";
#	print "rungatkApplyIndelRecal start:";
#	print scalar localtime(),qq(\n);
	print "rungatkApplyIndelRecal running...\n";

	my($workPath)="$projectFolder/variant_calling/";
	my($x);
	chdir($workPath);

	#qsub bam file to removed
	my(@attr);
	my($command);

	@attr=("-T ApplyRecalibration","-R $refSequence","-input $sample.SnpRecal.vcf","--ts_filter_level 99.0","-recalFile $sample.IndelRecal.model","-tranchesFile $sample.IndelTranches.model","-mode INDEL","-o $sample.recal.vcf");
	$command=join(" ","qsub -b y -wd \$PWD",$gatk,@attr);
#	print $command."\n";
	system("$command > runApplyIndelRecal.jobID"); 
	checkJobFinish("runApplyIndelRecal.jobID"); 


	#go back to project path
	chdir($projectFolder);
#	print "rungatkApplyIndelRecal end:";
#	print scalar localtime(),qq(\n);

}

sub runAnnovar{
	#enter the specific path
#	print "\n========================================================\n";
	my(@input)=@_;# runAnnovar(vcf_file,outputPath)
#	print "runAnnovar start:";
#	print scalar localtime(),qq(\n);
	print "runAnnovar running...\n";

	my($vcfFile)=$input[0];
	my($outputPath)=$input[1];
	my($x);
	chdir($outputPath);

	#qsub bam file to removed
	my(@attr);
	my($command);

	# Convert vcf to ann format
	unlink("$outputPath/AnnovarConvert.jobID") if (-e "$outputPath/AnnovarConvert.jobID"); # check if AnnovarConvert.jobID exist , revise every time
	@attr=($vcfFile,"-format vcf4","--includeinfo");
	$command=join(" ","qsub -b y -wd \$PWD -o $outputPath/$sample.ann -e $outputPath/$sample.ann.log","$AnnovarConvertAnn",@attr," >AnnovarConvert.jobID");
#	print $command."\n";
	system("$command"); 
	checkJobFinish("$outputPath/AnnovarConvert.jobID");


	# Annovar summary annotation
	unlink("$outputPath/AnnovarAnnotate.jobID") if (-e "$outputPath/AnnovarAnnotate.jobID"); # check if AnnovarAnnotate.jobID exist , revise every time


#############humandb database path

	@attr=("$outputPath/$sample.ann","--buildver hg19","--verdbsnp 135","--ver1000g 1000g2012apr","--veresp 6500","--genetype refgene","$AnnovarDb","--outfile $outputPath/$sample");


	$command=join(" ","qsub -b y -wd \$PWD -o $outputPath/$sample.annotate.out -e $outputPath/$sample.annotate.error","$Annovar_Annotation",@attr);
#	print $command."\n";
	system("$command > AnnovarAnnotate.jobID"); 
	#go back to project path

	my(@jobList)=extractJobIDList("AnnovarAnnotate.jobID");


	chdir($projectFolder);
#	print "runAnnovar end:";
#	print scalar localtime(),qq(\n);
	return(@jobList);

}

sub runFilter{
	#enter the specific path
	my(@input)=@_;# runFilter(vcf_file,outputPath)
#	print "\n========================================================\n";
#	print "runFilter start:";
#	print scalar localtime(),qq(\n);
	print "runFilter running...\n";

	my($vcfFile)=$input[0];
	my($outputPath)=$input[1];
	my($x);
	chdir($outputPath);

	my(@attr);
	my($command);

	# removeLowQual filter
	@attr=("$vcfFile","$outputPath/$sample.goodQual.vcf");
	$command=join(" ","qsub -b y -wd \$PWD",$removeLowQual,@attr);
	system("$command >removeLowQual.jobID");
	checkJobFinish("removeLowQual.jobID");


	# snpSift filter
	@attr=("filter","-f $sample.goodQual.vcf","\'(QUAL>=50)&(MQ>=20)&(GEN[0].AD[1]>=3)&(DP<201)\'",">$outputPath/$sample.goodQual.snpSift.vcf");
	$command=join(" ",$snpSift,@attr);
	open(OUT,">runsnpSift.sh");
	print OUT $command;
	close(OUT);
	$command=join(" ","qsub -wd \$PWD -e runsnpSift.error -o runsnpSift.out","runsnpSift.sh");
	system("$command > snpSift.jobID");
	checkJobFinish("snpSift.jobID");


	# my filter
	@attr=($sample,"$sample.goodQual.snpSift.vcf","$sample.filtered.vcf");
	$command=join(" ","qsub -b y -wd \$PWD",$myfilter,@attr);
	system("$command > myfilter.jobID");
	checkJobFinish("myfilter.jobID");


	#go back to project path
	chdir($projectFolder);
#	print "runFilter end:";
#	print scalar localtime(),qq(\n);

}

sub runvcf2bco{
	#enter the specific path
	my(@input)=@_;# runFilter(vcf_file,outputPath)
#	print "\n========================================================\n";
#	print "runvcf2bco start:";
#	print scalar localtime(),qq(\n);
	print "runvcf2bco running...\n";

	my($vcfFile)=$input[0];
	my($outputPath)=$input[1];
	my($x);
	chdir($outputPath);

	my(@attr);
	my($command);

	# run 
	@attr=($vcfFile,"$outputPath/$sample");
	$command=join(" ","qsub -b y -wd \$PWD -e runvcf2bco.error -o runvcf2bco.out",$vcf2bco,@attr);
	system("$command > runvcf2bco.jobID");

	#go back to project path
	my(@jobList)=extractJobIDList("runvcf2bco.jobID");

	chdir($projectFolder);
#	print "runvcf2bco end:";
#	print scalar localtime(),qq(\n);

	return(@jobList);
}

sub gatkEvaluation{ #output in the vcfFile path
	#enter the specific path
	my(@input)=@_;# gatkEvaluation(vcf_file)
#	print "\n========================================================\n";
#	print "gatkEvaluation start:";
#	print scalar localtime(),qq(\n);
	print "gatkEvaluation running...\n";

	my($vcfFile)=$input[0];
	my($x);

	my(@attr);
	my($command);

	# run 
	@attr=("-T VariantEval","-R $refSequence","--dbsnp $refDbsnp","--eval:set1 $vcfFile","-o $vcfFile.eval.gatkreport");
	$command=join(" ","qsub -b y -wd \$PWD -e $vcfFile.gatkEvaluation.error -o $vcfFile.gatkEvaluation.out",$gatk,@attr);
	system("$command > $vcfFile.gatkEvaluation.jobID");

	#go back to project path
	my(@jobList)=extractJobIDList("$vcfFile.gatkEvaluation.jobID");

	chdir($projectFolder);
#	print "gatkEvaluation end:";
#	print scalar localtime(),qq(\n);


	return(@jobList);
}

sub collectSummary{
	my($picardAlignmentSummary)="$projectFolder/bam/picard/$sample"."_collect_AlignmentSummary";
	my($picardDupSummary)="$projectFolder/bam/picard/$sample"."_duprmedMetrics.info";
	my($gatkDepthSummary)="$projectFolder/analysis/readDepth/$sample.metrics.sample_summary";
	my($gatkRawVcfSummary)="$projectFolder/variant_calling/project.$sample.raw.vcf.eval.gatkreport";
	my($gatkFilteredVcfSummary)="$projectFolder/filter/$sample.filtered.vcf.eval.gatkreport";

#	my($gatkRawVcfSummary)="/wa/ugoodlfy/test/sample_11.eval.gatkreport"; # test!!!!!!!!!!!!!!!!!!!!!!!!!!

	my(@info);
	my($tmp);

	$tmp=`head -n 10 $picardAlignmentSummary`;
	push(@info,$tmp);

	$tmp=`head -n 8 $picardDupSummary`;
	push(@info,$tmp);

	$tmp=`head -n 10 $gatkDepthSummary`;
	push(@info,$tmp);

	$tmp="=================== raw vcf summary ==========================";
	push(@info,$tmp);
	$tmp=`sed -n '9,14p' $gatkRawVcfSummary`;
	push(@info,$tmp);

	$tmp=`sed -n '94,99p' $gatkRawVcfSummary`;
	push(@info,$tmp);

	$tmp="=================== filtered vcf summary ==========================";
	push(@info,$tmp);
	$tmp=`sed -n '9,14p' $gatkFilteredVcfSummary`;
	push(@info,$tmp);

	$tmp=`sed -n '94,99p' $gatkFilteredVcfSummary`;
	push(@info,$tmp);
       
	chdir($projectFolder);
	open(OUT,">$sample.Summary");
       	foreach $tmp (@info){
       		print OUT $tmp."\n";
       	}
	close(OUT);

}
