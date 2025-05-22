#!/usr/bin/perl
#----------------------------------------------------------------------+
#                                                                      |
# FungiSNC -  small non-coding RNA calling and expression             |
#           analysis based on Small RNAs sequencing data               |
#                                                                      |
#----------------------------------------------------------------------+
#                                                                      |
#  AUTHOR: Qi Liu                                                      |
#  CONTACT: biolq668@gmail.com                                         |        
#                                                                      |
#  LICENSE:                                                            |
#  GNU General Public License, Version 3                               |
#  http://www.gnu.org/licenses/gpl.html                                |  
#                                                                      |
#  VERSION: V.1.1             				                           |
#                                                                      |
#----------------------------------------------------------------------+


use strict;
use Getopt::Long;
use FindBin '$Bin';

my $usage = <<USAGE;
 perl run_fungi.pl <parameters>
  -mismatch	<num>	The mismatch allowed in the tags mapping [1]
  -infile	<str>	The path of input file (fasta) [required *]
  -species	<str>	The reference species [ath]
  -minlen	<num>	The minimum value of length interval [16]
  -maxlen	<num>	The maximum value of length interval [32]
  -mindepth	<num>	The least number of alignments (Read Per Million) to call tRFs [20]
  -pvalue	<num>   P-value inferred based on Binomial method to distinguish random fragments [0.01]
  -percent	<num>	Minimum percentage of tRF tags in according tRF-producing region on tRNA [0.03]
  -ncpu	<num>   Number of threads used in the analysis [1]
  -outdir	<str>   The output directory of results [tRFtools_time]
USAGE


# ------------------------------------------------------------------
# Options
# ------------------------------------------------------------------
my $mismatch;
my $mindepth;
my $infile;
my $species;
my $minlen;
my $maxlen;
my $pvalue;
my $minpercent;
my $novelPreSoft;
my $ncpu;
my $outdir;
my $help;
my $email;
my $samplename;
my $jobid;
my $scorecut;
my $multimap;
my $paddistance;
my $phaslength;
my $phasmindepth;
my $intype;
my $mirAlignTool;
my $evalue;
my $bvalue;
my $vvalue;
my $miReapL;
my $miRDeepV;
my $tRFscorecut;
my $trnaMaxMap;
my $tRFpvalue;
my $alignType;
my $seedLength;
my $allMismatch;

GetOptions (
            "mismatch=i"    => \$mismatch,
			"infile=s"      => \$infile,
			"species=s"     => \$species,
			"minlen=i"      => \$minlen,
			"maxlen=i"      => \$maxlen,
			"mindepth=i"    => \$mindepth,
			"pvalue=f"      => \$pvalue,
			"percent=f"     => \$minpercent,
			"multimap=s"    => \$multimap,
			"paddistance=s" => \$paddistance,
			"phaslength=s"  => \$phaslength,
			"phasmindepth=s"=> \$phasmindepth,
			"newMirTool=s" 	=> \$novelPreSoft,
			"ncpu=i"        => \$ncpu,
			"outdir=s"      => \$outdir,
			"email=s"       => \$email,
			"jobid=s"       => \$jobid,
			"samplename=s"  => \$samplename,
			"intype=s"      => \$intype,
			"mirTool=s" 	=> \$mirAlignTool,
			"evalue=f"      => \$evalue,
			"bvalue=i"      => \$bvalue,
			"vvalue=i"      => \$vvalue,
			"miReapL=i"  	=> \$miReapL,
			"miRDeepV=f"    => \$miRDeepV,
			"tRFscorecut=i"	=> \$tRFscorecut,
			"trnaMaxMap=i"  => \$trnaMaxMap,
			"tRFpvalue=f"   => \$tRFpvalue,
			"alignType=s"	=> \$alignType,
			"seedLength=i"  => \$seedLength,
			"allMismatch=i" => \$allMismatch,
			'help!'         => \$help
			); 

#$species ||= "zma";
$mindepth ||= 20;
$minlen ||= 16;
$maxlen ||= 40;
$pvalue ||= 0.01;
$minpercent ||= 0.3;
$multimap ||= 50;
$paddistance ||= 150;
$phaslength ||= "21,24";
$novelPreSoft ||= "miRDeep2";
$ncpu ||= 20;
$email ||= "NA";
$jobid ||= "1234561";
$samplename ||= $jobid;
#$outdir ||= "./data/results/$jobid";
$scorecut ||= 10;
$intype ||= "fasta";
$mirAlignTool ||= "bowtie";
$evalue ||= 0.01;
$bvalue ||= 5;
$vvalue ||= 5;
$miReapL ||= 100;
$miRDeepV ||= 10;
$tRFscorecut ||= 20;
$trnaMaxMap ||= 30;
$tRFpvalue ||= 0.05;
$alignType ||= "seed";
$seedLength ||= 20;
$mismatch ||= 2;
$allMismatch ||= 2;

my $mismatchPara;
if($alignType eq "seed") {
	$mismatchPara = "-n $mismatch -l $seedLength";	
} else {
	$mismatchPara = "-v $allMismatch";	
}

if ($infile eq "" or $help) {
   print $usage,"\n";
   exit;
} elsif (! -e $infile) {
   print "\nInput file $infile is not exists!\n\n";
   exit;
}

### Test program ###
my $BWT=qx(which bowtie);
if ($BWT eq ""){
    print "\n\nPlease install bowtie before runing sRNAtools\n\n\n";
    exit;
}
my $SAMTOOLS=qx(which samtools);
if ($SAMTOOLS eq ""){
    print "\n\nPlease install samtools (>v1.2) before runing sRNAtools\n\n\n";
    exit;
}
my $BEDTOOLS=qx(which bedtools);
if ($BEDTOOLS eq ""){
    print "\n\nPlease install bedtools before runing sRNAtools\n\n\n";
    exit;
}


### Read database location configuration###
my $dbbasepath;
my $rfampath;
my $mirbasepath;
my $genomepath;
my $mrnapath;
my $lncRNApath;
my $circRNApath;
my $tRNApath;
my $repbasepath;
my $nt24path;
my $natrnapath;


if (! -f "DBCONFIG_Fungi.txt") {
    print "The databases config file 'DBCONFIG_Fungi.txt' is not exist in current directory!\n";
	exit;
}
open CONFIG,"DBCONFIG_Fungi.txt" or die "Can not open database config file 'DBCONFIG_Fungi.txt'";
while (<CONFIG>) {
   chomp;
   if (/dbbasepath\s*=\s*(\S+)/) {
       $dbbasepath = $1;
   } elsif (/genomepath\s*=\s*(\S+)/) {
       $genomepath = $1;
   } elsif (/mirbasepath\s*=\s*(\S+)/) {
       $mirbasepath = $1;
   } elsif (/rfampath\s*=\s*(\S+)/) {
       $rfampath = $1;
   } elsif (/mrnapath\s*=\s*(\S+)/) {
       $mrnapath = $1;
   } elsif (/lncRNApath\s*=\s*(\S+)/) {
       $lncRNApath = $1;
   } elsif (/circRNApath\s*=\s*(\S+)/) {
       $circRNApath = $1;
   } elsif (/tRNApath\s*=\s*(\S+)/) {
       $tRNApath = $1;
   } elsif (/24ntpath\s*=\s*(\S+)/) {
       $nt24path = $1;
   } elsif (/repbasepath\s*=\s*(\S+)/) {
       $repbasepath = $1;
   }  elsif (/natrnapath\s*=\s*(\S+)/) {
       $natrnapath = $1;
   } 
}						
close CONFIG;

my $genomeindex = "$genomepath/$species.fa";
#my $mirbasehairpin = "$mirbasepath/hairpins/$species.fa";
#my $mirbasemature=    "$mirbasepath/matures/$species.fa";
my $rfamindex = "$rfampath/$species.fa";
#my $rfam2kind = "$rfampath/rfam2kind.txt";
my $mRNAindex = "$mrnapath/$species.fa";
my $lncRNAindex = "$lncRNApath/$species"."_lncrnas.fa";
#my $circRNAindex = "$circRNApath/$species.fa";
my $tRNApri = "$tRNApath/$species.pretRNA.expand.fa";
my $tRNAmature = "$tRNApath/$species.mature.fa";
my $tRNAstructure = "$tRNApath/$species.structure.gff";
my $expandstufile = "$tRNApath/$species.expand.structure.txt";
#my $repbaseindex = "$repbasepath/$species.transposon.bed";
#my $repbasefamily = "$repbasepath/RepeatMasker.family.txt";
#my $locidb = "$nt24path/$species.24nt.cluster.bed";
my $natRNAindex = "$natrnapath/$species.nat.fa";

if(not -e $genomeindex) {
   print "[Error] $species.fa in not exists in directory $genomepath";
   exit;
#} elsif(not -e $mirbasehairpin) {
#   print "[Error] $species.fa in not exists in directory $mirbasepath";
#   exit;
#} elsif(not -e $mirbasemature) {
#   print "[Error] mature sequence $species.fa in not exists in directory $mirbasepath";
#   exit;
} elsif(not -e $rfamindex) {
   print "[Error] $species.fa in not exists in directory $rfampath";
   exit;
} elsif(not -e $mRNAindex) {
   print "[Error] $species.fa in not exists in directory $mrnapath";
   exit;
} elsif(not -e $tRNApri) {
   print "[Error] $species".".pretRNA.expand.fa in not exists in directory $tRNApath";
   exit;
} elsif(not -e $tRNAmature) {
   print "[Error] $species".".mature.fa in not exists in directory $tRNApath";
   exit;
} elsif(not -e $tRNAstructure) {
   print "[Error] $species.structure.gff in not exists in directory $tRNApath";
   exit;
} elsif(not -e $expandstufile) {
   print "[Error] $species.expand.structure.txt in not exists in directory $tRNApath";
   exit;
#} elsif(not -e $repbaseindex) {
#   print "[Error] Repeat/transposon annotation file $species".".transposon.bed in not exists in directory $repbasepath\n\n";
#   exit;
}


checkBowtieIndex($genomeindex);
#checkBowtieIndex($mirbasehairpin);
checkBowtieIndex($rfamindex);
checkBowtieIndex($mRNAindex);
checkBowtieIndex($tRNApri);
checkBowtieIndex($tRNAmature);


##creat result directory
my $resultdir = $outdir;
my $resultTemp = "$resultdir/temp";  ##put the temp results files
if (! -d "$resultTemp/temp") {
   system "mkdir -p $resultTemp";
} else {
   warn "The output directory of results '$resultdir' is already exist!\n";
   #system "mkdir $resultTemp";
   #exit;
}
#filter by length
#&get_log_status("Filter tags by length $minlen -- $maxlen","Tags mapping");
#system("perl ./program/filterByLen.pl $infile $minlen $maxlen > $infile.filter");

#my @filestat = stat ($infile);
#if($filestat[7]/(1024*1024) > 300) {
#	get_log_status("Adaptor information is wrong! Please check!","Error");
#}
open FILE, $infile;
my $totalFreq = 0;
while(<FILE>) {
	chomp $_;
	if ($_ =~ /^>(.+)_x(\d+)/) {
		my $freq = $2;
		$totalFreq += $freq;
	}
}
close FILE;
print $totalFreq."\n"; 
if($totalFreq > 1000000000000) {      
	get_log_status("To many sequene records in the file! The maximum of sequence records (total of frequence) is 50M.","Error");
	exit;
}


#mirtrace
&get_log_status("Running mirtrace","Tags mapping");
system("perl ./program/fa2fq.pl $infile > $resultTemp/converted.fq");
system("java -jar -Xms4G -Xmx4G program/mirtrace/mirtrace.jar trace $resultTemp/converted.fq --output-dir $resultdir/mirtrace");

##mapping to genome
&get_log_status("Mapping to genome","Tags mapping");
system("bowtie -f -p $ncpu $mismatchPara -S $genomeindex $infile | samtools view -bS - > $resultTemp/togenome.bam 2>>$resultdir/run.log");
system("bedtools bamtobed -i $resultTemp/togenome.bam > $resultTemp/togenome.bed");
system("samtools sort $resultTemp/togenome.bam -o $resultTemp/togenome.sort.bam");
system("samtools index $resultTemp/togenome.sort.bam");
system("samtools view -F4 $resultTemp/togenome.bam |awk '{print \$1}'|awk -F\"_[xX]\" '{unique+=1;total+=\$2}END{print total\"\\t\"unique}' > $resultTemp/all.mapping.stat.txt");
system("samtools view -F4 $resultTemp/togenome.bam |awk '{print \$1}'|uniq > $resultTemp/togenome.read.txt");

##miRNA
#&get_log_status("Mapping to miRNA","Tags mapping");
#system("bowtie -f -p $ncpu $mismatchPara --al $resultTemp/mirmapped.fa --un $resultTemp/nonmiR.fa -S $mirbasehairpin $infile | samtools view -bS - > $resultTemp/tomiRNA.bam 2>>$resultdir/run.log");
#system("samtools view -F4 $resultTemp/tomiRNA.bam |awk '{print \$1}'|uniq > $resultTemp/tomiRNA.read.txt");

##mapping siRNA to the tRNA sequences
&get_log_status("mapping to the tRNA","Tags mapping");
#system("bowtie -f -p $ncpu -a -m $trnaMaxMap $mismatchPara -S $tRNApri $resultTemp/nonmiR.fa |samtools view -bS - > $resultTemp/pritRNA.bam 2>>$resultdir/run.log");
system("bowtie -f -p $ncpu -a -m $trnaMaxMap $mismatchPara -S $tRNApri $infile |samtools view -bS - > $resultTemp/pritRNA.bam 2>>$resultdir/run.log");
system("samtools sort -o $resultdir/pritRNA.sort.bam $resultTemp/pritRNA.bam");
system("samtools index $resultdir/pritRNA.sort.bam");
system("samtools view -F4 $resultTemp/pritRNA.bam|awk '{print \$1}' > $resultTemp/pritRNA.read.txt");

#system("bowtie -f -p $ncpu -a -m $trnaMaxMap $mismatchPara --al $resultTemp/tRNAmapped.fa --un $resultTemp/nonmiRtRNA.fa  -S $tRNAmature $resultTemp/nonmiR.fa |samtools view -bS - > $resultTemp/tRNA.bam 2>>$resultdir/run.log");
system("bowtie -f -p $ncpu -a -m $trnaMaxMap $mismatchPara --al $resultTemp/tRNAmapped.fa --un $resultTemp/nonmiRtRNA.fa  -S $tRNAmature $infile |samtools view -bS - > $resultTemp/tRNA.bam 2>>$resultdir/run.log");
system("samtools sort -o $resultdir/tRNA.sort.bam $resultTemp/tRNA.bam");
system("samtools index $resultdir/tRNA.sort.bam");
system("samtools view -F4 $resultTemp/tRNA.bam|awk '{print \$1}' > $resultTemp/tRNA.read.txt");
system("cat $resultTemp/pritRNA.read.txt $resultTemp/tRNA.read.txt |sort|uniq > $resultTemp/totRNA.read.txt");

##other small non-coding RNA in Rfam
&get_log_status("Mapping to other non-coding sncRNA excluding tRNAs","Tags mapping");
system("bowtie -f -p $ncpu $mismatchPara --al $resultTemp/rfammapped.fa --un $resultTemp/nonmiRtRNARfam.fa $rfamindex $resultTemp/nonmiRtRNA.fa -S  |samtools view -bS - > $resultTemp/toRfam.bam 2>>$resultdir/run.log");
system("perl ./program/get_rfam_read.pl $resultTemp/toRfam.bam $resultTemp");
system("samtools view -F4 $resultTemp/toRfam.bam|awk '{print \$1}'|uniq > $resultTemp/toRfam.read.txt");


#novel miRNAs
&get_log_status("Novel miRNA detecting ","Tags mapping");
my $exlength ||= 100;
my $novelExpFile;
#system("cat $resultTemp/tomiRNA.read.txt $resultTemp/totRNA.read.txt $resultTemp/toRfam.read.txt > $resultTemp/classifed.read.txt");
system("cat $resultTemp/totRNA.read.txt $resultTemp/toRfam.read.txt > $resultTemp/classifed.read.txt");
if($novelPreSoft =~ /mireap/i){  
	system("perl ./program/get_unclassfied_reads_mapping_site_for_mireap.pl $resultTemp/classifed.read.txt $resultTemp/nonmiRtRNARfam.fa $resultTemp/togenome.bed $resultdir");
    print "mireap novel miRNA detecting ...\n";	
	system("perl ./program/mireap/bin/mireap.pl -i $resultdir/query_sequence_for_mireap.fa -m $resultdir/unclassfied_reads_mapping_site.txt -r $genomeindex -o $resultdir/ -A 18 -B 26 -a 20 -b 24 -e -18 -d 35 -p 14 -s 0 -f 10 -u 20 -v 0");
 	print "mireap novel miRNA complete. start expression analysis ...\n";	
	system("perl ./program/separate_novel_mature.pl $resultdir/mireap-xxx.aln $resultdir");
	system("perl ./program/novel_mirna_expression_mireap.pl $resultdir/mireap-xxx.aln $resultdir > $resultdir/novel_mirna_expression.list");
	system("perl ./program/parse_mireap_result_get_sequences.pl $resultdir");
	$novelExpFile = "$resultdir/novel_mirna_expression.list";
} elsif ($novelPreSoft =~ /miRDeep2/i) {
    system "perl ./program/excise_premirna.pl $genomeindex $resultTemp/togenome.bed $resultTemp/classifed.read.txt $exlength > $resultdir/precursors.fa";
    system "perl ./program/auto_mapping.pl $resultTemp/nonmiRtRNARfam.fa $resultdir/precursors.fa";
    system "perl ./program/miRDeep_v2.pl $resultdir/signatures $resultdir/precursors.fa > $resultdir/predictions";
	system "perl ./program/mirdeep_prediction_to_aln.pl $resultdir/predictions $resultTemp/nonmiRtRNARfam.fa $resultdir/precursors.fa $resultdir";
	system "perl ./program/separate_novel_mature.pl $resultdir/novel_aln_mirdeep.txt $resultdir";
	system "perl ./program/novel_mirna_expression_mirdeep.pl $resultdir/predictions $resultTemp/nonmiRtRNARfam.fa $resultdir/precursors.fa $resultdir > $resultdir/novel_mirna_expression.list"; 
	$novelExpFile = "$resultdir/novel_mirna_expression.list";
}
system("awk -F'\t' '{print \$8}' $novelExpFile > $resultTemp/tonovel_miRNA.read.txt");


#phasiRNA detecting
&get_log_status("phasiRNA detecting ","Tags mapping");
if (! -d "$resultdir/phasdetect") {
   system "mkdir -p $resultdir/phasdetect";
}
if (! -d "$resultdir/phasmerge") {
   system "mkdir -p $resultdir/phasmerge";
}
my @phasiLenTem = split(/,/,$phaslength);
foreach my $ele (@phasiLenTem) {	if (! -d "$resultdir/phastrigs_$ele") {
	   system "mkdir -p $resultdir/phastrigs_$ele";
	}
}
open CONFIG,">$resultdir/phasdetect/phasis.set" or die "Can not open database config file $resultdir/phasdetect/phasis.set";
print CONFIG "\@runType = G\n\@reference = $genomeindex\n\@userLibs = $resultTemp/nonmiRtRNARfam.fa\n\@libFormat = F\n\@phase = $phaslength\n\@index = $genomeindex\n\@minDepth = 3\n\@clustBuffer = 300\n\@mismat = 0\n";
close CONFIG;
system("python3 ./program/PHASIS/phasdetect.py -outdir $resultdir/phasdetect");
system("python3 ./program/PHASIS/phasmerge.py -mode merge -dir $resultdir/phasdetect -outdir $resultdir/phasmerge");
#foreach my $ele (@phasiLenTem) {
#	system("python3 ./program/PHASIS/phastrigs.py -mode auto -setdir $resultdir/phasdetect -dir $resultdir/phasmerge/$ele -mir $mirbasemature -outdir $resultdir/phastrigs_$ele") if -e "$resultdir/phasmerge/$ele";
#}
system("perl ./program/processPhasResult.pl $resultdir $phaslength $resultTemp/togenome.bed");

##filter
#system("perl ./program/filter.pl $infile $resultTemp/tomiRNA.read.txt $resultTemp/toRfam.read.txt $resultTemp/totRNA.read.txt $resultTemp/tonovel_miRNA.read.txt $resultTemp/tophasiRNA.read.txt $resultTemp/filter.fa");
system("perl ./program/filter.pl $infile $resultTemp/toRfam.read.txt $resultTemp/totRNA.read.txt $resultTemp/tonovel_miRNA.read.txt $resultTemp/tophasiRNA.read.txt $resultTemp/filter.fa");

##extract 24-nt siRNA and mapping 24-nt siRNAs to the reference genome
&get_log_status("extract 24-nt siRNA and mapping 24-nt siRNAs to the reference genome","Tags mapping");
system("perl ./program/get24ntSeq.pl $resultTemp/filter.fa $resultTemp");
print "\nmapping the filterd 24-nt siRNAs to the genome sequences ...\n";
system("bowtie -f -p $ncpu -a -m $multimap $mismatchPara --best --strata -S $genomeindex $resultTemp/24nt.filter.fa 2>>$resultTemp/run.log | samtools view -b -T $genomeindex - | samtools sort - > $resultTemp/24nt_sort.bam 2>>$resultdir/run.log");
system("samtools index $resultTemp/24nt_sort.bam");

##extract 21-24nt siRNA and mapping to the reference genome
&get_log_status("extract 21-24nt and mapping siRNAs to the reference genome","Tags mapping");
print "mapping the filterd 21_24-nt siRNAs to the genome sequences ...\n";
system("bowtie -f -p $ncpu -a -m $multimap -v $mismatch --best --strata -S $genomeindex $resultTemp/21_24nt.filter.fa 2>>$resultTemp/run.log | samtools view -b -T $genomeindex - | samtools sort - > $resultTemp/21_24nt_sort.bam 2>>$resultdir/run.log");
system("samtools index $resultTemp/21_24nt_sort.bam");

##mRNA derived siRNA
&get_log_status("Mapping to mRNA sequences","Tags mapping");
system("bowtie -f -p $ncpu $mismatchPara -S $mRNAindex --un $resultTemp/filter_no_mRNA.fa $resultTemp/filter.fa | samtools view -bS - > $resultTemp/tomRNA.bam 2>>$resultdir/run.log");
system("samtools view -F4 $resultTemp/tomRNA.bam|awk '{print \$1}'|uniq > $resultTemp/tomRNA.read.txt");

#system("cat $resultTemp/togenome.read.txt $resultTemp/tomiRNA.read.txt $resultTemp/toRfam.read.txt $resultTemp/totRNA.read.txt $resultTemp/tonovel_miRNA.read.txt $resultTemp/tophasiRNA.read.txt $resultTemp/tomRNA.read.txt > $resultTemp/all.mapped.read.txt");
system("cat $resultTemp/togenome.read.txt $resultTemp/toRfam.read.txt $resultTemp/totRNA.read.txt $resultTemp/tonovel_miRNA.read.txt $resultTemp/tophasiRNA.read.txt $resultTemp/tomRNA.read.txt > $resultTemp/all.mapped.read.txt");

my $filterSeq = "$resultTemp/filter_no_mRNA.fa";
## lncRNA derived siRNA
#if( -e $lncRNAindex) {
#	&get_log_status("Mapping to lncRNA sequences","Tags mapping");
#	system("bowtie -f -p $ncpu $mismatchPara -S $lncRNAindex --un $resultTemp/filter_no_mRNA_lnc.fa $filterSeq | samtools view -bS - > $resultTemp/lncRNA.bam 2>>$resultdir/run.log");
#	system("samtools view -F4 $resultTemp/lncRNA.bam|awk '{print \$1}'|uniq > $resultTemp/tolncRNA.read.txt");
#	$filterSeq = "$resultTemp/filter_no_mRNA_lnc.fa";
#	system("cat $resultTemp/tolncRNA.read.txt >> $resultTemp/all.mapped.read.txt");
#}

## circRNA derived siRNA
#if( -e $circRNAindex) {
#	&get_log_status("mmapping to circRNA sequences","Tags mapping");
#	system("bowtie -f -p $ncpu $mismatchPara -S $circRNAindex --un $resultTemp/filter_no_mRNA_lnc_circ.fa $filterSeq | samtools view -bS - > $resultTemp/circRNA.bam 2>>$resultdir/run.log");
#	system("samtools view -F4 $resultTemp/circRNA.bam|awk '{print \$1}'|uniq > $resultTemp/tocircRNA.read.txt");
#	$filterSeq = "$resultTemp/filter_no_mRNA_lnc_circ.fa";
#	system("cat $resultTemp/tocircRNA.read.txt >> $resultTemp/all.mapped.read.txt");
#}

## Natural Antisense Transcript derived siRNA
if( -e $natRNAindex) {
	&get_log_status("mmapping to NAT sequences","Tags mapping");
	system("bowtie -f -p $ncpu $mismatchPara -S $natRNAindex $filterSeq | samtools view -bS - > $resultTemp/natRNA.bam 2>>$resultdir/run.log");
	system("samtools view -F4 $resultTemp/natRNA.bam|awk '{print \$1}'|uniq > $resultTemp/toNAT.read.txt");
	system("cat $resultTemp/toNAT.read.txt >> $resultTemp/all.mapped.read.txt");
}


#detect tRF
&get_log_status("tRFs detecting","sncRNAs calling");
system("perl ./program/tRFs.pl --structurefile $tRNAstructure --expandStructurefile $expandstufile --pretRNAbamfile $resultdir/pritRNA.sort.bam --tRNAbamfile $resultdir/tRNA.sort.bam --genomebedfile $resultTemp/togenome.bed --mindepth $mindepth --pvalue $pvalue --percent 0.3 --outdir $resultdir");

#isomiR
&get_log_status("isoMir detecting","sncRNAs calling");

#system("perl ./program/isomiR2Function_rev.pl -t $resultTemp/mirmapped.fa -p $mirbasehairpin -a $mirbasemature -o $resultdir -g $resultTemp/togenome.bed");
#system("cat $resultdir/non-templated_raw.quant1 $resultdir/templated_raw.quant1|awk '{print \$1}' > $resultTemp/isomiR.read.txt");
#system("cat $resultTemp/isomiR.read.txt >> $resultTemp/all.mapped.read.txt");
#system("perl ./program/mapping_stat_isomir.pl $resultdir");


#miRNA and Rfam expression
&get_log_status("Statistics and known expression","sncRNAs calling");
if($mirAlignTool eq "bowtie") {
#	system("bowtie -f -p $ncpu -n 2 --al $resultTemp/mirmapped.fa -S $mirbasemature $infile | samtools view -bS - > $resultTemp/tomatureMiRNA.bam 2>>$resultdir/run.log");
	system("bowtie -f -p $ncpu -v 2 -a -m 50 $rfamindex $resultTemp/rfammapped.fa -S |samtools view -bS - > $resultTemp/toRfam.exp.bam 2>>$resultdir/run.log");
#	system("perl ./program/parse_bowtie_formatted_for_miRNA_exp.pl $resultTemp/tomatureMiRNA.bam $mirbasemature $resultTemp/mirmapped.fa $resultTemp/togenome.bed 18 30 $resultdir");
	system("perl ./program/parse_bowtie_formatted_for_rfam_exp.v2.pl $resultTemp/toRfam.exp.bam $rfamindex $resultTemp/rfammapped.fa $resultTemp/togenome.bed 18 45 $resultdir");
} else {
#	system("./program/external/megablast -i $resultTemp/mirmapped.fa -o $resultTemp/miresult.txt -d $mirbasemature -D 2 -W 12 -v $vvalue -b $bvalue -e $evalue -a $ncpu 2>> $resultTemp/mi.run.log");
	system("perl ./program/megablast_fliter.pl -i $resultTemp/miresult.txt > $resultTemp/miresult.txt.formated");
	system("./program/external/megablast -i $resultTemp/rfammapped.fa -o $resultTemp/rfamesult.txt -d $rfamindex -D 2 -W 12 -v $vvalue -b $bvalue -e $evalue -a $ncpu 2>> $resultTemp/rfam.run.log");
	system("perl ./program/megablast_fliter.pl -i $resultTemp/rfamesult.txt -d 0.9 -l 0 > $resultTemp/rfamesult.txt.formated");
	system("perl ./program/parse_megablast_formatted_for_rnacenter_exp.pl $resultTemp/rfamesult.txt.formated $rfamindex $resultTemp/rfammapped.fa $resultTemp/togenome.bed 18 45 $resultdir");
#	system("perl ./program/parse_megablast_formatted_for_miRNA_exp.pl $resultTemp/miresult.txt.formated $mirbasemature $resultTemp/mirmapped.fa $resultTemp/togenome.bed 18 30 $resultdir");
}

#system "perl ./program/mirna_aln.pl $resultdir/miRNA_result_download $inputfile $mirbaseindex $mirbasetxt >$resultdir/miRNA_aln.result";
#system "perl ./program/mirna_table.pl $resultdir/miRNA_aln.result $resultdir/sum_stat $resultdir > $resultdir/miRNA_stat_table_detail.result";
#system "awk 'BEGIN {OFS=\"\t\"} {print \$1,\$2,\$4,\$5,\$7,\$8,\$9,\$10,\$11}' $resultdir/miRNA_stat_table_detail.result > $resultdir/miRNA_stat_table.result"; 
#system "cp $goutputfile $resultdir/genome_result_download";


system("perl ./program/parse_bowtie_formatted_for_trna_exp.pl $resultTemp/tRNA.bam $tRNAmature $resultTemp/tRNAmapped.fa $resultTemp/togenome.bed 18 45 $resultdir");
system("/data1/users/liuqi/software/miniconda3/envs/R4.4.2/bin/Rscript ./program/tojson.R $resultdir");

##call 24nt clusters###
&get_log_status("24nt siRNA clusters calling","sncRNAs calling");
system("perl ./program/external/ShortStack --nohp --pad $paddistance --mincov $mindepth --bowtie_cores $ncpu --bamfile $resultTemp/24nt_sort.bam --genomefile $genomeindex --outdir $resultTemp/clusterCall 2>>$resultTemp/run.log\n");
system("perl ./program/parse_result_24nt.pl $resultTemp/clusterCall/Results.txt | bedtools sort > $resultTemp/24nt.clusters.sort.bed");


##loci expression### 
&get_log_status("24nt siRNA loci expression analysis ","sncRNAs calling");
system("bedtools bamtobed -i $resultTemp/24nt_sort.bam > $resultTemp/24nt_sort.bed");
##known loci expression###
#if(-e $locidb) {
#    print "24nt with known loci\n";
#	system("bedtools intersect -a $resultTemp/24nt_sort.bed -b $locidb -wo > $resultTemp/24nt.loci.intersact.bed");
#	system("bedtools intersect -a $locidb -b $repbaseindex -wo > $resultTemp/24nt.loci.intersact.repeat.bed");
#	system("perl ./program/getExp.pl $resultTemp/24nt.loci.intersact.bed $resultTemp/all.mapping.stat.txt $resultTemp/24nt.loci.intersact.repeat.bed > $resultdir/24nt.loci.exp.txt");
#} else {
    system("bedtools intersect -a $resultTemp/24nt_sort.bed -b $resultTemp/24nt.clusters.sort.bed -wo > $resultTemp/24nt.loci.intersact.bed");
#	system("bedtools intersect -a $resultTemp/24nt.clusters.sort.bed -b $repbaseindex -wo > $resultTemp/24nt.clusters.intersact.repeat.bed");
#	system("perl ./program/getExp.pl $resultTemp/24nt.loci.intersact.bed $resultTemp/all.mapping.stat.txt $resultTemp/24nt.clusters.intersact.repeat.bed > $resultdir/24nt.loci.exp.txt");
#}
system("perl ./program/getReadExp.pl $resultTemp/24nt.loci.intersact.bed $resultTemp/all.mapping.stat.txt $resultTemp/24nt.filter.fa > $resultdir/24nt.loci.exp.read.txt");
system("awk -F'\t' '{print \$1}' $resultdir/24nt.loci.exp.read.txt > $resultTemp/to24nt_siRNA.read.txt");
system("cat $resultTemp/to24nt_siRNA.read.txt >> $resultTemp/all.mapped.read.txt");
system("perl ./program/filter.pl $infile $resultTemp/all.mapped.read.txt $resultTemp/unmapped.fa");


##repeats/transposons tags statistics###
#&get_log_status("24-nt siRNA statistics on repeats/transposons","sncRNAs calling");
#system("bedtools intersect -a $resultTemp/24nt_sort.bed -b $repbaseindex -bed -wo > $resultTemp/repeats.intersact.bed");
#system("perl ./program/getRepeatMapping.pl $repbasefamily $resultTemp/repeats.intersact.bed > $resultdir/24nt.transposons.tags.stat.txt"); 

##basic mapping stat###
&get_log_status("Basic mapping statistics","sncRNAs calling");
system("samtools view -F 4 $resultTemp/24nt_sort.bam |awk '{print \$1}'|awk -F\"_x|_X\" '{unique+=1;total+=\$2} END{print total\"\t\"unique}' > $resultdir/24nt.mapped.stat.txt");
system("samtools view -f 4 $resultTemp/24nt_sort.bam |awk '{print \$1}'|awk -F\"_x|_X\" '{unique+=1;total+=\$2} END{print total\"\t\"unique}' > $resultdir/24nt.unmapped.stat.txt");
#system("perl ./program/mergeRepeatExp.pl 24nt $resultdir");
system("perl ./program/mapping_stat_chr.pl $resultTemp/24nt_sort.bed $resultdir/24nt.loci.exp.txt 24nt $resultdir");
#system("perl $Bin/program/mapping_stat_type.pl $resultTemp/tomiRNA.read.txt $resultTemp/toRfam.read.txt $resultdir/all.mapping.stat.txt  $resultTemp/mapped.stat.txt > $resultdir/types.stat.txt");


##call 21-24nt clusters###
&get_log_status("21-24nt siRNA clusters calling","sncRNAs calling");
system("perl ./program/external/ShortStack --nohp --pad $paddistance --mincov $mindepth --bowtie_cores $ncpu --bamfile $resultTemp/21_24nt_sort.bam --genomefile $genomeindex --outdir $resultTemp/clusterCall2 2>>$resultTemp/run.log\n");
system("perl ./program/parse_result_24nt.pl $resultTemp/clusterCall2/Results.txt | bedtools sort > $resultTemp/21_24nt.clusters.sort.bed");
##loci expression### 
&get_log_status("21-24nt siRNA loci expression analysis ","sncRNAs calling");
system("bedtools bamtobed -i $resultTemp/21_24nt_sort.bam > $resultTemp/21_24nt_sort.bed");
system("bedtools intersect -a $resultTemp/21_24nt_sort.bed -b $resultTemp/21_24nt.clusters.sort.bed -wo > $resultTemp/21_24nt.loci.intersact.bed");
#system("bedtools intersect -a $resultTemp/21_24nt.clusters.sort.bed -b $repbaseindex -wo > $resultTemp/21_24nt.clusters.intersact.repeat.bed");
#system("perl ./program/getExp_cluster.pl $resultTemp/21_24nt.loci.intersact.bed $resultTemp/all.mapping.stat.txt $resultTemp/21_24nt.clusters.intersact.repeat.bed $resultTemp/21_24nt.filter.fa > $resultdir/21_24nt.loci.exp.txt");
system("perl ./program_fungi/getExp_cluster.pl $resultTemp/21_24nt.loci.intersact.bed $resultTemp/all.mapping.stat.txt $resultTemp/21_24nt.filter.fa > $resultdir/21_24nt.loci.exp.txt");

##known loci expression###
system("perl ./program/getReadExp.pl $resultTemp/21_24nt.loci.intersact.bed $resultTemp/all.mapping.stat.txt $resultTemp/21_24nt.filter.fa > $resultdir/21_24nt.loci.exp.read.txt");
system("awk -F'\t' '{print \$1}' $resultdir/21_24nt.loci.exp.read.txt > $resultTemp/to21_24nt_siRNA.read.txt");
##repeats/transposons tags statistics###
#&get_log_status("21_24-nt siRNA statistics on repeats/transposons","sncRNAs calling");
#system("bedtools intersect -a $resultTemp/21_24nt_sort.bed -b $repbaseindex -bed -wo > $resultTemp/repeats.21_24nt.intersact.bed");
#system("perl ./program/getRepeatMapping.pl $repbasefamily $resultTemp/repeats.21_24nt.intersact.bed > $resultdir/21_24nt.transposons.tags.stat.txt"); 
&get_log_status("Basic mapping statistics","sncRNAs calling");
system("samtools view -F 4 $resultTemp/21_24nt_sort.bam |awk '{print \$1}'|awk -F\"_x|_X\" '{unique+=1;total+=\$2} END{print total\"\t\"unique}' > $resultdir/21_24nt.mapped.stat.txt");
system("samtools view -f 4 $resultTemp/21_24nt_sort.bam |awk '{print \$1}'|awk -F\"_x|_X\" '{unique+=1;total+=\$2} END{print total\"\t\"unique}' > $resultdir/21_24nt.unmapped.stat.txt");
system("perl ./program/mergeRepeatExp.pl 21_24nt $resultdir");
system("perl ./program/mapping_stat_chr.pl $resultTemp/21_24nt_sort.bed $resultdir/21_24nt.loci.exp.txt 21_24nt $resultdir");


##total mapping stat
system("perl ./program/mappingstat.pl $resultTemp $resultdir $infile");
system("perl ./program/mapping_stat_chr_total.pl $resultTemp $resultdir");
system("/data1/users/liuqi/software/miniconda3/envs/R4.4.2/bin/Rscript ./program/lengthStat.R $infile $jobid");
system("/data1/users/liuqi/software/miniconda3/envs/R4.4.2/bin/Rscript ./program/pdfHtmlMaker.R $jobid");
system("mv $resultTemp/24nt.clusters.sort.bed $resultTemp/all.mapping.stat.txt $resultTemp/24nt_sort.bed $resultdir");
#system("rm -fr $resultTemp");

if($email ne "NA") {
    system("python ./program/sendresultmail.py $email $jobid");
}

&get_log_status("Job completed! Thanks for using sRNAtools","Job complete");

sub checkBowtieIndex{
    my $file = shift;
    my $one = $file.".1.ebwt";
	my $two = $file.".2.ebwt";
	my $three = $file.".3.ebwt";
	my $four = $file.".4.ebwt";
	my $rev1 = $file.".rev.1.ebwt";
	my $rev2 = $file.".rev.2.ebwt";
	unless((-r $one) and (-r $two) and (-r $three) and (-r $four) and (-r $rev1) and (-r $rev2)) { 
	    print("Expected bowtie indices for the file not found. Trying to build them with bowtie-build ...\n");
	    system "bowtie-build $file $file > /dev/null";
	    unless((-r $one) and (-r $two) and (-r $three) and (-r $four) and (-r $rev1) and (-r $rev2)) {
		   print ("Build index for $file FAILED. Aborting.");
		   exit;
	    }
	}
}

sub get_log_status()
{ 
    my $info = $_[0];
	my $status = $_[1];
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime();
	$year += 1900;
	$mon +=1;
	my $logdate=sprintf("$info\t%04d-%02d-%02d %02d:%02d:%02d\n",$year,$mon,$mday,$hour,$min,$sec);
	print($logdate);
	system "echo '{\"type\":\"singlecase\",\"species\":\"$species\",\"intype\":\"$intype\",\"current\":\"$status\",\"message\":\"$info\"}' > $resultdir/status.json"; 
	system "echo '$logdate' >> $resultdir/log.txt";
	print $info."\n";
}
