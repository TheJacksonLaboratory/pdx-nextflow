#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

################Written By Anuj Srivastava#######################################

my $usage = <<'USAGE';

#####################EMASE_Analysis_PDX.pl#############################

        usage: EMASE_Analysis_PDX.pl [options]

                -i1= infile.fastq  [fastq file]
                -od= outDir
                -dd= dataDir
                -oa= outname
                -md= EMASE model
                -ei= EM iteration
                -rl= READ length
                -pc= pseudocount
                -te= tolerance
                -db= delete bam
                -dh= delete h5
                     
###################################################################


USAGE

my ( $infile1, $outDir, $dataDir, $outname, $model, $readLength, $EMi, $pseudoCount, $tolerance, $deleteBam, $deleteh5 );

$model = 3;
$readLength = 100;
$EMi = 999;
$pseudoCount = 0.0;
$tolerance = 0.0001;
$deleteBam = 'TRUE';
$deleteh5  = 'TRUE';

my $result = GetOptions(
    "i1=s"  => \$infile1,
    "od=s"  => \$outDir,
    "dd=s"  => \$dataDir,
    "oa=s"  => \$outname,
    "md=i"  => \$model,
    "rl=i"  => \$readLength,
    "ei=i"  => \$EMi,
    "pc=f"  => \$pseudoCount,
    "te=f"  => \$tolerance,
    "db=s"  => \$deleteBam,
    "dh=s"  => \$deleteh5
)
;


die $usage unless ($infile1);     ##### Mandatory arguments
die $usage unless ($outDir);      ##### Mandatory arguments
die $usage unless ($dataDir);     ##### Mandatory arguments
die $usage unless ($outname);     ##### Mandatory arguments
die $usage unless ($model);       ##### Mandatory arguments
die $usage unless ($readLength);  ##### Mandatory arguments
die $usage unless ($EMi);         ##### Mandatory arguments
die $usage unless ($pseudoCount); ##### Mandatory arguments
die $usage unless ($tolerance);   ##### Mandatory arguments
die $usage unless ($deleteBam);   ##### Mandatory arguments
die $usage unless ($deleteh5);    ##### Mandatory arguments



my $outname2    =`basename $outname`;

chomp $outname2;

&PDXReference();

if($deleteBam =~ /TRUE/ig)
  {
      system("rm $outDir/$outname2.bowtie.transcriptome.bam");
  }

if($deleteh5 =~ /TRUE/ig)
  {
      system("rm $outDir/$outname2.bowtie.transcriptome.h5");
  }



############################Subroutint to Run EMASE#######################

#####Single transcriptome run####

sub PDXReference
{

  if($infile1 =~ /^.*?\.bz2/)
  {
  system("bzcat $infile1 | bowtie -q -a --best --strata --sam -v 3 $dataDir/GRCh38_NOD_ShiLtJ_based_on_mm10.transcripts -    | samtools view -bS -F 4 - > $outDir/$outname2.bowtie.transcriptome.bam");
  }

 elsif($infile1 =~ /^.*?\.gz/)
  {
   system("zcat $infile1 | bowtie -q -a --best --strata --sam -v 3 $dataDir/GRCh38_NOD_ShiLtJ_based_on_mm10.transcripts -    | samtools view -bS -F 4 - > $outDir/$outname2.bowtie.transcriptome.bam");
  }

 else
  {
  system("bowtie -q -a --best --strata --sam -v 3 $dataDir/GRCh38_NOD_ShiLtJ_based_on_mm10.transcripts $infile1    | samtools view -bS -F 4 - > $outDir/$outname2.bowtie.transcriptome.bam"); 
  }

  system("bam-to-emase -a $outDir/$outname2.bowtie.transcriptome.bam  -i $dataDir/pooled.emase.transcriptome.info   -o $outDir/$outname2.bowtie.transcriptome.h5");

 system("run-emase -i $outDir/$outname2.bowtie.transcriptome.h5 -g $dataDir/pooled.emase.gene2transcript.tsv -L $dataDir/pooled.emase.transcriptome.info -M $model -o $outDir/$outname2.emase -m $EMi -r $readLength  -p $pseudoCount  -t $tolerance -c");

 return;

}

