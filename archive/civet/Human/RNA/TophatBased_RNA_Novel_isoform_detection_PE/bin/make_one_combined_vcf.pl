#!/usr/bin/perl
use strict;
use warnings;

open(FILEIN1, $ARGV[0]) || die "cannot open the $ARGV[0] file";####SNP   file
open(FILEIN2, $ARGV[1]) || die "cannot open the $ARGV[1] file";####Indel file

my $chrLine;

open(FILEOUT1, ">$ARGV[0].TMP") || die "cannot open the file";

my $ref      = $ARGV[2];
my $finalout = $ARGV[3];

####Getting CHROM line of indel file#####

while(my $readFile = <FILEIN2>)####Reading Indel file
{
  if(($readFile =~ /^\s*#.*$/) && ($readFile !~ /^\s*#CHROM.*$/))
  {
      next;
  }

  elsif($readFile =~ /^\s*#CHROM.*$/)
  {

      $chrLine = $readFile;
      next;
  }

  else
  {
     next;


  }

}

####Making a newSNP file with Indel CHROM LINE##### 


while(my $readFile1 = <FILEIN1>)####Reading Indel file
{
  if(($readFile1 =~ /^\s*#.*$/) && ($readFile1 !~ /^\s*#CHROM.*$/))
  {
      print FILEOUT1 "$readFile1";
      next;
  }

  elsif($readFile1 =~ /^\s*#CHROM.*$/)
  {

      print FILEOUT1 "$chrLine";
      next;
  }

  else
  {
     print FILEOUT1 "$readFile1";
     next;
 
  }

}

my $cmd1 = `java -Xmx2g -jar /opt/compsci/GATK/3.1-1/GenomeAnalysisTK.jar -R $ARGV[2]  -T CombineVariants --variant:SNP $ARGV[0].TMP --variant:INDEL $ARGV[1] -o $finalout `;

system($cmd1);

system("rm $ARGV[0].TMP $ARGV[0].TMP.idx");


