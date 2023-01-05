#!/usr/bin/perl
use strict;
use Getopt::Long;

#################Written By Anuj Srivastava#######################################

my $usage = <<'USAGE';

#########################GeneName_and_Normalization_TCGA.pl#######################

        usage: GeneName_and_Normalization_TCGA.pl [options]

               -i1 = infile1    [RSEM ENSEMBL .TCGA.genes.results]
               -i2 = infile2    [RSEM ENSEMBL .TCGA.isoforms.results]

##################################################################################

USAGE

my ($infile1, $infile2);

my $result = GetOptions("i1=s"  => \$infile1, "i2=s" => \$infile2);

die $usage unless ($infile1);            ##### Mandatory arguments
die $usage unless ($infile2);            ##### Mandatory arguments

###########################################################################################################################


my $seventyFifthQuartileGene  = quantileCal($infile1); 
my $seventyFifthQuartileIso   = quantileCal($infile2); 


open(FILEINGENE, $infile1) || die "cannot open the $infile1 file";####*genes.results 
open(FILEOUTGENENorm, ">$infile1.Normalized")   || die "cannot open the file";


my $flagGene = 0;
my $ExpectedcountG;

while(my $readFileGene = <FILEINGENE>)
{ 
  if($flagGene == 0 )
  {
    $flagGene = 1;
    print FILEOUTGENENorm "gene_id\ttranscript_id(s)\tnormalized_count\n";
    next;
  }

  if($flagGene == 1)
  {
    if($readFileGene =~ /^\s*(.*?)\s+(.*?)\s+(.*?)\s+(.*?)\s+(.*?)\s+(.*?)\s+(.*)$/)
    { 
       $ExpectedcountG = sprintf("%.2f",(($5/$seventyFifthQuartileGene)*1000));
       print FILEOUTGENENorm "$1\t$2\t$ExpectedcountG\n";      


    }
  }
}


##################################################################################################################



#######Isoforms.results modification#############


open(FILEINISO, $infile2) || die "cannot open the $infile2 file";####*genes.results 
open(FILEOUTISONorm, ">$infile2.Normalized")   || die "cannot open the file";

my $flagIso = 0;
my $ExpectedcountI;

while(my $readFileIso = <FILEINISO>)
{ 
  if($flagIso == 0 )
  {
    $flagIso = 1;
    print FILEOUTISONorm "transcript_id\tgene_id\tnormalized_count\n";
    next;
  }

  if($flagIso == 1)
  {
    if($readFileIso =~ /^\s*(.*?)\s+(.*?)\s+(.*?)\s+(.*?)\s+(.*?)\s+(.*?)\s+(.*?)\s+(.*)$/)
    { 
       $ExpectedcountI     = sprintf("%.2f",(($5/$seventyFifthQuartileIso)*300));
       print FILEOUTISONorm "$1\t$2\t$ExpectedcountI\n";

    }
  }
}



sub quantileCal
{

my $filename = shift;

open(FILEOUT1, ">$filename.TMP") || die "cannot open the file";

print FILEOUT1 "X=read.table(\"$filename\", header=T)\n";
print FILEOUT1 "head(X)\n";
print FILEOUT1 "Y=subset(X, X\$expected_count > 0)\n";
print FILEOUT1 "head(Y)\n";
print FILEOUT1 "Z=sort(Y[,5])\n";
print FILEOUT1 "quantile(Z)\n";


system("cat $filename.TMP | R --vanilla >$filename.TMP2 ");

open(FILEINTMP, "$filename.TMP2") || die "cannot open the file";

my $flag = 0;
my $value;

label:while(my $readFile = <FILEINTMP>)
{
 if($readFile =~ /quantile/)
  {
    $flag = 1; 
    next;
  }
  
 elsif($flag == 1)
  { 
    $flag = 2;
    next;
  }

 elsif($flag == 2)
  {
     if($readFile =~ /^\s*(.*?)\s+(.*?)\s+(.*?)\s+(.*?)\s+(.*)$/)
     {
       $value = $4;
       last label;
     }

  }
}

close(FILEOUT1);
close(FILEINTMP);

system("rm $filename.TMP $filename.TMP2 ");

return $value;

}

