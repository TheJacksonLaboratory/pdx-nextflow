#!/usr/bin/perl
use strict;
use warnings;

my $flag  = 0;
my @line  = ();
my @line2 = ();

open(FILEIN1, $ARGV[0]) || die "cannot open the file"; #####*_stat fiie
 
open(FILEIN2, $ARGV[1]) || die "cannot open the file"; #### *_aln_sum.txt

##############################################
my ($line1, $line2, $line3);

while(my $readFile1 = <FILEIN1>) ####reading filter/*_stat file
{
  if(($readFile1 =~ /^\s*Statistic.*Read.*1.*Read.*2.*$/) && ($flag == 0))
   {
        print "$readFile1";
        $flag = 1;
        next;
   }
   elsif(($flag == 1) && ($readFile1 !~ /Detailed QC statistics/))
  {
      $line1 = $readFile1;
      $flag = 2;
      next;   
  }
   elsif(($flag == 2) && ($readFile1 !~ /Detailed QC statistics/))
  {
      $line2 = $readFile1;
      $flag = 3;
      next;   
  }
   elsif(($flag == 3) && ($readFile1 !~ /Detailed QC statistics/))
  {
      $line3 = $readFile1;
      $flag = 4;
      next;   
  }
}


print "$line2";
print "$line3";
print "$line1";

$flag = 0;



while(my $readFile2 = <FILEIN2>) ####Reading aln_sum.txt file
{
 if($readFile2 =~ /^\s*## METRICS CLASS.*$/)
 {
    $flag = 1;
    next;
 }

 if($flag == 1)
 {
   chomp $readFile2;
   $flag = 2;
   @line = split("\t",$readFile2);
   next;
 }

 if($flag == 2)
 {
   chomp $readFile2;
   @line2 = split("\t",$readFile2);
   last;
 }

}

#print "@line\t@line2\n";

for (my $i=0; $i<=$#line; $i++)
{

  if(($i ==0 )  || ($i == 3))
  {
   print "$line[$i]\t$line2[$i]\n";
   next;
  }
  else
  {
   $line2[$i] = sprintf("%.4f", $line2[$i]);
    print "$line[$i]\t$line2[$i]\n"
  }

}
