#!/usr/bin/env perl
use strict;
use warnings;


open(FILEIN1, $ARGV[0]) || die "cannot open the $ARGV[0]";

my $flag = 0;
my $humanPercentage; 

while(my $readFile = <FILEIN1>)
{
 if(($readFile =~ /^\s*count.*?percent.*?class.*$/) && ($flag == 0))
 {
    $flag = 1;
    next;
 }
 elsif($flag == 1)
 {
   chomp $readFile;

   if($readFile =~ /^\s*.*?\s+(.*?)\s+.*$/)
   {
     $humanPercentage = $1;

   }
  last;
 }
}

if($humanPercentage < $ARGV[1])
{
  die "Human fraction of reads is less than $ARGV[1]% of total reads: The human percentage is $humanPercentage\n";
}
