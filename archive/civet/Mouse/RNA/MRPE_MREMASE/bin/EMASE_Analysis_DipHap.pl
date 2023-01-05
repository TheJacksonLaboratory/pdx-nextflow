#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

#######################################################

my $usage = <<'USAGE';

#####################EMASE_Analysis_DipHap.pl#############################

        usage: EMASE_Analysis_DipHap.pl [options]

                -gl= genome list   [genome list]
                -od= outDir
                -dd= dataDir
            
               
                     
###################################################################


USAGE

my ( $genomeList, $dataDir,$outDir  );


my $result = GetOptions(
    
    "gl=s"  => \$genomeList,
    "dd=s"  => \$dataDir,
    "od=s"  => \$outDir
    
)
;



die $usage unless ($genomeList);  ##### Mandatory arguments
die $usage unless ($dataDir);     ##### Mandatory arguments
die $usage unless ($outDir);      ##### Mandatory arguments



$genomeList     =`basename $genomeList`;


print"genomeList  is $genomeList\n";

chomp $genomeList;


###my $joinlist  = join(":",  $genomeList);  Preeti changed it because it was only taking the first reference 

my $joinlist  = $genomeList;
print"joinlist  is $joinlist\n";


##my @splitlist = split(":", $joinlist); Preeti changed it 

my @splitlist = split(":", $joinlist);




####################################



 &diploidOrpolyploid(\@splitlist);


###################################################
#####F1 or any number of combine reference#########
###################################################



sub diploidOrpolyploid
{

 my @nameStrains = @{$_[0]};
 my $allStrain = '';
 my $firstStrain = '';

 for(my $i=0; $i<=$#nameStrains; $i++)
 {

  my $modStrain = $nameStrains[$i];
     $modStrain =~ s/_//ig;
   
   print "Doing $nameStrains[$i]\n";
  `grep .  $dataDir/$nameStrains[$i]/$nameStrains[$i].transcripts.fa | awk '{ if(\$1 ~ />.*/) print \$1"_""$modStrain"; else print \$0}'                 >>$outDir/pooled.transcripts.fa`;

   `grep .  $dataDir/$nameStrains[$i]/$nameStrains[$i].emase.transcriptome.info  | awk '{ print \$1"_""$modStrain""\t"\$2}' >>$outDir/pooled.transcripts.info`;

   print "i is >>>> $i\n";
   if($i == 0)
   {
     $firstStrain = $modStrain;
	 print "allStrain at 0 is >> $allStrain\n";
     next;
   }

   elsif($i == 1)
   {
    $allStrain= "$firstStrain".","."$modStrain";
	print "allStrain at 1 is >> $allStrain\n";
    next;
   }

   elsif($i > 1)
   {
     $allStrain = "$allStrain".","."$modStrain";
	 print "allStrain at more than 1 is >> $allStrain\n";
     next;
   }
   
   print "allStrain $allStrain\n";
 }####End For
  
 
  #system("bowtie-build $outDir/pooled.transcripts.fa pooled.transcripts");
return;

}

