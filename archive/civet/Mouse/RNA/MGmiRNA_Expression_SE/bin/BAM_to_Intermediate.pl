#!/usr/bin/perl
use strict;

###################################################################################################################################
# Bam_to_Intermediate.pl
#
# This program uses samtools to read a BAM file and from this create a miRspring intermediate file. Please see comments below about 
#   setting up the default global variables to locate the default miRBase file. 
# For more information please visit http://mirspring.victorchang.edu.au.
#
# The following command will list the options available with this script:
#   perl Bam_to_Intermediate.pl
#
# Version 1.2
#	- Added: Non template addition feature (up to 2nt). This will continue to be worked up
#	-
#	- Fixed interrogating BAM file with just numerical chromosome ID (previously only work on records with a "chr" prefix
#	- Fixed samtools recognition
#	
# Version 1.3
#	- Fixed error detection in flanking sequence (was case sensitive in previous versions)
#	- Can use gff3 formatted files
#
#
#
# Author: David Thomas Humphreys
# Copyright (C) 2013, Victor Chang Cardiac Research Institute
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    see <http://www.gnu.org/licenses/> for GNU license.
####################################################################################################################################




####------SETUP A------------------------------------------------------------------------------------------------------------------------- 
#### Before running the software it is recommended to set the following SIX global variables. 

## The three variables below should contain the DIRECTORY of the miRbase files. 
my $Precursor_GFF_Dir  = "../miRbase_Files/";
my $Precursor_Seq_Dir  = "../miRbase_Files/";
my $Mature_Seq_Dir 	   = "../miRbase_Files/";
my $Samtools_Directory = '';	# ONLY Set this if the samtools directory if not set in the global path variable.

## This states the default mapping strategy you use to map data sets, either to a genome or a local version of miRbase. 
my $ReferenceFormat = 0;	# 0 = genome	1 = mirbase

## This states the default flanking sequence you provide in making a miRspring document. 
## If you are using unmodified miRbase files set this to 0 
my $Flank = 0;	# Flanking sequence on each side of precursor
####-------------------------------------------------------------------------------------------------------------------------------

my $parameters = {};


sub load_default_reference_files()
{
####------SETUP B------------------------------------------------------------------------------------------------------------------------- 
#### Before running the software it is recommended to set the default miRBase filenames


	$parameters->{gff} 			= "$Precursor_GFF_Dir/ TYPE FILE NAME HERE", if (! defined $parameters->{gff});
	$parameters->{precursor_seq_file} 	= "$Precursor_Seq_Dir/ TYPE FILE NAME HERE", if (! defined $parameters->{precursor_seq_file});
	$parameters->{mature_seq_file} 	= "$Precursor_Seq_Dir/ TYPE FILE NAME HERE", if (! defined $parameters->{precursor_seq_file});
####-------------------------------------------------------------------------------------------------------------------------------

	$parameters->{flank}     = $Flank, if (! defined $parameters->{flank});
	$parameters->{reference_format} = $ReferenceFormat, if (! defined $parameters->{reference_format});
	$parameters->{colour_mm} = 6, if (! defined $parameters->{colour_mm});
	$parameters->{mapQual}   = 0, if (! defined $parameters->{mapQual});
	$parameters->{nonTemplate} = 'OFF', if (! defined $parameters->{nonTemplate});
	$parameters->{mismatch} = 1, if (! defined $parameters->{mismatch});
	$parameters->{mismatch} = 0, if ($parameters->{nonTemplate} ne 'OFF');
	$parameters->{mismatch} = 2, if (($parameters->{nonTemplate} <0) && ($parameters->{nonTemplate} > 2)); # hopefully this will catch any errors
}


my $miRspring;			# Global hash to save final output
my $antisense_miRs;		# Global hash to save potential antisense or bidirectional miRs
my $ChromStats;         # Global hash to save chromosome numbers
my $Flanking_Seq;

# Test to see if samtools is installed correctly
{	my $Samtools_Version = '';
	my $cmd = "$Samtools_Directory"."samtools";
	my @SamtoolTest = `$cmd 2>&1`;		# The "2>&1" redirects STDERR to STDOUT. 

	#print "$cmd gave: @SamtoolTest";
	foreach (@SamtoolTest)
	{	
		if ($_ =~ m/^(Version: .*?)$/)
		{	$Samtools_Version = $1;	}
	}
	if ($Samtools_Version eq '')
	{ 	print "\nError with the following samtools path or installation: $cmd\nPlease check\n@SamtoolTest\n\n";	
		exit(0);
	}
	else {	print "\nDetected samtools installation, $Samtools_Version\n\n";		}
}	


# Test to see if enough parameters have been passed to script.
sub usage {
	print "\nUsage: $0 \n\n\t ";

    print "REQUIRED \n\t";
    print "-bam <input bam file> \n\t";
    print "-ml <minimum length> \n\t";
    print "-s <three letter species code>\n\t";

    print "-out <output file name with full path> \n\n\t";

    print "OPTIONAL \n\t";
    print "-mm <mismatches: '0' or '1', default is 1> \n\t";
	print "-nta <Non template addition: 1,2  or 'OFF', default is 'OFF'. When set -mm default is set to 0\n\t";
    print "-gff <Gene features file (gff), default is defined in script> \n\t";
    print "-mat <Mature sequence file, default is defined in script> \n\t";
    print "-pre <Precursor file, default is defined in script> \n\t";
    print "-flank <length of flanking sequence>\n\t";
	print "-ref <format of reference: \'0\' (genome) or \'1\' (custom), 0 (genome) is default\n\t";
	print "-cmm <color mismatch O(SOLiD data): 0 - 6, default is 6>\n\t";
	print "-mq <minimum mapping quality, default is 0>\n\n\n";

    exit(1);
}
if(scalar(@ARGV) == 0){
    usage();
}

# Parse the Command Line
&parse_command_line($parameters, @ARGV);

sub parse_command_line 
{
    my($parameters, @ARGV) = @_;

    my $next_arg;

    while(scalar @ARGV > 0)
	{	$next_arg = shift(@ARGV);
        if($next_arg eq "-s"){ $parameters->{species} = shift(@ARGV); }
        elsif($next_arg eq "-ml"){ $parameters->{min_length} = shift(@ARGV); }
        elsif($next_arg eq "-bam"){ $parameters->{bam} = shift(@ARGV); }
        elsif($next_arg eq "-mm"){ $parameters->{mismatch} = shift(@ARGV); }
		elsif($next_arg eq "-nta"){ $parameters->{nonTemplate} = shift(@ARGV); }
        elsif($next_arg eq "-gff"){ $parameters->{gff} = shift(@ARGV); }
        elsif($next_arg eq "-mat"){ $parameters->{mature_seq_file} = shift(@ARGV); }
        elsif($next_arg eq "-pre"){ $parameters->{precursor_seq_file} = shift(@ARGV); }
        elsif($next_arg eq "-out"){ $parameters->{output} = shift(@ARGV); }
        elsif($next_arg eq "-flank"){ $parameters->{flank} = shift(@ARGV); }
        elsif($next_arg eq "-ref"){ $parameters->{reference_format} = shift(@ARGV); }
        elsif ($next_arg eq "-cmm") {$parameters->{colour_mm} = shift(@ARGV); }
		elsif ($next_arg eq "-mq") {$parameters->{mapQual} = shift(@ARGV); }


        else { print "Invalid argument: $next_arg"; usage(); }
    }
	print "\nFlank is ".$parameters->{flank};

	# Check all essential options have been entered.
    my $Error_Log = '';
    $Error_Log .= "Species undefined (-s)\n\t", if (! defined $parameters->{species});
    $Error_Log .= "No BAM file defined (-bam)\n\t", if (! defined $parameters->{bam});
    $Error_Log .= "No output file defined (-out)\n\t", if (! defined $parameters->{output});
    $Error_Log .= "Minimum length not defined (-ml)\n\t", if (! defined $parameters->{min_length});
    if ($Error_Log ne '')
    {	print "\nThe following essential option(s) have not been defined:\n\t$Error_Log\n\n";
	usage();
}

	# Set the sequence files if user did not specify

	&load_default_reference_files($parameters->{species});


	print "\n\n**************    Parameters configures for BAM_to_miRspring script ********************************\n";
	while (my ($key, $value) = each %$parameters)
	{	print "$key\t$value\n";	}
	print "********************************************************************************************\n";

}

# Read header of BAM file to obtain reference layout and other information
my $HeaderFields;	# Will contain string of relevant header information
my $RefFields;	# A hash of all reference contig names
my $Headercmd = "$Samtools_Directory"."samtools view -H $parameters->{bam}";
my @HeaderInfo = `$Headercmd`;
($HeaderFields,$RefFields) = Extract_Header_From_Sam(@HeaderInfo);


# open mirbase precursor coordinate file
open PRECURSOR_GFF, $parameters->{gff} or die "Cannot open precursor file $parameters->{gff}";
#### Example of what a gff2 file looks like
# 1	.	miRNA	20669091	20669163	.	+	.	ACC="MI0000249"; ID="mmu-mir-206";
#### gff3 similar to gff2 but also has stem loop entries. 3rd column identified the entry:
# chr1	.	miRNA_primary_transcript	12425986	12426106	.	+	.	ID=MI0021869;Alias=MI0021869;Name=mmu-mir-6341
# chr1	.	miRNA	12426016	12426038	.	+	.	ID=MIMAT0025084;Alias=MIMAT0025084;Name=mmu-miR-6341;Derives_from=MI0021869

# Identify if file is gff2 or gff3 from file extension:
my $GFF_file_format = 2;	# Default is gff2 file type.
$GFF_file_format = 3, if ($parameters->{gff} =~ m/\.gff3$/);

my (@miR_Name, @miR_Start, @miR_End, @miR_Strand); 	# Storage for miR coordinates
my ($i, $last_chr, $Index) = (0,0,0);
my $Unique_miR_ID;					# Hash that stores miR name index ID
while (<PRECURSOR_GFF>)
{	chomp();	# Remove return character
	if ($_ =~ /^chr(.*?)$/) # Strip of "chr" in first column if it exists in gff file
	{	$_ = $1;	}
	my @Record = split("\t",$_);
	next, if (scalar(@Record) < 7);	# miRNA records should have 8 fields\
	next, if (($GFF_file_format == 3) && ($Record[2] eq 'miRNA'));	# Only want stem loop entries, not mature miRNA entries.

	if ($RefFields->{$Record[0]} == '')
	{	# Occasionally there is no chrY in bam file! This causes an error later on in the scripts, therefore make the entry here
		# Note that this "fix" will then result in an error from samtools, and can be ignored.
		$RefFields->{$RefFields->{"Max"}} = $Record[0];
		$RefFields->{$Record[0]} = $RefFields->{"Max"}++;
	}
		
	$Record[0] = $RefFields->{$Record[0]}, if (($Record[0] eq 'X') || ($Record[0] eq 'Y') || ($Record[0] eq 'M'));
			
	if ($Record[0] > $RefFields->{"Max"})
	{	print "\nSuspect miRNA genome coordinate from $parameters->{species} gff file : $Record[0], max should be ".$RefFields->{"Max"};
		next;
	}

	my $Stem_Loop_ID='';
	if ($GFF_file_format == 2)
	{	$Record[8] =~ m/\"(\w\w\w-\w\w\w-.*?)\"/;	
	#	$Record[8] = m/($parameters->{species}-.*?)\"/;		# This is potentially more accurate but may be (rare) cases where it could hinder.
		$Stem_Loop_ID = $1;
	}
	elsif ($GFF_file_format == 3)
	{	$Record[8] =~ m/.*?\=(\w\w\w-\w\w\w-.*?)$/;	
		$Stem_Loop_ID = $1;
	}

	$i = 0, if ($last_chr ne $Record[0]);
	$miR_Name[$Record[0]][$i] = $Stem_Loop_ID;
	$Unique_miR_ID->{$Stem_Loop_ID} = $Index;
	$miR_Start[$Record[0]][$i]= $Record[3];
	$miR_End[$Record[0]][$i]= $Record[4];
	$miR_Strand[$Record[0]][$i] = $Record[6];
	$i++;
	$last_chr = $Record[0];
	$Index++;
}
close PRECURSOR_GFF;
print "\n Uploaded miR coordinates from $last_chr (max is ".$RefFields->{"Max"}.")chromosomes/contigs/miR precursors.\nIndexed a total of $Index miR precurors\nLargest non numeric chromosome is $RefFields->{LargestNonNumeric}";

# Record all antisense pairs
my $Antisense_pairs;
for(my $i = 1; $i <= $last_chr; $i++)
{
	for(my $j = 0; $miR_Start[$i][$j] > 0; $j++)
	{	# I could assume the miRbase file is in chronological order and only check the neighbour. 
		# Instead will check all entries on current chromosome to see if they overlap at all.
		for (my $k = $j+1; $miR_Start[$i][$k] > 0; $k++)
		{
			if ((abs ($miR_Start[$i][$j] -  $miR_Start[$i][$k]) < 20) || (abs ($miR_End[$i][$j] -  $miR_End[$i][$k]) < 20))
			{	$Antisense_pairs->{$miR_Name[$i][$j]} = $miR_Name[$i][$k];
				$Antisense_pairs->{$miR_Name[$i][$k]} = $miR_Name[$i][$j];
			}
		}
	}
}

my $Precursor_Seqs;
open PRECURSOR_SEQ, $parameters->{precursor_seq_file} or die "Cannot open precursor file $parameters->{precursor_seq_file}";
while(<PRECURSOR_SEQ>)
{	if ($_ =~ m/^>(.*?),/)
	{       $Precursor_Seqs->{$1} = <PRECURSOR_SEQ>;	}
}
close PRECURSOR_SEQ;
print "\nUpload precursor sequences from $last_chr chromosomes";

if ($parameters->{reference_format} == 0)
{	# samtools -view SAMFILE:chr x-y
	for(my $i = 1; $i <= $last_chr; $i++)
	{       my $j = 0;
		my $Chrom = $i;
		if ($i > $RefFields->{'LargestNonNumeric'})
		{	$Chrom  = $RefFields->{$i};	
			
		}
		print "\nChr $RefFields->{$i} ";
		while($miR_Start[$i][$j] > 0)
		{	my $PreLabel = $RefFields->{"prelabel"};
			my $Command = "$Samtools_Directory"."samtools view $parameters->{bam} $PreLabel$Chrom:$miR_Start[$i][$j]-$miR_End[$i][$j]";
			print ".";
			my $Start = $miR_Start[$i][$j] - $parameters->{flank};
			my $End = $miR_End[$i][$j] + $parameters->{flank};
			my @sam = `$Command`;
			Extract_miRD_From_Sam($i, $Start, $End, $j, @sam);
			$j++;
		}
	}
}
else	# mirbase reference
{	my $Headercmd = "$Samtools_Directory"."samtools view -H $parameters->{bam}";
	my @HeaderInfo = `$Headercmd`;
	my $RefIndex;
	($HeaderFields,$RefIndex) = Extract_Header_From_Sam(@HeaderInfo);
	print "\nmiRBase reference with ".$RefIndex->{Max}. " entries\n";

	for(my $i = 1; $i < $RefIndex->{Max}; $i++)
	{	my $samcmd = "$Samtools_Directory"."samtools view $parameters->{bam} $RefIndex->{$i}";
		my @sam = `$samcmd`;
		print ".";
		Extract_From_miRBase_Reference($i, @sam);
	}
}


# Create the intermediate file!
open OUTPUT, ">$parameters->{output}", or die "cannot create output file: $parameters->{output}; $!";
print OUTPUT "$HeaderFields\n";
print OUTPUT "miRNA idex\tStart pos\tmiRNA name\tFreq\tLength\tERROR";		# Header
while(my ($key, $value) = each %$miRspring)
{	my @Display =  split ("\t",$key );
	print OUTPUT "\n$Display[0]\t$Display[1]\t$Display[2]\t$value\t$Display[3]\t$Display[4]\t";
}
print OUTPUT "\n";
close OUTPUT;

exit(0);

# Antisense data!
open OUTPUT, ">$parameters->{output}"."_antisense", or die "cannot create output file: antisense $parameters->{output}; $!";
print OUTPUT "$HeaderFields\n";
print OUTPUT "miRNA idex\tStart pos\tmiRNA name\tFreq\tLength\tERROR";		# Header
while(my ($key, $value) = each %$antisense_miRs)
{	my @Display =  split ("\t",$key );
	print OUTPUT "\n$Display[0]\t$Display[1]\t$Display[2]\t$value\t$Display[3]\t$Display[4]\t";
}
print OUTPUT "\n";
close OUTPUT;

exit(0);

sub Extract_miRD_From_Sam
{   	my $Chr = shift;
	my $Start = shift;
	my $End = shift;
	my $Index = shift;
	my @Sam_Input = @_;

	my $Name = $miR_Name[$Chr][$Index];
	my $Strand = $miR_Strand[$Chr][$Index];
	my $DebugCount = 0;
	my $checkthis = 0;
	my $AntisenseCount = 0;

#	print "\n$Name has a total counts of ".scalar(@Sam_Input);

	foreach(@Sam_Input)
	{	# 667_1107_1401   0       chr2    180133127       100     18M17H  *       0 0       TCAGCTGTTNNNNNNNNN      :75434111%%%%%))##      RG:Z:Library4_4 NH:i:1 CM:i:1  NM:i:8  CQ:Z:@@@@@@@@@*@@@@/@@@@-@@@@-@@>@*@@@@-        CS:Z:T02123211030231023023302010303131121
		# 1396_563_435_F3	0	hsa-mir-671	29	100	23M	*	0	23	AGGAAGCCCTGGAGGGGCTGGAG	26>8?&&@:?2*<>7/)5?9/-=	NM:i:0  CM:i:2	CS:Z:T32020230021222002321022220302000303	XI:Z:hsa-miR-671-5p|{hsa-miR-671-5p}|29_52|	CQ:Z:26>8?&&@:?2*<>7/)5?9/-=/%9%8<<%)927  IH:i:1

       	if ($_ =~ m/CM:i:(\d)/)
		{	next, if ($1 > $parameters->{colour_mm});	}

		my @Record = split("\t",$_);

		next, if ($Record[4] < $parameters->{mapQual});	# Check read has requested mapping quality.
		
		
		$Record[2] =~ m/$RefFields->{"prelabel"}(.*?)$/;
		my $Chr = $1;
		$Chr = $ChromStats->{$parameters->{species}}->{"M"}, if ($1 eq 'M');
		$Chr = $ChromStats->{$parameters->{species}}->{"X"}, if ($1 eq 'X');
		$Chr = $ChromStats->{$parameters->{species}}->{"Y"}, if ($1 eq 'Y');

		if ($Record[5] =~ m/(\d+)M.*?/)
		{	my $length = $1;
			my $Precursor_Pos = $Record[3] - $Start;
			my $Name_ID = $Name;

			my $RefSeq = $Record[9];
			my $AntisenseFlag = 0;

			if ($Record[1] == 16)
			{
				$Precursor_Pos = $End - $Record[3] - $length+1; 	# Reverse orientation
				$RefSeq = reverse($Record[9]);
				$RefSeq =~ tr/ATCGatcg/TAGCtagc/;
				if ($Strand eq '+')
				{
					$Name_ID = "$Name_ID"."AntiSense";
					$AntisenseFlag = 1;
				}
			}
			elsif ($Strand eq '-')
			{	$Name_ID = "$Name_ID"."AntiSense";
				$AntisenseFlag = 1;
			}

			my $Error;
			if ($_ =~ m/NM:i:(\d+)/)
			{   # miRspring only set up to accept 1 internal mismatch (editing event) or 2 mismatches at end of read (Non template addition events).
				if ($parameters->{nonTemplate} ne "OFF" )
				{	if (($1 > 0) && ($AntisenseFlag == 0))
					{	# Need to grab reference and then compare to the Bam field.
						my $StartPos = $Precursor_Pos;
						my $Reference_Seq = substr($Precursor_Seqs->{$Name},$StartPos,$length);
						$Error = IsMatch($Reference_Seq,$RefSeq);

						my @Error_Position_List = split(' ',$Error); # 10GA = nucleotide 10 changes from G to A
						
						

						# Now check to see if first "mismatch" is within a defined distance from end of tag.
						my $NTAflag = 0;
						foreach my $ErrorCode(@Error_Position_List)
						# if (($length - $Error_Position_List[0]) <= $parameters->{nonTemplate})
						{	
							$ErrorCode =~ /^(\d+)\w(\w)/;
							my $ErrorPosition = $1;
							my $Replace = lc($2);

							if (($length - $ErrorPosition) < $parameters->{nonTemplate})
							{	# we have a NTA. Consider all remaining sequences a NTA event
								if ($NTAflag == 0)
								{	$Error = "_$Replace";	}
								else	# only allowing a 2nt NTA, therefore if enter this part it must be last base
								{	$Error .= "$Replace";
									$NTAflag = 0;
								}
								$NTAflag = $ErrorPosition, if ($NTAflag == 0);
							}
							else
							{	$Error = '';	}
						}
						# else we ignore all errors and report 						
						if (($NTAflag !=0) && ($NTAflag - $length < 0))
						{   # The last base matched reference, convert it to lower case
							my $NTAseq = lc(substr($Reference_Seq, $length-1));
							$Error .= "$NTAseq";
						#	print "\n$Error from $Reference_Seq vs $RefSeq with length of $length @Error_Position_List", if (length($Error) > 3);
						}
						#print "\n$Error";
					}
				}
				elsif ($parameters->{mismatch} ne "OFF")
				{	$checkthis++;
					if (($1 > 1) || ($1 > $parameters->{mismatch}))
					{	next;	}
					elsif (($1 > 0) && ($AntisenseFlag == 0))
					{	# Need to grab reference and then compare to the Bam field.


						my $StartPos = $Precursor_Pos;
						my $Reference_Seq = substr($Precursor_Seqs->{$Name},$StartPos,$length);
						$Error = IsMatch($Reference_Seq,$RefSeq);
						my @Error_Position_List = split(' ',$Error);	# Should check what is returned. If people use 

						# Check to make sure number of errors is sensible. If greater than 5 then something very wrong... 
						# perhaps not using correct genome coordinates 
						# Perhaps not using correct flank file
						my $TotalErrorDetected = scalar(@Error_Position_List);
#						print "\nWarning large number of errors ($TotalErrorDetected) detected for $Name_ID\n$Reference_Seq\n$RefSeq", if ($TotalErrorDetected > 2);


						next, if (scalar @Error_Position_List > 1);   # This is a check as some aligners actually have more mismatches than they report in the NM:i field!?!
					}
					elsif (($1 > 0) && ($AntisenseFlag == 1))
					{	$Error = "1mm";	$AntisenseCount++; }	# Will need to implement a detailed error report
					elsif ($AntisenseFlag == 1)
					{	$AntisenseCount++;	}
				}
			}
			else
			{	print "\nCannot detect mismatch count";	}

			
			
			# At the moment the mismatch information is only used for sense transcripts. Stay tuned for antisense data.
			if ($AntisenseFlag == 0)
			{	$miRspring->{"$Unique_miR_ID->{$Name}\t$Precursor_Pos\t$Name_ID\t$length\t$Error"}++;
				$DebugCount++;	
			}
			elsif ($Antisense_pairs->{$Name} eq '')
			{	$antisense_miRs->{"$Unique_miR_ID->{$Name}\t$Precursor_Pos\t$Name_ID\t$length\t$Error"}++;	}

	    }
		else
		{	print "\nCannot obtain $Name length from $Record[5] in $_";	}
		
	}
#	print "\t\t........Detected $DebugCount counts for this hairpin. Check trace = $checkthis + $AntisenseCount";
}


sub IsAntisenseNeighbour()
{      my $Chr = shift;
	my $Index = shift;
	my $AntisenseNeighbour = 0;

	if ( (abs($miR_Start[$Chr][$Index] - $miR_Start[$Chr][$Index-1]) < 20) && (abs($miR_End[$Chr][$Index] - $miR_End[$Chr][$Index-1]) < 20) )
	{	$AntisenseNeighbour++; }
	if ( (abs($miR_Start[$Chr][$Index] - $miR_Start[$Chr][$Index+1]) < 20) && (abs($miR_End[$Chr][$Index] - $miR_End[$Chr][$Index+1]) < 20) )
       {       $AntisenseNeighbour++; }

	return($AntisenseNeighbour);
}
# multimapper->{Chr}->{Start}->{Length} = "ChrA:Start-Length,ChrB:Start-Length";

sub Extract_From_miRBase_Reference()
{
	my $Index = shift;
	my @Sam_Input = @_;
	my $Name;


	foreach(@Sam_Input)
	{	# QName			Flag	Rname		pos	mapQ CIGAR	RNext TLen
		# 1396_563_435_F3	0	hsa-mir-671	29	100	23M	*	0	23	AGGAAGCCCTGGAGGGGCTGGAG	26>8?&&@:?2*<>7/)5?9/-=	NM:i:0  CM:i:2	CS:Z:T32020230021222002321022220302000303	XI:Z:hsa-miR-671-5p|{hsa-miR-671-5p}|29_52|	CQ:Z:26>8?&&@:?2*<>7/)5?9/-=/%9%8<<%)927  IH:i:1
		my @Record = split("\t",$_);
		next, if ($Record[4] < $parameters->{mapQual});	# Check read has requested mapping quality.
		
		if ($Record[5] =~ m/(\d+)M.*?/)
		{	my $length = $1;
			$Name = $Record[2];

			my $StartPos = $Record[3]+$parameters->{flank}-1;

			my $Error;
			if ($Record[11] =~ m/NM:i:(\d+)/)
            {   #print " extracted $1";
                if (($1 > 1) || ($1 > $parameters->{mismatch}))
				{	next;	}
				elsif ($1 > 0)
                {	# Need to grab reference and then compare to the Bam field.
					my $Reference_Seq = substr($Precursor_Seqs->{$Name},$StartPos,$length);
					$Error = IsMatch($Reference_Seq,$Record[9]);
                }
            }
			$miRspring->{"$Unique_miR_ID->{$Name}\t$StartPos\t$Name\t$length\t$Error"}++;
	    }
		else
		{	print "\nCannot obtain $Name length from $Record[5] in $_";	}
	}
}

sub Extract_Header_From_Sam
{  	my @Sam_Input = @_;
	my $KeepThis = extract_miRspring_Settings();
	my $Ref_Index;	# Hash lookup for all references, whether they be chr or miRs.
	my $Ref_Count = 1;	# Keeps track of chromosome number (or miRbase entries)
	my $Obscure_Refs;	# Hash to record non-numerical indexes
	my $Chr_prelabel = 'chr';	# To help descriminate between bam files that use chr1 or just 1

	foreach(@Sam_Input)
	{
		if ($_ =~ m/^\@RG/) # @RG     ID:Library6_6   PL:SOLID        PU:2_4  LB:Library6     PI:0    DS:50F  DT:2012-06-01T15:42:33+1000     SM:Traude_Nico_Bene_3
		{ 	$KeepThis .= $_;	}
		elsif ($_ =~ m/^\@PG/) # @PG     ID:Victor_2012-06-01T15:42:33+1000_2_4  PN:LifeScope    VN:2.5
		{	$KeepThis .= $_;	}
		elsif ($_ =~ m/^\@CO/) #@CO     RG:Library6_6   IX:6    II:6    LD:Bene B       LT:Fragment     AT:Small RNA    DE:2012-06-04T16:10:19+1000     CU:18743080     CT:18743080     SP:RNA-Seq SD:      BX:0    TN:25   TX:50   EC:0    ER:0    CO:Administrator        PJ:Hentze Lab   SO:Traude Bene Nico     PN:Victor
		{	$KeepThis .= $_;	}
		elsif ($_ =~ m/^\@SQ\tSN:(.*?)\t.*?/)
		{	my $Chr_Description = $1;	
			if (($ReferenceFormat==0) && ($Chr_Description =~ m/chr(.*?)$/))
			{	my $Contig = $1;
				if ($Contig =~ /^\d+$/)
				{
					$Ref_Index->{$Contig} = $Contig;
					$Ref_Count++;	
				}
				else
				{	$Obscure_Refs->{$Contig} = 1;	}	# Save the name of non integer contigs		
			}
			else # miRbase mapping OR chromosomes without a 'chr' prelabel
			{	$Ref_Index->{$Ref_Count++} = $Chr_Description;
				$Chr_prelabel = '';
			}
		}
	}
	$KeepThis =~ s/"/\\"/g;
	$Ref_Index->{"Max"} = $Ref_Count;
	return ($KeepThis, $Ref_Index), if ($ReferenceFormat==1);

	$Ref_Index->{'LargestNonNumeric'} = ($Ref_Count -1);
	if (defined $Obscure_Refs->{'X'})
	{	$Ref_Index->{$Ref_Count} = 'X';
		$Ref_Index->{'X'} = $Ref_Count++;
	}
	if (defined $Obscure_Refs->{'Y'})
       {	$Ref_Index->{$Ref_Count} = 'Y'; 
		$Ref_Index->{'Y'} = $Ref_Count++;
	}
	if (defined $Obscure_Refs->{'M'}) # What about reference annotated at MT??
       {	$Ref_Index->{$Ref_Count} = 'M';
		$Ref_Index->{'M'} = $Ref_Count++;
	}


	while(my ($keys, $values) = each %$Obscure_Refs)
	{
		next, if (($keys eq 'X') || ($keys eq 'Y') || ($keys eq 'M'));
		$Ref_Index->{$Index} = $keys;
		$Ref_Index->{$keys} = $Ref_Count++;
	}
	$Ref_Index->{"Max"} = $Ref_Count;
	$Ref_Index->{"prelabel"} = $Chr_prelabel;
	return ($KeepThis, $Ref_Index);
}

sub extract_miRspring_Settings()
{	my $miRspringHeader ="\@MIRSPRING BAM_to_Intermediate.pl settings<br>";
        $miRspringHeader   .= "Species (-s): ".$parameters->{species}."<br>";
        $miRspringHeader   .= "Minimum length (-ml): ".$parameters->{min_length}."<br>";
        $miRspringHeader   .= "BAM filename (-bam): ".$parameters->{bam}."<br>";
        $miRspringHeader   .= "Number of mismatches (-mm): ".$parameters->{mismatch}."<br>";
	$miRspringHeader   .= "Non template addition (-nta): ".$parameters->{nonTemplate}."<br>";
        $miRspringHeader   .= "Number of flanking nucleotides (-flank): ".$parameters->{flank}."<br>";
        $miRspringHeader   .= "Reference format (-ref): ".$parameters->{reference_format}."<br>";
        $miRspringHeader   .= "Color mismatches (-cmm): ".$parameters->{colour_mm}."<br>";
		$miRspringHeader   .= "Minimum mapping quality (-mq): ".$parameters->{mapQual}."<br>";
		# Next three entries have a unique mixture of characters before the <br> so that the Intermediate_to_miRspring.pl script can read these parameters directly
        $miRspringHeader   .= "GFF coordinate file (-gff): ".$parameters->{gff}."\t <br>";
        $miRspringHeader   .= "Mature sequence file (-mat): ".$parameters->{mature_seq_file}."\t<br>";
        $miRspringHeader   .= "Hairpin precursor file (-pre): ".$parameters->{precursor_seq_file}.".<br>";
        $miRspringHeader   .= "Intermediate file (-out): ".$parameters->{output}."<br>\n";
		
		return $miRspringHeader;
}

sub IsMatch()
{	my $str_source = uc(shift);
	my $str_base = uc(shift);
        my $mismatch_pos;

	my @source = split //,$str_source;  #split first rather than substr
	my @base = split //, $str_base;

	for (my $i=0; $i < length($str_source); $i++)
        {       $mismatch_pos .= ($i+1)."$source[$i]$base[$i] ", if ($source[$i] ne $base[$i]);	}

	return $mismatch_pos;
}
