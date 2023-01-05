#!/usr/bin/perl

###downloading the CosmicCoding and CosmicNonCoding vcf file from Sanger database#####


system ("wget ftp://ngs.sanger.ac.uk/production/cosmic/Cosmic*.vcf.gz");
system ("gunzip *.gz");
system ("cat CosmicCoding*  CosmicNonCoding*");
system ("rm  CosmicCoding*  CosmicNonCoding*");
