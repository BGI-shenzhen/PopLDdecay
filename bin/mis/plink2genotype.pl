#!/usr/bin/perl-w
use strict;

use Data::Dumper;
use Getopt::Long;

#############Befor  Start  , open the files ####################

sub usage
{
	print STDERR <<USAGE;
	2016-04-22       hewm\@genomics.cn

	Usage:		perl $0  -inPED in.ped -inMAP in.map  -outGenotype out.genotype

	Options
	-inPED        <s> :  Input the ped file path
	-inMAP        <s> :  Input the map file path
	-outGenotype  <s> :  Output the genotype path

	-help             :  show this help

USAGE
}

my ($help,$inPED,$inMAP,$outGenotype);


GetOptions(
	"help"=>\$help,
	"inPED:s"=>\$inPED,
	"inMAP:s"=>\$inMAP,
	"outGenotype:s"=>\$outGenotype,
);

if( defined($help) || !defined($outGenotype)   || !defined($inPED)   || !defined($inMAP)  )
{
	usage;
	exit ;
}


open (IA,"$inPED") || die "input file can't open $!";
open (IB,"$inMAP") || die "input file can't open $!";
open (OA,">$outGenotype") || die "input file can't open $!";

my @Data;

my $sampleNum=0;
while(<IA>)
{
	chomp ;
	my @inf=split ;
	push @Data,\@inf;
	$sampleNum++;
}
close IA;
$sampleNum--;
my $Flag=6;
my %number2base=();
$number2base{"0"}="-";
$number2base{"1"}="A";
$number2base{"2"}="C";
$number2base{"3"}="G";
$number2base{"4"}="T";


while(<IB>)
{
	chomp ;
	my @inf=split ;
	my $chr=$inf[0];
	my $site=abs($inf[-1]);
	my $Genotype="";
	my $SedFlag=$Flag+1;
	   my $A=$Data[0]->[$Flag];
	   if (exists $number2base{$A})
	   {
		   $Genotype=$number2base{$A};
	   }
	   else
	   {
		    $Genotype=$A;
	   }
	    $A=$Data[0]->[$SedFlag];
	   if (exists $number2base{$A})
	   {
		   $Genotype=$Genotype." ".$number2base{$A};
	   }
	   else
	   {
		    $Genotype=$Genotype." ".$A;
	   }

	foreach my $samID (1..$sampleNum)
	{
	    $A=$Data[$samID]->[$Flag];
	   if (exists $number2base{$A})
	   {
		   $Genotype==$Genotype." ".$number2base{$A};
	   }
	   else
	   {
		    $Genotype=$Genotype." ".$A;
	   }
	    $A=$Data[$samID]->[$SedFlag];
	   if (exists $number2base{$A})
	   {
		   $Genotype=$Genotype." ".$number2base{$A};
	   }
	   else
	   {
		    $Genotype=$Genotype." ".$A;
	   }
	}

	$Flag++;
	$Flag++;
	print OA "$chr\t$site\t$Genotype\n";
}
close IB ;
close OA ;

######################swiming in the sky and flying in the sea #############################
