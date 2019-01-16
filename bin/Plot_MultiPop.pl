#!/usr/bin/perl -w
use strict;

use Data::Dumper;
use Getopt::Long;

#############Befor  Start  , open the files ####################

sub usage
{
	print STDERR <<USAGE;
	2016-07-04       hewm\@genomics.cn

	Usage:  perl $0  -inList  LDdecayResult.list  -output OUT

	Options
	-inList    <s> :  Muti LDDecat Stat File of muti Pop.
	                  Format:[FilePath  PopID]
	-output    <s> :  Output Figure File Prefix

	-bin1      <n> :  the size bin for mean r^2/D' of Short Dist[10]
	-bin2      <n> :  the size bin for mean r^2/D' of Long Dist [100]
	-break     <n> :  break point to distinguish Short or Long Dist[100] 
	-maxX      <n> :  max X coordinate Dist to plot LDdecay[kb] [maxDist]
	-measure   <s> :  use the [r2/D/both] to measure LD decay [r2]
	-method    <s> :  Plot methold (MeanBin/HW/MedianBin/PercentileBin)[MeanBin]
	-percent   <f> :  percent ratio(0-1) for PercentileBin method [0.5]
	-keepR         :  keep the R script for draw the LDdecay Fig

	-help          :  show this help

USAGE
}

my ($help,$inList,$output,$bin1,$bin2,$break,$keepR,$maxX,$measure,$method,$percent);


GetOptions(
	"help"=>\$help,
	"inList:s"=>\$inList,
	"output:s"=>\$output,
	"bin1:s"=>\$bin1,
	"bin2:s"=>\$bin2,
	"keepR"=>\$keepR,
	"break:s"=>\$break,
	"maxX:s"=>\$maxX,
	"measure:s"=>\$measure,
	"method:s"=>\$method,
	"percent:s"=>\$percent,
);

if( defined($help) || (!defined($output))   ||  (!defined($inList)) )
{
	usage;
	exit ;
}


$bin1||=10;
$bin2||=100;
$break||=100;
$method||="MeanBin";
$percent||=0.5;
$measure||="r2"; my $RowPlot=2;
if ($measure eq  "D"  ||  $measure eq  "D'" )  { $RowPlot=3 ;}
elsif ($measure eq  "b"  ||  $measure eq  "both"  ||   $measure eq  "Both" ||   $measure eq  "B")  { $RowPlot=4 ;}

if ( $method eq "PercentileBin" &&  $percent==0.5) {$method="MedianBin" ;}


################ Do what you want to do #######################

my $bin=$bin2 ;
my $Small_bin=$bin1 ;

my $maxX_tmp=0;
my $maxY_tmp=0;

my @PopIDAryy=();

my $Rshell="";




########################################################################################################################

my @ColArry ; 
$ColArry[0]="maroon";
$ColArry[1]="black";
$ColArry[2]="Darkblue";
$ColArry[3]="Purple";
$ColArry[4]="DarkGreen";
$ColArry[5]="DarkOrange3";
$ColArry[6]="DimGrey";
$ColArry[7]="Brown";
$ColArry[8]="Orange";
$ColArry[9]="Cyan";
$ColArry[10]="Grey";
$ColArry[11]="LightSkyBlue";
$ColArry[12]="Gold";
$ColArry[13]="IndianRed2";
$ColArry[14]="DeepSkyBlue2";
$ColArry[15]="Chartreuse2";
$ColArry[16]="Orchid";
$ColArry[17]="greenyellow";
$ColArry[18]="SlateGrey";
$ColArry[19]="LightSlateGray";
$ColArry[20]="Gold1";
$ColArry[21]="SaddleBrown";
$ColArry[22]="MidnightBlue";
$ColArry[23]="NavyBlue";
$ColArry[24]="YellowGreen";
$ColArry[25]="CornflowerBlue";
$ColArry[26]="MediumBlue";
$ColArry[27]="Blue";
$ColArry[28]="RoyalBlue";
$ColArry[29]="DeepPink";
$ColArry[30]="DeepSkyBlue1";
$ColArry[31]="MediumTurquoise";
$ColArry[32]="Turquoise";
$ColArry[33]="SpringGreen3";
$ColArry[34]="CadetBlue";
$ColArry[35]="DodgerBlue";
$ColArry[36]="MediumAquamarine";
$ColArry[37]="Aquamarine";
$ColArry[38]="Green";
$ColArry[39]="DarkOliveGreen";
$ColArry[40]="DarkSeaGreen";
$ColArry[41]="SeaGreen";
$ColArry[42]="MediumSeaGreen";
$ColArry[43]="LightSeaGreen";
$ColArry[44]="PaleGreen";
$ColArry[45]="SpringGreen";
$ColArry[46]="LawnGreen";
$ColArry[47]="Chartreuse";
$ColArry[48]="MedSpringGreen";
$ColArry[49]="GreenYellow";
$ColArry[50]="LimeGreen";
$ColArry[51]="OrangeRed2";
$ColArry[52]="DeepSkyBlue";
$ColArry[53]="Magenta3";
$ColArry[54]="SteelBlue";
$ColArry[55]="LightSteelBlue";
$ColArry[56]="LightBlue";
$ColArry[57]="PowderBlue";
$ColArry[58]="PaleTurquoise";
$ColArry[59]="DarkTurquoise";
$ColArry[60]="ForestGreen";
$ColArry[61]="OliveDrab";
$ColArry[62]="DarkKhaki";
$ColArry[63]="PaleGoldenrod";
$ColArry[64]="LtGoldenrodYello";
$ColArry[65]="Pink";
$ColArry[66]="LightYellow";
$ColArry[67]="LightGoldenrod";
$ColArry[68]="goldenrod";
$ColArry[69]="DarkGoldenrod";
$ColArry[70]="RosyBrown";
$ColArry[71]="IndianRed";
$ColArry[72]="Sienna";
$ColArry[73]="Peru";
$ColArry[74]="DarkSlateGray";
$ColArry[75]="Burlywood";
$ColArry[76]="SkyBlue";
$ColArry[77]="Beige";
$ColArry[78]="Wheat";
$ColArry[79]="SandyBrown";
$ColArry[80]="Tan";
$ColArry[81]="Chocolate";
$ColArry[82]="Firebrick";
$ColArry[83]="DarkSalmon";
$ColArry[84]="Salmon";
$ColArry[85]="LightSalmon";
$ColArry[86]="Coral";
$ColArry[87]="LightCoral";
$ColArry[88]="Tomato";
$ColArry[89]="OrangeRed";
$ColArry[90]="HotPink";
$ColArry[91]="LightPink";
$ColArry[92]="PaleVioletRed";
$ColArry[93]="Maroon";
$ColArry[94]="MediumVioletRed";
$ColArry[95]="VioletRed";
$ColArry[96]="Violet";
$ColArry[97]="Plum";
$ColArry[98]="DarkSlateBlue";
$ColArry[99]="SlateBlue";
$ColArry[100]="MediumSlateBlue";
$ColArry[101]="LightSlateBlue";
$ColArry[102]="MediumOrchid";
$ColArry[103]="DarkOrchid";
$ColArry[104]="DarkViolet";
$ColArry[105]="BlueViolet";
$ColArry[106]="MediumPurple";
$ColArry[107]="Thistle";
$ColArry[108]="Snow1";
$ColArry[109]="Snow2";
$ColArry[110]="Snow3";
$ColArry[111]="Snow4";
$ColArry[112]="Seashell1";
$ColArry[113]="Seashell2";
$ColArry[114]="Seashell3";
$ColArry[115]="Seashell4";
$ColArry[116]="AntiqueWhite1";
$ColArry[117]="AntiqueWhite2";
$ColArry[118]="AntiqueWhite3";
$ColArry[119]="AntiqueWhite4";
$ColArry[120]="Bisque1";


########################################################################################################################


###########################################  MeanBin  ##################################################################

if ($method  eq "MeanBin" || $method  eq "Mean" ||  $method  eq "mean" )  
{

	if  ($inList =~s/\.gz$/\.gz/)
	{
		open LIST,"gzip -cd  $inList | "  || die "input file can't open $!" ;
	}
	else
	{
		open LIST,"$inList"  || die "input file can't open $!" ;
	}

	while(<LIST>)
	{
		my $FileThis=$_;
		chomp $FileThis ;
		my @SplltPath=split /\s+/,$FileThis ;
		if ( $#SplltPath<1 )
		{
			print "FileList Format wrong, should be (Two columns) :\n";
			print "Stat.FilePath.stat.gz   PopulationIDA\n";
			exit(1);
		}

		my %hash_cout=();
		my %hash_coutV2=();
		my %hash_RR_sum=();
		my %hash_D_sum=();

		my %Small_cout=();
		my %Small_coutV2=();
		my %Small_RR_sum=();
		my %Small_D_sum=();

		$FileThis=$SplltPath[0];
		my $PopID=$SplltPath[1];
		push  @PopIDAryy , $PopID;


		&ReadFile2HashV2($FileThis,$bin1,$bin2,\$maxX_tmp ,\%hash_cout, \%hash_coutV2, \%hash_D_sum,\%hash_RR_sum, \%Small_cout,\%Small_coutV2,\%Small_D_sum, \%Small_RR_sum);

		open OA,">$output.$PopID" || die "output file can't open $!" ;

		print OA "#Dist\tMean_r^2\tMean_D'\tSum_r^2\tSum_D'\tNumberPairs\n";

		my $Cal_D=0;


		foreach my $k (sort {$a<=>$b} keys %Small_cout)
		{
			my $mean_R=$Small_RR_sum{$k}/$Small_cout{$k};
			if  ($mean_R> $maxY_tmp) {	$maxY_tmp=$mean_R ;	}
			if (exists $Small_coutV2{$k})
			{
				my $mean_D=$Small_D_sum{$k}/$Small_coutV2{$k};
				print OA ($k+1)*$Small_bin,"\t$mean_R\t$mean_D\t$Small_RR_sum{$k}\t$Small_D_sum{$k}\t$Small_cout{$k}\n";
				$Cal_D=1;
			}
			else
			{
				print OA ($k+1)*$Small_bin,"\t$mean_R\tNA\t$Small_RR_sum{$k}\tNA\t$Small_cout{$k}\n";
			}
		}

		foreach my $k (sort {$a<=>$b} keys %hash_cout)
		{
			my $mean_R=$hash_RR_sum{$k}/$hash_cout{$k};
			if  ($mean_R> $maxY_tmp) {	$maxY_tmp=$mean_R ;	}
			if  ( exists  $hash_coutV2{$k} )
			{
				my $mean_D=$hash_D_sum{$k}/$hash_coutV2{$k};
				print OA ($k+1)*$bin ,"\t$mean_R\t$mean_D\t$hash_RR_sum{$k}\t$hash_D_sum{$k}\t$hash_cout{$k}\n";
				$Cal_D=1;
			}
			else
			{
				print OA ($k+1)*$bin ,"\t$mean_R\tNA\t$hash_RR_sum{$k}\tNA\t$hash_cout{$k}\n";
			}
		}

		close OA ;

		if  ( $RowPlot>2  &&  $Cal_D==0 )
		{
			print "The $PopID input result only with measure [r2], can't find the D', To calculate the D' result , you should run the core program [PopLDdecay] with parameter [ -OutType 2]\n";
			print "To run the measure r2 LD decay, please set the parameter [ -measure r2 ] \n";
			exit ;	
		}

	}

	close LIST ;


	$maxX_tmp=int($maxX_tmp/1000);
	$maxX||=$maxX_tmp;
	if  ($maxY_tmp>0.96)
	{
		$maxY_tmp=1;
	}


	my $PopNumber=$#PopIDAryy ;

	my $PopName=$PopIDAryy[0];
	my $Plot=<<LOVE;
read.table("$output.$PopName")->E$PopName;
plot(E$PopName\[,1\]/1000,E$PopName\[,2\],type="l",col="$ColArry[0]",main="LD decay",xlab="Distance(Kb)",xlim=c(0,$maxX),ylim=c(0,$maxY_tmp),ylab=expression(r^{2}),bty="n")
LOVE

	my $legendCol="\"$ColArry[0]\"";
	my $legendName="\"$PopName\"";
	my $legendlty=1;

	for (my $IDE=1; $IDE<= $PopNumber; $IDE++)
	{

		$PopName=$PopIDAryy[$IDE];
		$Plot.=<<LOVE;
read.table("$output.$PopName")->E$PopName;
lines(E$PopName\[,1\]/1000,E$PopName\[,2\],col="$ColArry[$IDE]")
LOVE

		$legendCol.=",\"$ColArry[$IDE]\"";
		$legendName.=",\"$PopName\"";
		$legendlty.=",1";
	}


	$Rshell=<<HEWM;
pdf("$output.pdf")
$Plot
legend("topright",c($legendName),col=c($legendCol),cex=1,lty=c($legendlty),bty="n");
dev.off()
png("$output.png")
$Plot
legend("topright",c($legendName),col=c($legendCol),cex=1,lty=c($legendlty),bty="n");
dev.off()

HEWM


if ( $RowPlot==3 )
{

	$PopName=$PopIDAryy[0];
	$Plot=<<LOVE;
read.table("$output.$PopName")->E$PopName;
plot(E$PopName\[,1\]/1000,E$PopName\[,3\],type="l",col="$ColArry[0]",main="LD decay",xlab="Distance(Kb)",xlim=c(0,$maxX),ylim=c(0,1.0),ylab="D'",bty="n")
LOVE

	for (my $IDE=1; $IDE<= $PopNumber; $IDE++)
	{
		$PopName=$PopIDAryy[$IDE];
		$Plot.=<<LOVE;
read.table("$output.$PopName")->E$PopName;
lines(E$PopName\[,1\]/1000,E$PopName\[,3\],col="$ColArry[$IDE]")
LOVE
	}

	$Rshell=<<HEWM;
pdf("$output.pdf")
$Plot
legend("topright",c($legendName),col=c($legendCol),cex=1,lty=c($legendlty),bty="n");
dev.off()
png("$output.png")
$Plot
legend("topright",c($legendName),col=c($legendCol),cex=1,lty=c($legendlty),bty="n");
dev.off()

HEWM
	}
	elsif ( $RowPlot==4 )
	{

		$PopName=$PopIDAryy[0];
		my $Plot2=<<LOVE;
read.table("$output.$PopName")->E$PopName;
plot(E$PopName\[,1\]/1000,E$PopName\[,3\],type="l",col="$ColArry[0]",xlab="",xlim=c(0,$maxX),ylim=c(0,1.0),ylab="",bty="n",lty=2,yaxt="n")
axis(4,col="red",col.ticks="red",col.axis="red",lty=2)
mtext("D'",side=4,line=-1,col="red")
LOVE

		for (my $IDE=1; $IDE<= $PopNumber; $IDE++)
		{
			$PopName=$PopIDAryy[$IDE];
			$Plot2.=<<LOVE;
read.table("$output.$PopName")->E$PopName;
lines(E$PopName\[,1\]/1000,E$PopName\[,3\],col="$ColArry[$IDE]", lty=2)
LOVE
		}

		$Rshell=<<HEWM;
pdf("$output.pdf")
$Plot
legend("topright",c($legendName),col=c($legendCol),cex=1,lty=c($legendlty),bty="n");
par(new=T)
$Plot2
dev.off()
png("$output.png")
$Plot
legend("topright",c($legendName),col=c($legendCol),cex=1,lty=c($legendlty),bty="n");
par(new=T)
$Plot2
dev.off()

HEWM


	}

}


###########################################    PercentileBin  HW #############################################################

elsif ( $method eq "HW" || $method eq "MedianBin" || $method eq "PercentileBin" || $method eq "Percentile")
{

	if  ($inList =~s/\.gz$/\.gz/)
	{
		open LIST,"gzip -cd  $inList | "  || die "input file can't open $!" ;
	}
	else
	{
		open LIST,"$inList"  || die "input file can't open $!" ;
	}

	my $SamplesSize=0;

	while(<LIST>)
	{
		my $FileThis=$_;
		chomp $FileThis ;
		my @SplltPath=split /\s+/,$FileThis ;
		if ( $#SplltPath<1 )
		{
			print "FileList Format wrong, should be (Two columns) :\n";
			print "Stat.FilePath.stat.gz   PopulationIDA\n";
			exit(1);
		}

		my %hashRR=(); 
		my %hashDD=();
		my %hashTT=();
		my $FileType=0;

		$FileThis=$SplltPath[0];
		my $PopID=$SplltPath[1];
		push  @PopIDAryy , $PopID;

		&ReadFile2HashV1( $FileThis,\$maxX_tmp,\$SamplesSize,\$FileType,\%hashDD,\%hashRR,\%hashTT);


		if  ( $RowPlot>2  &&  $FileType< 2 )
		{
			print "The $PopID input result only with measure [r2], can't find the D', To calculate the D' result , you should run the core program [PopLDdecay] with parameter [-OutType 4]\n";
			print "To run the measure r2 LD decay, please set the parameter [ -measure r2 ] \n";
			exit ;	
		}



		if  ( $method eq "HW" )
		{

			open OD,">$output.$PopID.rdtmp" || die "output file can't open $!";
			print OD "#Dist\tR^2\tR^2_count\tD'\tD_count\n";
			my  $ALL_Pairwise=0;
			foreach my $k ( sort {$a<=>$b} keys %hashTT)
			{
				my $cc=$hashTT{$k};
				foreach my $jj ( sort {$a<=>$b} keys  %$cc )
				{
					my $DD=$hashDD{$k}{$jj} ;  $DD||=0;
					my $RR=$hashRR{$k}{$jj} ;  $RR||=0;
					print OA $k,"\t$jj\t$RR\t$jj\t$DD\n";
					$ALL_Pairwise+=$RR;
				}
			}
			close OD;

			if  ( $ALL_Pairwise>100000 || $ALL_Pairwise< 0 )
			{
				print "\t\tHW methold is only suitable for small data, we recommend the total Pairwise Number <100000\n";
				print "\t\tHowever Now all Pairwise Number is $ALL_Pairwise ($PopID)\t; we suggest you may use the [MeanBin] to plot\n";
				exit (1);
			}
			%hashTT=(); undef  %hashTT ;

		}
		else
		{
###########################################    PercentileBin  HW #############################################################
			my %newRRS=();
			%hashTT=();
			foreach my $k ( keys %hashRR)
			{  
				my $cc=$hashRR{$k};
				if  ($k<$break)
				{
					my $key_bin=int($k/$bin1);
					foreach my $jj (keys %$cc)
					{
						$newRRS{$key_bin}{$jj}+=$hashRR{$k}{$jj};
					}
				}
				else
				{
					my $key_bin=int($k/$bin2);
					foreach my $jj (keys %$cc)
					{
						$hashTT{$key_bin}{$jj}+=$hashRR{$k}{$jj};
					}
				}
			}
			%hashRR=(); undef %hashRR;
			open ORR,">$output.$PopID.rrtmp" || die "output file can't open $!";

			foreach  my $k ( sort {$a<=>$b}  keys %newRRS)
			{
				my $cc=$newRRS{$k};
				my $Total=0;
				foreach my $jj (keys %$cc)
				{
					$Total+=$newRRS{$k}{$jj};
				}
				my $Cutoff=int($Total* $percent);
				$Total=0;
				foreach my $jj (  sort {$a<=>$b} keys %$cc)
				{
					$Total+=$newRRS{$k}{$jj};
					if  ($Total>=$Cutoff)
					{
						print ORR  $k*$bin1,"\t$jj\n";
						last ;
					}
				}
			}

			%newRRS=(); undef %newRRS;
			foreach  my $k ( sort {$a<=>$b}  keys %hashTT)
			{
				my $cc=$hashTT{$k};
				my $Total=0;
				foreach my $jj (keys %$cc)
				{
					$Total+=$hashTT{$k}{$jj};
				}
				my $Cutoff=int($Total* $percent);
				$Total=0;
				foreach my $jj (  sort {$a<=>$b} keys %$cc)
				{
					$Total+=$hashTT{$k}{$jj};
					if  ($Total>=$Cutoff)
					{
						print ORR  $k*$bin2,"\t$jj\n";
						last ;
					}
				}
			}

			close ORR;

			%hashTT=(); undef %hashTT;
#			my %newDDL=();			my %newDDS=();

			foreach my $k ( keys %hashDD)
			{  
				my $cc=$hashDD{$k};
				if  ($k<$break)
				{
					my $key_bin=int($k/$bin1);
					foreach my $jj (keys %$cc)
					{
						$hashRR{$key_bin}{$jj}+=$hashDD{$k}{$jj};
					}
				}
				else
				{
					my $key_bin=int($k/$bin2);
					foreach my $jj (keys %$cc)
					{
						$hashTT{$key_bin}{$jj}+=$hashDD{$k}{$jj};
					}
				}
			}

			%hashDD=(); undef %hashDD;
			open ODD,">$output.$PopID.ddtmp" || die "output file can't open $!";
			foreach  my $k ( sort {$a<=>$b}  keys %hashRR)
			{
				my $cc=$hashRR{$k};
				my $Total=0;
				foreach my $jj (keys %$cc)
				{
					$Total+=$hashRR{$k}{$jj};
				}
				my $Cutoff=int($Total* $percent);
				$Total=0;
				foreach my $jj (  sort {$a<=>$b} keys %$cc)
				{
					$Total+=$hashRR{$k}{$jj};
					if  ($Total>=$Cutoff)
					{
						print ODD  $k*$bin1,"\t$jj\n";
						last;
					}
				}
			}

			%hashRR=(); undef %hashRR;

			foreach  my $k ( sort {$a<=>$b}  keys %hashTT)
			{
				my $cc=$hashTT{$k};
				my $Total=0;
				foreach my $jj (keys %$cc)
				{
					$Total+=$hashTT{$k}{$jj};
				}
				my $Cutoff=int($Total* $percent);
				$Total=0;
				foreach my $jj (  sort {$a<=>$b} keys %$cc)
				{
					$Total+=$hashTT{$k}{$jj};
					if  ($Total>=$Cutoff)
					{
						print ODD  $k*$bin2,"\t$jj\n";
						last ;
					}
				}
			}
			close ODD;

			%hashTT=(); undef %hashTT;

		}  

	} 	###### end while

	$maxX_tmp=int($maxX_tmp/1000);
	$maxX||=$maxX_tmp;
	my $PopNumber=$#PopIDAryy ;
	my $PopName=$PopIDAryy[0];


	if ($method eq "MedianBin" || $method eq "PercentileBin" || $method eq "Percentile" )
	{

#############################################
		$percent=$percent*100;
		$percent="$percent%";
		$method=~s/Bin$//g;

		$PopName=$PopIDAryy[0];

		my $Plot=<<LOVE;
read.table("$output.$PopName.rrtmp")->E$PopName;
plot(E$PopName\[,1\]/1000,E$PopName\[,2\],type="l",col="$ColArry[0]",main="LD decay $method ($percent) ",xlab="Distance(Kb)",xlim=c(0,$maxX),ylim=c(0,1),ylab=expression(r^{2}),bty="n")
LOVE

		my $legendCol="\"$ColArry[0]\"";
		my $legendName="\"$PopName\"";
		my $legendlty=1;

		for (my $IDE=1; $IDE<= $PopNumber; $IDE++)
		{

			$PopName=$PopIDAryy[$IDE];
			$Plot.=<<LOVE;
read.table("$output.$PopName.rrtmp")->E$PopName;
lines(E$PopName\[,1\]/1000,E$PopName\[,2\],col="$ColArry[$IDE]")
LOVE

			$legendCol.=",\"$ColArry[$IDE]\"";
			$legendName.=",\"$PopName\"";
			$legendlty.=",1";
		}




		$PopName=$PopIDAryy[0];
		my $Plot2=<<LOVE;
	read.table("$output.$PopName.ddtmp")->E$PopName;
	plot(E$PopName\[,1\]/1000,E$PopName\[,2\],type="l",col="$ColArry[0]",main="LD decay $method ($percent) ",xlab="Distance(Kb)",xlim=c(0,$maxX),ylim=c(0,1.0),ylab="D'",bty="n")
LOVE
		my $Plot3="";
		for (my $IDE=1; $IDE<= $PopNumber; $IDE++)
		{

			$PopName=$PopIDAryy[$IDE];
			$Plot3.=<<LOVE;
read.table("$output.$PopName.ddtmp")->E$PopName;
lines(E$PopName\[,1\]/1000,E$PopName\[,2\],col="$ColArry[$IDE]")
LOVE
		}



		if ( $RowPlot==2)
		{

			$Rshell=<<HEWM;
pdf("$output.pdf")
$Plot
legend("topright",c($legendName),col=c($legendCol),cex=1,lty=c($legendlty),bty="n");
dev.off()
png("$output.png")
$Plot
legend("topright",c($legendName),col=c($legendCol),cex=1,lty=c($legendlty),bty="n");
dev.off()

HEWM
		}
		elsif ($RowPlot==3)
		{

			$Rshell=<<HEWM;
pdf("$output.pdf")
$Plot2
$Plot3
legend("topright",c($legendName),col=c($legendCol),cex=1,lty=c($legendlty),bty="n");
dev.off()
png("$output.png")
$Plot2
$Plot3
legend("topright",c($legendName),col=c($legendCol),cex=1,lty=c($legendlty),bty="n");
dev.off()

HEWM
		}
		elsif ($RowPlot==4)
		{
			$Plot3=~s/col=/lty=2,col=/g;
			$PopName=$PopIDAryy[0];
			$Plot2=<<LOVE;
	read.table("$output.$PopName.ddtmp")->E$PopName;
	plot(E$PopName\[,1\]/1000,E$PopName\[,2\],type="l",col="$ColArry[0]",xlab="",xlim=c(0,$maxX),ylim=c(0,1.0),ylab="",bty="n",yaxt="n",lty=2)
	axis(4,col="red",col.ticks="red",col.axis="red",lty=2)
	mtext("D'",side=4,line=-1,col="red")
LOVE

			$Rshell=<<HEWM;
pdf("$output.pdf")
$Plot
legend("topright",c($legendName),col=c($legendCol),cex=1,lty=c($legendlty),bty="n");
par(new=T)
$Plot2
$Plot3
dev.off()
png("$output.png")
$Plot
legend("topright",c($legendName),col=c($legendCol),cex=1,lty=c($legendlty),bty="n");
par(new=T)
$Plot2
$Plot3
dev.off()
HEWM

		}
	}
	else
	{
		######################################   HW    ####################################################

		$PopName=$PopIDAryy[0];
		my $Plot=<<LOVE;
n=$SamplesSize;
read.table("$output.$PopName.rdtmp")->data;
sumCount<-sum(data[,3])
Row<-nrow(data);
dist  <- vector(length = sumCount)
rsq  <- vector(length = sumCount)
count=1;
for(i in 1:Row)
{
  for (j in 1:data[i,3])
  {
	 dist[count]<-data[i,1]
	 rsq[count]<-data[i,2]
	 count <- count + 1
  }
}
file <- data.frame(dist,rsq);
Cstart <- c(C=0.1)
modelC <- nls(rsq ~ ((10+C*dist)/((2+C*dist)*(11+C*dist)))*(1+((3+C*dist)*(12+12*C*dist+(C*dist)^2))/(n*(2+C*dist)*(11+C*dist))),data=file, start=Cstart, control=nls.control(maxiter=100))
rho <- summary(modelC)\$parameters[1]
newrsq <- ((10+rho*file\$dist)/((2+rho*file\$dist)*(11+rho*file\$dist)))*(1+((3+rho*file\$dist)*(12+12*rho*file\$dist+(rho*file\$dist)^2))/(n*(2+rho*file\$dist)*(11+rho*file\$dist)))

newfile <- data.frame(file\$dist, newrsq)
newfile <- newfile[order(newfile\$file.dist),]
newfile$PopName  <- base::unique(newfile)

plot(newfile$PopName\$dist/1000, newfile$PopName\$rsq, type="l",col="$ColArry[0]",main="LD decay", xlab="Distance (Kb)", ylab=expression(r^{2}),bty="n",xlim=c(0,$maxX))

LOVE



		my $legendCol="\"$ColArry[0]\"";
		my $legendName="\"$PopName\"";
		my $legendlty=1;


		for (my $IDE=1; $IDE<= $PopNumber; $IDE++)
		{

			$PopName=$PopIDAryy[$IDE];
			$Plot.=<<LOVE;
read.table("$output.$PopName.rdtmp")->data;
sumCount<-sum(data[,3])
Row<-nrow(data);
dist  <- vector(length = sumCount)
rsq  <- vector(length = sumCount)
count=1;
for(i in 1:Row)
{
  for (j in 1:data[i,3])
  {
	 dist[count]<-data[i,1]
	 rsq[count]<-data[i,2]
	 count <- count + 1
  }
}
file <- data.frame(dist,rsq);
Cstart <- c(C=0.1)
modelC <- nls(rsq ~ ((10+C*dist)/((2+C*dist)*(11+C*dist)))*(1+((3+C*dist)*(12+12*C*dist+(C*dist)^2))/(n*(2+C*dist)*(11+C*dist))),data=file, start=Cstart, control=nls.control(maxiter=100))
rho <- summary(modelC)\$parameters[1]
newrsq <- ((10+rho*file\$dist)/((2+rho*file\$dist)*(11+rho*file\$dist)))*(1+((3+rho*file\$dist)*(12+12*rho*file\$dist+(rho*file\$dist)^2))/(n*(2+rho*file\$dist)*(11+rho*file\$dist)))

newfile <- data.frame(file\$dist, newrsq)
newfile <- newfile[order(newfile\$file.dist),]
newfile$PopName  <- base::unique(newfile)
lines(newfile$PopName\$dist/1000, newfile$PopName\$rsq, type="l",col="$ColArry[$IDE]",main="LD decay", xlab="Distance (Kb)", ylab=expression(r^{2}),bty="n",xlim=c(0,$maxX))

LOVE


			$legendCol.=",\"$ColArry[$IDE]\"";
			$legendName.=",\"$PopName\"";
			$legendlty.=",1";
		}





		my $Plot2=<<LOVE;
n=$SamplesSize;
read.table("$output.$PopName.rdtmp")->data;
sumCount<-sum(data[,5])
Row<-nrow(data);
dist  <- vector(length = sumCount)
rsq  <- vector(length = sumCount)
count=1;
for(i in 1:Row)
{
  for (j in 1:data[i,5])
  {
	 dist[count]<-data[i,1]
	 rsq[count]<-data[i,4]
	 count <- count + 1
  }
}
file <- data.frame(dist,rsq);
Cstart <- c(C=0.1)
modelC <- nls(rsq ~ ((10+C*dist)/((2+C*dist)*(11+C*dist)))*(1+((3+C*dist)*(12+12*C*dist+(C*dist)^2))/(n*(2+C*dist)*(11+C*dist))),data=file, start=Cstart, control=nls.control(maxiter=100))
rho <- summary(modelC)\$parameters[1]
newrsq <- ((10+rho*file\$dist)/((2+rho*file\$dist)*(11+rho*file\$dist)))*(1+((3+rho*file\$dist)*(12+12*rho*file\$dist+(rho*file\$dist)^2))/(n*(2+rho*file\$dist)*(11+rho*file\$dist)))

newfile <- data.frame(file\$dist, newrsq)
newfile <- newfile[order(newfile\$file.dist),]
newfile$PopName  <- base::unique(newfile)
LOVE



		if ($RowPlot==3)
		{
			$Plot2.=<<LOVE;
	plot(newfile$PopName\$dist/1000, newfile$PopName\$rsq, type="l",col="$ColArry[0]",main="LD decay", xlab="Distance (Kb)", ylab="D'",bty="n",xlim=c(0,$maxX))
LOVE
		}
		elsif ($RowPlot==4)
		{
			$Plot2.=<<LOVE;
	plot(newfile$PopName\$dist/1000, newfile$PopName\$rsq, type="l",col="$ColArry[0]", xlab="", ylab="",bty="n",xlim=c(0,$maxX),yaxt="n",lty=2)
	axis(4,col="red",col.ticks="red",col.axis="red",lty=2)
	mtext("D'",side=4,line=-1,col="red") 
LOVE
		}




		for (my $IDE=1; $IDE<= $PopNumber; $IDE++)
		{

			$PopName=$PopIDAryy[$IDE];
			$Plot2.=<<LOVE;
read.table("$output.$PopName.rdtmp")->data;
sumCount<-sum(data[,5])
Row<-nrow(data);
dist  <- vector(length = sumCount)
rsq  <- vector(length = sumCount)
count=1;
for(i in 1:Row)
{
  for (j in 1:data[i,5])
  {
	 dist[count]<-data[i,1]
	 rsq[count]<-data[i,4]
	 count <- count + 1
  }
}
file <- data.frame(dist,rsq);
Cstart <- c(C=0.1)
modelC <- nls(rsq ~ ((10+C*dist)/((2+C*dist)*(11+C*dist)))*(1+((3+C*dist)*(12+12*C*dist+(C*dist)^2))/(n*(2+C*dist)*(11+C*dist))),data=file, start=Cstart, control=nls.control(maxiter=100))
rho <- summary(modelC)\$parameters[1]
newrsq <- ((10+rho*file\$dist)/((2+rho*file\$dist)*(11+rho*file\$dist)))*(1+((3+rho*file\$dist)*(12+12*rho*file\$dist+(rho*file\$dist)^2))/(n*(2+rho*file\$dist)*(11+rho*file\$dist)))

newfile <- data.frame(file\$dist, newrsq)
newfile <- newfile[order(newfile\$file.dist),]
newfile$PopName  <- base::unique(newfile)
lines(newfile$PopName\$dist/1000, newfile$PopName\$rsq, type="l",col="$ColArry[$IDE]",main="LD decay", xlab="Distance (Kb)", ylab=expression(r^{2}),bty="n",xlim=c(0,$maxX))

LOVE

		}


		if ($RowPlot==2)
		{
			$Rshell=<<HEWM;
		pdf("$output.pdf")
		$Plot
		legend("topright",c($legendName),col=c($legendCol),cex=1,lty=c($legendlty),bty="n");
		dev.off()
HEWM
		}
		elsif ($RowPlot==3)
		{
			$Rshell=<<HEWM;
		pdf("$output.pdf")
		$Plot2
		legend("topright",c($legendName),col=c($legendCol),cex=1,lty=c($legendlty),bty="n");
		dev.off()
HEWM
		}
		elsif ($RowPlot==4)
		{
			$Plot2=~s/col=/lty=2,col=/g;
			$Rshell=<<HEWM;
		pdf("$output.pdf")
		$Plot
		legend("topright",c($legendName),col=c($legendCol),cex=1,lty=c($legendlty),bty="n");
		par(new=T)
		$Plot2
		dev.off()
HEWM
		}

	}
}
else
{
	print "\t\tPara [-method] wrong, Please set it one of [MeanBin/HW/MedianBin/PercentileBin]\n";
	exit(1);
}




my $R="/usr/bin/Rscript";
if  ( !(-e $R) )
{
	$R="/ifshk4/BC_PUB/biosoft/newblc/03.Soft_ALL/R-3.4.1/bin/Rscript";
}


if  ( !(-e $R) )
{
	$R=`which Rscript`;chomp $R;
}


open TMP,">$output.r" || die "output file can't open $!" ;

print TMP $Rshell ;
close TMP;


if  ( !(-e $R) )
{
	print "Can't find the [ Rscript ] bin, You shoud install the R First,then:\t\t";
	print " Rscript  $output.r  \n";
	exit(1);
}


system (" $R  $output.r  ");

if  (  defined($keepR) )
{
	system ("echo  $R  $output.r  ");
}
else
{
	system ("rm -rf   $output.r    $output.*tmp  ");
}




######################swiming in the sky and flying in the sea #############################




sub ReadFile2HashV1
{
	my $inFile=shift ;
	my $maxX_tmp=shift;
	my $SamplesSize=shift;
	my $FileType=shift;
	my $hash_D_sum=shift ;
	my $hash_RR_sum=shift ;
	my $hash_TT_sum=shift ;



	if  ($inFile =~s/\.gz$/\.gz/)
	{
		open AA,"gzip -cd  $inFile | "  || die "input file can't open $!" ;
	}
	else
	{
		open AA,"$inFile"  || die "input file can't open $!" ;
	}

	$_=<AA> ; chomp $_ ;
	if  ($_=~s/#Dist/#Dist/)
	{
		my @infSplit=split ;
		if ($#infSplit==5 )  { ${$FileType}=0;  print "Only Methold [MeanBin] can support this File format\n#Dist\tMean_r^2\tMean_D'\tSum_r^2\tSum_D'\tNumberPairs\nWrong File format at $inFile\nyou may try the [MeanBin] plot,or try core program [PopLDdecay] with parameter [ -OutType 4/5]\n"; exit(1);}		
		elsif ($#infSplit==2){ ${$FileType}=1;}
		elsif ($#infSplit==4){ ${$FileType}=2;}

		@infSplit=();
		$_=<AA>;  chomp ;  @infSplit=split ;
		if  (${$SamplesSize}==0) {${$SamplesSize}=$infSplit[-1];}
		elsif ( (${$SamplesSize})!= $infSplit[-1]  )
		{
			print "sample Num is diff with other Files at $inFile\n"; exit(1);
		}

	}
	else
	{
		print "Input File $inFile Format wrong,Can't found [#Dist] at first line\n";
		exit(1);
	}


	if  (${$FileType} == 1 )
	{
		while (<AA>)
		{
			chomp  ;
			my @inf=split ;
			next if ( $_=~s/#/#/);
			if  ($inf[0]>${$maxX_tmp})
			{
				${$maxX_tmp}=$inf[0];
			}
			($hash_RR_sum->{$inf[0]}{$inf[1]})+=($inf[2]);
			($hash_TT_sum->{$inf[0]}{$inf[1]})+=($inf[2]);
		}
	}
	elsif  (${$FileType} == 2  )
	{
		while (<AA>)
		{
			chomp  ;
			my @inf=split ;
			next if ( $_=~s/#/#/);
			if  ($inf[0]>${$maxX_tmp})
			{
				${$maxX_tmp}=$inf[0];
			}
			($hash_RR_sum->{$inf[0]}{$inf[1]})+=($inf[2]);
			($hash_TT_sum->{$inf[0]}{$inf[1]})+=($inf[2]);
			($hash_D_sum->{$inf[0]}{$inf[3]})+=($inf[4]);
			($hash_TT_sum->{$inf[0]}{$inf[3]})+=($inf[4]);
		}
	}
	close AA;
}





sub ReadFile2HashV2
{

	my $inFile=shift ;
	my $Small_bin=shift; 
	my $bin=shift;
	my $maxX_tmp=shift;
	my $hash_cout=shift ;
	my $hash_coutV2=shift ;
	my $hash_D_sum=shift ;
	my $hash_RR_sum=shift ;
	my $Small_cout=shift ;
	my $Small_coutV2=shift ;
	my $Small_D_sum=shift ;
	my $Small_RR_sum=shift ;



	my $FileType=0;

	if  ($inFile =~s/\.gz$/\.gz/)
	{
		open AA,"gzip -cd  $inFile | "  || die "input file can't open $!" ;
	}
	else
	{
		open AA,"$inFile"  || die "input file can't open $!" ;
	}

	$_=<AA> ; chomp $_ ;
	if  ($_=~s/#Dist/#Dist/)
	{
		my @infSplit=split ;
		if ($#infSplit==5){ $FileType=0;}
		elsif ($#infSplit==2){ $FileType=1;}
		elsif ($#infSplit==4){ $FileType=2;}
	}
	else
	{
		print "Input File $inFile Format wrong,Can't found [#Dist] at first line\n";
		exit(1);
	}

	if  ($FileType == 0 )
	{
		while (<AA>)
		{
			chomp  ;
			my @inf=split ;
			next if ($_=~s/#/#/);
			if  ($inf[0]>${$maxX_tmp})
			{
				${$maxX_tmp}=$inf[0];
			}
			if ($inf[0]>=$break)
			{
				$inf[0]=int($inf[0]/$bin);
				($hash_cout->{$inf[0]})+=$inf[-1] ;
				($hash_coutV2->{$inf[0]})+=$inf[-1] if  ($inf[-2] ne  "NA" );
				($hash_D_sum->{$inf[0]})+=$inf[-2]  if  ($inf[-2] ne  "NA" );
				($hash_RR_sum->{$inf[0]})+=$inf[-3];
			}
			else
			{
				$inf[0]=int($inf[0]/$Small_bin);
				($Small_cout->{$inf[0]})+=$inf[-1] ;
				($Small_coutV2->{$inf[0]})+=$inf[-1]  if  ($inf[-2] ne  "NA" );
				($Small_D_sum->{$inf[0]})+=$inf[-2]   if  ($inf[-2] ne  "NA" );
				($Small_RR_sum->{$inf[0]})+=$inf[-3];
			}
		}
	}
	elsif  ($FileType == 1 )
	{
		while (<AA>)
		{
			chomp  ;
			my @inf=split ;
			next if ( $_=~s/#/#/);
			if  ($inf[0]>${$maxX_tmp})
			{
				${$maxX_tmp}=$inf[0];
			}
			if ($inf[0]>=$break)
			{
				$inf[0]=int($inf[0]/$bin);
				($hash_cout->{$inf[0]})+=$inf[2] ;
				($hash_RR_sum->{$inf[0]})+=($inf[1]*$inf[2]);
			}
			else
			{
				$inf[0]=int($inf[0]/$Small_bin);
				($Small_cout->{$inf[0]})+=$inf[2] ;
				($Small_RR_sum->{$inf[0]})+=($inf[1]*$inf[2]);
			}
		}
	}
	elsif  ($FileType == 2  )
	{
		while (<AA>)
		{
			chomp  ;
			my @inf=split ;
			next if ($_=~s/#/#/);
			if  ($inf[0]>${$maxX_tmp})
			{
				${$maxX_tmp}=$inf[0];
			}
			if ($inf[0]>=$break)
			{
				$inf[0]=int($inf[0]/$bin);
				($hash_cout->{$inf[0]})+=$inf[2] ;
				($hash_RR_sum->{$inf[0]})+=($inf[1]*$inf[2]);
				($hash_coutV2->{$inf[0]})+=$inf[4];
				($hash_D_sum->{$inf[0]})+=($inf[4]*$inf[3]);
			}
			else
			{
				$inf[0]=int($inf[0]/$Small_bin);
				($Small_cout->{$inf[0]})+=$inf[2] ;
				($Small_RR_sum->{$inf[0]})+=($inf[1]*$inf[2]);				
				($Small_coutV2->{$inf[0]})+=$inf[4] ;
				($Small_D_sum->{$inf[0]})+=($inf[4]*$inf[3]);
			}
		}
	}
	close AA ;
}


######################swiming in the sky and flying in the sea #############################



