#!/usr/bin/perl -w
use strict;

use Data::Dumper;
use Getopt::Long;

#############Befor  Start  , open the files ####################

sub usage
{
	print STDERR <<USAGE;
	2016-04-22       hewm\@genomics.cn

	Usage:    perl $0  -inFile  LDdecay.stat.gz  -output OUT

	Options
	-inFile    <s> :  Input PopLDDecay OutPut Stat File
	-inList    <s> :  Input FileList if multi-File of PopLDDecay OutPut Stat
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

my ($help,$inFile,$inList,$output,$bin1,$bin2,$break,$keepR,$maxX,$measure,$method,$percent);


GetOptions(
	"help"=>\$help,
	"inFile:s"=>\$inFile,
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

if( defined($help) || !defined($output))
{
	usage;
	exit ;
}
if( (!defined($inFile)) && (!defined($inList)))
{
	usage;
	exit ;
}

$bin1||=10;
$bin2||=100;
$method||="MeanBin";
$break||=100;  $percent||=0.5;
$measure||="r2";    my $RowPlot=2;  
if ($measure eq  "D"  ||  $measure eq  "D'" )  { $RowPlot=3 ;} 
elsif ($measure eq  "b"  ||  $measure eq  "both"  ||   $measure eq  "Both"  ||   $measure eq  "B")  { $RowPlot=4 ;}

my $Rshell="";
my $maxX_tmp=0;

if ( $method eq "PercentileBin" &&  $percent==0.5) {$method="MedianBin" ;}

# MeanBin/HW/MedianBin/PercentileBin
if ( $method eq "HW" || $method eq "MedianBin" || $method eq "PercentileBin" || $method eq "Percentile" ) # HW =="Hill and Weir"
{
	my %hashRR=(); 
	my %hashDD=();
	my %hashTT=();
	my $SamplesSize=0;
	my $FileType=0;
	if ( defined($inFile) )
	{
		&ReadFile2HashV1( $inFile,\$maxX_tmp,\$SamplesSize,\$FileType,\%hashDD,\%hashRR,\%hashTT);
	}

	if ( defined($inList) )
	{
		if  ($inList =~s/\.gz$/\.gz/)
		{
			open LIST,"gzip -cd  $inList | "  || die "input file can't open $!" ;
		}
		else
		{
			open LIST,"$inList"  || die "input file can't open $!" ;
		}
		while($_=<LIST>)
		{
			my $FileThis=$_; chomp $FileThis;
			&ReadFile2HashV1($FileThis,\$maxX_tmp,\$SamplesSize,\$FileType,\%hashDD,\%hashRR);
		}
		close LIST ;
	}

	if ( $RowPlot  >2 && $FileType< 2 )
	{
		print "The input result only with measure [r2], can't find the D', To calculate the D' result , you should run the core program [PopLDdecay] with parameter [ -OutType 4]\n";
		print "To run the measure r2 LD decay, please set the parameter [ -measure r2 ] \n";
		exit(1);
	}

	open OA,">$output.bin" || die "output file can't open $!" ;

	print OA"#Dist\tR^2\tR^2_count\n";
	my $ALL_Pairwise=0;
	foreach my $k ( sort {$a<=>$b} keys %hashRR)
	{
		my $cc=$hashRR{$k};
		foreach my $jj (sort {$a<=>$b} keys %$cc)
		{
			print OA $k,"\t$jj\t$hashRR{$k}{$jj}\n";
			$ALL_Pairwise+=$hashRR{$k}{$jj};
		}
	}

	close OA;

	if ( $FileType == 2  && $RowPlot  >2 )
	{
		open  OD,">$output.bin" || die "output file can't open $!" ;
		print OD "#Dist\tR^2\tR^2_count\tD'\tD_count\n";
		foreach my $k ( sort {$a<=>$b} keys %hashTT)
		{
			my $cc=$hashTT{$k};
			foreach my $jj ( sort {$a<=>$b} keys  %$cc )
			{
				my $DD=$hashDD{$k}{$jj} ;  $DD||=0;
				my $RR=$hashRR{$k}{$jj} ;  $RR||=0;
				print OA $k,"\t$jj\t$RR\t$jj\t$DD\n";
			}
		}
		close OD;
	}

	%hashTT=();  undef %hashTT;

	$maxX_tmp=int($maxX_tmp/1000);
	$maxX||=$maxX_tmp;



	if  ( $method eq "HW" )
	{
		if  ( $ALL_Pairwise>100000 || $ALL_Pairwise< 0 )
		{
			print "\t\tHW methold is only suitable for small data, we recommend the total Pairwise Number <100000\n";
			print "\t\tHowever Now all Pairwise Number is $ALL_Pairwise\n; we suggest you may use the [MeanBin] to plot\n";
			exit (1);
		}
		if ( $RowPlot  == 2 )
		{
			$Rshell=<<LOVE;
# see web : https://jujumaan.com/2017/07/15/linkage-disequilibrium-decay-plot
read.table("$output.bin")->data
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
n=$SamplesSize ;
Cstart <- c(C=0.1)
modelC <- nls(rsq ~ ((10+C*dist)/((2+C*dist)*(11+C*dist)))*(1+((3+C*dist)*(12+12*C*dist+(C*dist)^2))/(n*(2+C*dist)*(11+C*dist))),data=file, start=Cstart, control=nls.control(maxiter=100))
rho <- summary(modelC)\$parameters[1]
newrsq <- ((10+rho*file\$dist)/((2+rho*file\$dist)*(11+rho*file\$dist)))*(1+((3+rho*file\$dist)*(12+12*rho*file\$dist+(rho*file\$dist)^2))/(n*(2+rho*file\$dist)*(11+rho*file\$dist)))

newfile <- data.frame(file\$dist, newrsq)
newfile <- newfile[order(newfile\$file.dist),]
newfile  <- base::unique(newfile)

pdf("$output.pdf")
plot(newfile\$dist/1000, newfile\$rsq, type="l",col="blue",main="LD decay", xlab="Distance (Kb)", ylab=expression(r^{2}),bty="n",xlim=c(0,$maxX))
dev.off()
png("$output.png")
plot(newfile\$dist/1000, newfile\$rsq, type="l",col="blue",main="LD decay", xlab="Distance (Kb)", ylab=expression(r^{2}),bty="n",xlim=c(0,$maxX))
dev.off()

LOVE

		}
		elsif( $RowPlot ==3)
		{


			$Rshell=<<LOVE;
read.table("$output.bin")->data

sumCount<-sum(data[,5])
Row<-nrow(data);
dist  <- vector(length = sumCount)
D  <- vector(length = sumCount)
count=1;
for(i in 1:Row)
{
  for (j in 1:data[i,5])
  {
	 dist[count]<-data[i,1]
	 D[count]<-data[i,4]
	 count <- count + 1
  }
}


file <- data.frame(dist,D);
n=$SamplesSize ;
Cstart <- c(C=0.1)
modelC <- nls(D ~ ((10+C*dist)/((2+C*dist)*(11+C*dist)))*(1+((3+C*dist)*(12+12*C*dist+(C*dist)^2))/(n*(2+C*dist)*(11+C*dist))),data=file, start=Cstart, control=nls.control(maxiter=100))
rho <- summary(modelC)\$parameters[1]
newD <- ((10+rho*file\$dist)/((2+rho*file\$dist)*(11+rho*file\$dist)))*(1+((3+rho*file\$dist)*(12+12*rho*file\$dist+(rho*file\$dist)^2))/(n*(2+rho*file\$dist)*(11+rho*file\$dist)))

newfile <- data.frame(file\$dist, newD)
newfile <- newfile[order(newfile\$file.dist),]
newfile  <- base::unique(newfile)

pdf("$output.pdf")
plot(newfile\$dist/1000, newfile\$D, type="l",col="blue",main="LD decay", xlab="Distance (Kb)", ylab="D'",bty="n",xlim=c(0,$maxX))
dev.off()
png("$output.png")
plot(newfile\$dist/1000, newfile\$D, type="l",col="blue",main="LD decay", xlab="Distance (Kb)", ylab="D'",bty="n",xlim=c(0,$maxX))
dev.off()

LOVE


		}
		elsif( $RowPlot ==4)
		{

			$Rshell=<<LOVE;
read.table("$output.bin")->data
Row<-nrow(data);

sumCountD<-sum(data[,5])
distD  <- vector(length = sumCountD)
DD  <- vector(length = sumCountD)
countD=1;

sumCountR<-sum(data[,3])
dist  <- vector(length = sumCountR)
rsq  <- vector(length = sumCountR)
countR=1;


for(i in 1:Row)
{
  for (j in 1:data[i,5])
  {
	 distD[countD]<-data[i,1]
	 DD[countD]<-data[i,4]
	 countD <- countD + 1
  }

  for (j in 1:data[i,3])
  {
	 dist[countR]<-data[i,1]
	 rsq[countR]<-data[i,2]
	 countR <- countR + 1
  }

}




n=$SamplesSize ;
file <- data.frame(dist,rsq);
Cstart <- c(C=0.1)
modelC <- nls(rsq ~ ((10+C*dist)/((2+C*dist)*(11+C*dist)))*(1+((3+C*dist)*(12+12*C*dist+(C*dist)^2))/(n*(2+C*dist)*(11+C*dist))),data=file, start=Cstart, control=nls.control(maxiter=100))
rho <- summary(modelC)\$parameters[1]
newrsq <- ((10+rho*file\$dist)/((2+rho*file\$dist)*(11+rho*file\$dist)))*(1+((3+rho*file\$dist)*(12+12*rho*file\$dist+(rho*file\$dist)^2))/(n*(2+rho*file\$dist)*(11+rho*file\$dist)))

newfile <- data.frame(file\$dist, newrsq)
newfile <- newfile[order(newfile\$file.dist),]
newfile  <- base::unique(newfile)


fileEE <- data.frame(distD,DD);
Cstart <- c(C=0.1)
modelC <- nls(DD ~ ((10+C*distD)/((2+C*distD)*(11+C*distD)))*(1+((3+C*distD)*(12+12*C*distD+(C*distD)^2))/(n*(2+C*distD)*(11+C*distD))),data=fileEE, start=Cstart, control=nls.control(maxiter=100))
rhoDD <- summary(modelC)\$parameters[1]
newDD <- ((10+rhoDD*fileEE\$distD)/((2+rhoDD*fileEE\$distD)*(11+rhoDD*fileEE\$distD)))*(1+((3+rhoDD*fileEE\$distD)*(12+12*rhoDD*fileEE\$distD+(rhoDD*fileEE\$distD)^2))/(n*(2+rhoDD*fileEE\$distD)*(11+rhoDD*fileEE\$distD)))

newfileEE <- data.frame(fileEE\$distD, newDD)
newfileEE <- newfileEE[order(newfileEE\$file.distD),]
newfileEE  <- base::unique(newfileEE)


pdf("$output.pdf")
plot(newfile\$dist/1000, newfile\$rsq, type="l",col="blue",main="LD decay", xlab="Distance (Kb)", ylab=expression(r^{2}),bty="n",ylim=c(0,1),xlim=c(0,$maxX))
par(new=T)
plot(newfileEE\$distD/1000, newfileEE\$DD, type="l",col="red",xlab="", ylab="",bty="n",xlim=c(0,$maxX),yaxt="n")
axis(4,col="red",col.ticks="red",col.axis="red")
mtext("D'",side=4,line=-1,col="red")
dev.off()

png("$output.png")
plot(newfile\$dist/1000, newfile\$rsq, type="l",col="blue",main="LD decay", xlab="Distance (Kb)", ylab=expression(r^{2}),bty="n",ylim=c(0,1),xlim=c(0,$maxX))
par(new=T)
plot(newfileEE\$distD/1000, newfileEE\$DD, type="l",col="red",xlab="", ylab="",bty="n",xlim=c(0,$maxX),yaxt="n")
axis(4,col="red",col.ticks="red",col.axis="red")
mtext("D'",side=4,line=-1,col="red")

dev.off()

LOVE

		}

	}

	##################

	else 
	{
		#	my %newRRL=();
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
		 %hashRR=();		undef  %hashRR ; 
		open ORR,">$output.rrtmp" || die "output file can't open $!";
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

		 %hashDD=();		undef %hashDD;
		open ODD,">$output.ddtmp" || die "output file can't open $!";
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

		%hashRR=();undef %hashRR;

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

		%hashTT=();undef %hashTT;


		$percent=$percent*100;
		$percent="$percent%";
		$method=~s/Bin$//g;

		if ( $RowPlot  == 2 )
		{

			$Rshell=<<LOVE;
read.table("$output.rrtmp")->data
pdf("$output.pdf")
plot(data[,1]/1000,data[,2],type="l",col="blue",main="LD decay $method ($percent)",xlab="Distance(Kb)",xlim=c(0,$maxX),ylab=expression(r^{2}),bty="n")
dev.off()
png("$output.png")
plot(data[,1]/1000,data[,2],type="l",col="blue",main="LD decay $method ($percent)",xlab="Distance(Kb)",xlim=c(0,$maxX),ylab=expression(r^{2}),bty="n")
dev.off()

LOVE
		}
		elsif ( $RowPlot==3 )
		{
			$Rshell=<<LOVE;
read.table("$output.ddtmp")->data
pdf("$output.pdf")
plot(data[,1]/1000,data[,2],type="l",col="blue",main="LD decay $method ($percent)",xlab="Distance(Kb)",xlim=c(0,$maxX),ylab="D'",bty="n")
dev.off()
png("$output.png")
plot(data[,1]/1000,data[,2],type="l",col="blue",main="LD decay $method ($percent)",xlab="Distance(Kb)",xlim=c(0,$maxX),ylab="D'",bty="n")
dev.off()

LOVE
		}
		elsif (  $RowPlot==4 )
		{
			$Rshell=<<LOVE;
read.table("$output.rrtmp")->data
read.table("$output.ddtmp")->dataDD
pdf("$output.pdf")
plot(data[,1]/1000,data[,2],type="l",col="blue",main="LD decay $method ($percent)",xlab="Distance(Kb)",xlim=c(0,$maxX),ylab=expression(r^{2}),bty="n")
par(new=T)
plot(dataDD[,1]/1000,dataDD[,2],type="l",col="red",xlab="",ylab="",xlim=c(0,$maxX),bty="n",yaxt="n")
axis(4,col="red",col.ticks="red",col.axis="red")
mtext("D'",side=4,line=-1,col="red")
dev.off()

png("$output.png")
plot(data[,1]/1000,data[,2],type="l",col="blue",main="LD decay $method ($percent)",xlab="Distance(Kb)",xlim=c(0,$maxX),ylab=expression(r^{2}),bty="n")
par(new=T)
plot(dataDD[,1]/1000,dataDD[,2],type="l",col="red",xlab="",ylab="",xlim=c(0,$maxX),bty="n",yaxt="n")
axis(4,col="red",col.ticks="red",col.axis="red")
mtext("D'",side=4,line=-1,col="red")
dev.off()


LOVE
		}




	}

}
elsif ($method eq "MeanBin" ) 
{

	open OA,">$output.bin" || die "output file can't open $!" ;

	my %hash_cout=();
	my %hash_coutV2=();
	my %hash_RR_sum=();
	my %hash_D_sum=();

	my %Small_cout=();
	my %Small_coutV2=();
	my %Small_RR_sum=();
	my %Small_D_sum=();


	if ( defined($inFile) )
	{
		&ReadFile2HashV2($inFile,$bin1,$bin2,\$maxX_tmp ,\%hash_cout, \%hash_coutV2, \%hash_D_sum,\%hash_RR_sum, \%Small_cout,\%Small_coutV2,\%Small_D_sum, \%Small_RR_sum);
	}

	if ( defined($inList) )
	{

		if  ($inList =~s/\.gz$/\.gz/)
		{
			open LIST,"gzip -cd  $inList | "  || die "input file can't open $!" ;
		}
		else
		{
			open LIST,"$inList"  || die "input file can't open $!" ;
		}

		while($_=<LIST>)
		{
			my $FileThis=$_; chomp $FileThis;
			&ReadFile2HashV2($FileThis,$bin1,$bin2,\$maxX_tmp ,\%hash_cout, \%hash_coutV2, \%hash_D_sum,\%hash_RR_sum, \%Small_cout,\%Small_coutV2,\%Small_D_sum, \%Small_RR_sum);
		}
		close LIST ;
	}

	$maxX_tmp=int($maxX_tmp/1000);

	$maxX||=$maxX_tmp;

	print OA "#Dist\tMean_r^2\tMean_D'\tSum_r^2\tSum_D'\tNumberPairs\n";

	my $Cal_D=0;

	foreach my $k (sort {$a<=>$b} keys %Small_cout)
	{
		my $mean_R=$Small_RR_sum{$k}/$Small_cout{$k};

		if  (exists $Small_coutV2{$k})
		{
			my $mean_D=$Small_D_sum{$k}/$Small_coutV2{$k};
			print OA ($k+1)*$bin1,"\t$mean_R\t$mean_D\t$Small_RR_sum{$k}\t$Small_D_sum{$k}\t$Small_cout{$k}\n";
			$Cal_D=1;
		}
		else
		{
			print OA ($k+1)*$bin1,"\t$mean_R\tNA\t$Small_RR_sum{$k}\tNA\t$Small_cout{$k}\n";
		}
	}

	foreach my $k (sort {$a<=>$b} keys %hash_cout)
	{
		my $mean_R=$hash_RR_sum{$k}/$hash_cout{$k};
		if  (exists $hash_coutV2{$k})
		{
			my $mean_D=$hash_D_sum{$k}/$hash_coutV2{$k};
			print OA ($k+1)*$bin2 ,"\t$mean_R\t$mean_D\t$hash_RR_sum{$k}\t$hash_D_sum{$k}\t$hash_cout{$k}\n";
			$Cal_D=1;
		}
		else
		{
			print OA ($k+1)*$bin2 ,"\t$mean_R\tNA\t$hash_RR_sum{$k}\tNA\t$hash_cout{$k}\n";
		}
	}
	close OA ;

	if  ( $RowPlot>2  &&  $Cal_D==0 )
	{
		print "The input result only with measure [r2], can't find the D', To calculate the D' result , you should run the core program [PopLDdecay] with parameter [ -OutType 2]\n";
		print "To run the measure r2 LD decay, please set the parameter [ -measure r2 ] \n";
		exit ;
	}


	$Rshell=<<LOVE;
read.table("$output.bin")->data
pdf("$output.pdf")
plot(data[,1]/1000,data[,2],type="l",col="blue",main="LD decay",xlab="Distance(Kb)",xlim=c(0,$maxX),ylab=expression(r^{2}),bty="n")
dev.off()
png("$output.png")
plot(data[,1]/1000,data[,2],type="l",col="blue",main="LD decay",xlab="Distance(Kb)",xlim=c(0,$maxX),ylab=expression(r^{2}),bty="n")
dev.off()

LOVE

	if ( $RowPlot==3 )
	{
		$Rshell=<<LOVE;
read.table("$output.bin")->data
pdf("$output.pdf")
plot(data[,1]/1000,data[,3],type="l",col="blue",main="LD decay",xlab="Distance(Kb)",xlim=c(0,$maxX),ylab="D'",bty="n")
dev.off()
png("$output.png")
plot(data[,1]/1000,data[,3],type="l",col="blue",main="LD decay",xlab="Distance(Kb)",xlim=c(0,$maxX),ylab="D'",bty="n")
dev.off()

LOVE

	}
	elsif (  $RowPlot==4 )
	{
		$Rshell=<<LOVE;
read.table("$output.bin")->data
pdf("$output.pdf")
plot(data[,1]/1000,data[,2],type="l",col="blue",main="LD decay",xlab="Distance(Kb)",xlim=c(0,$maxX),ylab=expression(r^{2}),bty="n")
par(new=T)
plot(data[,1]/1000,data[,3],type="l",col="red",xlab="",ylab="",xlim=c(0,$maxX),bty="n",yaxt="n")
axis(4,col="red",col.ticks="red",col.axis="red")
mtext("D'",side=4,line=-1,col="red")
dev.off()

png("$output.png")
plot(data[,1]/1000,data[,2],type="l",col="blue",main="LD decay",xlab="Distance(Kb)",xlim=c(0,$maxX),ylab=expression(r^{2}),bty="n")
par(new=T)
plot(data[,1]/1000,data[,3],type="l",col="red",xlab="",ylab="",xlim=c(0,$maxX),bty="n",yaxt="n")
axis(4,col="red",col.ticks="red",col.axis="red")
mtext("D'",side=4,line=-1,col="red")
dev.off()

LOVE
	}


}
else
{
	print "\t\tPara [-method] wrong, Please set it one of [MeanBin/HW/MedianBin/PercentileBin]\n";
	exit(1);
}









#############################   R  plot ###############################

open TMP,">$output.r" || die "output file can't open $!" ;
print TMP $Rshell ;
close TMP;

my $R="/usr/bin/Rscript";

if  ( !(-e $R) )
{
	$R=`which Rscript`;chomp $R;
}



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
	system ("rm -rf   $output.r     $output.*tmp  ");
	if  ( -e "$output.bin" )
	{
	system ("gzip  $output.bin ");
	}
}

################## sub funtion #########


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
		if ($#infSplit==5 )  { ${$FileType}=0;  print "Only Methold [MeanBin] can support this File format\n#Dist\tMean_r^2\tMean_D'\tSum_r^2\tSum_D'\tNumberPairs\nWrong File format at $inFile\n,you may try the [MeanBin] plot,or try core program [PopLDdecay] with parameter [ -OutType 4/5]\n"; exit(1);}
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
######################swiming in the sky and flying in the sea #############################


