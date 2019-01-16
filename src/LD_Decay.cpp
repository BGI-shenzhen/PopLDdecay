#ifndef LDDecay_H_
#define LDDecay_H_

#include "HeadIN.h"
#include "FileDeal.h"
#include "Calculate.h"
#include "ProMethod1.h"
#include "ProMethod2.h"
#include "PariWiseCal.h"
#include "EHHDeal.h"

using namespace std;

void  LDdecaySNP_help()
{
	cout <<""
		"\n"
		"\tUsage: PopLDDecay -InVCF  <in.vcf.gz>  -OutStat <out.stat>\n"
		"\n"
		"\t\t-InVCF         <str>     Input SNP VCF Format\n"
		"\t\t-InGenotype    <str>     Input SNP Genotype Format\n"
		"\t\t-OutStat       <str>     OutPut Stat Dist ~ r^2/D' File\n"
		"\n"
		"\t\t-SubPop        <str>     SubGroup Sample File List[ALLsample]\n"
		"\t\t-MaxDist       <int>     Max Distance (kb) between two SNP [300]\n"
		"\t\t-MAF           <float>   Min minor allele frequency filter [0.005]\n"
		"\t\t-Het           <float>   Max ratio of het allele filter [0.88]\n"
		"\t\t-Miss          <float>   Max ratio of miss allele filter [0.25]\n"
		"\t\t-EHH           <str>     To Run EHH Region decay set StartSite [NA]\n"
		"\t\t-OutFilterSNP            OutPut the final SNP to calculate\n"
		"\t\t-OutType       <int>     1: R^2 result 2: R^2 & D' result 3:PairWise LD Out[1]\n"
		"\t\t                         See the Help for more OutType [1-8] details\n"
		//		"\t\t-Method        <int>     Select the Cal agorithm [1]\n"
		//		"\t\t                         1:Low MEM  2 May Big MEM\n"
		"\t\t\n"
		"\t\t-help                    Show more help [hewm2008 v3.40]\n"
		"\n";
}

void More_HelpLDDecay()
{
	cout<<""
		"\n"
		"\t\t More Help document please see the Manual.pdf file\n"
		" Para [-OutType] can be [1-8],(1-3) are recommended\n"
		"        [-OutType 1] is the fastest for only cal (Dist ~ R^2) for MeanBin method plot\n"
		"        [-OutType 2] will OutPut the Stat (Dist ~ r^2 & D') result for R^2 & D' MeanBin method plot\n"
		"        [-OutType 3] will OutPut one more result of PairWise LD compaire result(with Dist~r^2)\n"
		"        [-OutType 4] will OutPut the Stat (Dist ~ r^2 & D' ~ Number) result for R^2 & D' MeanBin/HW/MedianBin/PercentileBin plot\n"
		"        [-OutType 5] will OutPut the Stat (Dist ~ r^2 ~ Number) result for R^2 MeanBin/HW/MedianBin/PercentileBin plot\n"
		"        [-OutType 6] will OutPut one more result of PairWise LD compaire result(with Dist~r^2/D')\n"
		"        [-OutType 7] will OutPut one more result of PairWise LD compaire result(with Dist~r^2/D'/LOD)\n"
		"        [-OutType 8] will OutPut one more result of PairWise LD compaire result(with Dist~r^2/D'/LOD/CIlow/CIhi)\n"
		" Para [-Method] can be 1 or 2, default agorithm [1],at most time agorithm 1 is faster than agorithm 2, and agorithm 2 may will cost Big MEM\n"
		" Para [-EHH] format should be chr:site, such like -EHH chr1:5000000 will give out the EHH Decay of this site nearby distance\n"
        "\n"
		"\n"
		" Here, we provide four classic cases to demonstrate the application of this software, four situation will be show how to follow to get the LD decay figure out.\n"
		"\n"
		"1). One population\n"
		"    This situation (one population with all chromosomes together) is encountered by most users, and this situation is the simplest to carry out.\n"
		"\n"
		"              ./bin/PopLDdecay -InVCF  ALLchr.vcf.gz  -OutStat  LDDecay.stat.gz\n"
		"              perl bin/Plot_OnePop.pl  -inFile  LDDecay.stat.gz -output  Out.Prefix\n"
		"\n"
		"    Note:\n"
		"      This will generate the two finale figures named 'Out.Prefix.png' and 'Out.Prefix.pdf'\n"
		"\n"
		"2). Muti population\n"
		"    This is common situation in the LD decay analysis. For example, if there are 50 samples (wild1, wild2, wild3...wild25, cul1, cul2, cul3...cul25) in the VCF file,\n"
		"    To compare the LD decay of these two groups (wild vs cultivation), first of all, put their sample names into own file list for each group, column or row is ok.\n"
		"\n"
		"             ./bin/PopLDdecay -InVCF  In.vcf.gz  -OutStat  wild.stat.gz  -SubPop wildName.list\n"
		"             ./bin/PopLDdecay -InVCF  In.vcf.gz  -OutStat   cul.stat.gz  -SubPop culName.list\n"
		"             #   created manually  muti.list by yourself\n"
		"             perl bin/Plot_MutiPop.pl -inList  muti.list  -output  OutputPrefix\n"
		"\n"
		"    Note:\n"
		"      A. The <wildName.list> can list as follow(column or row is ok):\n"
		"                      wild1\n"
		"                      wild2\n"
		"                      ...\n"
		"                      wild25\n"
		"      B. The format of <muti.list> had two columns, the file path of population result and the population flag, such as:\n"
		"                      /ifshk7/BC_PS/Lddecay/wild.stat.gz   wild\n"
		"                      /ifshk7/BC_PS/Lddecay/cul.stat.gz    cultivation\n"
		"\n"
		"3). One population with multi-chr\n"
		"    One population with multiple chromosome VCF files. For example, if there are 3 chromosomes VCF files (Chr1, Chr2 and Chr3) as the input.\n"
		"\n"
		"            ./bin/PopLDdecay -InVCF  Chr1.vcf.gz  -OutStat  Chr1.stat.gz\n"
		"            ./bin/PopLDdecay -InVCF  Chr2.vcf.gz  -OutStat  Chr2.stat.gz\n"
		"            ./bin/PopLDdecay -InVCF  Chr3.vcf.gz  -OutStat  Chr3.stat.gz\n"
		"            ls  `pwd`/Chr*.stat.gz   > chr.list\n"
		"            perl bin/Plot_OnePop.pl -inList  chr.list  -output  OutputPrefix\n"
		"\n"
		"   Note:\n"
		"     A. It can run in parallel when calculating the chromosomes' statistics files.\n"
		"     B. The files list only store the file path, which is diff with the multi-population list\n"
		"     C. It will generate the file 'OutputPrefix.bin' is the summary statistics file of all chromosomes, and same format with the chromosomes' statistics files.\n"
		"     D. The <chr.list> format can be generated by as above command 'ls Chr*.stat.gz   > chr.list'\n"
		"\n"
		"4). Muti population with multi-chr\n"
		"    Muti population with multiple chromosome VCF files. For example, if there are 2 chromosomes VCF files (Chr1, Chr2) as the input.\n"
		"\n"
		"             ./bin/PopLDdecay -InVCF  Chr1.vcf.gz  -OutStat  W.Chr1.stat.gz -SubPop wildName.list\n"
		"             ./bin/PopLDdecay -InVCF  Chr2.vcf.gz  -OutStat  W.Chr2.stat.gz -SubPop wildName.list\n"
		"             ./bin/PopLDdecay -InVCF  Chr1.vcf.gz  -OutStat  C.Chr1.stat.gz -SubPop culName.list\n"
		"             ./bin/PopLDdecay -InVCF  Chr2.vcf.gz  -OutStat  C.Chr2.stat.gz -SubPop culName.list\n"
		"             ls  `pwd`/W.Chr*.stat.gz   > W.chr.list\n"
		"             perl bin/Plot_OnePop.pl -inList  W.chr.list  -output  Wild.cat\n"
		"             ls  `pwd`/C.Chr*.stat.gz   > C.chr.list\n"
		"             perl bin/Plot_OnePop.pl -inList  C.chr.list  -output  Cul.cat\n"
		"             perl bin/Plot_MutiPop.pl -inList  muti.list  -output  OutputPrefix\n"
		"\n"
		"    Note:\n"
		"     A. The format of <muti.list> had two columns , the file path of population result and the population flag, such as:\n"
		"                      /ifshk7/BC_PS/Lddecay/Wild.cat.bin    wild\n"
		"                      /ifshk7/BC_PS/Lddecay/Cul.cat.bin     cultivation\n"
		"\n";
}

int LDdecay_help01(int argc, char **argv , In3str1v * paraFA04, Para_18 * para_18)
{
	if (argc <2 ) {LDdecaySNP_help();return 0;}

	for(int i = 1; i < argc ; i++)
	{
		if(argv[i][0] != '-')
		{
			cerr << "command option error! please check." << endl;
			return 0;
		}
		string flag=argv[i] ;
		flag=replace_all(flag,"-","");

		if (flag  == "InVCF" ||  flag  == "i")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFA04->InStr1=argv[i];
		}
		else if (flag  == "InGenotype")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFA04->InStr3=argv[i];
		}
		else if (flag  == "SubPop"  ||  flag  == "s")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFA04->SubPop=argv[i];
		}
		else if (flag  ==  "OutStat" ||  flag  == "o")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFA04->InStr2=argv[i];
		}
		else if (flag  == "Het" )
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			para_18->Het=atof(argv[i]);
		}
		else if (flag == "MAF")
		{
			if(i + 1== argc) {LogLackArg(flag);return 0;}
			i++;
			para_18->MAF=atof(argv[i]);
		}
		else if (flag == "Miss")
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0;}
			i++;
			para_18->Miss=atof(argv[i]);
		}
		else if (flag == "MaxDist" )
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0;}
			i++;
			paraFA04->InInt=atoi(argv[i]);
		}
		else if (flag == "EHH" )
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0;}
			i++;
			paraFA04->EHH=argv[i];
			if ((paraFA04->EHH).find(':') ==string::npos)
			{
				cerr<<"\tPara [-EHH] should be [chr:Site],such [chr1:5000]"<<endl;
				return 0;
			}
		}
		else if (flag == "OutPairLD" || flag == "OutType")
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0;}
			i++;
			paraFA04->TF=atoi(argv[i]);
		}
		else if (flag == "Method")
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0;}
			i++;
			paraFA04->Method=atoi(argv[i]);
		}
		else if (flag == "OutFilterSNP")
		{
			paraFA04->TF2=false;
		}
		else if (flag == "help"  ||  flag == "h")
		{
			More_HelpLDDecay();return 0;
		}
		else
		{
			cerr << "UnKnow argument -"<<flag<<endl;
			return 0;
		}
	}

	if  ( (paraFA04->InStr2).empty() ||( (paraFA04->InStr1).empty() &&  (paraFA04->InStr3).empty()) )
	{
		cerr<< "lack argument for the must"<<endl ;
		return 0;
	}

	if ( (!(paraFA04->EHH).empty())  &&   ((para_18->MAF)<0.05) )
	{
		(para_18->MAF)=0.05;
		cerr<<"\t\t -MAF for -EHH should >=0.05 ; so -MAF  default 0.05 Now"<<endl;
	}

	if ( (paraFA04->TF) >8  || (paraFA04->TF) <0 )
	{
		cerr<<"\t\t-OutType   should be  [0-8]"<<endl;
		return 0;
	}


	string Stat=(paraFA04->InStr2);
	string ext =Stat.substr(Stat.rfind('.') ==string::npos ? Stat.length() : Stat.rfind('.') + 1);
	if (ext == "gz")
	{
		(paraFA04->InStr2)=(paraFA04->InStr2).substr(0,(paraFA04->InStr2).length()-3);
	}

	Stat=(paraFA04->InStr2);
	ext =Stat.substr(Stat.rfind('/') ==string::npos ? Stat.length() : Stat.rfind('/') + 1);

	if (ext != "stat")
	{
		ext =Stat.substr(Stat.rfind('.') ==string::npos ? Stat.length() : Stat.rfind('.') + 1);
		if (ext == "stat")
		{
			(paraFA04->InStr2)=(paraFA04->InStr2).substr(0,(paraFA04->InStr2).length()-5);
		}
	}


	return 1 ;
}


int main(int argc, char *argv[])
{
	In3str1v *paraFA04 = new In3str1v;
	Para_18 * para_18 = new Para_18 ;
	if ( (LDdecay_help01(argc, argv, paraFA04, para_18)==0))
	{
		delete paraFA04 ;
		delete para_18 ;
		return 1 ;
	}

	(paraFA04->InInt)=(paraFA04->InInt)*1000;

	//	char buf[1024*1024]; 	setbuf(stdout, buf);
	//*///////////////////////////Test  Out File is OK ///////////////////////////////*//
	string Stat=(paraFA04->InStr2);
	Stat=(paraFA04->InStr2)+".stat.gz";	
	ogzstream OUTTest ((Stat).c_str());
	if((!OUTTest.good()))
	{
		cerr << "open OUT File error: "<<Stat<<endl;
		delete para_18;
		delete  paraFA04 ; return 1;
	}
	OUTTest.close();


	map <string,map <llong, vector <BaseType> > >  SNPList ;
	int Flag_for_pro=0;

	/////*               VCF   IN Deal //////////////////////*////
	if (!(paraFA04->InStr1).empty())
	{
		bool PhaseVCF=DeleVCFPhase(paraFA04->InStr1);
		if  ( PhaseVCF )
		{
			if ((paraFA04->SubPop).empty())
			{
				Read_VCF_IN_Phase( paraFA04, para_18 , SNPList, Flag_for_pro);
			}
			else
			{
				Read_SubPopVCF_IN_Phase( paraFA04, para_18 , SNPList, Flag_for_pro);
			}
		}
		else
		{
			if ((paraFA04->SubPop).empty())
			{
				Read_VCF_IN( paraFA04, para_18 , SNPList, Flag_for_pro);
			}
			else
			{
				Read_SubPopVCF_IN( paraFA04, para_18 , SNPList, Flag_for_pro);
			}
		}
	}

	/////*               Genotype IN Deal //////////////////////*////
	if (!(paraFA04->InStr3).empty())
	{
		if  ((paraFA04->SubPop).empty())
		{
			Read_Genotype_IN(paraFA04, para_18 ,SNPList,Flag_for_pro);
		}
		else
		{
			Read_SubPopGenotype_IN(paraFA04, para_18 ,  SNPList , Flag_for_pro );
		}
	}



	if (!(paraFA04->EHH).empty())
	{
		EHH_Region_LDDecay(paraFA04 , para_18, SNPList, Flag_for_pro );
		delete paraFA04 ;
		delete para_18 ;
		return 0 ;		
	}



	//*///////////////////////////PairWise Compare//////////////////////////////////*//

	if ( ((paraFA04->TF)==4) ||  ((paraFA04->TF)==5)  )
	{
		PairWiseComNewOUT_A(paraFA04, para_18 , SNPList, Flag_for_pro );
	}
	else
	{
		StarRsult *All_Stat = new StarRsult [((paraFA04->InInt)+1)];

		if ((paraFA04->Method)==2)
		{
			PairWiseComV2(paraFA04, para_18 , SNPList ,  All_Stat , Flag_for_pro );
		}
		else
		{
			PairWiseComV1(paraFA04, para_18 , SNPList ,  All_Stat , Flag_for_pro );
		}
		////* OUT Stat  File /*/////
		OUTStatFile( paraFA04, para_18 , All_Stat );
		delete [] All_Stat ;
	}

	cout<<"Used [perl  ../bin/Plot_XX.pl ] to Plot the LDdecay"<<endl;

	delete para_18 ;
	delete paraFA04 ;
	return 0;

}


#endif // LDDecay_H_ //
///////// swimming in the sky and flying in the sea ////////////


