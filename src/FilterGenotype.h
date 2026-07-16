#ifndef FilterGenotype_H_
#define FilterGenotype_H_

using namespace std ;

int parse_cmd_18(int argc, char **argv, Para_18 * para_18  )
{
//	if (argc <=2 ) {print_usage_18();return 0;}
	int err_flag = 0;

	for(int i = 1; i < argc || err_flag; i++)
	{
		if(argv[i][0] != '-')
		{
			cerr << "command option error! please check." << endl;
			return 0;
		}
		string flag=argv[i] ;
		flag=replace_all(flag,"-","");

		if (flag  == "InPut" )
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			para_18->input=argv[i];
		}
		else if (flag  == "Het" )
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			para_18->Het=atof(argv[i]);
		}
		else if (flag  == "MAF" )
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0;}
			i++;
			para_18->MAF=atof(argv[i]);
		}
		else if (flag  == "Miss" )
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0;}
			i++;
			para_18->Miss=atof(argv[i]);
		}
		else if (flag  ==  "OutPut")
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0;}
			i++;
			para_18->output=argv[i];
		}
		else if (flag  == "Cut3base")
		{
			para_18->Cut3base=1;
		}
	}
	return 1 ;
}


int Filter_genotype_main(int argc, char *argv[])
{    
	Para_18 * para_18 = new Para_18;
	if( parse_cmd_18(argc, argv, para_18 )==0)
	{
		delete  para_18 ;
		return 0 ;
	}

	igzstream IN (para_18->input.c_str(),ifstream::in); // ifstream  + gz 
	ogzstream OUT (para_18->output.c_str());

	if(!IN.good())
	{
		cerr << "open IN File error: "<<para_18->input<<endl;
		return 1;
	}

	if(!OUT.good())
	{
		cerr << "open OUT File error: "<<para_18->output<<endl;
		return 1;
	}

	map <string,string > SNP_back_Allele ;
	SNP_back_Allele["M"]="AC";SNP_back_Allele["K"]="GT";SNP_back_Allele["Y"]="CT";
	SNP_back_Allele["R"]="AG";SNP_back_Allele["W"]="AT";SNP_back_Allele["S"]="CG";
	SNP_back_Allele["C"]="CC";SNP_back_Allele["G"]="GG";SNP_back_Allele["T"]="TT";
	SNP_back_Allele["A"]="AA";
	SNP_back_Allele["-"]="NN"; SNP_back_Allele["N"]="NN";

	vector<string> inff ;

	while(!IN.eof())
	{
		string  line ;
		getline(IN,line);
		if (line.length()<=0  || line[0] == '#')  { continue ; }
		inff.clear();
		split(line,inff," \t");
		map <char , int>  Count ;
		int AA=inff.size();
		int het_count=0;
		for (int ii=2 ; ii<AA ; ii++)
		{
			string A_tmp=SNP_back_Allele[inff[ii]];
			Count[A_tmp[0]]++;
			Count[A_tmp[1]]++;
			if ( A_tmp[0]  != A_tmp[1])
			{
				het_count++;
			}
		}

		int sample=(AA-2);
		int miss=0 ;
		char best_base='N';
		char sed_base='N';
		int Max=0;
		int SeD=0;
		map <char,int>  :: iterator it=Count.begin();
		int Base_Number=0;

		for ( ; it!=Count.end(); it++ )
		{
			if ( (it->first ) == 'N' )
			{
				miss=it->second;
				continue ;
			}
			else if ((it->second)  > Max )
			{
				SeD=Max;
				sed_base=best_base;
				Max=(it->second);
				best_base=it->first;
			}
			else if ( (it->second)  >= SeD )
			{
				SeD=(it->second);
				sed_base=it->first;
			}
			Base_Number++;
		}

		if ( ( (miss*1.0)/(sample*2.0) )>(para_18->Miss) )
		{
			continue ;
		}
		else if ( (het_count*1.0/sample)>(para_18->Het) ) 
		{
			continue ;
		}
		else if ( ((SeD*1.0)/(sample*2.0))< (para_18->MAF))
		{
			continue ;
		}
		else if (Base_Number<2)
		{
			continue ;
		}		
		else if ( ((para_18->Cut3base)==1)   && ( Base_Number>2 ) )
		{
			continue ;
		}
		OUT<<line<<"\n";
	}

	OUT.close();
	IN.close();
	delete para_18 ;
	return 0;


}

#endif //FilterGenotype_H_
///////// swimming in the sky and flying in the sea ////////////
