#ifndef EHH_H_
#define EHH_H_

using namespace std;


string int2str( unsigned short int val )
{
	ostringstream out;
	out<<val;
	return out.str();
}


int EHH_Region_LDDecay(In3str1v * paraFA04, Para_18 * para_18 , map <string,map <llong, vector <BaseType> > >  & SNPList, int & Flag_for_pro)
{

	string Stat=(paraFA04->InStr2)+".ehh.gz";
	ogzstream OUT ((Stat).c_str());
	if((!OUT.good()))
	{
		cerr << "open OUT File error: "<<Stat<<endl;
		return 0;
	}


	vector<string> inf;
	split((paraFA04->EHH),inf,":");
	llong StartSite=atoi(inf[1].c_str());
	string chrName=inf[0];
	unsigned short int ValueE=2;
	map <string,map <llong, vector <BaseType> > >  :: iterator itSNP=SNPList.find(chrName);
	if (itSNP == SNPList.end())
	{
		cerr<<"\t\tInPut Para -EHH  chromosome ["<<chrName<<"]  can't be found in the SNP dataset\n";		
		return 0;
	}  

	llong BeginRegion=StartSite-(paraFA04->InInt);   if (BeginRegion<0) { BeginRegion=0;}	
	llong EndRegion=StartSite+(paraFA04->InInt);

	map<llong, vector <BaseType> >  :: iterator key2=(itSNP->second).begin(); ;

	map <llong, vector <BaseType> >  NewSNPList ;

	int count =0 ;
	int TotalSNP=0;

	statementVar Var;
	Var.Asize= (key2->second).size();

	for(   ; key2!=(itSNP->second).end(); key2++)
	{
		if ( ( key2->first )<  BeginRegion)
		{
			continue;
		}
		if  ( ( key2->first )> EndRegion)
		{
			break;
		}
		TotalSNP++;
		bool tr=true  ;
		for (Var.i=0 ;Var.i<Var.Asize ; (Var.i)++)
		{
			if ((key2->second)[Var.i].Value == ValueE ) 
			{
				tr=false ;
				break;
			}
		}

		if  (tr)
		{
			count++;
			NewSNPList.insert(map <llong, vector <BaseType> >  :: value_type(key2->first,(key2->second))) ;
		}
	}

	cout<<"##Start EHH region :"<<chrName<<" "<<BeginRegion<<" "<<EndRegion<<"; In This Region TotalSNP Number is "<<TotalSNP<<",No Missing SNP Site Number is "<<count<<endl;
	if (count < 8)
	{
		cerr<<"EHH should be No Missing SNP Number  > 8 \n";
		return 0;
	}
	else if ( count  > 168888)
	{
		cerr<<"Warning: EHH Region SNP Number too much,you may use the small region or more stringent conditions to filter the SNP\n";
	}

	OUT<<"#Chr\tSite\tDist\tEHH_all\tEHH_0\tEHH_1\n";
	OUT<<"#SNP_Number\t"<<count<<"\n";

	cout<<"##Begin Cal EHH...\n";
	int fengmu=Var.Asize*(Var.Asize-1);

	map<llong, vector <BaseType> > :: iterator key3=NewSNPList.end();
	for ( key2 = NewSNPList.begin() ;  key2!=NewSNPList.end()   ; key2++)
	{
		if (( key2 ->first ) <=  StartSite)
		{
			key3=key2 ;
		}
		else
		{
			break ;
		}
	}

	if (key3==NewSNPList.begin()) {key3++;}
	if (key3==NewSNPList.end()) {key3--; key3--;}

	map <llong ,PairInfoV1 >   ResultEHH;
	vector <string> Haplotype(Var.Asize ,"") ;

	PairInfoV1 EHH_Value; 
	int CalTmp=0;

	for (Var.i=0 ;Var.i<Var.Asize ; (Var.i)++)
	{
		Haplotype[Var.i]=int2str(((key3->second)[Var.i]).Value);
		if ((((key3->second)[Var.i]).Value) == 0)
		{
			CalTmp++;
		}
	}

	EHH_Value.RR=(CalTmp*(CalTmp-1)*1.0)/fengmu;
	CalTmp=Var.i-CalTmp;
	EHH_Value.D=((CalTmp*(CalTmp-1)*1.0)/fengmu)+EHH_Value.RR;
	ResultEHH.insert(map <llong ,PairInfoV1> :: value_type(key3->first,EHH_Value));

	key2 = key3 ;	key2 --;

	for (   ;  ; key2--)
	{
		for (Var.i=0 ;Var.i<Var.Asize ; (Var.i)++)
		{
			Haplotype[Var.i]=Haplotype[Var.i]+int2str(((key2->second)[Var.i]).Value);
		}
		if (Haplotype[0].empty() || (Haplotype[0].size()<2)) { continue ;}
		map <string , int > HapNum;
		HapNum[Haplotype[0]]=1;
		map <string , int >  :: iterator HapNum_key;

		for (Var.i=1 ;Var.i<Var.Asize ; (Var.i)++)
		{
			HapNum_key=HapNum.find(Haplotype[Var.i]);
			if (HapNum_key==HapNum.end())
			{
				HapNum.insert(map <string , int > :: value_type( Haplotype[Var.i], 1) );
			}
			else
			{
				(HapNum_key->second)++;
			}
		}
		HapNum_key=HapNum.begin();
		int fengzi=0;
		int fengzi_A=0;
		for (  ;  HapNum_key!=HapNum.end();  HapNum_key++ )
		{
			CalTmp=((HapNum_key->second)*((HapNum_key->second)-1));
			fengzi+=CalTmp;
			if ( (HapNum_key->first)[0] == '0' )
			{
				fengzi_A+=CalTmp;
			}
		}

		EHH_Value.D=fengzi*1.0/fengmu;
		EHH_Value.RR=fengzi_A*1.0/fengmu;

		ResultEHH.insert(map <llong ,PairInfoV1> :: value_type(key2->first,EHH_Value));	
		if (EHH_Value.D < 0.088) {break;}
		if  (key2 == NewSNPList.begin())
		{
			break;
		}
	}




	key2=key3;
	for (Var.i=0 ;Var.i<Var.Asize ; (Var.i)++)
	{
		Haplotype[Var.i]=int2str(((key3->second)[Var.i]).Value);
	}

	for ( key2++ ;   key2!=NewSNPList.end() ; key2++)
	{
		for (Var.i=0 ;Var.i<Var.Asize ; (Var.i)++)
		{
			Haplotype[Var.i]=Haplotype[Var.i]+int2str(((key2->second)[Var.i]).Value);
		}

		if  ( Haplotype[0].empty() || Haplotype[0].size()<2) { continue ;}
		map <string , int > HapNum;
		HapNum[Haplotype[0]]=1;
		map <string , int >  :: iterator HapNum_key;
		for (Var.i=1 ;Var.i<Var.Asize ; (Var.i)++)
		{
			HapNum_key=HapNum.find(Haplotype[Var.i]);
			if (HapNum_key==HapNum.end())
			{
				HapNum.insert(map <string , int > :: value_type( Haplotype[Var.i], 1) );
			}
			else
			{
				(HapNum_key->second)++;
			}
		}

		HapNum_key=HapNum.begin();
		int fengzi=0;
		int fengzi_A=0;
		for (  ;  HapNum_key!=HapNum.end();  HapNum_key++ )
		{
			CalTmp=((HapNum_key->second)*((HapNum_key->second)-1));
			fengzi+=CalTmp;
			if ( (HapNum_key->first)[0] == '0' )
			{
				fengzi_A+=CalTmp;
			}
		}

		EHH_Value.D=fengzi*1.0/fengmu;
		EHH_Value.RR=fengzi_A*1.0/fengmu;

		ResultEHH.insert(map <llong ,PairInfoV1> :: value_type(key2->first,EHH_Value));
		if (EHH_Value.D < 0.088) {break;}
		if  (key2 == NewSNPList.begin())
		{
			break;
		}



	}

	map <llong ,PairInfoV1> :: iterator MapPaRsult ;
	for ( key2 = NewSNPList.begin() ;  key2!=NewSNPList.end()   ; key2++)
	{
		int Dis=(key2->first)-StartSite;
		MapPaRsult=ResultEHH.find(key2->first);
		if (MapPaRsult==ResultEHH.end())
		{
		OUT<<chrName<<"\t"<<key2->first<<"\t"<<Dis<<"\t0.0000\t0.0000\t0.0000\n";
		}
		else
		{
		OUT<<chrName<<"\t"<<key2->first<<"\t"<<Dis<<"\t"<<setprecision(4)<<setiosflags(ios::right)<<setiosflags(ios::fixed)<<(MapPaRsult->second).D<<"\t"<<(MapPaRsult->second).RR<<"\t"<<((MapPaRsult->second).D-(MapPaRsult->second).RR)<<"\n";
		}

	}
	OUT.close();







	/////////////////  plot EHH ////////////////////// 

	char   buf[2048]={'\0'};
	string cc="which  Rscript  2> /dev/null ";
	memset( buf, '\0', sizeof(buf) );
	FILE   *stream ;
	stream=popen(cc.c_str(),"r") ;
	fread( buf, sizeof(char), sizeof(buf), stream);
	string binPath=buf;
	binPath=binPath.substr(0,binPath.length()-1);

	string OutPlotr=(paraFA04->InStr2)+".tmp.r";
	ofstream  OUTR (OutPlotr.c_str());
	OUTR<<""
		"\n"
		"read.table(\""<<Stat<<"\")->r;\n"
		"pdf(\""<<(paraFA04->InStr2)<<".ehh.pdf\");\n"
		"plot(r[,3]/1000,r[,4],col=\"blue\",type=\"l\",ylab=\"EHH\",main=\"EHH decay\",,bty=\"n\",xlab=\"Distance from core region "<<chrName<<":"<<StartSite<<" (Kb)\")\n"
		"dev.off();\n"
		"png(\""<<(paraFA04->InStr2)<<".ehh.png\");\n"
		"plot(r[,3]/1000,r[,4],col=\"blue\",type=\"l\",ylab=\"EHH\",main=\"EHH decay\",,bty=\"n\",xlab=\"Distance from core region "<<chrName<<":"<<StartSite<<" (Kb)\")\n"
		"dev.off();\n"
		"\n"
		<<endl;
	OUTR.close();

	if (binPath == "")
	{
		cout <<"\twarning: can't find the [Rscript] in your $PATH ; no png Figure Out"<<endl;
		cout <<"\t\tRscript "<<OutPlotr<<endl;
	}
	else
	{
		cc=binPath+"\t"+OutPlotr ;
		if (paraFA04->TF)
		{
			cc=binPath+"\t"+OutPlotr+"  ; rm  -rf  "+ OutPlotr  ;
		}
		std::system(cc.c_str()) ;
	}

	cc="rm -rf "+(paraFA04->InStr2)+".stat.gz";
	std::system(cc.c_str());

	return 1;
}
#endif
