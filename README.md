

# PopLDdecay
<b>[PopLDdecay: a fast and effective tool for linkage disequilibrium decay analysis based on variant call format files](https://doi.org/10.1093/bioinformatics/bty875)</b>

The <i><b>[PopLDDecay article](https://doi.org/10.1093/bioinformatics/bty875)</b></i> has been published in <b> [Bioinformatics](https://doi.org/10.1093/bioinformatics/bty875) magazine</b>, please cited this article if possible

PMID: [30321304](https://www.ncbi.nlm.nih.gov/pubmed/30321304)   &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;            DOI:[10.1093/bioinformatics/bty875](https://doi.org/10.1093/bioinformatics/bty875)

###  1) Install
------------

<b> [Download](https://github.com/BGI-shenzhen/PopLDdecay/archive/v3.40.tar.gz) </b>

</br>
Method1 For <b>linux/Unix</b> and <b> macOS </b>
<pre>
        git clone https://github.com/BGI-shenzhen/PopLDdecay.git 
        cd PopLDdecay; chmod 755 configure; ./configure;
        make;
        mv PopLDdecay  bin/;    #     [rm *.o]
</pre>

**Note:** If fail to link,try to <b>re-install</b> the libraries [**_zlib_**](https://zlib.net/)


Method2 For <b>linux/Unix</b> and <b> macOS </b>
<pre>
        tar -zxvf  PopLDdecayXXX.tar.gz
        cd PopLDdecayXXX;
        cd src;
        make ; make clean                            # or [sh make.sh]
        ../bin/PopLDdecay
</pre>
**Note:** If fail to link,try to <b>re-install</b> the libraries [**_zlib_**](https://zlib.net/)


###  2) Example
------------

see more detailed Usage in the <b>[Documentation](https://github.com/BGI-shenzhen/PopLDdecay/blob/master/Manual.pdf)</b>


* 1) Calculate LD decay 
<pre>
      # 1)  For gatk VCF file deal , run PopLDdecay  direct
            ./bin/PopLDdecay    -InVCF  SNP.vcf.gz  -OutStat LDdecay   
      # 2)  For plink [.ped .map], chang plink 2 genotype first  2) run PopLDdecay  
            perl bin/mis/plink2genotype.pl    -inPED in.ped -inMAP in.map  -outGenotype out.genotype ;      ./bin/PopLDdecay        -InGenotype out.genotype -OutStat LDdecay 
      # 3)  To Calculate the <b>subgroup GroupA LDdecay</b> in VCF Files   # put GroupA sample name into GroupA_sample.list
            ./bin/PopLDdecay   -InVCF  <in.vcf.gz>  -OutStat <out.stat>   -SubPop    GroupA_sample.list
</pre>

* 2) draw the Figure
```
        #    2.1  For one Population
        perl  bin/Plot_OnePop.pl  -inFile   LDdecay.stat.gz  -output  Fig
        #    2.2  For one Population  muti chr          # List Format [chrResultPathWay]
        perl  bin/Plot_OnePop.pl  -inList   Chr.ReslutPath.List  -output Fig
        #    2.3  For muti Population                   #  List Format :[Pop.ResultPath  PopID ]
        perl  bin/Plot_MutiPop.pl  -inList  Pop.ReslutPath.list  -output Fig
```
* 3) see the result  [LDdecay.stat.gz] and [Fig.png Fig.pdf]

###  3) Introduction
------------
Linkage disequilibrium (LD) decay[1] is the most important and most common analysis in the population resequencing[2]. Special in the self-pollinated crops, the LD decay may not only reveal much about domestication and breed history[3], but also can reveal gene flow phenomenon, selection regions[1].However, to measure the LD decay, it takes too much resources and time by using currently existent software and tools. The LD decay studies also generate extraordinarily large amounts of data to temporary storage when you using the mainstream software "Haploview"[4], the classical LD processing tools. Effective use and analysis to get the LD decay result remains a difficult task for individual researchers. Here, we introduce PopLDdecay, a simple- efficient software for LD decay analysis, which processes the Variant Call Format (VCF)[5] file to produce the LD decay statistics results and plot the LD decay graphs. PopLDdecay is designed to use compressed data files as input or output to save storage space and it facilitates faster and more computationally efficient than the currently existent softwares. This software makes the LD decay pipeline significantly
* <b> Parameter description</b>
```php
	Usage: PopLDDecay -InVCF  <in.vcf.gz>  -OutStat <out.stat>

		-InVCF       <str>    Input SNP VCF Format
		-InGenotype  <str>    Input SNP Genotype Format
		-OutStat     <str>    OutPut Stat Dist ~ r^2 File

		-SubPop      <str>    SubGroup SampleList of VCFFile [ALLsample]
		-MaxDist     <int>    Max Distance (kb) between two SNP [300]
		-MAF         <float>  Min minor allele frequency filter [0.005]
		-Het         <float>  Max ratio of het allele filter [0.88]
		-Miss        <float>  Max ratio of miss allele filter [0.25]
		-EHH         <str>    To Run EHH Region decay set StartSite [NA]
		-OutFilterSNP         OutPut the final SNP to calculate
		-OutType     <int>    1: R^2 result 2: R^2 & D' result 3:PairWise LD Out[1]
		                      See the Help for more OutType [1-8] details
		
		-help                 Show more help [hewm2008 v3.40]

```


###  4) Results
------------
some LD decay images which I draw in the paper before.

* [50 Rices NBT](http://www.nature.com/nbt/journal/v30/n1/images/nbt.2050-F2.jpg)
* [31 soybeans  NG]( http://www.nature.com/ng/journal/v42/n12/images/ng.715-F1.jpg)

###  5) Discussing
------------
- [:email:](https://github.com/BGI-shenzhen/PopLDdecay) hewm2008@gmail.com / hewm2008@qq.com
- join the<b><i> QQ Group : 125293663</b></i>
- <i><b>[PopLDDecay article](https://doi.org/10.1093/bioinformatics/bty875)</b></i> has been published in <b> [Bioinformatics](https://doi.org/10.1093/bioinformatics/bty875) magazine</b>, please cited this article if possible
- PMID: [30321304](https://www.ncbi.nlm.nih.gov/pubmed/30321304)   &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;            DOI:[10.1093/bioinformatics/bty875](https://doi.org/10.1093/bioinformatics/bty875)

######################swimming in the sky and flying in the sea #############################
