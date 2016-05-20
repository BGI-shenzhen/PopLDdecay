# PopLDdecay
PopLDdecay: A new simple and efficient software for Linkage Disequilibrium Decay analysis based Variant Call Format



1) Install
        tar -zxvf  PopLDdecayXXX.tar.gz
        cd PopLDdecayXXX;
        cd src;
        sh  make.sh
        ../bin/PopLDdecay
  
2) Example
        1) Calculate LD decay
        ./bin/PopLDdecay    -InVCF  SNP.vcf.gz  -OutStat LDdecay
        2) draw the Figure
        perl  bin/Plot_LDDecay.pl  -inFile   LDeecay.stat.gz  -output  Fig
        3) you can see the  result  [LDeecay.stat.gz]  and   [Fig.png Fig.pdf]



3) Introduction
Linkage disequilibrium (LD) decay[1] is the most important and most common
analysis in the population resequencing[2]. Special in the self-pollinated
crops, the LD decay may not only reveal much about domestication and breed
history[3], but also can reveal gene flow phenomenon, selection regions[1].
However, to measure the LD decay, it takes too much resources and time by
using currently existent software and tools. The LD decay studies also
generate extraordinarily large amounts of data to temporary storage when you
using the mainstream software "Haploview"[4], the classical LD
processing tools. Effective use and analysis to get the LD decay result
remains a difficult task for individual researchers. Here, we introduce
PopLDdecay, a simple- efficient software for LD decay analysis, which
processes the Variant Call Format (VCF)[5] file to produce the LD decay
statistics results and plot the LD decay graphs. PopLDdecay is designed to use
compressed data files as input or output to save storage space and it
facilitates faster and more computationally efficient than the currently
existent softwares. This software makes the LD decay pipeline significantly
simplifies and user-friendly for the users.
