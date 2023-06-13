#The function below will return graph contain segnificant GO terms (Benjamini<0.05), Fold enrichment,and how many genes are found in set of GO terms.

#segnificant GO terms (Benjamini<0.05), Fold enrichment,and how many genes are found,
#for shared GO term from two sets of GO terms that went through anti-join (withdrow) of a diffrent set of GO terms.
#(EXPLANATION BY EXEMPLE: There are two set of genes, each set we devide to two lists (list 1 & 2). GO terms for each list is generated.
#We want to get the GO terms that are shown though list 1 and NOT in list 2, let's call it list 3 (notice: list 3 contain GO term and not genes).
#We want to have at the end of the process the SHARED GO terms for the two list 3 we generated (Remember that we have two sets of genes).

#NOTICE!:
#file 1 = txt file contain GO terms you are intrested in
#file 2 = diffrent txt file contain GO terms you are intrested in
#file 3 = txt file contain GO terms you want to withdraw form file 1
#file 4 = txt file contain GO terms you want to withdraw form file 2

library('tidyr')
library('reshape2')
library('stringr')
library('dplyr')
library('ggplot2')
library('ggpubr')

Shared_GO_terms_2files_with_anti_join<-function(file1,file2,file3,file4){
  file1.sig<-subset(file1,Benjamini<0.05)
  file2.sig<-subset(file2,Benjamini<0.05)
  file3.sig<-subset(file3,Benjamini<0.05)
  file4.sig<-subset(file4,Benjamini<0.05)
  
  
  file1.sig<-file1.sig %>% 
    mutate(Category = substr(Category, 8,9),
           Term = stringr::str_remove(Term, "(.*?)~"))
  
  file1.sig$Term<-paste(file1.sig$Term,file1.sig$Category)
  
  file1.sig1<-file1.sig[,c('Term','Count','Fold.Enrichment','Benjamini')]
  
  
  file2.sig<-file2.sig %>% 
    mutate(Category = substr(Category, 8,9),
           Term = stringr::str_remove(Term, "(.*?)~"))
  
  file2.sig$Term<-paste(file2.sig$Term,file2.sig$Category)
  
  file2.sig1<-file2.sig[,c('Term','Count','Fold.Enrichment','Benjamini')]
  
  
  
  file3.sig<-file3.sig %>% 
    mutate(Category = substr(Category, 8,9),
           Term = stringr::str_remove(Term, "(.*?)~"))
  
  file3.sig$Term<-paste(file3.sig$Term,file3.sig$Category)
  
  file3.sig1<-file3.sig[,c('Term','Count','Fold.Enrichment','Benjamini')]
  
  
  
  file4.sig<-file4.sig %>% 
    mutate(Category = substr(Category, 8,9),
           Term = stringr::str_remove(Term, "(.*?)~"))
  
  file4.sig$Term<-paste(file4.sig$Term,file4.sig$Category)
  
  file4.sig1<-file4.sig[,c('Term','Count','Fold.Enrichment','Benjamini')]
  
  
  file1.sig1<-anti_join(file1.sig1,file3.sig1,by='Term')
  file2.sig1<-anti_join(file2.sig1,file4.sig1,by='Term')
  
  
  file1.file2.val<-inner_join(file1.sig1,file2.sig1,by='Term')
  
  file1.file2.val[,ncol(file1.file2.val)+1]<--log10(file1.file2.val$Benjamini.x)
  file1.file2.val[,ncol(file1.file2.val)+1]<--log10(file1.file2.val$Benjamini.y)
  
  
  colnames(file1.file2.val)<-c('Term','Count Genes file1',
                               'Fold.ERC file1','Benj file1','Count Genes file2',
                               'Fold.ERC file2','Benj file2',
                               '-log10(Benj file1)','-log10(Benj file2)')
  
  #Export the data
  #write_xlsx(file1.file2.val,"C:/..../similar_term_file1_file2.xlsx")
  
  ###CONSTRUCTING THE DATA###
  file1.file2.val2<-file1.file2.val[,c(1,8,9)]
  file1.file2.val2<-melt(file1.file2.val2,id='Term')
  colnames(file1.file2.val2)<-c('Term','var.-log10(Benj)','val.-log10(Benj)')
  file1.file2.val2$num<-seq.int(nrow(file1.file2.val2))
  
  file1.file2.val3<-file1.file2.val[,c(1,2,5)]
  file1.file2.val3<-melt(file1.file2.val3,id='Term')
  colnames(file1.file2.val3)<-c('Term','var.Count','val.Count')
  file1.file2.val3$num<-seq.int(nrow(file1.file2.val3))
  
  
  file1.file2.val1<-file1.file2.val[,c(1,3,6)]
  file1.file2.val1<-melt(file1.file2.val1,id='Term')
  colnames(file1.file2.val1)<-c('Term','var.ERC','val.ERC')
  file1.file2.val1$num<-seq.int(nrow(file1.file2.val1))
  
  
  d<-inner_join(file1.file2.val2,file1.file2.val3,by='num')
  dd<-inner_join(d,file1.file2.val1,by='num')
  dd<-dd[,-c(1,5)]
  
  #Visualize the data
  Benj<-ggplot(dd,aes(reorder(Term,-`val.-log10(Benj)`),`val.-log10(Benj)`,fill=`var.-log10(Benj)`)) +
    geom_bar(position="dodge", stat="identity") +
    scale_x_discrete(labels=function(x) str_wrap(x,width = 6)) +
    ylab('-log10(Benjamini)')+
    scale_fill_manual(values = c("#33bd68", "#539ed6"))
  
  ERC<-ggplot(dd,aes(reorder(Term,-`val.-log10(Benj)`),`val.ERC`,fill=`var.ERC`)) +
    geom_bar(position="dodge", stat="identity") +
    scale_x_discrete(labels=function(x) str_wrap(x,width = 7)) +
    ylab('Fold Enrichment')+
    theme(axis.text.x=element_blank(), 
          axis.ticks.x=element_blank() 
    ) +xlab("")+
    scale_fill_manual(values = c("#33bd68", "#539ed6"))
  
  Count<-ggplot(dd,aes(reorder(Term,-`val.-log10(Benj)`),`val.Count`,fill=`var.Count`)) +
    geom_bar(position="dodge", stat="identity") +
    scale_x_discrete(labels=function(x) str_wrap(x,width = 7)) +
    ylab('Gene count')+
    theme(axis.text.x=element_blank(), 
          axis.ticks.x=element_blank() 
    ) +xlab("")+
    scale_fill_manual(values = c("#33bd68", "#539ed6"))
  
  figure<-ggarrange(ERC, Count,Benj,
                    ncol = 1, nrow = 3)
  return(figure)
}

file1<-read.delim('gene_list1.txt')
file2<-read.delim('gene_list2.txt')
file3<-read.delim('gene_list3.txt')
file4<-read.delim('gene_list4.txt')

Shared_GO_terms_2files_with_anti_join(file1,file2,file3,file4)