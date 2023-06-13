#1 Upload set of genes to DAVID, after 'submit list' click 'Functional Annotation Chart'. then download the data to .txt file.
# (Notice: The data should containing at least Category,Term,Fold enrichment and Count columns).
#2 Repeat the process for another set of genes you are intrested to compare with the first set.

#The function below will return graph contain segnificant GO terms (Benjamini<0.05), Fold enrichment,
#and how many genes are found for each GO term, That are in file1 AND NOT in file2.

library("dplyr")
library('ggplot2')
library('egg')

fun_in1_not_in2<-function(file1,file2){
  
  #Cleaning the data#
  file1.only<-anti_join(file1,file2,by='Term')
  file1.only<-file1.only[file1.only$Benjamini<0.05,]
  
  file1.only<-file1.only %>% 
    mutate(Category = substr(Category, 8,9),
           Term = stringr::str_remove(Term, "(.*?)~"))
  
  file1.only$Term<-paste(file1.only$Term,file1.only$Category)
  
  file1.only1<-file1.only[,c("Term","Count","Fold.Enrichment","Benjamini")]
  file1.only1$'-log10(Benjamini)'<--log10(file1.only1$Benjamini)
  
  #Visualizing the data#
  names_Benj1<- ggplot(file1.only1, aes(x=reorder(Term,-log10(Benjamini)),
                                        y=-log10(Benjamini)))+
    geom_bar(position="dodge", stat="identity")+
    coord_flip()
  
  Fold_En1<-ggplot(file1.only1, aes(x=reorder(Term,-log10(Benjamini)),
                                    y=Fold.Enrichment))+
    geom_bar(position="dodge", stat="identity")+
    coord_flip()+
    theme(axis.text.y=element_blank(),
          axis.ticks.y=element_blank())+
    xlab("")
  
  Count1<-ggplot(file1.only1, aes(x=reorder(Term,-log10(Benjamini)),
                                  y=Count))+
    geom_bar(position="dodge", stat="identity")+
    coord_flip()+
    theme(axis.text.y=element_blank(),
          axis.ticks.y=element_blank())+
    xlab("")
  
  figure<-ggarrange(names_Benj1, Fold_En1, Count1,
                    ncol = 3, nrow = 1)
  return(figure)  
}

file1<-read.delim('gene_list1.txt')
file2<-read.delim('gene_list2.txt')

dafun<-fun_in1_not_in2(file1,file2)


