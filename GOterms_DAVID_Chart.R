#1 Upload set of genes to DAVID, after 'submit list' click 'Functional Annotation Chart'. then download the data to .txt file.
# (Notice: The data should containing at least Category,Term,Fold enrichment and Count columns).

#The function below will return graph contain segnificant GO terms (Benjamini<0.05), Fold enrichment,
#and how many genes are found for each GO term for the set of genes that were analyiezd by DAVID.
#(Notice: An option to add title to the plot at line 61 of the script).

library('dplyr')
library('ggplot2')
library('egg')

Significant_GO_terms<-function(file){
  
  #Cleaning data#
  file <-file %>% 
    mutate(Category = substr(Category, 8,9),
           Term = stringr::str_remove(Term, "(.*?)~"))
  
  file$Term<-paste(file$Term,file$Category)
  
  file<-file[file$Benjamini<0.05,]
  
  file$'-log10(Benj)'<--log10(file$Benjamini)
  
  file<-file[,c("Term","Count","Fold.Enrichment","-log10(Benj)")]
  
  colnames(file)<-c('Term','Count','Fold.ERC','-log10(Benj)')
  
  #Visulaize the data#
  names_Benj<- ggplot(file, aes(x=reorder(Term,`-log10(Benj)`),
                                y=`-log10(Benj)`))+
    geom_bar(position="dodge", stat="identity")+
    coord_flip()+
    xlab('Significant GO terms')
  
  Fold_En<-ggplot(file, aes(x=reorder(Term,`-log10(Benj)`),
                            y=Fold.ERC))+
    geom_bar(position="dodge", stat="identity")+
    coord_flip()+
    theme(axis.text.y=element_blank(),
          axis.ticks.y=element_blank())+
    xlab("")
  
  Count<-ggplot(file, aes(x=reorder(Term,`-log10(Benj)`),
                          y=Count))+
    geom_bar(position="dodge", stat="identity")+
    coord_flip()+
    theme(axis.text.y=element_blank(),
          axis.ticks.y=element_blank())+
    xlab("")
  
  figure2<- ggarrange(names_Benj, Fold_En, Count,
                      ncol = 3, nrow = 1)
  return(figure2)
}

file<-read.delim('GenesToTry.txt')
Significant_GO_terms(file)

