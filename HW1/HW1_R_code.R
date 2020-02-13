#1

install.packages(c('tidyverse','nycflights13','ggplot2','gridExtra','grid','reshape','stringr'))
library(tidyverse)
library(nycflights13)
library(ggplot2)
library(gridExtra)
library(grid)
library(reshape)
library(stringr)

##(a)
count = as_tibble(read_csv("count.csv"))
mrna = as_tibble(read_csv("mRNA.csv"))
run_info = as_tibble(read_csv("Run_Info.csv"))
treatment_info = as_tibble(read_csv("Treatment_Info.csv"))

#show that these keys uniquely identifies observations in these datasets:
count %>% count(ensgene) %>% filter(n > 1)
mrna %>% count(x) %>% filter(n > 1)
run_info %>% count(Run) %>% filter(n > 1)
treatment_info %>% count(Sample) %>% filter(n > 1)

##(b)
colnames(treatment_info)[1] = 'sample'
metadata = as_tibble(full_join(run_info, treatment_info, by = 'sample'))
metadata

##(c)
metadata %>% 
  group_by(Treatment) %>%
  summarise(mean = mean(avgLength, na.rm = TRUE))
metadata %>% 
  group_by(Cell_type) %>%
  summarise(mean = mean(avgLength, na.rm = TRUE))
#Plot the average length by treatment type or average length by cell type
options(repr.plot.width=7,repr.plot.height=3)
theme_update(plot.title = element_text(hjust = 0.5,size=8))
p1 = ggplot(data = metadata, aes(x = Treatment , y = avgLength))+geom_jitter(width = .5,height = 0,size = 1,aes(colour = Treatment))+ggtitle('Average length by treatment type')
p2 = ggplot(data = metadata, aes(x = Cell_type , y = avgLength))+geom_jitter(width = .5,height = 0,size = 1,aes(colour = Cell_type))+ggtitle('Average length by cell type')
grid.arrange(p1, p2, nrow = 1, top = textGrob('Average length by treatment type and cell type',gp=gpar(fontsize=20,font=1)))

##(d)
mnc <- count %>%
  pivot_longer(
    cols = -ensgene,
    names_to = 'Run'
  ) %>%
  left_join(metadata, by = 'Run')
dim(mnc)
mnc[1:10,]

##(e)
bp1 = ggplot(data = mnc)
bp1 + geom_boxplot(aes(y = value, x = Treatment, fill = Treatment))
ndf<-data.frame(mnc,log_count=log2(mnc$value+1))
bp1_l = ggplot(data = ndf)
bp1_l + geom_boxplot(aes(y = log_count, x = Treatment, fill = Treatment))

##(f)
bp2 = ggplot(data = mnc)
bp2 + geom_boxplot(aes(y = value, x = Cell_type, fill = Treatment))
bp2_l = ggplot(data = ndf)
bp2_l + geom_boxplot(aes(y = log_count, x = Cell_type, fill = Cell_type))

##(g)
bp3 = ggplot(data = mnc)
bp3 + geom_boxplot(aes(y = value, x = Cell_type, fill = Treatment))
bp3_l = ggplot(data = ndf)
bp3_l + geom_boxplot(aes(y = log_count, x = Cell_type, fill = Treatment))

##(h)
ext_gene = ndf[which(ndf$ensgene=='ENSG00000064607'),]
ext_gene
options(repr.plot.width=12,repr.plot.height=3)
theme_update(plot.title = element_text(hjust = 0.5,size=12))
ggplot(data = ext_gene,aes(x = sample, y = 0, col = Treatment))+geom_jitter(size = 5,width = 0, height = 0)+ggtitle('Expression profile of gene ENSG00000064607  against the samples')+theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())

#2

##(a)
mRNA_first<-read.csv("mRNA.csv")
(mRNA<-toString(mRNA_first[[1]]))
str_locate_all(mRNA,"AUG")

##(b)
str_locate_all(mRNA,c('UGA','UAA','UAG'))

##(c)
# first
sf = substr(mRNA,(str_locate_all(mRNA,"AUG")[[1]][1]),(str_locate_all(mRNA,"AUG")[[1]][2])-1)
print(c(nchar(sf),sf))
#between
for (i in 2:(length(str_locate_all(mRNA,"AUG")[[1]])/2-1)){
  s = substr(mRNA,(str_locate_all(mRNA,"AUG")[[1]][i]),(str_locate_all(mRNA,"AUG")[[1]][i+1]))
  print(c(nchar(s),s))
}
#last
sl = substr(mRNA,(str_locate_all(mRNA,"AUG")[[1]][(length(str_locate_all(mRNA,"AUG")[[1]])/2)]),length(str_split(mRNA,'')[[1]]))
print(c(nchar(sl),sl))

##(d)
genetic_code = function(x){
  s0 = ''
  for (i in 1:(nchar(x)/3)){
    if ((substr(x,3*i-2,3*i)=='UUU')|(substr(x,3*i-2,3*i)=='UUC')){
      s0 = str_c(s0,'-Phe')
    }
    if ((substr(x,3*i-2,3*i)=='UUA')|(substr(x,3*i-2,3*i)=='UUG')|(substr(x,3*i-2,3*i)=='CUU')|(substr(x,3*i-2,3*i)=='CUC')|(substr(x,3*i-2,3*i)=='CUA')|(substr(x,3*i-2,3*i)=='CUG')){
      s0 = str_c(s0,'-Leu')
    }
    if ((substr(x,3*i-2,3*i)=='UCU')|(substr(x,3*i-2,3*i)=='UCC')|(substr(x,3*i-2,3*i)=='UCA')|(substr(x,3*i-2,3*i)=='UCG')|(substr(x,3*i-2,3*i)=='AGU')|(substr(x,3*i-2,3*i)=='AGC')){
      s0 = str_c(s0,'-Ser')
    }
    if ((substr(x,3*i-2,3*i)=='UAU')|(substr(x,3*i-2,3*i)=='UAC')){
      s0 = str_c(s0,'-Tyr')
    }
    if ((substr(x,3*i-2,3*i)=='UGU')|(substr(x,3*i-2,3*i)=='UGC')){
      s0 = str_c(s0,'-Cys')
    }
    if ((substr(x,3*i-2,3*i)=='CCU')|(substr(x,3*i-2,3*i)=='CCC')|(substr(x,3*i-2,3*i)=='CCA')|(substr(x,3*i-2,3*i)=='CCG')){
      s0 = str_c(s0,'-Pro')
    }
    if ((substr(x,3*i-2,3*i)=='CAU')|(substr(x,3*i-2,3*i)=='CAC')){
      s0 = str_c(s0,'-His')
    }
    if ((substr(x,3*i-2,3*i)=='UGG')){
      s0 = str_c(s0,'-Trp')
    }
    if ((substr(x,3*i-2,3*i)=='CAA')|(substr(x,3*i-2,3*i)=='CAG')){
      s0 = str_c(s0,'-Gln')
    }
    if ((substr(x,3*i-2,3*i)=='CGU')|(substr(x,3*i-2,3*i)=='CGC')|(substr(x,3*i-2,3*i)=='CGA')|(substr(x,3*i-2,3*i)=='CGG')|(substr(x,3*i-2,3*i)=='AGA')|(substr(x,3*i-2,3*i)=='AGG')){
      s0 = str_c(s0,'-Arg')
    }
    if ((substr(x,3*i-2,3*i)=='AUU')|(substr(x,3*i-2,3*i)=='AUC')|(substr(x,3*i-2,3*i)=='AUA')){
      s0 = str_c(s0,'-Ile')
    }
    if ((substr(x,3*i-2,3*i)=='AUG')){
      s0 = str_c(s0,'-Met')
    }
    if ((substr(x,3*i-2,3*i)=='ACU')|(substr(x,3*i-2,3*i)=='ACC')|(substr(x,3*i-2,3*i)=='ACA')|(substr(x,3*i-2,3*i)=='ACG')){
      s0 = str_c(s0,'-Thr')
    }
    if ((substr(x,3*i-2,3*i)=='AAU')|(substr(x,3*i-2,3*i)=='AAC')){
      s0 = str_c(s0,'-Asn')
    }
    if ((substr(x,3*i-2,3*i)=='AAA')|(substr(x,3*i-2,3*i)=='AAG')){
      s0 = str_c(s0,'-Lys')
    }
    if ((substr(x,3*i-2,3*i)=='GUU')|(substr(x,3*i-2,3*i)=='GUC')|(substr(x,3*i-2,3*i)=='GUA')|(substr(x,3*i-2,3*i)=='GUG')){
      s0 = str_c(s0,'-Val')
    }
    if ((substr(x,3*i-2,3*i)=='GCU')|(substr(x,3*i-2,3*i)=='GCC')|(substr(x,3*i-2,3*i)=='GCA')|(substr(x,3*i-2,3*i)=='GCG')){
      s0 = str_c(s0,'-Ala')
    }
    if ((substr(x,3*i-2,3*i)=='GGU')|(substr(x,3*i-2,3*i)=='GGC')|(substr(x,3*i-2,3*i)=='GGA')|(substr(x,3*i-2,3*i)=='GGG')){
      s0 = str_c(s0,'-Gly')
    }
    if ((substr(x,3*i-2,3*i)=='GAU')|(substr(x,3*i-2,3*i)=='GAC')){
      s0 = str_c(s0,'-Asp')
    }
    if ((substr(x,3*i-2,3*i)=='GAA')|(substr(x,3*i-2,3*i)=='GAC')){
      s0 = str_c(s0,'-Glu')
    }
    if ((substr(x,3*i-2,3*i)=='UAA')|(substr(x,3*i-2,3*i)=='UAG')|(substr(x,3*i-2,3*i)=='UGA')){
      break
    }
  }
  return(substring(s0,2))
}

##(e)
cat(genetic_code(sf),sep="\n")
for (i in 2:(length(str_locate_all(mRNA,"AUG")[[1]])/2-1)){
  s = substr(mRNA,(str_locate_all(mRNA,"AUG")[[1]][i]),(str_locate_all(mRNA,"AUG")[[1]][i+1]))
  cat(genetic_code(s),sep="\n")
}
cat(genetic_code(sl),sep="\n")

##(f)
c = (nchar(genetic_code(sf))+1)/4
print(str_c("First one has ", c))
print("For Second to the one before last one")
for (i in 2:(length(str_locate_all(mRNA,"AUG")[[1]])/2-1)){
  s = substr(mRNA,(str_locate_all(mRNA,"AUG")[[1]][i]),(str_locate_all(mRNA,"AUG")[[1]][i+1]))
  c = (nchar(genetic_code(s))+1)/4
  print(str_c(i,": ",c))
}
c = (nchar(genetic_code(sl))+1)/4
print(str_c("Last one has ", c))