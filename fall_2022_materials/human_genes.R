library(dplyr)
library(ggplot2)
library(gridExtra)
library(glue)

human_genes<-readr::read_csv("human_genes.csv")


human_genes<-human_genes %>% select(name, size, description)


human_genes %>% dplyr::sample_n(20)


# define a custom function:
my_func <- function(df) {
  df %>% 
    dplyr::summarise(MeanLength = mean(size),
                     SD_Length=sd(size))
} 

# sample repeatedly (100 times) 100 rows of df and store in a list
# apply the custom function to each sample in the list,
# bind rows together and create an index column, all in a "pipe":


human_genes_15<-filter(human_genes, size<=15000)

all_summ<-human_genes %>% summarise(MeanLength=mean(size), SDLength=sd(size))

p0<-ggplot(human_genes_15, aes(x=size)) + 
  geom_histogram(aes(y = after_stat(count / sum(count))), bins = 25,color="lightgray", fill="black")+
  ylab("Population Relative Frequency") +
  xlab("Gene Length (in base pairs") +
  xlim(0,15000)+
  ylim(0,0.20)+
  ggtitle("Histogram of human gene lengths (in basepairs)")+
  geom_vline(all_summ, mapping=aes(xintercept=MeanLength), color="orange")

set.seed(1)
my_tab10000<-replicate(10000, sample_n(human_genes, 5), simplify = FALSE) %>%
  lapply(., my_func) %>% 
  bind_rows %>%
  mutate(replicate = 1:n()) 
set.seed(1)
my_tab20000<-replicate(20000, sample_n(human_genes, 5), simplify = FALSE) %>%
  lapply(., my_func) %>% 
  bind_rows %>%
  mutate(replicate = 1:n()) 

set.seed(1)
my_tab1000<-replicate(1000, sample_n(human_genes, 5), simplify = FALSE) %>%
  lapply(., my_func) %>% 
  bind_rows %>%
  mutate(replicate = 1:n()) 

set.seed(1)
my_tab100<-replicate(100, sample_n(human_genes, 5), simplify = FALSE) %>%
  lapply(., my_func) %>% 
  bind_rows %>%
  mutate(replicate = 1:n()) 

set.seed(1)
my_tab10<-replicate(10, sample_n(human_genes, 5), simplify = FALSE) %>%
  lapply(., my_func) %>% 
  bind_rows %>%
  mutate(replicate = 1:n())

set.seed(1)
my_tab2100<-replicate(100, sample_n(human_genes, 100), simplify = FALSE) %>%
  lapply(., my_func) %>% 
  bind_rows %>%
  mutate(replicate = 1:n())

set.seed(1)
my_tab21000<-replicate(1000, sample_n(human_genes, 100), simplify = FALSE) %>%
  lapply(., my_func) %>% 
  bind_rows %>%
  mutate(replicate = 1:n())

set.seed(1)
samp1<-sample_n(human_genes,100)
set.seed(1)
replicate(10000, sample_n(samp1, 99, replace=T), simplify=F) %>%
  lapply(., my_func) %>% 
  bind_rows %>%
  mutate(replicate = 1:n()) %>% summarise(SSAMP=sd(MeanLength))

my_summ10<-my_tab10 %>% summarise(GrandMean=mean(MeanLength), MeanSD=mean(SD_Length))
my_summ100<-my_tab100 %>% summarise(GrandMean=mean(MeanLength), MeanSD=mean(SD_Length))

my_summ1000<-my_tab1000 %>% summarise(GrandMean=mean(MeanLength), MeanSD=mean(SD_Length))
my_summ10000<-my_tab10000 %>% summarise(GrandMean=mean(MeanLength), MeanSD=mean(SD_Length))
my_summ20000<-my_tab20000 %>% summarise(GrandMean=mean(MeanLength), MeanSD=mean(SD_Length))
my_summ2100<-my_tab2100 %>% summarise(GrandMean=mean(MeanLength), MeanSD=mean(SD_Length))
my_summ21000<-my_tab21000 %>% summarise(GrandMean=mean(MeanLength), MeanSD=mean(SD_Length))

p1<-ggplot(my_tab100, aes(x=MeanLength)) + 
  geom_histogram(aes(y = after_stat(count / sum(count))), color="lightblue", fill="black")+
  ylab("Proportion of samples") +
  xlab("Gene Length (in basepairs)") +
  xlim(0,15000) +
  ylim(0,0.20) +
  geom_vline(data=all_summ,mapping=aes(xintercept=MeanLength), color="orange")+
  geom_vline(data=my_summ100,mapping=aes(xintercept=GrandMean), color="cornflowerblue")+
  ggtitle("Sampling distribution of mean gene length (n=5),100 reps") 

p1b<-ggplot(my_tab20000, aes(x=MeanLength)) + 
  geom_histogram(aes(y = after_stat(count / sum(count))), color="lightblue", fill="black")+
  ylab("Proportion of samples") +
  xlab("Gene Length (in basepairs)") +
  xlim(0,15000) +
  ylim(0,0.20) +
  geom_vline(data=all_summ,mapping=aes(xintercept=MeanLength), color="orange")+
  geom_vline(data=my_summ20000,mapping=aes(xintercept=GrandMean), color="cornflowerblue")+
  ggtitle("Sampling distribution of mean gene length (n=5),20,000 reps") 


p2<-ggplot(my_tab1000, aes(x=MeanLength)) + 
  geom_histogram(aes(y = after_stat(count / sum(count))), color="lightblue", fill="black")+
  ylab("Proportion of samples") +
  xlab("Gene Length (in basepairs)") +
  xlim(0,15000) +
  ylim(0,0.20) +
  geom_vline(data=all_summ, mapping=aes(xintercept=MeanLength), color="orange")+
  geom_vline(data=my_summ1000, mapping=aes(xintercept=GrandMean), color="cornflowerblue")+
  ggtitle("Sampling distribution of mean gene length (n=5),1000 reps") 

p3<-ggplot(my_tab2100, aes(x=MeanLength)) + 
  geom_histogram(aes(y = after_stat(count / sum(count))), color="lightblue", fill="black")+
  ylab("Proportion of samples") +
  xlab("Gene Length (in basepairs)") +
  xlim(0,15000) +
  ylim(0,0.20) +
  geom_vline(data=all_summ, mapping=aes(xintercept=MeanLength), color="orange")+
  geom_vline(data=my_summ2100, mapping=aes(xintercept=GrandMean), color="cornflowerblue")

p4<-ggplot(my_tab21000, aes(x=MeanLength)) + 
  geom_histogram(aes(y = after_stat(count / sum(count))), color="lightblue", fill="black")+
  ylab("Proportion of samples") +
  xlab("Gene Length (in basepairs)") +
  geom_vline(data=all_summ, mapping=aes(xintercept=MeanLength), color="orange")+
  geom_vline(data=my_summ21000, mapping=aes(xintercept=GrandMean), color="cornflowerblue")+
  ggtitle("Sampling distribution of mean gene length (n=100),1000 reps") 


library(cowplot)


rbind(my_summ10,my_summ100, my_summ1000) %>% mutate(Number_replicates=c(10,100,1000)) %>% grid.table
all_summ %>% grid.table

my_tab21000<-replicate(1000, sample_n(human_genes, 100), simplify = FALSE) %>%
  lapply(., my_func) %>% 
  bind_rows %>%
  mutate(replicate = 1:n())

my_summ21000<-my_tab21000 %>% summarise(GrandMean=mean(MeanLength), MeanSD=mean(SD_Length))


my_tab21000<-my_tab21000 %>% mutate(SEM=sd(MeanLength)/sqrt(100))


set.seed(1)
sample_n(human_genes,5) %>% summarise(mean=mean(size),sd=sd(size))

sample_n(human_genes,100) %>% summarise(mean=mean(size),sd=sd(size))

p6<-ggplot(my_tab21000, aes(x=SEM)) + 
  geom_histogram(color="lightblue", fill="black")+
  ylab("Proportion of samples") +
  xlab("Estimate of the SEM (n=100)") 
 # xlim(0,15000) +
 # ylim(0,0.20) +
 # geom_vline(data=all_summ, mapping=aes(xintercept=MeanLength), color="orange")
#  geom_vline(data=my_summ1000, mapping=aes(xintercept=GrandMean), color="cornflowerblue")+
  #ggtitle("Sampling distribution of mean gene length (n=100),1000 reps") 


my_tab20000<-my_tab20000 %>%mutate(SEM = SD_Length/sqrt(5))
my_tab21000<-my_tab21000 %>%mutate(SEM = SD_Length/sqrt(100))
p7<-ggplot(my_tab20000, aes(x=SEM)) + 
  geom_histogram(aes(y = after_stat(count / sum(count))), color="lightblue", fill="black")+
  ylab("Proportion of samples") +
  xlab("SEM (in basepairs)") +
  #geom_vline(data=all_summ, mapping=aes(xintercept=MeanLength), color="orange")+
 # geom_vline(data=my_summ20000, mapping=aes(xintercept=GrandMean), color="cornflowerblue")+
  ggtitle("Sampling distribution of mean gene length (n=5),20000 reps") 

my_tab20000
set.seed(1)
my_tab_2_20000<-replicate(20000, sample_n(human_genes, 100), simplify = FALSE) %>%
  lapply(., my_func) %>% 
  bind_rows %>%
  mutate(replicate = 1:n()) 


my_CI<-function(df, ci=0.95){
  meanx=mean(df$size)
  sdx=sd(df$size)
  ns<-nrow(df)
  Z<-ifelse(ci==0.95, 1.96, ifelse(ci==0.90, 1.6449, ifelse(ci==0.99, 2.5758, NA)))
  CI_x_upper=meanx + Z*(sdx/sqrt(ns))
  CI_x_lower=meanx - Z*(sdx/sqrt(ns))
  data.frame(Mean=meanx,CI_upper=CI_x_upper, CI_lower=CI_x_lower, Sample_size=ns, CI=ci)
  
}

set.seed(1)
ci_10<-bind_rows(replicate(10, sample_n(human_genes, 100, replace=T), simplify=F) %>%
  lapply(., my_CI) %>% 
  bind_rows %>%
  mutate(replicate = 1:n()),
  replicate(10, sample_n(human_genes, 100, replace=T), simplify=F) %>%
    lapply(., my_CI, ci=0.90) %>% 
    bind_rows %>%
    mutate(replicate = 1:n())
)
  
ci_10<-ci_10 %>% mutate(width=CI_upper-CI_lower)


set.seed(10)
ci_100<-bind_rows(replicate(100, sample_n(human_genes, 100, replace=T), simplify=F) %>%
                   lapply(., my_CI) %>% 
                   bind_rows %>%
                   mutate(replicate = 1:n()),
                 replicate(100, sample_n(human_genes, 100, replace=T), simplify=F) %>%
                   lapply(., my_CI, ci=0.90) %>% 
                   bind_rows %>%
                   mutate(replicate = 1:n())
)

ci_100<-ci_100 %>% mutate(width=CI_upper-CI_lower)
set.seed(1)
ci_100_50<-replicate(50, sample_n(human_genes, 100, replace=T), simplify=F) %>%
  lapply(., my_CI) %>% 
  bind_rows %>%
  mutate(replicate = 1:n()) 
ci_100_50<-ci_100_50 %>% mutate(width=CI_upper-CI_lower)
ci_500_50<-replicate(50, sample_n(human_genes, 500, replace=T), simplify=F) %>%
  lapply(., my_CI) %>% 
  bind_rows %>%
  mutate(replicate = 1:n()) 

##
set.seed(1234)
ci_5_20_100_1000<-bind_rows(replicate(10, sample_n(human_genes, 5, replace=T), simplify=F) %>%
                   lapply(., my_CI) %>% 
                   bind_rows %>%
                   mutate(replicate = 1:n()),
                 replicate(10, sample_n(human_genes, 20, replace=T), simplify=F) %>%
                   lapply(., my_CI, ci=0.95) %>% 
                   bind_rows %>%
                   mutate(replicate = 1:n()),
                 replicate(10, sample_n(human_genes, 100, replace=T), simplify=F) %>%
                   lapply(., my_CI, ci=0.95) %>% 
                   bind_rows %>%
                   mutate(replicate = 1:n()),
                 replicate(10, sample_n(human_genes, 1000, replace=T), simplify=F) %>%
                   lapply(., my_CI, ci=0.95) %>% 
                   bind_rows %>%
                   mutate(replicate = 1:n())
)

ci_5_20_100_1000<-ci_5_20_100_1000 %>% mutate(width=CI_upper-CI_lower)


ggplot(ci_5_20_100_1000,aes(y=factor(replicate), x=Mean)) +
  geom_point() +
  geom_pointrange(aes(xmin=CI_lower, xmax=CI_upper), width=0.1, size=0.2) + 
  facet_wrap(vars(Sample_size))+
  ylab("Replicate") + xlab("Gene Length (in basepairs)") +
  geom_vline(data=all_summ, mapping=aes(xintercept=MeanLength), color="orange") +
  ggtitle("95% CI estimated from sample of genes") +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank() 
  )
  

ggplot(ci_100_50,aes(y=factor(replicate), x=Mean)) +
  geom_point() +
  geom_pointrange(aes(xmin=CI_lower, xmax=CI_upper), width=0.1, size=0.2) + 
  ylab("Replicate") + xlab("Gene Length (in basepairs)") +
  geom_vline(data=all_summ, mapping=aes(xintercept=MeanLength), color="orange") +
  ggtitle("95% CI estimated from sample of genes (n=100)")+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank() 
  )


ggplot(ci_100_100,aes(y=factor(replicate), x=Mean)) +
  geom_point() +
  geom_pointrange(aes(xmin=CI_lower, xmax=CI_upper), width=0.1, size=0.2) + 
  ylab("Replicate") + xlab("Gene Length (in basepairs)") +
  geom_vline(data=all_summ, mapping=aes(xintercept=MeanLength), color="orange") +
  ggtitle("95% CI estimated from sample of genes (n=100)")+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank() 
  )


toplot<-rbind(ci_50, ci_100_50, ci_500_50)
onesamp<-sample_n(ci_50,5)
ggplot(onesamp,aes(y=factor(replicate), x=Mean)) +
  geom_point() +
  geom_pointrange(aes(xmin=CI_lower, xmax=CI_upper), width=0.1, size=0.2) + 
  ylab("Replicate") + xlab("Gene Length (in basepairs)") +
  geom_vline(data=all_summ, mapping=aes(xintercept=MeanLength), color="orange") +
  ggtitle("95% CI estimated from sample of genes (n=100)")+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        element_blank())
ggplot(onesamp, aes(x=Mean))+geom_histogram() + 
  xlab("Mean gene length")+
  facet_wrap(vars(Sample_size),nrow = 3)


#
library(dplyr)
library(rsample)
library(purrr)
library(boot)
human_genes2<-select(human_genes,size)
my_CI2<-function(data, i){
  data<-as.data.frame(data[i,])
  ci=0.95
  meanx=mean(data[i,]$size)
  sdx=sd(data[i,]$size)
  ns<-nrow(data[i,])
  #Z<-ifelse(ci==0.95, 1.96, ifelse(ci==0.90, 1.6449, NA))
  Z=1.96
  CI_x_upper=meanx + Z*(sdx/sqrt(ns))
  CI_x_lower=meanx - Z*(sdx/sqrt(ns))
  (c(meanx,sdx,CI_x_upper, CI_x_lower))
  
}

#set.seed(2014)
res10<-boot(human_genes[sample(nrow(human_genes), size=100),], my_CI2,R=10)
res100<-boot(human_genes[sample(nrow(human_genes), size=100),], my_CI2,R=100)
res1000<-boot(human_genes[sample(nrow(human_genes), size=100),], my_CI2,R=1000)
res10000<-boot(human_genes[sample(nrow(human_genes), size=100),], my_CI2,R=10000)

set.seed(123)
res<-replicate(3, sample_n(select(human_genes,size), 100, replace=T), simplify=F)
res[[1]]<-res[[1]] %>% mutate(Rep=1)
res[[2]]<-res[[2]] %>% mutate(Rep=2)
res[[3]]<-res[[3]] %>% mutate(Rep=3)
res<-do.call(rbind, res)

library(plotrix)

my_CI<-function(x=size){
  n=length(x)
  s=sd(x)
  xbar=mean(x)
  margin <- qt(0.975,df=n-1)*s/sqrt(n)
  low <- xbar - margin
  high <- xbar + margin
  
  c(low, high)
}
  
res2<-res %>% 
  group_by(Rep) %>% 
  summarise(Mean=mean(size),SEM=std.error(size), SD=sd(size),CI=my_CI(size)[2]-my_CI(size)[1]) %>%
  ungroup

res3<-data.frame(Rep=1, 
                 size=3387, 
                 Stat=c("Data","SEM","SD","CI_95"), 
                 lower=c(NA,3387-c(255,2551,(1012/2))),
                 upper=c(NA,3387+c(255,2551,(1012/2)))
                 )
res3$Stat<-factor(res3$Stat, ordered=T, levels=c("Data","SEM","CI_95","SD"))
res4<-as.data.frame(filter(res, Rep==1) %>% 
                      mutate(Stat="Data",lower=3397, upper=3387)) %>%
                      select(Rep, size, Stat, lower, upper)
p<-ggplot(NULL, aes(x=Stat,y=size)) +
  geom_errorbar(data=res3,aes(ymin=lower, ymax=upper),width=0.1)+
  geom_jitter(data=res4, aes(x=Stat, y=size), col="cornflowerblue") +
  theme_bw()

  
