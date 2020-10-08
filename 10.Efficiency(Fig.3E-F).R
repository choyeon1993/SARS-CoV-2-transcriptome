#scatter plot for process efficiency
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(ggsci)

setwd('G:/01.projects/COVID-19/scripts/Datasets/RibozeroDedupJunctions/')
outdir = 'G:/01.projects/COVID-19/scripts/Outputs'
dir.create(outdir)
outdir = 'G:/01.projects/COVID-19/scripts/Outputs/Efficiency'
dir.create(outdir)

VB10j = read.table('VB10_summary.csv', sep = ',', header = T)
VC10j = read.table('VC10_summary.csv', sep = ',', header = T)

VB12j = read.table('VB12_summary.csv', sep = ',', header = T)
VC12j = read.table('VC12_summary.csv', sep = ',', header = T)

totalcov = read.table('StrandCoverages.csv', header = T, sep = ',')
rownames(totalcov) = totalcov$X
totalcov = totalcov[,-1]

VB10j$CPM_totalViral = VB10j$count/totalcov$Virus.Rep1[3]*1e6
VB10j$CPM_totalJunc = VB10j$count/sum(VB10j$count) * 1e6
VC10j$CPM_totalViral = VC10j$count/totalcov$Virus.Rep2[3]*1e6
VC10j$CPM_totalJunc = VC10j$count/sum(VC10j$count) * 1e6
VB12j$CPM_totalViral = VB12j$count/totalcov$Virus.RDV.Rep1[3]*1e6
VB12j$CPM_totalJunc = VB12j$count/sum(VB12j$count) * 1e6
VC12j$CPM_totalViral = VC12j$count/totalcov$Virus.RDV.Rep2[3]*1e6
VC12j$CPM_totalJunc = VC12j$count/sum(VC12j$count) * 1e6

bothevents1 = intersect(VB10j[VB10j$strand == '+','event'], VB10j[VB10j$strand == '-','event'])
bothevents2 = intersect(VC10j[VC10j$strand == '+','event'], VC10j[VC10j$strand == '-','event'])
bothevents_v10 = intersect(bothevents1, bothevents2)

bothevents1 = intersect(VB12j[VB12j$strand == '+','event'], VB12j[VB12j$strand == '-','event'])
bothevents2 = intersect(VC12j[VC12j$strand == '+','event'], VC12j[VC12j$strand == '-','event'])
bothevents_v12 = intersect(bothevents1, bothevents2)

bothevents = intersect(bothevents_v10,bothevents_v12)


# intersect(VB10j$event[1:10], VB12j$event[1:10])
# intersect(VB10j$event[1:10], bothevents1)
# intersect(VB12j$event[1:10], bothevents1)



aggplots = function(rawdata, title, outdir){
  
  rawdata$count = rawdata$CPM_totalViral
  data = merge(rawdata[rawdata$strand == '+',c( 'count','event')], rawdata[rawdata$strand == '-',], by.x = 'event', by.y = 'event', all = T)
  
  data$eff = data$count.x/data$count.y
  data$eff[(!is.na(data$count.y)) &(is.na(data$count.x))] = -100
  
  
  subdata = data[(!is.na(data$eff)) ,]
  p=ggplot(subdata) + geom_point(aes(count.y, eff, color = name)) + geom_text_repel(aes(count.y, eff, label = name)) +
    xlab('Anti-Sense Strand Count') + ylab("Transcription Efficiency") + labs(title = title)
  ggsave(file.path(outdir, paste0(title, '_aggname_scatterplot.pdf')),p, units = 'in', width = 6,height = 4)
  ggplot(subdata) + geom_col(aes(name, eff,fill =name)) + xlab('Pattern') + ylab('Transcription Efficiency')+
    theme_bw()+scale_fill_npg()+ labs(title = title)
  ggsave(file.path(outdir, paste0(title, '_aggname_barplot.pdf')), units = 'in', width = 5,height = 3)
  
    write.csv(subdata, file.path(outdir, paste0(title, '_efficiencies.csv')), row.names = F, quote = F)
  p=ggplot(subdata) + geom_point(aes(count.y, eff, color = name)) + 
    xlab('Percentage') + ylab("Transcription Efficiency") + labs(title = title)
  ggsave(file.path(outdir, paste0(title, '_aggname_count_scatterplot.pdf')),p, units = 'in', width = 5,height = 3)
  ggplot(subdata) + geom_col(aes(name, eff,fill =name)) + theme_bw()+
    theme_bw()+scale_fill_npg()+xlab('Pattern') + ylab('Transcription Efficiency') + labs(title = title)
  ggsave(file.path(outdir, paste0(title, '_aggname_count_barplot.pdf')), units = 'in', width = 5,height = 3)
  
  return(subdata)
}


subVB10 = VB10j[VB10j$event %in% bothevents,]

vb10 = aggplots(VB10j[VB10j$event %in% bothevents,], 'VB10', outdir)
vc10 = aggplots(VC10j[VC10j$event %in% bothevents,], 'VC10', outdir)
vb12 = aggplots(VB12j[VB12j$event %in% bothevents,], 'VB12', outdir)
vc12 = aggplots(VC12j[VC12j$event %in% bothevents,], 'VC12', outdir)
all = rbind(vb10,vc10,vb12,vc12)
all$Group = 'NA'
#all$group[(all$sample == 'HB6') | (all$sample == 'HC6')] = 'H6'
all$Group[(all$sample == 'VB10') | (all$sample == 'VC10')] = 'Virus'
all$Group[(all$sample == 'VB12') | (all$sample == 'VC12')] = 'Virus + RDV'
write.csv(all, file.path(outdir, 'all_aggPatterns_efficiency.csv'))


library(tidyr)
library(dplyr)
library(ggpubr)
data = all[,c('name','Group','eff')]
colnames(data)[1] = 'sgRNAs'
df_stat = tbl_df(data) %>%
  #gather(sample, eff, -name) %>%  # 将"宽数据"转化为"长数据"
  group_by(sgRNAs, Group) %>%         # 将数据分组
  summarise(mean=mean(eff, na.rm=T), sd=sd(eff, na.rm=T)) %>% # 计算每组数据的mean和sd
  ungroup()
str(df_stat)
dodge <- position_dodge(width=.9)
df_stat = df_stat[rank(df_stat$mean),]
df_stat$sgRNAs = factor(df_stat$sgRNAs, levels =  unique(df_stat$sgRNAs))
#compaired = list(c('Virus','Virus + RDV'))
#compaired = list(as.character(df_stat$sgRNAs))
sampletest = function(df){
  df = as.data.frame(df)
  return(t.test(df$eff[df$Group == 'Virus'], df$eff[df$Group == 'Virus + RDV'], alternative = 'less')$p.value)
}

# test = aggregate(data, by=list(data$sgRNAs), FUN = sampletest)
ttest = data.frame(sgRNAs = unique(data$sgRNAs), ttest = 0, pstar = 'ns')
for (gene in unique(data$sgRNAs)){
  p= sampletest(data[data$sgRNAs == gene,])
  ttest$ttest[ttest$sgRNAs == gene] = p
  if (p <= 0.05){
    ttest$pstar = '*'
  }
  if (p <= 0.01){
    ttest$pstar = '**'
  }
  if (p <= 0.001){
    ttest$pstar = '***'
  }
  
}

data$sgRNAs = factor(data$sgRNAs, levels = c('S','ORF3a', 'E','M','ORF6','ORF7a','ORF8','N'))
ggbarplot(data, x = "sgRNAs", y = "eff", add = "mean_se",
          color = "Group", #palette = "jco",
          position = position_dodge(0.8)) +
  scale_fill_npg()+
  stat_compare_means(aes(group = Group), , method = 't.test', method.args = list(alternative = 'greater'),
                     label = "p.signif", label.y = 100)+
  labs(x = 'sgRNAs', y='Efficiency') 
ggsave(file.path(outdir, 'EfficiencyOfGenebody_bar.pdf'), dpi = 600, units = 'in', width = 5,height = 3)
#aggall = aggregate(all[,c('name','sample','eff')],by = list(all$group,all$name), 'mean')

data = all[,c('name','Group','eff','count.y')]
colnames(data)[1] = 'sgRNAs'
colnames(data)[4] = 'NegCounts'
aggdata = aggregate(data[c('eff', 'NegCounts')],by = list(data$sgRNAs, data$Group),'mean')
colnames(aggdata)[1] = 'sgRNAs'
colnames(aggdata)[2] = 'Samples'
ggplot(aggdata,aes(NegCounts, eff,color = Samples)) + geom_point(aes(shape = Samples)) + labs(x = 'Negative Strand Counts (CPM)',y= 'Transcription Efficiency') + 
  geom_text_repel(label = aggdata$sgRNAs) + scale_fill_npg() + theme_bw()
ggsave(file.path(outdir, 'EfficiencyOfGenebody_scatter.pdf'), dpi=600, units = 'in',width = 4.5,height = 3)




##########################Efficiencies of TRS-L and TRSs #################
aggplots(VB10[VB10$event %in% bothevents_v10,], 'VB10_TRSL&TRSs.', outdir)
aggplots(VC10[VC10$event %in% bothevents_v10,], 'VC10_TRSL&TRSs.', outdir)


