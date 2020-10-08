#scatter plot for process efficiency
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(ggsci)

plots = function(rawdata, title, outdir){
  
  
  data = merge(rawdata[rawdata$strand == '+',c('count','event')], rawdata[rawdata$strand == '-',], by.x = 'event', by.y = 'event', all = T)
  
  data$eff = data$count.x/data$count.y
  data$eff[(!is.na(data$count.y)) &(is.na(data$count.x))] = -100
  subdata = data[(!is.na(data$eff)) & (nchar(as.character(data$name))>0) &(data$pattern == 'TRSL-can'),]
  subdata = subdata[order(subdata$name),]
  p=ggplot(subdata) + geom_point(aes(count.y, eff, color = name)) + geom_text_repel(aes(count.y, eff, label = name)) +
    xlab('Anti-Sense Strand Count') + ylab("Transcription Efficiency") + labs(title = title)
  ggsave(file.path(outdir, paste0(title, '_TRSL_scatterplot.pdf')),p, units = 'in', width = 6,height = 4)
  ggplot(subdata) + geom_col(aes(name,eff, fill =name)) + xlab('Gene Body') + ylab('Transcription Efficiency') + labs(title = title)
  ggsave(file.path(outdir, paste0(title, '_TRSL_barplot.pdf')), units = 'in', width = 6,height = 4)
  ggplot(subdata) + geom_col(aes(eff,event, fill =name)) + ylab('Junction Events') + xlab('Transcription Efficiency') + labs(title = title)
  ggsave(file.path(outdir, paste0(title, '_TRSL_event_barplot.pdf')), units = 'in', width = 6,height = 8)
  
  return(subdata[,c('event', 'sample', 'eff','name')])
  
  
  # subdata = data[(!is.na(data$eff)) ,]
  # p=ggplot(subdata) + geom_point(aes(log2(count.y), eff, color = pattern)) +
  #   xlab('Anti-Sense Strand Count') + ylab("Transcription Efficiency") + labs(title = title)
  # ggsave(file.path(outdir, paste0(title, '_pattern_scatterplot.pdf')),p, units = 'in', width = 6,height = 4)
  # ggplot(subdata) + geom_col(aes(eff,pattern, fill =pattern)) + ylab('Pattern') + xlab('Transcription Efficiency') + labs(title = title)
  # ggsave(file.path(outdir, paste0(title, '_pattern_barplot.pdf')), units = 'in', width = 8,height = 4)
  # 
}

#setwd('G:/01.projects/COVID-19/Virus_Junction/Ribozero/jumps_bystrand/summary')
#outdir = 'G:/01.projects/COVID-19/Virus_Junction/Ribozero/jumps_bystrand/summary'
setwd('G:/01.projects/COVID-19/Virus_Junction/Ribozero/jumps_bystrand/dedup/summary')
outdir = 'G:/01.projects/COVID-19/Virus_Junction/Ribozero/jumps_bystrand/dedup/summary/Vero'
dir.create(outdir)


#HB6 = read.table('HB6_summary.csv',sep = ',', header = T)
#HC6 = read.table('HC6_summary.csv',sep = ',', header = T)
VB10 = read.table('VB10_summary.csv', sep = ',', header = T)
VC10 = read.table('VC10_summary.csv', sep = ',', header = T)

VB12 = read.table('VB12_summary.csv', sep = ',', header = T)
VC12 = read.table('VC12_summary.csv', sep = ',', header = T)


#hb = plots(HB6[HB6$count > 5,], 'HB6', outdir)
#hc = plots(HC6[HC6$count > 5,], 'HC6', outdir)
# vb10 = plots(VB10[VB10$count > 5,], 'VB10', outdir)
# vc10 = plots(VC10[VC10$count > 5,], 'VC10', outdir)
# vb12 = plots(VB12[VB12$count > 5,], 'VB12', outdir)
# vc12 = plots(VC12[VC12$count > 5,], 'VC12', outdir)
vb10 = plots(VB10, 'VB10', outdir)
vc10 = plots(VC10, 'VC10', outdir)
vb12 = plots(VB12, 'VB12', outdir)
vc12 = plots(VC12, 'VC12', outdir)


#all = rbind(hb,hc,vb10,vc10,vb12,vc12)
all = rbind(vb10,vc10,vb12,vc12)
all$group = 'NA'
#all$group[(all$sample == 'HB6') | (all$sample == 'HC6')] = 'H6'
all$group[(all$sample == 'VB10') | (all$sample == 'VC10')] = 'V10'
all$group[(all$sample == 'VB12') | (all$sample == 'VC12')] = 'V12'
write.csv(all, file.path(outdir, 'Vero_efficiency.csv'))

############################################################################
bothevents1 = intersect(VB10[VB10$strand == '+','event'], VB10[VB10$strand == '-','event'])
bothevents2 = intersect(VC10[VC10$strand == '+','event'], VC10[VC10$strand == '-','event'])
bothevents_v10 = intersect(bothevents1, bothevents2)

bothevents1 = intersect(VB12[VB12$strand == '+','event'], VB12[VB12$strand == '-','event'])
bothevents2 = intersect(VC12[VC12$strand == '+','event'], VC12[VC12$strand == '-','event'])
bothevents_v12 = intersect(bothevents1, bothevents2)

bothevents = intersect(bothevents_v10,bothevents_v12)

bothevents = bothevents2
subVB10 = VC10[VC10$event %in% bothevents,]
ggplot() + geom_point(aes(subVB10[subVB10$strand == '+','count'], subVB10[subVB10$strand == '-', 'count']))
ggplot() + geom_point(aes(log10(subVB10[subVB10$strand == '+','count']), log10(subVB10[subVB10$strand == '-', 'count'])))

VC10$group = 'Sequenced'
psudocounts = VC10[!(VC10$event %in% bothevents) & (VC10$strand == '+'), ]
psudocounts['count'] = 0.5
psudocounts['strand'] = '-'
psudocounts$group = '(+)-Specific'

psudocounts2 = VC10[!(VC10$event %in% bothevents) & (VC10$strand == '-'), ]
psudocounts2['count'] = 0.25
psudocounts2['strand'] = '+'
psudocounts2$group = '(-)-Specific'

subVB10 = rbind(VC10, psudocounts, psudocounts2)
subVB10 = subVB10[order(subVB10$event),]
ggplot() + geom_bin2d(aes(subVB10[subVB10$strand == '+','count'], subVB10[subVB10$strand == '-', 'count']), bins =100)+
  scale_fill_continuous(type = "viridis") + theme_bw() + labs(x = 'Sense Strand Counts', y = 'Anti-Sense Strand Counts') +
  theme(axis.text.x=element_text(size=12), axis.text.y=element_text(size=12),
        plot.title = element_text(lineheight=.8, face="bold", hjust=0.1, size =16),
        panel.background = element_rect(fill="white", colour='gray'),
        panel.grid.major = element_line(linetype = 'dotted', colour = "gray"),
        panel.border=element_rect(fill='transparent',
                                  color='black'))

ggsave(file.path(outdir, 'Density_counts_VC10.pdf'), units = 'in', width = 6, height=5)
ggplot() + geom_point(aes(log2(subVB10[subVB10$strand == '+','count']), log2(subVB10[subVB10$strand == '-', 'count']))) + 
  labs(x = 'Sense Strand Counts', y = 'Anti-Sense Strand Counts') + 
  theme(axis.text.x=element_text(size=12), axis.text.y=element_text(size=12),
        plot.title = element_text(lineheight=.8, face="bold", hjust=0.1, size =16),
        panel.background = element_rect(fill="white", colour='gray'),
        panel.grid.major = element_line(linetype = 'dotted', colour = "gray"),
        panel.border=element_rect(fill='transparent',
                                  color='black')) + geom_abline(slope = 1)
ggsave(file.path(outdir, 'logcounts_VC10.pdf'), units = 'in', width = 6, height=5)

ggplot() + geom_bin2d(aes(log2(subVB10[subVB10$strand == '+','count']), log2(subVB10[subVB10$strand == '-', 'count'])), bins = 60) +
  scale_fill_continuous(type = "viridis") + theme_bw() + labs(x = 'Log2(Sense Strand Counts)', y = 'Log2(Anti-Sense Strand Counts)') + 
  theme(axis.text.x=element_text(size=12), axis.text.y=element_text(size=12),
        plot.title = element_text(lineheight=.8, face="bold", hjust=0.1, size =16),
        panel.background = element_rect(fill="white", colour='gray'),
        panel.grid.major = element_line(linetype = 'dotted', colour = "gray"),
        panel.border=element_rect(fill='transparent',
                                  color='black')) + geom_abline(slope = 1)
ggsave(file.path(outdir, 'Density_log2counts_VC10.pdf'), units = 'in', width = 6, height=5)

#########################################################################
############################################################################




aggplots = function(rawdata, title, outdir){
  
  
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
  
  
  # data = merge(rawdata[rawdata$strand == '+',c('event', 'count')], rawdata[rawdata$strand == '-',], by.x = 'event', by.y = 'event', all = T)
  # 
  # data$eff = data$count.x*1e6/(data$count.y*1e6)
  # data$eff[(!is.na(data$count.y)) &(is.na(data$count.x))] = -100
  # 
  # 
  # subdata = data[(!is.na(data$eff)) ,]
  write.csv(subdata, file.path(outdir, paste0(title, '_efficiencies.csv')), row.names = F, quote = F)
  p=ggplot(subdata) + geom_point(aes(count.y, eff, color = name)) + 
    xlab('Percentage') + ylab("Transcription Efficiency") + labs(title = title)
  ggsave(file.path(outdir, paste0(title, '_aggname_count_scatterplot.pdf')),p, units = 'in', width = 5,height = 3)
  ggplot(subdata) + geom_col(aes(name, eff,fill =name)) + theme_bw()+
    theme_bw()+scale_fill_npg()+xlab('Pattern') + ylab('Transcription Efficiency') + labs(title = title)
  ggsave(file.path(outdir, paste0(title, '_aggname_count_barplot.pdf')), units = 'in', width = 5,height = 3)
  
  return(subdata)
}



# HB6 = read.csv('HB6_pct.csv',sep = ',', header = T)
# HC6 = read.csv('HC6_pct.csv',sep = ',', header = T)
# VB10 = read.csv('VB10_pct.csv', sep = ',', header = T)
# VC10 = read.csv('VC10_pct.csv', sep = ',', header = T)
# 
# VB12 = read.csv('VB12_pct.csv', sep = ',', header = T)
# VC12 = read.csv('VC12_pct.csv', sep = ',', header = T)
# 
VB10 = read.table('VB10_summary.csv', sep = ',', header = T)
VB10$count = VB10$count * 1e6/sum(VB10$count)
VC10 = read.table('VC10_summary.csv', sep = ',', header = T)
VC10$count = VC10$count * 1e6/sum(VC10$count)
VB12 = read.table('VB12_summary.csv', sep = ',', header = T)
VB12$count = VB12$count * 1e6/sum(VB12$count)
VC12 = read.table('VC12_summary.csv', sep = ',', header = T)
VC12$count = VC12$count * 1e6/sum(VC12$count)

bothevents = intersect(bothevents_v10,bothevents_v12)
vb10 = aggplots(VB10[VB10$event %in% bothevents,], 'VB10', outdir)
vc10 = aggplots(VC10[VC10$event %in% bothevents,], 'VC10', outdir)
vb12 = aggplots(VB12[VB12$event %in% bothevents,], 'VB12', outdir)
vc12 = aggplots(VC12[VC12$event %in% bothevents,], 'VC12', outdir)

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

#p.values <- sapply(split(data, data$sgRNAs), function(x){t.test(x$eff[x$Group == 'Virus'], x$eff[x$Group == 'Virus + RDV'], alternative = 'less')$p.value})
#p.values <- sapply(split(data, data$sgRNAs), function(x){t.test(eff~Group, x, alternative = 'less')$p.value})

#labels <- symnum(p.values, corr = FALSE, cutpoints = c(0,  .001,.01,.05, 1), symbols = c("***","**","*","n.s."))
#y.values <- sapply(split(data, data$sgRNAs), function(x){max(sapply(split(x, x$Group), function(xx){barplot(x$eff, plot=F)$stats[5, ]}))})+2

# df_stat %>% ggplot() +
#   geom_bar(aes(x=sgRNAs, y=mean, fill=Group),
#            stat="identity", position=dodge) +
#   geom_errorbar(aes(x=sgRNAs, ymin=mean-sd, ymax=mean+sd,  group = Group),
#                 stat="identity", position=dodge, width=.3)+
#   #scale_x_continuous(breaks = c(2,6,10,14,18,22,26,30)) +
#   labs(x = 'sgRNAs', y='Efficiency') +
#   theme_bw()+scale_fill_npg() +
#   geom_signif(annotations = ttest$pstar[1:2],  xmin = c(1,3), xmax = c(2,4))
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

