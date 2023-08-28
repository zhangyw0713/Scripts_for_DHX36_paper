
read.table("/lustre/zhangyw/myProject/Dhx36/new_RNAseq/00_deseq2/DHX36_wt_vs_ko_raw.tsv",sep="\t",header=T)->df1
-1*df1$log2FoldChange->df1$log2Fc_ko_wt
df1[which(df1$padj<0.05 & df1$log2Fc_ko_wt>=log2(1.5)),]->up
df1[which(df1$padj<0.05 & df1$log2Fc_ko_wt<=log2(1/1.5)),]->down

setwd("/lustre/zhangyw/myProject/Dhx36/new_RNAseq/chrRNAseq/03_dif")
read.table("chrRNAseq_wt_vs_ko_raw.tsv",sep="\t",header=T)->df
-1*df$log2FoldChange->df$log2Fc_ko_wt

merge(df1,df,by="row.names")->upc
upc[,c(1,8,15)]->upcd
colnames(upcd)<-c("gene","total","chr")
2^upcd$total->upcd$totalfc
2^upcd$chr->upcd$chrfc
upcd$fcfc=upcd$totalfc/upcd$chrfc
upcd[which(upcd$fcfc>1.5),]->testfc
testfc[,1]->post_up_gene
#str_split_fixed(testfc[,1],"\\.",n=2)[,1]->post_up_gene
length(post_up_gene)
upcd[which(upcd$fcfc<1/1.5),]->testfc
#str_split_fixed(testfc[,1],"\\.",n=2)[,1]->post_down_gene
testfc[,1]->post_down_gene
length(post_down_gene)

ggplot(data=upcd, aes(x=total, y=chr)) +
  geom_point(data=subset(upcd,upcd$fcfc <= 1.5 & upcd$fcfc >=1/1.5),color="gray",alpha=0.6) +
  geom_point(data=subset(upcd,upcd$fcfc > 1.5),color="#fc8d62",alpha=0.9) +
  geom_point(data=subset(upcd,upcd$fcfc < 1/1.5),color="#94C47D",alpha=0.9) +
  xlab("FC in whole cell (KO/WT)") + ylab("FC in chromatin fraction (KO/WT)") +   
  theme_classic(base_size=15)+   
  theme(panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))+ 
  theme(plot.title = element_text(size=15,hjust = 0.5),legend.position='none')+xlim(-4,4)+ylim(-4,4)+
  annotate(geom="text", x=2.5, y=-3, label="Post-transcriptionally upregulated genes",color="#fc8d62")+annotate(geom="text", x=-2, y=3, label="Post-transcriptionally downregulated genes",color="#94C47D")
 ggsave("fc_choose_ko_wt.pdf",width=8.3,height=8)
 
 