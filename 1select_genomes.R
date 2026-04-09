#!/usr/bin/env Rscript
#===============================================================================
# GTDB release 226 genome quality filtering and representative genome selection
# GTDB release 226 基因组质量过滤与代表基因组挑选
#
# Workflow:
# 1. Load bacterial and archaeal metadata
# 2. CheckM2 quality filtering (completeness >= 50%, contamination < 10%)
# 3. Archaeal superphylum classification
# 4. GUNC chimera removal
# 5. Phylum-level filtering:
#    a. Each phylum must have at least 2 distinct classes
#    b. Each phylum must have at least 25 genomes
# 6. Select representative genomes at order level and genus level
# 7. Output filtered genome lists
#
# 流程:
# 1. 加载细菌和古菌元数据
# 2. CheckM2质量过滤 (completeness >= 50%, contamination < 10%)
# 3. 古菌超门分类
# 4. GUNC嵌合体去除
# 5. 门水平过滤:
#    a. 每个门至少包含2个不同的纲
#    b. 每个门至少包含25个基因组
# 6. 挑选目水平和属水平代表基因组
# 7. 输出过滤后的基因组列表
#
# Representative genome selection strategy:
#
# Order-level representatives:
# - One representative genome per order
# - Application 1: Machine learning model training (SAA frequency prediction)
# - Application 2: Tree of Life construction (phylogenetic tree)
# - Application 3: Molecular clock and reconciliation analysis
#
# Genus-level representatives:
# - One representative genome per genus
# - Application: HMM profile construction for sulfonate metabolism
#
# Quality score for representative selection:
#   score = completeness - 5 * contamination
#
# 代表基因组挑选策略:
#
# 目水平代表基因组:
# - 每个目挑选一个代表基因组
# - 应用1: 机器学习模型训练 (SAA频率预测)
# - 应用2: 生命树构建 (系统发育树)
# - 应用3: 分子钟和调和分析
#
# 属水平代表基因组:
# - 每个属挑选一个代表基因组
# - 应用: 磺酸盐代谢的HMM模型构建
#
# 代表基因组挑选的质量评分:
#   评分 = 完整度 - 5 × 污染度
#===============================================================================

###质控后的基因组每个属挑选一个作为基础数据集
setwd("D:/databases/gtdb/226")
library("tidyverse")

bac <- data.table::fread("bac120_metadata_r226.tsv")
bac <- bac[which(bac$gtdb_representative == "t"),]
bac <- bac[,c(1,3,4,20,59,62,82,83)]
bac <- bac[which(bac$checkm2_completeness >= 50 & bac$checkm2_contamination < 10),]
bac <- separate(bac,col = "gtdb_taxonomy",into = c("d","p","c","o","f","g","s"),sep = ";")

ar <- data.table::fread("ar53_metadata_r226.tsv")
ar <- ar[which(ar$gtdb_representative == "t"),]
ar <- ar[,c(1,3,4,20,59,62,82,83)]
ar <- ar[which(ar$checkm2_completeness >= 50 & ar$checkm2_contamination < 10),]
ar <- separate(ar,col = "gtdb_taxonomy",into = c("d","p","c","o","f","g","s"),sep = ";")
ar$d[grep("Thermoproteota",ar$p)] <- "TACK"
ar$d[grep("Korarchaeota",ar$p)] <- "TACK"
ar$d[grep("Halobacteriota",ar$p)] <- "Euryarchaeota"
ar$d[grep("Thermoplasmatota",ar$p)] <- "Euryarchaeota"
ar$d[grep("Methanobacteriota",ar$p)] <- "Euryarchaeota"
ar$d[grep("Hydrothermarchaeota",ar$p)] <- "Euryarchaeota"
ar$d[grep("Nanobdellota",ar$p)] <- "DPANN"
ar$d[grep("Micrarchaeota",ar$p)] <- "DPANN"
ar$d[grep("Aenigmatarchaeota",ar$p)] <- "DPANN"
ar$d[grep("Iainarchaeota",ar$p)] <- "DPANN"
ar$d[grep("Nanohalarchaeota",ar$p)] <- "DPANN"
ar$d[grep("Altiarchaeota",ar$p)] <- "DPANN"
ar$d[grep("Asgardarchaeota",ar$p)] <- "Asgard"

superphylum <- c("DPANN", "Asgard", "TACK", "Euryarchaeota")
ar <- ar %>% filter(d %in% superphylum)
table(ar$d)
all <- rbind(bac,ar)

gunc <- data.table::fread("GUNC_226.tsv")
gunc$genome <- gsub("_protein","",gunc$genome)
all <- merge(all,gunc[,c(1,13)],by.x = "accession",by.y = "genome",all.x = T)
all <- all[all$pass.GUNC == "TRUE",]

##规范门名
all[grep("p__Bacteroidota_",all$p),5] <- "p__Bacteroidota"
all[grep("p__Bdellovibrionota_",all$p),5] <- "p__Bdellovibrionota"
all[grep("p__Desulfobacterota_",all$p),5] <- "p__Desulfobacterota"
all[grep("p__Methanobacteriota_",all$p),5] <- "p__Methanobacteriota"
all[grep("p__Myxococcota_",all$p),5] <- "p__Myxococcota"
all[grep("p__Nitrospirota_",all$p),5] <- "p__Nitrospirota"
t <- all %>% group_by(p) %>%
  summarise(n_class = n_distinct(c)) %>%
  arrange(desc(n_class))
all <- all %>%  group_by(p) %>% mutate(n_order = n_distinct(o)) %>%  
  ungroup() %>% filter(n_order >= 2) %>%  select(-n_order)
all <- all %>% group_by(p) %>% mutate(p_count = n()) %>% ungroup() %>% 
  filter(p_count >= 25) %>% select(-p_count)
t1 <- table(all$p) %>% as.data.frame()

all$d <- gsub("d__","",all$d)
all$p <- gsub("p__","",all$p)
all$c <- gsub("c__","",all$c)
all$o <- gsub("o__","",all$o)
all$f <- gsub("f__","",all$f)
all$g <- gsub("g__","",all$g)
all$s <- gsub("s__","",all$s)

# 每个属保留一个最佳:完整度-污染度*5 最大
genome <- all %>%
  group_by(g) %>%
  slice_max(checkm2_completeness- 5 *checkm2_contamination, n = 1, with_ties = FALSE) %>%
  ungroup()

# 每个目保留一个最佳:完整度-污染度*5 最大
genomes <- all %>%
  group_by(o) %>%
  slice_max(checkm2_completeness- 5 *checkm2_contamination, n = 1, with_ties = FALSE) %>%
  ungroup()

write.table(genome,"reanalysis/genome_genus.tsv",col.names = T,row.names = F,quote = F,sep = "\t")
write.table(genomes,"reanalysis/genome_order.tsv",col.names = T,row.names = F,quote = F,sep = "\t")