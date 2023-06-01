########################################################
#             Prj: EMMAX Downstream analysis
#             Assignment: haplotype Heatmap
#             Author: Shawn Wang
#             Date: Tue Apr 18
#             Location: HENU
#########################################################
#

# Prepare  ----------------------------------------------------------
TEST = "FALSE"
options(stringsAsFactors = F)
options(warn = -1)
suppressMessages(if (!require('getopt')) install.packages('getopt'))
suppressMessages(if (!require('crayon')) install.packages('crayon'))
# Args ------------------------------------------------------------------
##> crayon
msg_yes = green$bold$italic
msg_no = red$bold$italic
msg_run = blue$bold$italic$underline
msg_warning = yellow$bold$italic
message(msg_yes(
  paste0("\n#=========================================================#\n",
         "#        Prj: EMMAX-mate ver 0.0.1\n",
         "#        Assignment: Haplotype Heatmap\n",
         "#        Author: Shawn Wang <shawnwang2016@126.com>\n",
         "#        Date: Wed 19 Apr, 2023\n",
         "#=========================================================#\n\n"
  )
))
##> getopt
command=matrix(c(
  'help', 'h', 0, 'logic', 'help information',
  'genofile', 'g', 1, 'character', 'Emmax-down file: xxx.gwas. ',
  'famfile', 'f', 1, 'character','Plink file: xxx.fam',
  'regions', 'r', 2, 'character','Key chr regions: chr:start-end',
  'phenotype', 'p', 2, 'character','phenotype file',
  'tree', 't', 2, 'character', 'NJ tree, .nwk file',
  'colorkey', 'c', 2, 'character', 'color of snp type. default is: Missing:Heterozygous:Minor:Major = grey:blue:orange:cyan',
  'output', 'o', 2, 'character', 'output prefix'
),byrow = T, ncol = 5)
args = getopt(command)

##> help information
if (!is.null(args$help)) {
  message(msg_yes(paste(getopt(command, usage = T), "\n")))
  q(status=1)
}

##> default value
if (is.null(args$genofile)){
  message(msg_no("Need .gwas file! please read help carefully!"))
  message(msg_yes(paste(getopt(command, usage = T), "\n")))
  q(status = 1)
}

if (is.null(args$famfile)){
  message(msg_no("Need .fam file! please read help carefully!"))
  message(msg_yes(paste(getopt(command, usage = T), "\n")))
  q(status = 1)
}

if (is.null(args$tree)){
  message(msg_no("Need .nwk file! please read help carefully!"))
  message(msg_yes(paste(getopt(command, usage = T), "\n")))
  q(status = 1)
}

if (is.null(args$regions)){
  message(msg_no("Need chr region! in format: Chr:Start-End,  please read help carefully!"))
  message(msg_yes(paste(getopt(command, usage = T), "\n")))
  q(status = 1)
}

if (is.null(args$phenotype)){
  message(msg_no("Need phenotype file!"))
  message(msg_yes(paste(getopt(command, usage = T), "\n")))
  q(status = 1)
}

if (is.null(args$colorkey)){
  message(msg_warning("The colorkey is not set, and the current time is used by default"))
  args$colorkey = c("grey","blue","orange","cyan")
}

if (is.null(args$output)){
  message(msg_warning("The output file name is not set, and the current time is used by default"))
  args$output = paste0('haplotype.',gsub(pattern = ' ',replacement = '_',Sys.time()))
}
# library for data cleaning and visualize ---------------------------------
suppressMessages(if (!require('tidyverse')) install.packages('tidyverse'))
suppressMessages(if (!require('devtools')) install.packages('devtools'))
suppressMessages(if (!require("treeio")) devtools::install_github("YuLab-SMU/treeio"))
suppressMessages(if (!require("ggtree")) devtools::install_github("YuLab-SMU/ggtree"))
suppressMessages(if (!require("ComplexHeatmap")) devtools::install_github("jokergoo/ComplexHeatmap"))
suppressMessages(if (!require('ape')) install.packages('ape'))
suppressMessages(if (!require('phylogram')) install.packages('phylogram'))
suppressMessages(if (!require('writexl')) install.packages('writexl'))
suppressMessages(if (!require('magick')) install.packages('magick'))
# test --------------------------------------------------------------------
if(TEST == "TRUE") {
  genopath = "/Users/shawn/GeekSpace/棉花类黄酮课题/21.Gh_383_GWAS/02.Progress/01.Raw_genotype/average/petal/work_dir/Morin_Flower_miss_result/LDFilebyChr/A11_Gh_383_miss.geno"
  fampath = "/Users/shawn/GeekSpace/棉花类黄酮课题/21.Gh_383_GWAS/01.Data/01.raw_g/Gh_383.maf0.05.int0.8.fam"
  treepath = "/Users/shawn/GeekSpace/棉花类黄酮课题/21.Gh_383_GWAS/01.Data/01.raw_g/Morin_A11_91k.nwk"
  phenotype = "/Users/shawn/GeekSpace/棉花类黄酮课题/21.Gh_383_GWAS/01.Data/01.raw_g/phenotype.txt"
  regions = "A11:91108756-93315562"
  output = "test"
  colorkey = 'grey:blue:orange:cyan'
  } else {
  genopath = args$genofile
  fampath = args$famfile
  treepath = args$tree
  phenotypepath = args$phenotype
  output = args$output
  regions = args$regions
  colorkey = args$colorkey
  }


# import file and set variables -------------------------------------------

message(msg_run(paste0(
  "Start running ..."
)))

message(msg_run("Step1. Import data ..."))
message(msg_yes("\n-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-\n"))
genofile <- read.delim(genopath,header = F)

sample_order <- read.delim(fampath,header = F)

tree <- read.newick(treepath)

phenotype <- read.delim(phenotypepath,header = T)

colnames(phenotype)[1] = "Sample"

regions2 = str_split(string = regions,pattern = ":",n = 2,simplify = T)[,2]

start = str_split(string = regions2,pattern = "-",n = 2,simplify = T)[,1] %>% as.numeric()

end = str_split(string = regions2,pattern = "-",n = 2,simplify = T)[,2] %>% as.numeric()

colorkey = str_split(colorkey,":",4,T) %>% as.character()
# data transformation -----------------------------------------------------

message(msg_run("Step2. Transformed genotype data to numeric format ..."))
message(msg_yes("\n-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-\n"))

##> add sample name and convert to long table
geno_clean <-
  genofile %>% 
  dplyr::filter(
    V2 >= start & V2 <= end ## filter snps
  ) %>% tidyr::separate(col = V3,into = c(sample_order$V2),sep = ' ',convert = T) %>% 
  pivot_longer(-c(V1,V2),names_to = 'sample',values_to = 'geno') %>% 
  set_names('chr','pos','sample','geno')%>% 
  mutate(
    tag1 = case_when(
      geno == "N" ~ "Missing",
      geno != "N" & geno != "A" & geno != "T" & geno != "G" & geno != "C" ~ "Heterozygous",
      geno == "A" | geno == "T" | geno == "G" | geno == "C" ~ "Homogeneous"
    )
  ) 
##> find major and minor allel.

geno_final_long <-
  geno_clean %>% 
  filter(tag1 == "Homogeneous") %>% 
  group_by(chr,pos,tag1,geno) %>% 
  summarise(n = n()) %>% 
  group_by(chr,pos) %>% 
  mutate(
    tag2 = case_when(
      n ==  max(n) ~ 'Major',
      n == min(n) ~ 'Minor'
    )
  ) %>% 
  ungroup() %>% 
  select(-n) %>% 
  right_join(geno_clean) %>% 
  mutate(
    tag = case_when(
      is.na(tag2) ~ tag1,
      TRUE ~ tag2
    ),
    tagNum = case_when(
      tag == 'Missing' ~ 0,
      tag == 'Heterozygous' ~ 1,
      tag == 'Minor' ~ 2,
      tag == 'Major' ~ 3
    )
  ) %>% 
  select(-tag1,tag2)

##> convert numeric genodata and tag genodata

geno_num <- 
  geno_final_long %>%
  mutate(SNP = paste0(chr,':',pos)) %>% 
  select(SNP,sample,tagNum) %>% 
  pivot_wider(names_from = sample, values_from = tagNum)

##> convert order plhlo tree and convert to dendrogram

tree_laddered <- ladderize(tree)
tree_dend <- as.dendrogram.phylo(tree_laddered)
tree_clado <- as.cladogram(tree_dend)

# ComplexHeatmap ----------------------------------------------------------

message(msg_run("Step3. Ordered Haplotype heatmap by laddered NJ tree ..."))
message(msg_yes("\n-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-\n"))
##> ordered numeric table by original tree.
tbl_ht <- geno_num %>% 
  column_to_rownames("SNP") %>% 
  select(tree$tip.label) %>% 
  t()
message(msg_run("Step4. Haplotype heatmap using all SNPs"))
message(msg_yes("\n-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-\n"))
##> draw heatmap
ht_num_raster <-
  Heatmap(
    matrix = tbl_ht,
    col = circlize::colorRamp2(breaks = c(0,1,2,3),colors = colorkey),
    show_row_names = F,show_column_names = F,
    cluster_columns = F,cluster_rows = tree_clado,use_raster = T,
    heatmap_legend_param = list(
      color_bar = "discrete",
      labels = c("Missing","Heterozygous","Minor","Major"),
      title = "Genotype"
    ),
    column_title = paste0("SNP\n(Chr:",paste(regions,collapse = "_"),")"),
    row_title = paste0("Sample\n(n = ",nrow(tbl_ht),")"),border = T,border_gp = gpar(size = 1,color = 'black')
  )

##> export data
pdf(file = paste0(output,"_raster.pdf"),width = 10,height = 10)
ht_out <- draw(ht_num_raster)
dev.off()

##> high res plot
ht_num <-
  Heatmap(
    matrix = tbl_ht,
    col = circlize::colorRamp2(breaks = c(0,1,2,3),colors = colorkey),
    show_row_names = F,show_column_names = F,
    cluster_columns = F,cluster_rows = tree_clado,use_raster = F,
    heatmap_legend_param = list(
      color_bar = "discrete",
      labels = c("Missing","Heterozygous","Minor","Major"),
      title = "Genotype"
    ),
    column_title = paste0("SNP\n(Chr:",paste(regions,collapse = "_"),")"),
    row_title = paste0("Sample\n(n = ",nrow(tbl_ht),")"),border = T,border_gp = gpar(size = 1,color = 'black')
  )
pdf(file = paste0(output,"_high_res.pdf"),width = 10,height = 10)
ht_out <- draw(ht_num)
dev.off()

message(msg_run("Step5. Haplotype heatmap using top 1000 (low missing rate) SNPs ..."))
message(msg_yes("\n-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-\n"))

##> thin genotype data, slice the top 1000 snps with lower missing rate.
tmp_tbl =   tbl_ht %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("SNP")
tbl_ht_order_thin_1000 <- 
 tmp_tbl%>% 
  mutate(
    zerocount = rowSums(.[,-1] != 0)
  ) %>% 
  arrange(zerocount) %>% 
  slice(1:1000)%>% 
  select(-zerocount) %>% 
  right_join(tmp_tbl %>% select("SNP"),by="SNP") %>% 
  drop_na() %>% 
  distinct() %>% 
  column_to_rownames("SNP") %>% 
  t()

ht_num_thin_high_res <-
  Heatmap(
    matrix = tbl_ht_order_thin_1000,
    col = circlize::colorRamp2(breaks = c(0,1,2,3),colors = colorkey),
    show_row_names = F,show_column_names = F,
    cluster_columns = F,cluster_rows = tree_clado,use_raster = F,
    heatmap_legend_param = list(
      color_bar = "discrete",
      labels = c("Missing","Heterozygous","Minor","Major"),
      title = "Genotype"
    ),
    column_title = paste0("SNP\n(Chr:",paste(regions,collapse = "_"),")"),
    row_title = paste0("Sample\n(n = ",nrow(tbl_ht),")"),border = T,border_gp = gpar(size = 1,color = 'black')
  )

pdf(file = paste0(output,"_thin_high_res.pdf"),width = 10,height = 10)
ht_out_thin =draw(ht_num_thin_high_res)
dev.off()

##> arrange with hclust
ht_hclust <- 
  Heatmap(
    matrix = tbl_ht_order_thin_1000,
    col = circlize::colorRamp2(breaks = c(0,1,2,3),colors = colorkey),
    show_row_names = F,show_column_names = F,
    cluster_columns = F,cluster_rows = T,use_raster = F,
    heatmap_legend_param = list(
      color_bar = "discrete",
      labels = c("Missing","Heterozygous","Minor","Major"),
      title = "Genotype"
    ),
    column_title = paste0("SNP\n(Chr:",paste(regions,collapse = "_"),")"),
    row_title = paste0("Sample\n(n = ",nrow(tbl_ht),")"),border = T,border_gp = gpar(size = 1,color = 'black')
  )

pdf(file = paste0(output,"_hclust.pdf"),width = 10,height = 10)
ht_out_thin =draw(ht_hclust)
dev.off()

message(msg_run("Step5. Export results ..."))
message(msg_yes("\n-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-\n"))
##> export ordered genotype file

sample_order <- row_order(ht_out)

tbl_ht_order <- 
  tbl_ht %>% 
    t() %>% 
    as.data.frame() %>% 
    select(all_of(sample_order)) %>% 
    rownames_to_column("SNP")

tbl_ht_order_thin <- 
  tmp_tbl%>% 
  mutate(
    zerocount = rowSums(.[,-1] != 0)
  ) %>% 
  arrange(zerocount) %>% 
  slice(1:1000)%>% 
  select(-zerocount) %>% 
  right_join(tmp_tbl %>% select("SNP"),by="SNP") %>% 
  drop_na() %>% 
  distinct() %>% 
  select(SNP,all_of(sample_order))
out = list(
  all_snp = tbl_ht_order,
  top1000 = tbl_ht_order_thin
)
write_xlsx(x = out,path = paste0(output,"_ordered_geno_num.xlsx"))

##> export phenotype in order

phenotype_order <- 
  left_join(
    data.frame(Sample = colnames(tbl_ht_order)[-1]),
    phenotype,by = "Sample"
  )
writexl::write_xlsx(phenotype_order,path = paste0(output,"_order_phenotype.xlsx"))
##> export ComplexHeatmap for InteractiveComplexHeatmap

save(tree_clado,tbl_ht,ht_num_thin_high_res,ht_hclust,ht_num,file = paste0(output,"_ComplexHeatmap.rds"))

message(msg_yes("All finish!"))
