########################################################
#             Prj: EMMAX Downstream analysis
#             Assignment: visulization
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
         "#        Assignment: Visualization and sig SNP extract\n",
         "#        Author: Shawn Wang <shawnwang2016@126.com>\n",
         "#        Date: Wed 19 Apr, 2023\n",
         "#=========================================================#\n\n"
  )
))
##> getopt
command=matrix(c(
  'help', 'h', 0, 'logic', 'help information',
  'emmax_ps', 'e', 1, 'character', 'Input emmax result file: emmax output, xxx.ps. ',
  'trait', 't', 1, 'character','Trait name',
  'annotation', 'a', 2, 'character','SNP annotation',
  'out_filename', 'o', 2, 'character', 'input metadata file: compound metadata, include compound annotations',
  'cutoff', 'f', 2, 'integer', '-log10P cutoff,default = 4',
  'img_type', 'i', 2, 'character', 'output file type. jpg, pdf or ',
  'point_size', 's', 2, 'character', 'The point size of 0-4-6 ',
  'point_color', 'w', 2, 'character', 'The point color of 4-6  '
),byrow = T, ncol = 5)
args = getopt(command)

##> help information
if (!is.null(args$help)) {
  message(msg_yes(paste(getopt(command, usage = T), "\n")))
  q(status=1)
}

##> default value
if (is.null(args$emmax_ps)){
  message(msg_no("Need emmax.ps file! please read help carefully!"))
  message(msg_yes(paste(getopt(command, usage = T), "\n")))
  q(status = 1)
}
if (is.null(args$trait)){
  message(msg_no("Need trait tag! please read help carefully!"))
  message(msg_yes(paste(getopt(command, usage = T), "\n")))
  q(status = 1)
}
if (is.null(args$out_filename)){
  message(msg_warning("The output file name is not set, and the current time is used by default"))
  args$out_filename = paste0('trait.',gsub(pattern = ' ',replacement = '_',Sys.time()))
}
if (is.null(args$cutoff)){
  message(msg_warning("No cutoffing parameter is set, the default is - log10(P) > 4"))
  args$cutoff = 4
}
if (is.null(args$annotation)){
  message(msg_warning("Missing SNP annotation file, the step of adding annotations will be skipped."))
  args$annotation = NULL
}
if (is.null(args$img_type)){
  message(msg_warning("default img file type: jpg."))
  args$img_type = "jpg"
}
if (is.null(args$point_size)){
  message(msg_warning("default point size: 0.3_0.5_0.5."))
  args$point_size = "0.3_0.5_0.5"
}
if (is.null(args$point_color)){
  message(msg_warning("default point color: red_green."))
  args$point_color = "red_green"
}
# library for data cleaning and visualize ---------------------------------
suppressMessages(if (!require('tidyverse')) install.packages('tidyverse'))
suppressMessages(if (!require('CMplot')) install.packages('CMplot'))
# test --------------------------------------------------------------------
if(TEST == "TRUE") {
  emmax_ps = "../03.Progress/Emmax/Morin_emmax_QK_log2.ps"
  trait = "Morin"
  out_filename = "Morin_QK"
  cutoff = 4
  annotation = ""
  img_type = "jpg"
  point_size = "0.3_0.5_0.5."
  point_color = "red_green"
} else {
  emmax_ps = args$emmax_ps
  trait = args$trait
  out_filename = args$out_filename
  cutoff = args$cutoff
  annotation = args$annotation
  img_type = args$img_type
  point_size = args$point_size %>% stringr::str_split(.,"_",3,T) %>% as.numeric()
  point_color = args$point_color %>% stringr::str_split(.,"_",2,T) 
}

# run ---------------------------------------------------------------------
##> Step1. vis by CMplot

message(msg_run(paste0(
  "Start running ..."
)))
message(msg_yes("\n-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-\n"))
dir.create(out_filename,showWarnings = F,recursive = T)
message(msg_yes(paste0("Create a new result folder: ",out_filename,"/")))
message(msg_yes("\n-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-\n"))
setwd(out_filename)
raw <- read.table(paste0("../",emmax_ps),header = F,sep = "\t")
anno <- read.table(paste0("../",annotation),header = T,sep = "\t")
message(msg_run("Step1. Start drawing QQ plot and Manhattan plot ..."))
message(msg_yes("\n-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-\n"))
clean <-
  raw %>%
  setNames(c("SNP","SE","pval")) %>%
  mutate(Chromosome = stringr::str_split(SNP,":",2,T)[,1],
         Position = stringr::str_split(SNP,":",2,T)[,2] %>% as.numeric()) %>%
  select(SNP,Chromosome,Position,pval) %>%
  filter(pval > 0 & pval <= 1) %>%
  set_names("SNP","Chromosome","Position",trait)
anno <-
  anno %>%
  select(Uploaded_variation,Feature_type,Gene,Consequence,Amino_acids,Codons,IMPACT) %>% 
  setNames(
    c("SNP","Feature_type","Gene","Consequence","aa_variation","Codons","Impact")
  ) %>%
  mutate(
    Chromosome = stringr::str_split(SNP,":",2,T)[,1],
    Position = stringr::str_split(SNP,":",2,T)[,2] %>% as.numeric()
  ) %>%
  relocate(Chromosome,.after = SNP) %>%
  relocate(Position,.after = Chromosome)
message(msg_run("QQ plot ..."))
CMplot(clean,plot.type="q",conf.int.col=NULL,box=TRUE,file=img_type,memo="",dpi=300,cex = .3,
       file.output=TRUE,verbose=TRUE)
message(msg_run("Manhattan plot ..."))
CMplot(clean, plot.type="m", LOG10=TRUE, ylim=NULL, threshold=c(1e-6,1e-4),threshold.lty=c(1,2),
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"),signal.col=c(point_color[2],point_color[1]),signal.cex=c(point_size[3],point_size[1]),
       signal.pch=c(19,19),file=img_type,memo="",dpi=300,file.output=TRUE,verbose=TRUE,cex = point_size[1],
       width=22,height=8)
message(msg_yes("Step1 finish!"))

message(msg_yes("\n-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-\n"))
##> Step2 run snp cutoff.

snp_cutoffing <-
  clean %>%
  rename('p' = trait) %>%
  mutate(
    nlogp = -log10(p)
  ) %>%
  filter(nlogp >= cutoff)

snp_cutoff_igv <-
  snp_cutoffing %>%
  select(Chromosome,Position,SNP,p) %>%
  setNames(c("CHR","BP","SNP","P"))
write.table(snp_cutoff_igv,file = paste0(out_filename,'_filter.gwas'),sep = "\t",quote = F,row.names = F,col.names = T)

message(msg_run(paste0("You can view significant SNP (-log10(P) > 4) by import: ",out_filename,'_filter.gwas to IGV')))
message(msg_yes("\n-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-\n"))

if(is.null(anno)) {
  message(msg_no("No annotation files. Skip step2."))
  message(msg_yes("\n-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-\n"))
} else {
  message(msg_run("Step2. SNP filtering and SNP annotation"))
  message(msg_yes("\n-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-\n"))
  snp_cutoffing <-
    snp_cutoffing %>%
    left_join(.,anno)
}


writexl::write_xlsx(snp_cutoffing,paste0(out_filename,"_filter.xlsx"))

message(msg_run(paste0("Key SNP information: ",out_filename,"_filter.xlsx")))
message(msg_yes("\n-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-\n"))

message(msg_yes("Step2 finish!"))

message(msg_yes("\n-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-\n"))

##> Step3 export genotype file
message(msg_run("Step3.Export genotype file for LDBlockShow"))
message(msg_yes("\n-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-\n"))

clean %>%
  select(Chromosome,Position,trait) %>%
  write.table(.,paste0(out_filename,"_genotype_all.txt"),
              row.names = F,col.names = F,quote = F,sep = '\t')

message(msg_run(paste0("LDBlockShow genotype file (all snp): ",out_filename,"_genotype_all.txt")))
message(msg_yes("\n-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-\n"))
dir.create(path = "LDFilebyChr",showWarnings = F,recursive = T)
tmp_sum =
  snp_cutoffing %>%
  group_by(Chromosome) %>%
  summarise(n = n()) %>%
  slice_max(order_by = n,n = 5) %>%
  filter(n >= 10) %>% pull(Chromosome)
##> write
if(length(tmp_sum)==0) {
  message(msg_no("The number of Significant SNP in same chromosome small than 10. Skip this step ..."))
  message(msg_yes("\n-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-\n"))
} else {
  message(msg_yes("Significant SNPs (Top5 chromosomes) ..."))
  message(msg_yes("\n-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-\n"))
  tmp =
    clean %>%
    select(Chromosome,Position,trait)
  map(.x = tmp_sum,.f = function(.x) {
    tmp %>%
      filter(Chromosome == .x) %>%
      write.table(
        .,paste0("LDFilebyChr/",.x,"_genotype.txt"),
        row.names = F,col.names = F,quote = F,sep = '\t'
      )
  })
  message(msg_yes("Step3 finish!"))
  message(msg_yes("\n-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-\n"))
}

message(msg_run("After defining the observation range in IGV, LD block division can be performed through the LDblockShow file. The recommended code is:"))
message(msg_yes("\n-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-\n"))
message(msg_yes("\nLDBlockShow \\ \n    -InGenotype genotype.txt \\ \n    -OutPut MATE_all \\ \n    -InGWAS gwasfile.txt \\ \n    -InGFF  Gh.gff \\ \n    -Region  A06:4142987:4305257 \\ \n    -Cutline 6 \\ \n    -OutPdf \\ \n    -SeleVar 2 \\ \n    -TopSite &"))
message(msg_yes("\n-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-\n"))
message(msg_yes("All finish!"))
