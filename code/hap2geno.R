###########################################################
#       Prj: EMMAX-mate
#       Assignment: convert hmp to LDBlockShow genotype
#       Author: Shawn Wang
#       Date: Web 19 Apr, 2023
#       Location: HENU
###########################################################

# prepare -----------------------------------------------------------------

TEST = "FALSE"
options(stringsAsFactors = F)
options(warn = -1)
options(readr.show_progress = TRUE)
suppressMessages(if (!require('getopt')) install.packages('getopt'))
suppressMessages(if (!require('crayon')) install.packages('crayon'))
suppressMessages(if (!require('readr')) install.packages('readr'))

# Args ------------------------------------------------------------------
##> crayon
msg_yes = green$bold$italic
msg_no = red$bold$italic
msg_run = blue$bold$italic$underline
msg_warning = yellow$bold$italic
message(msg_yes(
  paste0("\n#=========================================================#\n",
         "#        Prj: EMMAX-mate ver 0.0.1\n",
         "#        Assignment: Convert .hmp to LDBlockShow genotype\n",
         "#        Author: Shawn Wang <shawnwang2016@126.com>\n",
         "#        Date: Wed 19 Apr, 2023\n",
         "#=========================================================#\n\n"
  )
))
##> getopt
command=matrix(c(
  'help', 'h', 0, 'logic', 'help information',
  'hmp', 'e', 1, 'character', 'genotype in hapmap format, DO NOT accept numeric or Diploid.',
  'split_tag', 's', 1, 'numeric', 'Should the results be eported split by chromosome?',
  'rename_chr', 'r', 1, 'character', 'Should chromosomes be renamed, new_name\told_name',
  'threads', 't', 1, 'numeric', 'Number of threads used for file read and write',
  'prefix', 'p', 2, 'character','output file prifix, eg: Gh_383, DO NOT contains Special Characters such as(" ","/","\","-","$",...), ONLY ("." or "_") accepted!'
),byrow = T, ncol = 5)
args = getopt(command)

##> help information
if (!is.null(args$help)) {
  message(msg_yes(paste(getopt(command, usage = T), "\n")))
  q(status=1)
}

##> default value
if (is.null(args$hmp)){
  message(msg_no("Need hapmap file! please read help carefully!"))
  message(msg_yes(paste(getopt(command, usage = T), "\n")))
  q(status = 1)
}
if (is.null(args$prefix)){
  message(msg_warning(
    paste0("The output file prefix is not set,Use current time:\n",
           paste0('genotype.',gsub(pattern = ' ',replacement = '_',Sys.time())),
           "\nas prefix"))
  )
  args$prefix = paste0('genotype.',gsub(pattern = ' ',replacement = '_',Sys.time()))
}
if (is.null(args$rename_chr)){
  message(msg_warning(
    "For polyploid genomes, EMMAX will require changing chromosome names to numbers. \nIf you have made changes, please add the `-r` parameter to change the chromosome names back to ensure consistency with the chromosome numbering in the .gff and .gwas files.\n will skip rename step"
  )
  )
  args$prefix = NULL
}
if (is.null(args$split)){
  args$split = 0
}
if (is.null(args$threads)){
  args$threads = 8
}
# library for data cleaning and visualize ---------------------------------
suppressMessages(if (!require('tidyverse')) install.packages('tidyverse'))
# test --------------------------------------------------------------------
if(TEST == "TRUE") {
  hmp = "r_1000.txt"
  prefix = "Gh_383"
  split_tag = 1
  rename_chr = "chr_cn2.txt"
  n_threads = 8
} else {
  hmp = args$hmp
  prefix = args$prefix
  split_tag = args$split_tag
  rename_chr = args$rename_chr
  n_threads = args$threads
}


# runing ------------------------------------------------------------------

message(msg_run(paste0(
  "Start running ..."
)))
message(msg_yes("\n-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-\n"))
message(msg_yes("Import hapmap file ...."))
geno_hmp <- readr::read_delim(file = hmp,delim = "\t",progress = show_progress(),num_threads = n_threads)

message(msg_yes("\n-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-\n"))

message(msg_run("Convert to LDBlockShow required genotype format"))

geno_hmp <- geno_hmp[,-c(1,2,5:11)]

geno_hmp <-
  geno_hmp %>%
  unite("SNP",-c(1,2),sep = " ")

if (is.null(rename_chr)) {
  message(msg_yes("\n-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-\n"))
  message(msg_warning("not detected rename chr file, skip this step."))
} else {
  message(msg_yes("\n-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-\n"))
  message(msg_run("Start rename chromosome"))
  chr_rename <- read.delim(rename_chr,header = T,sep = "\t") %>% setNames(c("new","chrom"))
  geno_hmp <-
    geno_hmp %>%
    left_join(.,chr_rename,by = "chrom") %>%
    mutate(chrom = new) %>%
    select(-new)
}

message(msg_yes("\n-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-\n"))

message(msg_yes("Convert finish! "))

if(split_tag == 1) {
  colnames(geno_hmp)[1] = "chrom"
  chr = geno_hmp %>%
    pull(chrom) %>% unique()
  map(.x = chr,.f = function(.x) {
    tmp.x =
    geno_hmp %>%
      filter(.[,1] == .x)
  #  write.table(x = tmp.x,file = paste0(.x,'_',prefix,".geno"),row.names = F,sep = "\t",quote = F,col.names = F)
    write_delim(x = tmp.x,
                file = paste0(.x,'_',prefix,".geno"),
                delim = "\t",
                quote = 'none',
                col_names = F,
                num_threads = n_threads,
                progress = show_progress())
  })
  message(msg_yes("\n-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-\n"))

  message(msg_yes(
    paste0("Please check your file at:\n",getwd(),paste0("/chr_",prefix,".geno"))
  ))
  message(msg_yes("\n-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-\n"))
} else {
 # write.table(geno_hmp,paste0(prefix,".geno"),row.names = F,sep = "\t",quote = F,col.names = F)
  write_delim(x = geno_hmp,
              file = paste0(prefix,".geno"),
              delim = "\t",
              quote = 'none',
              col_names = F,
              num_threads = n_threads,
              progress = show_progress())
  message(msg_yes("\n-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-\n"))

  message(msg_yes(
    paste0("Please check your file at:\n",getwd(),"/",paste0(prefix,".geno"))
  ))
  message(msg_yes("\n-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-\n"))
}

