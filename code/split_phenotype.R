########################################################
#             Prj: EMMAX Downstream analysis
#             Assignment: split phenotype
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
         "#        Assignment: split phenotype\n",
         "#        Author: Shawn Wang <shawnwang2016@126.com>\n",
         "#        Date: Wed 19 Apr, 2023\n",
         "#=========================================================#\n\n"
  )
))
##> getopt
command=matrix(c(
  'help', 'h', 0, 'logic', 'help information',
  'tfam', 't', 1, 'character', 'transposed plink file: .tfam, xxx.tfam. ',
  'phenotype', 'p', 1, 'character','phenotype file: samples in 1st column, traits value in following columns'
),byrow = T, ncol = 5)
args = getopt(command)

##> help information
if (!is.null(args$help)) {
  message(msg_yes(paste(getopt(command, usage = T), "\n")))
  q(status=1)
}

##> default value
if (is.null(args$tfam)){
  message(msg_no("Need emmax.tfam file! please read help carefully!"))
  message(msg_yes(paste(getopt(command, usage = T), "\n")))
  q(status = 1)
}
if (is.null(args$phenotype)){
  message(msg_no("Need phenotype file! please read help carefully!"))
  message(msg_yes(paste(getopt(command, usage = T), "\n")))
  q(status = 1)
}

# library for data cleaning and visualize ---------------------------------
suppressMessages(if (!require('tidyverse')) install.packages('tidyverse'))
suppressMessages(if (!require('CMplot')) install.packages('CMplot'))
# test --------------------------------------------------------------------
if(TEST == "TRUE") {

} else {
  tfam = args$tfam
  phenotype = args$phenotype
}

# run ---------------------------------------------------------------------
##> Step1. vis by CMplot

message(msg_run(paste0(
  "Start running ..."
)))
message(msg_yes("\n-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-\n"))

a <- read.table(tfam,header = F,sep = " ") %>% select('V1','V2') %>% setNames(c('ID',"NAME"))
b <- read.table(phenotype,header = T,sep = "\t")
colnames(b)[1] = "ID"

c = left_join(a,b,by = "ID")

for (i in 3:ncol(c)) {
	  out.name = colnames(c)[i]
  d <-
    c %>%
    select(ID,NAME,all_of(out.name))
  message(msg_yes(paste0(
    "export ", out.name," (",(i-2),"/",(ncol(c)-2),")"
    )))
  write.table(d,paste0(out.name,'.txt'),row.names = F,col.names = F,sep = "\t",quote = F)
}

message(msg_yes("\n-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-\n"))
message(msg_run("finish!"))
