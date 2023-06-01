########################################################
#             Prj: EMMAX Downstream analysis
#             Assignment: make emmax cov
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
		  paste0("#=========================================================#\n",
			     "#        Prj: EMMAX-mate ver 0.0.1\n",
				 "#        Assignment: make emmax cov format\n",
				 "#        Author: Shawn Wang <shawnwang2016@126.com>\n",
				 "#        Date: Wed 19 Apr, 2023\n",
				 "#=========================================================#\n\n"
						      )
		  ))
##> getopt
command=matrix(c(
		   'help', 'h', 0, 'logic', 'help information',
		    'nosex', 'n', 1, 'character', 'plink file, xxx.nosex ',
		    'Q', 'q', 1, 'character','admixture file, xxx.Q',
		    'outfile', 'o', 2, 'character','outformat file name',
		    'num', 'k', 2, 'numeric','num of used k or eigenvec'
		     ),byrow = T, ncol = 5)
args = getopt(command)

##> help information
if (!is.null(args$help)) {
	  message(msg_yes(paste(getopt(command, usage = T), "\n")))
  q(status=1)
}

##> default value
if (is.null(args$nosex)){
	  message(msg_no("Need plink file xxx.nosex! please read help carefully!"))
  message(msg_yes(paste(getopt(command, usage = T), "\n")))
    q(status = 1)
}
if (is.null(args$Q)){
	  message(msg_no("Need admixture cov file, xxx.Q! please read help carefully!"))
  message(msg_yes(paste(getopt(command, usage = T), "\n")))
    q(status = 1)
}
if (is.null(args$outfile)){
	  message(msg_warning("The output file name is not set, and the current time is used by default"))
  args$outfile = paste0('cov.',gsub(pattern = ' ',replacement = '_',Sys.time()))
}
if (is.null(args$num)){
	  message(msg_warning("The parameter n is not set, use default 3"))
  args$num = 3
}
# library for data cleaning and visualize ---------------------------------
suppressMessages(if (!require('tidyverse')) install.packages('tidyverse'))
# test --------------------------------------------------------------------
if(TEST == "TRUE") {

} else {
	nosex = args$nosex
    Q = args$Q
    outfile = args$outfile
	num = args$num-1
}
print(num)
message(msg_yes("\n-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-\n"))
message(msg_yes("Import data ..."))
nosex.f <- read.table(nosex,header = F,sep = "\t")
Q.f <- read.table(Q,header = F,sep = " ")
message(msg_yes("\n-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-\n"))
message(msg_yes("Data format ..."))
out <-
	  Q.f %>%
	  select(1:num) %>% 
	  mutate(s1 = nosex.f$V1,
		     s2 = nosex.f$V2,
			 tag = 1) %>%
			 relocate(s1,.before = V1) %>%
			 relocate(s2,.after = s1) %>%
			 relocate(tag,.after = s2)
##> path
get_abs_path <- function(x) {
	current_dir <- getwd()
    file_name <- x
    file_path <- file.path(current_dir, file_name)
    absolute_path <- normalizePath(file_path)
    return(absolute_path)
}
outfile_path = suppressWarnings(get_abs_path(paste0(outfile,"_emmax.cov")))
message(msg_yes("\n-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-\n"))
message(msg_yes(
		  paste0("Export Covariate ...\n","The result file: ",outfile_path)
		  ))

write.table(out,file = outfile_path,row.names = F,col.names = F,sep = "\t",quote = F)

message(msg_yes("\n-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-\n"))
message(msg_yes("All finish!"))

