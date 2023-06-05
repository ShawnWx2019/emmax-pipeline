library(httr)
library(jsonlite)
library(xml2)
library(conflicted)
conflict_prefer(name = 'rename',winner = 'dplyr')
get_GO_anno <- function(TERMID) {
  if(length(TERMID) == 1) {
    TERMID.url = str_replace(TERMID,pattern = "\\:",replacement = "%3A")
  } else {
    TERMID.url = str_replace(TERMID,pattern = "\\:",replacement = "%3A") %>%
      paste(.,collapse = "%2C")
  }
  url = "https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/"
  requestURL <- paste0(url,TERMID.url)
  r <- GET(requestURL, accept("application/json"))

  stop_for_status(r)

  json <- toJSON(content(r))
  x_list <-fromJSON(json)
  out = data.frame(
    TERM = TERMID,
    NAME = "",
    ONT = ""
  )
  tryCatch({
    out = data.frame(
      TERM = x_list$results$id %>% unlist,
      NAME = x_list$results$name %>% unlist,
      ONT = x_list$results$aspect %>% unlist
    )
  },error = function(e){
    message(paste0(TERMID," no match result"))
  })
  return(out)
}
load("/Volumes/1T/06.博士杂优项目/HeterosisPrj/rawdata.Rdata")
xx =
CRIanno %>%
  mutate(
    GOID = str_extract_all(GO_Term,"(?<=\\()G.*?(?=\\))",simplify =TRUE)
    )

##> 1. 基因功能注释
##> 2. gene2go
##> 3. gene2kegg

##> add cottonFGD
filePath <- dir(path = "../demo/",pattern = ".xlsx")
CottonFGD = map_dfr(.x = filePath,.f = function(.x){
  readxl::read_xlsx(paste0("../demo/",.x))
}) %>%
  setNames(c("GeneID","GeneName","Description"))
##> merge annotation
gene_anno <-
  CRIanno %>%
  select(transcript_id,Chr,start,end,strand,Function1,Function2,Function3,UniProt) %>%
  rename('GeneID' = "transcript_id") %>%
  left_join(CottonFGD) %>%
  relocate(GeneName,.before = Function1) %>%
  relocate(Description,.after = GeneName) %>%
  mutate(
    start = str_remove_all(start,"_|,") %>% as.numeric(),
    end = str_remove_all(end,"_|,") %>% as.numeric(),
    GeneName = case_when(
      GeneName == "NA" ~ GeneID,
      TRUE ~ GeneName
    ),
    Description = case_when(
      Description == "NA" ~ "-",
      TRUE ~ Description
    ),
    Chr = case_when(
      str_detect(Chr,"Contig") ~ "Contig",
      TRUE ~ Chr
    )
  ) %>%
  group_by(GeneID) %>%
  slice_head(n = 1) %>% ungroup()
writexl::write_xlsx(x = gene_anno,path = "../demo/CRI_v2_anno.xlsx")

##> t2g t2n

t2g.go <- read.delim("../demo/CRI.TM1.V2.TBtools.gene2go",header = F) %>%
  select(V2,V1) %>% set_names("TERM","GENE") %>% unique()

t2n.go <- get_GO_anno(TERMID = unique(t2g.go$TERM))


##>
t2g.kegg <- read.delim("../demo/CRI.TM1.V2TBtools.gene2ko",header = F)
t2n.kegg <- read.delim("../demo/TBtools.Plant.KEGG.backEnd.20200724",header = F,quote = "")

t2n.kegg.1 <-
  t2n.kegg %>%
  select(V1,V3) %>%
  mutate(term = str_extract(V3,"\\d++(?= )"),
         name = str_remove(V3,"\\d++(?= )")) %>%
  rename("gene" = "V1") %>%
  select(gene,term,name) %>%
  distinct()
t2n2 <-
  t2n.kegg %>%
  select(V1,V4) %>%
  mutate(V4 = str_remove(V4,"B  ")) %>%
  mutate(term = str_extract(V4,"\\d++(?= )"),
         name = str_remove(V4,"\\d++(?= )")) %>%
  rename("gene" = "V1") %>%
  select(gene,term,name) %>%
  distinct()
t2n_merge <- rbind(t2n.kegg.1,t2n2)

t2g.kegg = left_join(t2g.kegg,t2n_merge,c("V2" = "gene")) %>%
  select(term,V1) %>%
  distinct() %>%
  set_names("TERM","GENE")
t2n.kegg <- t2n_merge %>%
  select(term,name) %>%
  set_names('TERM',"NAME")

enrich_db_out <- list(
  t2g.go = t2g.go,
  t2n.go = t2n.go,
  t2g.kegg = t2g.kegg,
  t2n.kegg = t2n.kegg
)

writexl::write_xlsx(x = enrich_db_out,path = "../demo/Enrichmentdb.xlsx")
