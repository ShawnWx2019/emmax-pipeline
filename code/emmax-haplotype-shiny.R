#################################################################
#	                 prj: shiny app
#	                 Assignment: haplotye
#	                 Author: Shawn Wang
#	                 Date: May 23, 2023
#                  Version: V0.0.1
#################################################################
suppressMessages(if (!require('tidyverse')) install.packages('tidyverse'))
suppressMessages(if (!require('devtools')) install.packages('devtools'))
suppressMessages(if (!require("treeio")) devtools::install_github("YuLab-SMU/treeio"))
suppressMessages(if (!require("ggtree")) devtools::install_github("YuLab-SMU/ggtree"))
suppressMessages(if (!require("ComplexHeatmap")) devtools::install_github("jokergoo/ComplexHeatmap"))
suppressMessages(if (!require("InteractiveComplexHeatmap")) devtools::install_github("jokergoo/InteractiveComplexHeatmap"))
suppressMessages(if (!require('ape')) install.packages('ape'))
suppressMessages(if (!require('phylogram')) install.packages('phylogram'))
suppressMessages(if (!require('writexl')) install.packages('writexl'))
suppressMessages(if (!require('magick')) install.packages('magick'))
suppressMessages(if (!require('colourpicker')) install.packages('colourpicker'))
suppressMessages(if (!require('shiny')) install.packages('shiny'))
suppressMessages(if (!require('shinyjs')) install.packages('shinyjs'))
suppressMessages(if (!require('dashboardthemes')) install.packages('dashboardthemes'))
suppressMessages(if (!require('shinydashboard')) install.packages('shinydashboard'))
suppressMessages(if (!require('DT')) install.packages('DT'))
suppressMessages(if (!require('conflicted')) BiocManager::install('conflicted',update = FALSE))
suppressMessages(if (!require('shinythemes')) install.packages('shinythemes'))
suppressMessages(if (!require('shinyjqui')) install.packages('shinyjqui'))
suppressMessages(if (!require('excelR')) install.packages('excelR'))
suppressMessages(if (!require('shinyWidgets')) install.packages('shinyWidgets'))
suppressMessages(if (!require('ggsignif')) install.packages('ggsignif'))
suppressMessages(if (!require('ggbeeswarm')) install.packages('ggbeeswarm'))
options(shiny.maxRequestSize = 300*1024^2)
options(scipen = 6)

# conflict ----------------------------------------------------------------

conflict_prefer("select","dplyr")
conflict_prefer("filter","dplyr")
conflict_prefer("rename","dplyr")
conflict_prefer("desc","dplyr")
conflict_prefer("cor","stats")
conflict_prefer("colourInput","colourpicker")


# functions ---------------------------------------------------------------

convert_geno2num <- function(genofile,start,end,sample_order) {
  ##> add sample name and convert to long table
  colnames(genofile) = c("V1","V2","V3")
  geno_clean2 <-
    genofile %>%
    dplyr::filter(
      V2 >= start & V2 <= end ## filter snps
    ) %>% tidyr::separate(col = V3,into = c(sample_order$V2),sep = ' ',convert = T) %>%
    mutate(SNP = paste(V1,V2,sep = ":")) %>%
    select(-V1,-V2) %>%
    column_to_rownames("SNP") %>%
    t() %>% as.data.frame()
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
    select(-tag1,-tag2)

  ##> convert numeric genodata and tag genodata

  geno_num <-
    geno_final_long %>%
    mutate(SNP = paste0(chr,':',pos)) %>%
    select(SNP,sample,tagNum) %>%
    pivot_wider(names_from = sample, values_from = tagNum)

  geno_out =
    list(
      geno_num = geno_num,
      geno_clean2 = geno_clean2
    )

  return(geno_out)

}

convert_tree <- function(tree) {
  ##> convert order plhlo tree and convert to dendrogram
  tree_laddered <- ladderize(tree)
  tree_dend <- as.dendrogram.phylo(tree_laddered)
  tree_clado <- as.cladogram(tree_dend)
  return(tree_clado)
}


get_snps <- function(geno_tbl,geno_num,pheno,snps,heter = TRUE){
  geno_snps_raw <-
    geno_num %>%
    column_to_rownames("SNP") %>%
    t() %>% as.data.frame() %>%
    select(all_of(snps)) %>%
    rownames_to_column("Sample") %>%
    pivot_longer(!Sample,values_to = "value",names_to = "snpsite") %>%
    filter(value != 0)
  if(isFALSE(heter)) {
    geno_snps_raw <-
      geno_snps_raw %>%
      filter(value != 1)
  }
  geno_snps_sample_ID <-
    geno_snps_raw %>%
    pivot_wider(names_from = "snpsite",values_from = "value",values_fill = NA) %>%
    drop_na() %>%
    select(Sample)

  tmp.p = pheno %>%
    select(Sample,trait,group)

  geno_tbl2 <-
    geno_tbl %>%
    select(all_of(snps)) %>%
    rownames_to_column("Sample") %>%
    inner_join(geno_snps_sample_ID,by = "Sample") %>%
    left_join(tmp.p,by = "Sample")
  haplotype_tmp <-
    geno_tbl2 %>%
    unite(haplotype,all_of(snps),sep = "+") %>%
    select(Sample,haplotype)
  geno_tbl2 <- left_join(
    geno_tbl2,haplotype_tmp,by = "Sample"
  ) %>%
    filter(!is.na(trait))
  geno_tbl2_clean <-
    geno_tbl2 %>%
    select(Sample,haplotype,group,trait)
  out <- list(
    geno_tbl2 = geno_tbl2,
    clean = geno_tbl2_clean
  )
  return(out)
}
# plot function -----------------------------------------------------------

reform_region_haptbl <- function(x,tbl.group) {
  tbl.group <-
    tbl.group %>%
    setNames(c("Sample","group"))
  tbl.out <- left_join(x,tbl.group,by = "Sample") %>%
    relocate(group,.after = haplotype)
  return(tbl.out)
}

##> for snps input data


##> DA function
run_da <- function(hap1,hap2,mat,method) {
  mat1 <- mat %>%
    filter(haplotype == hap1) %>%
    pull(trait)
  mat2 <- mat %>%
    filter(haplotype == hap2) %>%
    pull(trait)
  if (method == "t.test") {
    res = t.test(mat1,mat2)
    pval = res$p.value
  } else if (method == "wilcox.test") {
    res = wilcox.test(mat1,mat2)
    pval = res$p.value
  }
  out =data.frame(
    left = hap1,
    right = hap2,
    pval = pval
  )
  return(out)
}

##> 生成pairwise比对组合
get_hap_pairwise = function(n,hap_index) {
  if(n == 2) {
    out_tbl = data.frame(
      left = "hap1",
      right = "hap2"
    )
  } else {
    left <- c()
    right <- c()
    for (i in n:1) {
      left = c(left,rep(i,i-1))
    }
    for (i in c((n-1):2)) {
      right = c(right,c(i:1))
    }
    out_tbl = data.frame(
      left = paste0("hap",left),
      right = paste0("hap",c(right,1))
    )
  }
  k_v_tbl <-
    data.frame(hap = hap_index,
               tag = paste0("hap",seq(1:n)))
  out_tbl2 =
    out_tbl %>%
    left_join(
      .,k_v_tbl,by = c('left'='tag')
    ) %>% mutate(left = hap) %>%
    select(-hap) %>%
    left_join(
      .,k_v_tbl,by = c('right'='tag')
    ) %>% mutate(right = hap)%>%
    select(-hap)
  return(out_tbl2)
}


get_da_result <- function(input_tbl,DA_method,fig_type,sig_tag) {
  if (sig_tag == 'numeric') {
    sig_judge = FALSE
  } else {
    sig_judge = TRUE
  }
  input_tbl <-
    input_tbl %>%
    mutate(trait = as.numeric(trait)) %>%
    drop_na()
  ##> Step1 get haplotype index
  hap_index = input_tbl$haplotype %>% unique()
  ##> 1.1 pairwise-index
  c_list <- get_hap_pairwise(n = length(hap_index),hap_index = hap_index)
  ##> 1.2 pairwise DA
  diff_res.test <-
    map2_dfr(.x = c_list$left,.y = c_list$right,.f = function(.x,.y) {
      tryCatch({
        run_da(hap1 = .x,hap2 = .y,mat = input_tbl,method = DA_method)
      }, error = function(e) {
        cat (paste0("Warning: ",.x,"_vs_",.y," failed\n"))
      })
    })
  ##> Step2. rename haplotype and get significant result.
  ##> 2.1 add number of haplotype
  tag <-
    input_tbl %>%
    group_by(haplotype) %>%
    summarise(n=n()) %>%
    mutate(tag = paste0(haplotype,"\n(n=",n,")")) %>%
    ungroup()
  ##> 2.2 re-union data (haplotype with number)
  if("group" %in% colnames(input_tbl)) {
    clean_tbl <-
      tag %>% as.data.frame() %>%
      right_join(.,input_tbl,by = 'haplotype',multiple = "all") %>%
      mutate(haplotype = tag) %>%
      select(Sample,haplotype,group,trait) %>% ungroup()
  } else {
    clean_tbl <-
      tag %>% as.data.frame() %>%
      right_join(.,input_tbl,by = 'haplotype',multiple = "all") %>%
      mutate(haplotype = tag) %>%
      select(Sample,haplotype,trait) %>% ungroup()
  }

  ##> 2.3 pick out pval < 0.05
  sig_res <-
    diff_res.test %>%
    filter(pval < 0.05) %>%
    left_join(tag,c('left' = 'haplotype')) %>%
    mutate(left = tag) %>%
    select(-tag) %>%
    left_join(tag,c('right' = 'haplotype')) %>%
    mutate(right = tag) %>%
    select(-tag)
  ##> 2.4 significant different pair for ggsignif
  c_list = map2(.x = sig_res$left,.y = sig_res$right,.f = function(.x,.y) {
    c(.x,.y)
  })
  ##> Step3. generated raw fig
  plt <-
    ggplot(clean_tbl %>% filter(haplotype > 1),mapping = aes(x = haplotype,y = trait,fill = haplotype)) +
    stat_boxplot(geom = "errorbar",width = 0.3,position = position_dodge(width = 1),color='black')+
    geom_boxplot(outlier.colour = 'black',
                 outlier.size = 2,
                 outlier.shape = 17,
                 alpha = .5,width =.5)+
    guides(color=guide_legend("haplotype"),fill = "none")+
    xlab("")+
    ylab(paste0("Content (ng/g)"))+
    theme_bw()+
    theme(
      axis.text = element_text(size = 12,color = 'black'),
      axis.title = element_text(size = 13,color = "black"),
      panel.border = element_rect(linewidth = 1),
      plot.title = element_text(hjust = .5,size = 13),
      legend.position = 'right'
    )
  if(length(c_list) > 0) {
    plt = plt +
      ggsignif::geom_signif(
        comparisons = c_list,
        test = DA_method,map_signif_level = sig_judge
      )
  }

  if("group" %in% colnames(clean_tbl)) {
    if (fig_type == "boxplot") {
      plt_out = plt
    } else if (fig_type == "boxplot+beeswarm") {
      plt_out = plt +
        geom_beeswarm(aes(x = haplotype,y = trait, color = group),
                      size = 1,alpha = .8)+
        guides(color = guide_legend(override.aes = list(size = 3)))+
        labs(color = "Group")
    }
  } else {
    if (fig_type == "boxplot") {
      plt_out = plt
    } else if (fig_type == "boxplot+beeswarm") {
      plt_out = plt +
        geom_beeswarm(aes(x = haplotype,y = trait),
                      size = 1,alpha = .8)+
        guides(color = guide_legend(override.aes = list(size = 3)))+
        labs(color = "Group")
    }
  }
  data_out <-
    list(
      DA_result = diff_res.test,
      tbl = clean_tbl,
      plt = plt_out
    )
  return(data_out)

}





# 01. UI =========================
## logo
customLogo <- shinyDashboardLogoDIY(

  boldText = "ShawnLearnBioinfo"
  ,mainText = "Emmax-haplotype"
  ,textSize = 14
  ,badgeText = "V0.0.6.230526"
  ,badgeTextColor = "white"
  ,badgeTextSize = 2
  ,badgeBackColor = "orange"
  ,badgeBorderRadius = 3
)


##> tag
hr_main = tags$hr(style = "border-top: 6px double #008080; border-bottom: 3px solid #008080;")
hr_bar = tags$hr(style = "border-top: 3px dotted #008080;")
hr_head = tags$hr(
  style = "border: 0;
    padding-top: 1.5px;
    background: linear-gradient(to right, transparent, #008080, transparent);"
)
ui <- shinyUI(
  navbarPage(
    theme = shinytheme("spacelab"),
    customLogo,
    tabPanel(
      useShinyjs(),
      title = "haplotype partitioning",
      icon = icon("bars"),
      sidebarLayout(
        div(id = "Sidebar",
            sidebarPanel(
              width = 2,
              tags$h3("Data input",style = 'color: #008080'),
              hr_bar,
              fileInput(
                inputId = "genopath",
                label = "genotype file (require)",
                accept = c(".txt",".geno")
              ),
              p("Translation: A 3-column table with tab as the delimiter. The first column is Chr, the second column is Pos, and the third column represents the genotype of the SNP in the population (separated by spaces).",
                style = "color: #7a8788;font-size: 12px; font-style:Italic"),
              fileInput(
                inputId = "fampath",
                label = "fam file (require)",
                accept = c(".fam")
              ),
              p(".fam file generated by plink",
                style = "color: #7a8788;font-size: 12px; font-style:Italic"),
              fileInput(
                inputId = "treepath",
                label = "tree file (option)",
                accept = c(".nwk")
              ),
              p(".nwk tree file (NJ tree or ML tree).",
                style = "color: #7a8788;font-size: 12px; font-style:Italic"),
              fileInput(
                inputId = "phenopath",
                label = "phenotype file  (require)",
                accept = c(".txt",".xls",".xlsx")
              ),
              p("Sample ID in 1st column, traits in following column..",
                style = "color: #7a8788;font-size: 12px; font-style:Italic"),
              hr_bar,
              tags$h3("Regions",style = 'color: #008080'),
              hr_bar,
              textInput(
                inputId = "chr",
                label = "chr",
                value = "A01"
              ),
              textInput(
                inputId = "start",
                label = "start",
                value = "100000"
              ),
              textInput(
                inputId = "end",
                label = "end",
                value = "200000"
              ),
              actionButton(inputId = "confirm_input",label = "Confirm input",icon = icon("check"))
            )),
        mainPanel(
          fluidPage(
            actionButton('toggleSidebar',"Toggle sidebar"),
            tabsetPanel(
              tabPanel(
                title = "Heatmap",height = "500px",width = "100%",
                icon = icon('fire'),
                div(id = "parameters",
                    style = "position: fixed;
                           top: 120px; right: -300px; width: 300px; height: 80%;
                           background-color: #f7f7f7; border: 1.5px solid #008080; overflow-y: scroll;
                           ; border-top-left-radius: 10px; border-bottom-left-radius: 10px; padding: 20px;",
                    h3("Allele color",style = "color:#darkgreen"),
                    hr_main,
                    colourInput(
                      inputId = "Missing",
                      label = "Missing",
                      value = "grey"
                    ),
                    hr_head,
                    colourInput(
                      inputId = "Heterozygous",
                      label = "Heterozygous",
                      value = "blue"
                    ),
                    hr_head,
                    colourInput(
                      inputId = "Minor",
                      label = "Minor",
                      value = "orange"
                    ),
                    hr_head,
                    colourInput(
                      inputId = "Major",
                      label = "Major",
                      value = "cyan"
                    ),
                    h3("Other Parameters",style = "color:#darkgreen"),
                    hr_main,
                    selectInput(
                      inputId = "treeMethod",
                      label = "cluster method",
                      choices = c("select a method","hclust","..."),
                      multiple = F,
                      selected = 'select a method'
                    ),
                    hr_head,
                    textInput(
                      inputId = "SNPnum",
                      label = "SNP filter (top n)",
                      value = "use all snps"
                    ),
                    hr_head,
                    radioButtons(
                      inputId = "use_raster",
                      label = "use raster",
                      choices = c("TRUE","FALSE"),
                      selected = 'FALSE')
                ),
                tags$h3("Parameter settings",style = 'color: #008080'),
                hr_main,
                div(actionButton("para_select1", "Show | Hide",icon = icon(name = 'option-vertical',lib = 'glyphicon')),
                    style = "margin-bottom: 15px;"),
                actionButton(inputId = "ht_para",label = "Confirm heatmap parameter."),
                tags$h3("Heatmap",style = 'color: #008080'),
                hr_main,
                actionButton(inputId = "show_ht",label = "Show heatmap",icon = icon('play')),
                htmlOutput("heatmap_output"),
                tags$h3("Add haplotype information",style = 'color: #008080'),
                hr_main,
                actionButton(inputId = "update_table",label = "update phenotype data",icon = icon("bell")),
                excelOutput("phenotable"),
                tags$h3("Confirm haplotype table",style = 'color: #008080'),
                hr_main,
                DT::dataTableOutput(outputId = "check_haplotype"),
                downloadButton(
                  outputId = "down_tbl_00",
                  label = 'Download',
                  icon = icon('download')
                ),
                tags$script("
                        function toggleParameters() {
                          var div = document.getElementById('parameters');
                          var right = div.style.right;
                          if (right === '-300px') {
                            div.style.right = '0px';
                            div.style.transition = 'right 0.8s ease-in-out';
                            return;
                          } else if (right === '0px') {
                            div.style.right = '-300px';
                            div.style.transition = 'right 0.8s ease-in-out';
                            return;
                          }
                        }
                      ")
              )
            )
          )
        )
      )
    ),
    tabPanel(
      useShinyjs(),
      title = "Phenotypic Differences between Haplotypes",
      icon = icon("hippo"),
      sidebarLayout(
        div(id = "Sidebar2",
            sidebarPanel(
              width = 2,
              tags$h3("DA parameters",style = 'color: #008080'),
              hr_bar,
              fileInput(
                inputId = "sample_group",
                label = "Sample Group (with header option)",
                accept = c(".txt",".xls")
              ),
              p("Classification of samples can be based on geographical distribution or grouping according to population structure, or any other classification with biological significance.",
                style = "color: #7a8788;font-size: 12px; font-style:Italic"),
              tags$h3("DA parameters",style = 'color: #008080'),
              hr_bar,
              selectInput(
                inputId = "M_test",
                label = "Methods of Hypothesis Testing",
                choices = c("t.test","wilcox.test"),
                selected = "t.test",
                multiple = F
              ),
              selectInput(
                inputId = "sig_tag",
                label = "Significant tag",
                choices = c("numeric","*"),
                selected = "numeric",
                multiple = F
              ),
              selectInput(
                inputId = "fig_type",
                label = "Figure type",
                choices = c("boxplot","boxplot+beeswarm","beeswarm"),
                selected = 'boxplot+beeswarm',
                multiple = F
              ),
              selectInput(
                inputId = "trait",
                label = "trait for analysis",
                choices = c("..."),
                selected = '...',
                multiple = F
              ),
              actionButton(inputId = "confirm2",label = "confirm parameters",icon = icon("check"))
            )
        ),
        mainPanel(
          fluidPage(
            actionButton("toggleSidebar2",
                         "Toggle sidebar"),
            tabsetPanel(
              tabPanel(
                title = "haplotype",height = "500px",width ="100%",
                icon = icon("image"),
                tags$h3("Start",style = 'color: #008080'),
                hr_main,
                p("We will remove haplotypes that have only 1 sample.",
                  style = "color: blue;font-size: 13px; font-style:Italic"),
                actionButton(inputId = "showfig",label = "Show figure"),
                tags$h3("Different analysis result",style = 'color: #008080'),
                hr_main,
                DT::dataTableOutput("DA_tbl"),
                downloadButton(
                  outputId = "down_tbl_01",
                  label = "Download",
                  icon = icon('download')
                ),
                tags$h3("Plot",style = 'color: #008080'),
                hr_main,
                jqui_resizable(
                  plotOutput(outputId = "DA_plt")
                ),
                textInput(inputId = "width1",
                          label = "width",
                          value = 10),
                textInput(inputId = "height1",
                          label = "height",
                          value = 10),
                actionButton(
                  inputId = 'adjust1',
                  label = "Set fig size",
                  icon = icon('check')
                ),
                downloadButton(
                  outputId = "down_plt_01",
                  label = "Download",
                  icon = icon('image')
                ),
                tags$h3("Plot table",style = 'color: #008080'),
                hr_main,
                DT::dataTableOutput("DA_tblp"),
                downloadButton(
                  outputId = "down_tbl_02",
                  label = "Download",
                  icon = icon('download')
                ),
              ),
              tabPanel(
                title = "Selected SNPs",height = "500px",width ="100%",
                icon = icon("paperclip"),
                tags$h3("Choose snps",style = 'color: #008080'),
                hr_main,
                selectInput(
                  inputId = "SNPs",
                  label = "Select SNPs",
                  choices = c("..."),
                  selected = "...",
                  multiple = T
                ),
                radioButtons(
                  inputId = "heter",
                  label = "Contain heter snps",
                  choices = c("FALSE","TRUE"),
                  selected = "FALSE"
                ),
                actionButton(inputId = "get_snp",label = "confirm selection",icon = icon("check")),
                tags$h3("Check SNPs",style = 'color: #008080'),
                hr_main,
                DT::dataTableOutput(outputId = "checkSNPs"),
                downloadButton(
                  outputId = "down_tbl_03",
                  label = "Download",
                  icon = icon('download')
                ),
                tags$h3("Start",style = 'color: #008080'),
                p("We will remove haplotypes that have only 1 sample.",
                  style = "color: blue;font-size: 13px; font-style:Italic"),
                hr_main,
                actionButton(inputId = "showfig2",label = "Show figure"),
                tags$h3("Different analysis result",style = 'color: #008080'),
                hr_main,
                DT::dataTableOutput("DA_tbl2"),
                downloadButton(
                  outputId = "down_tbl_04",
                  label = "Download",
                  icon = icon('download')
                ),
                tags$h3("Plot",style = 'color: #008080'),
                hr_main,
                jqui_resizable(
                  plotOutput("DA_plt2")
                ),
                textInput(inputId = "width2",
                          label = "width",
                          value = 10),
                textInput(inputId = "height2",
                          label = "height",
                          value = 10),
                actionButton(
                  inputId = 'adjust2',
                  label = "Set fig size",
                  icon = icon('check')
                ),
                downloadButton(
                  outputId = "down_plt_02",
                  label = "Download",
                  icon = icon('image')
                )
              )
            )
          )
        )
      )
    ),
    tabPanel(
      useShinyjs(),
      title = "Other evidence",
      icon = icon(name = "pizza-slice",lib = "font-awesome"),
      sidebarLayout(
        div(id = "Sidebar3",
            sidebarPanel(
              width = 2,
              tags$h3("Import database",style = 'color: #008080'),
              hr_bar,
              fileInput(
                inputId = "ExpMat",
                label = "Expression matrix \n(tab-delimed table with header, option)",
                accept = c(".txt",".xls")
              ),
              p("Expression matrix, genes in 1st column, samples in 1 row.",
                style = "color: #7a8788;font-size: 12px; font-style:Italic"),
              fileInput(
                inputId = "SnpAnno",
                label = "Snp Annotation \n(tab-delimed table with header, option)",
                accept = c(".txt",".xls")
              ),
              p("The first 3 columns is SNP.ID, Chr and Position. The column name of the column where the gene is located must be 「Gene」. Besides, the SNP id needs to be completely consistent with the genotype file.",
                style = "color: #7a8788;font-size: 12px; font-style:Italic"),
              textInput(inputId = "regions",label = "regions or SNP id",value = "11:100000-200000"),
              actionButton(inputId = "status_check",label = "confirm parameters",icon = icon("check"))
            )
        ),
        mainPanel(
          fluidPage(
            actionButton("toggleSidebar3",
                         "Toggle sidebar"),
            tabsetPanel(
              tabPanel(
                title = "Check your input",height = "500px",width ="100%",
                icon = icon("box"),
                tags$h3("Status",style = 'color: #008080'),
                hr_main,
                htmlOutput("ExpmatCheck"),
                htmlOutput("SNPannoCheck"),
                tags$h3("Expression matrix",style = 'color: #008080'),
                hr_main,
                DT::dataTableOutput("expmat_tbl"),
                tags$h3("SNP annotation",style = 'color: #008080'),
                hr_main,
                DT::dataTableOutput("snpAnno_tbl")
              ),
              tabPanel(
                title = "Candidate gene expression profile",height = "500px",width ="100%",
                icon = icon("image"),
                div(id = "parameters2",
                    style = "position: fixed;
                           top: 120px; right: -300px; width: 300px; height: 80%;
                           background-color: #f7f7f7; border: 1.5px solid #008080; overflow-y: scroll;
                           ; border-top-left-radius: 10px; border-bottom-left-radius: 10px; padding: 20px;",
                    h3("Color key",style = "color:#darkgreen"),
                    hr_main,
                    colourInput(
                      inputId = "col_min",
                      label = "min",
                      value = "blue"
                    ),
                    hr_head,
                    colourInput(
                      inputId = "col_mid",
                      label = "mid",
                      value = "white"
                    ),
                    hr_head,
                    colourInput(
                      inputId = "col_max",
                      label = "max",
                      value = "red"
                    ),
                    hr_head,
                    textInput(
                      inputId = "col_break",
                      label = 'break',
                      value = "-2,0,2"
                    ),
                    h3("Show names",style = "color:#darkgreen"),
                    hr_main,
                    radioButtons(
                      inputId = 'show_colname',
                      label = "show column names",
                      choices = c("TRUE","FALSE"),
                      selected = "TRUE"
                    ),
                    hr_head,
                    radioButtons(
                      inputId = 'show_rownames',
                      label = "show row names",
                      choices = c("TRUE","FALSE"),
                      selected = "TRUE"
                    ),
                    hr_head,
                    radioButtons(
                      inputId = "cluster_col",
                      label = "cluster column",
                      choices = c("TRUE","FALSE"),
                      selected = 'FALSE'),
                    hr_head,
                    radioButtons(
                      inputId = 'cluster_col',
                      label = "cluster row",
                      choices = c("TRUE","FALSE"),
                      selected = "TRUE"
                    )
                ),
                div(actionButton("para_select2", "Show | Hide",icon = icon(name = 'option-vertical',lib = 'glyphicon')),
                    style = "margin-bottom: 15px;"),
                actionButton(inputId = "ht_para2",label = "Confirm heatmap parameter."),
                tags$h3("Show heatmap",style = 'color: #008080'),
                hr_main,
                jqui_resizable(
                  plotOutput("expmat_ht")
                ),
                tags$script("
                        function toggleParameters2() {
                          var div = document.getElementById('parameters2');
                          var right = div.style.right;
                          if (right === '-300px') {
                            div.style.right = '0px';
                            div.style.transition = 'right 0.8s ease-in-out';
                            return;
                          } else if (right === '0px') {
                            div.style.right = '-300px';
                            div.style.transition = 'right 0.8s ease-in-out';
                            return;
                          }
                        }
                      ")
              )
            )
          )
        )
      )
    )
  )

)


server <- function(input,output,session) {
  #> toggle side bar
  observeEvent(input$toggleSidebar, {
    shinyjs::toggle(id = "Sidebar")
  })

  observeEvent(input$para_select1,{
    shinyjs::runjs('toggleParameters();')
  })

  observeEvent(input$para_select2,{
    shinyjs::runjs('toggleParameters2();')
  })
  genofile <- reactive({
    file1 <- input$genopath
    if(is.null(file1)){return()}
    read.delim(file = file1$datapath,header = F)
  })
  tree <- reactive({
    file2 <- input$treepath
    if(is.null(file2)){return()}
    read.tree(file = file2$datapath)
  })

  order_method = reactive({
    if(is.null(tree())) {"hclust"} else {c("tree","hclust")}
  })

  observe({
    updateSelectInput(session,inputId = "treeMethod",choices = order_method())
  })
  famfile <- reactive({
    file3 <- input$fampath
    if(is.null(file3)){return()}
    read.delim(file = file3$datapath,header = F)
  })
  phenfile <- reactive({
    file4 <- input$phenopath
    if(is.null(file4)){return()}
    read.delim(file = file4$datapath,header = T)
  })
  #> set reactiveValue
  ht_obj <- reactiveValues(data=NULL)
  downloads <- reactiveValues(data = NULL)
  #> Set parameter
  observeEvent(
    input$confirm_input,
    {
      if(is.null(genofile)){return()}
      if(is.null(famfile)){return()}
      if(is.null(phenfile)){return()}
      #> change colnames
      colnames(phenfile())[1] == "Sample"
      ht_obj$chr = as.character(input$chr)
      ht_obj$start = as.numeric(input$start)
      ht_obj$end = as.numeric(input$end)
      ht_obj$sample_order = famfile()
      ht_obj$regions = paste0(ht_obj$chr,":",ht_obj$start,"-",ht_obj$end)
      ht_obj$geno_num = list()
      ht_obj$tbl_ht = list()
      ht_obj$tree = list()
      ht_obj$tree_clado = list()
      ht_obj$topn = list()
      ht_obj$use_raster = list()
      progress_confirm = c(
        "Step1. Convert phenotype data into numeric",
        "Step2. Set cluster method",
        "Step3. Set SNP filter method"
      )
      withProgress(message = "confirm input data and parameters",value = 0,
                   expr = {
                     for (i in 1:3) {
                       incProgress(1/3,detail = progress_confirm[i])
                       if(i == 1) {
                         ht_obj$geno_reform = convert_geno2num(
                           genofile = genofile(),
                           start = ht_obj$start,
                           end = ht_obj$end,
                           sample_order = ht_obj$sample_order
                         )
                         ht_obj$geno_num = ht_obj$geno_reform$geno_num
                       } else if(i == 2) {
                         if(is.null(tree())) {
                           ht_obj$tbl_ht <-
                             ht_obj$geno_num %>%
                             column_to_rownames("SNP") %>%
                             t()
                         } else {
                           ht_obj$tree = tree()
                           ht_obj$tree_clado = as.cladogram(
                             as.dendrogram.phylo(
                               ladderize(ht_obj$tree)
                             )
                           )
                           ht_obj$tbl_ht <-
                             ht_obj$geno_num %>%
                             column_to_rownames("SNP") %>%
                             select(ht_obj$tree$tip.label) %>%
                             t()
                         }
                       } else if (i == 3) {
                         ht_obj$topn = ncol(ht_obj$tbl_ht)
                         ht_obj$use_raster = as.logical(input$use_raster)
                       } else {
                         return()
                       }
                     }
                   })
      observe({
        updateSelectInput(session,inputId = "SNPnum",choices = ht_obj$topn )
      })
    }
  )
  observeEvent(
    input$ht_para,
    {
      if(is.null(genofile)){return()}
      if(is.null(famfile)){return()}
      if(is.null(phenfile)){return()}
      ht_obj$col_miss = as.character(input$Missing)
      ht_obj$col_hete = as.character(input$Heterozygous)
      ht_obj$col_minor = as.character(input$Minor)
      ht_obj$col_major = as.character(input$Major)
      ht_obj$order_method = as.character(input$treeMethod)
      ht_obj$colorkey = c(ht_obj$col_miss,ht_obj$col_hete,ht_obj$col_minor,ht_obj$col_major)
      ht_obj$tmp_tbl = list()
      ht_obj$final_tbl = list()
      ht_obj$ht = list()
      ht_obj$geno_num = list()
      ht_obj$ordered_pheno = list()
      progress_ht = c(
        "Step1. Reform data",
        "Step2. Draw heatmap",
        "Step3. Ordered phenotype data by Heatmap cluster."
      )
      withProgress(message = "Draw heatmap",value = 0,
                   expr = {
                     for (i in 1:3) {
                       incProgress(1/3,detail = progress_ht[i])
                       if(i == 1) {
                         ht_obj$tmp_tbl =
                           ht_obj$tbl_ht %>%
                           t() %>% as.data.frame() %>% rownames_to_column("SNP")
                         ht_obj$final_tbl =
                           ht_obj$tmp_tbl %>%
                           mutate(
                             zerocount = rowSums(.[,-1] != 0)
                           ) %>%
                           arrange(zerocount) %>%
                           slice(1:ht_obj$topn)%>%
                           select(-zerocount) %>%
                           right_join(ht_obj$tmp_tbl %>% select("SNP"),by="SNP") %>%
                           drop_na() %>%
                           distinct() %>%
                           column_to_rownames("SNP") %>%
                           t()
                       } else if (i == 2) {
                         if(ht_obj$order_method == "hclust") {
                           ht_obj$ht =
                             Heatmap(
                               matrix = ht_obj$final_tbl,
                               col = circlize::colorRamp2(breaks = c(0,1,2,3),colors = ht_obj$colorkey),
                               show_row_names = F,show_column_names = F,
                               cluster_columns = F,cluster_rows = T,
                               use_raster = ht_obj$use_raster,
                               heatmap_legend_param = list(
                                 color_bar = "discrete",
                                 labels = c("Missing","Heterozygous","Minor","Major"),
                                 title = "Genotype"
                               ),
                               column_title = paste0("SNP\n(Chr:",paste(ht_obj$regions,collapse = "_"),")"),
                               row_title = paste0("Sample\n(n = ",nrow(ht_obj$final_tbl),")"),border = T,border_gp = gpar(size = 1,color = 'black')
                             )
                         } else {
                           ht_obj$ht =
                             Heatmap(
                               matrix = ht_obj$final_tbl,
                               col = circlize::colorRamp2(breaks = c(0,1,2,3),colors = ht_obj$colorkey),
                               show_row_names = F,show_column_names = F,
                               cluster_columns = F,cluster_rows = ht_obj$tree_clado,
                               use_raster = ht_obj$use_raster,
                               heatmap_legend_param = list(
                                 color_bar = "discrete",
                                 labels = c("Missing","Heterozygous","Minor","Major"),
                                 title = "Genotype"
                               ),
                               column_title = paste0("SNP\n(Chr:",paste(ht_obj$regions,collapse = "_"),")"),
                               row_title = paste0("Sample\n(n = ",nrow(ht_obj$final_tbl),")"),border = T,border_gp = gpar(size = 1,color = 'black')
                             )
                         }
                       } else if (i == 3) {
                         ht_obj$geno_num =
                           ht_obj$final_tbl %>%
                           t() %>%
                           as.data.frame() %>%
                           select(all_of(row_order(ht_obj$ht))) %>%
                           rownames_to_column("SNP")
                         ht_obj$ordered_pheno =
                           left_join(data.frame(Sample = colnames(ht_obj$geno_num)[-1]),phenfile(),by = "Sample") %>%
                           mutate(haplotype = "") %>%
                           relocate(haplotype,.after = Sample)
                       } else {
                         return()
                       }
                     }
                   })
    }
  )

  observeEvent(
    input$show_ht,
    {
      shinyjs::toggle(id = "Sidebar")
      InteractiveComplexHeatmapWidget(input, output, session, ht_obj$ht,output_id = "heatmap_output")
      output$phenotable <- renderExcel({
        if(is.null(ht_obj$ht)) {return()}
        if(is.null(phenfile())) {return()}
        excelTable(data = ht_obj$ordered_pheno,showToolbar = T,search = T,autoFill = T,)
      })
    }
  )

  observeEvent(
    input$update_table,
    {
      ht_obj$phenotable2= list()
      progress_update = c(
        "Step1. update information",
        "Finish!"
      )
      withProgress(message = "Update",value = 0,
                   expr = {
                     for (i in 1:2) {
                       incProgress(1/2,progress_update[i])
                       if(i == 1) {
                         ht_obj$phenotable2 = excel_to_R(input$phenotable)
                       } else {return()}
                     }
                   }
      )
    }
  )

  output$check_haplotype = DT::renderDataTable(
    DT::datatable(
      {
        input$update_table
        if(is.null(ht_obj$ht)) {return()}
        if(is.null(phenfile())) {return()}
        ht_obj$phenotable2
      },
      extensions = 'Buttons',
      options = list(
        autoWidth = T,
        dom = 'Bfrtip',
        scrollX = T,
        lengthMenu = list(c(5, 15, -1), c('5', '15', 'All')),
        pageLength = 15,
        scrollY = "400px",
        buttons = list(
          list(
            extend = "collection",
            text = 'Show All',
            action = DT::JS("function ( e, dt, node, config ) {
                                    dt.page.len(-1);
                                    dt.ajax.reload();
                                }")),
          list(
            extend = "collection",
            text = 'Show Less',
            action = DT::JS("function ( e, dt, node, config ) {
                              dt.page.len(15);
                              dt.ajax.reload();}")
          )
        )
      )
    )
  )


  ##> DA
  DA_obj <- reactiveValues(data=NULL)
  groupFile <- reactive({
    file5 <- input$sample_group
    if(is.null(file5)){return()}
    read.delim(file = file5$datapath,header = T)
  })
  observe({
    updateSelectInput(session,inputId = "trait",choices = colnames(ht_obj$phenotable2)[-c(1,2)])
  })
  observe({
    updateSelectInput(session,inputId = "SNPs",choices = ht_obj$geno_num$SNP)
  })
  observeEvent(
    input$confirm2,
    {
      DA_obj$group_table= list()
      DA_obj$g_tbl = list()
      DA_obj$da_test = list()
      DA_obj$sig_tag = list()
      DA_obj$fig_type = list()
      DA_obj$trait = list()
      DA_obj$SNPs = list()
      progress_update2 = c(
        "Step1. import group data",
        "Step2. confirm parameters!"
      )
      withProgress(message = "Update",value = 0,
                   expr = {
                     for (i in 1:2) {
                       incProgress(1/2,progress_update2[i])
                       if(i == 1) {
                         if(is.null(groupFile())) {return()}
                         DA_obj$g_tbl <- groupFile() %>% as.data.frame()
                       } else if ( i == 2) {
                         DA_obj$da_test = as.character(input$M_test)
                         DA_obj$sig_tag = as.character(input$sig_tag)
                         DA_obj$fig_type = as.character(input$fig_type)
                         DA_obj$trait = as.character(input$trait)
                         DA_obj$snps = as.character(input$snps)
                       }
                     }
                   }
      )
    }
  )

  observeEvent(
    input$showfig,
    {
      DA_obj$hap_res= list()
      progress_update3 = c(
        "Step1. Start",
        "Step2. Finish!"
      )
      withProgress(message = "showFig",value = 0,
                   expr = {
                     for (i in 1:2) {
                       incProgress(1/2,progress_update3[i])
                       if(i == 1) {
                         if(is.null(groupFile())) {
                           ht_obj$input_tmp = ht_obj$phenotable2 %>% select(Sample,haplotype,DA_obj$trait) %>%
                             setNames(c("Sample","haplotype","trait"))
                           DA_obj$hap_res <- get_da_result(
                             input_tbl = ht_obj$input_tmp,
                             DA_method = DA_obj$da_test,
                             fig_type = DA_obj$fig_type,
                             sig_tag = DA_obj$sig_tag
                           )
                         } else {
                           ht_obj$input_tmp = ht_obj$phenotable2 %>% select(Sample,haplotype,DA_obj$trait) %>%
                             setNames(c("Sample","haplotype","trait"))
                           ht_obj$input_tbl = reform_region_haptbl(x = ht_obj$input_tmp,tbl.group = DA_obj$g_tbl)
                           DA_obj$hap_res <- get_da_result(
                             input_tbl = ht_obj$input_tbl,
                             DA_method = DA_obj$da_test,
                             fig_type = DA_obj$fig_type,
                             sig_tag = DA_obj$sig_tag
                           )
                         }
                       }
                     }
                   }
      )
      DA_obj$plt_s <- DA_obj$hap_res$plt
    }
  )
  output$DA_tbl = DT::renderDataTable({
    input$showfig
    if(is.null(DA_obj$hap_res)){return()}
    DA_obj$hap_res$DA_result %>% as.data.frame()
  })
  output$DA_tblp = DT::renderDataTable({
    input$showfig
    if(is.null(DA_obj$hap_res)){return()}
    DA_obj$hap_res$tbl %>% as.data.frame()
  })
  output$DA_plt = renderPlot({
    input$showfig
    if(is.null(DA_obj$hap_res)){return()}
    DA_obj$plt_s
  })

  ##> SNPs
  snp_obj <- reactiveValues(data=NULL)
  observeEvent(
    input$get_snp,
    {
      if(is.null(DA_obj$hap_res)) {return()}
      snp_obj$SNPs  = as.character(input$SNPs)
      snp_obj$heter = as.logical(input$heter)
      snp_obj$snps_tbl <- get_snps(
        geno_tbl = ht_obj$geno_reform$geno_clean2,
        geno_num = ht_obj$geno_reform$geno_num,
        pheno = DA_obj$hap_res$tbl,
        snps = snp_obj$SNPs,
        heter = snp_obj$heter
      )
    }
  )

  observeEvent(
    input$showfig2,
    {
      snp_obj$da_res = get_da_result(
        input_tbl = snp_obj$snps_tbl$clean,
        DA_method = DA_obj$da_test,
        fig_type = DA_obj$fig_type,
        sig_tag = DA_obj$sig_tag
      )
    }
  )

  output$checkSNPs = DT::renderDataTable({
    input$get_snps
    if(is.null(DA_obj$hap_res)) {return()}
    snp_obj$snps_tbl$geno_tbl2
  })

  output$DA_tbl2 = DT::renderDataTable({
    input$showfig2
    if(is.null(snp_obj$da_res)) {return()}
    snp_obj$da_res$DA_result
  })
  output$DA_plt2 = renderPlot({
    input$showfig2
    if(is.null(snp_obj$da_res)) {return()}
    snp_obj$da_res$plt
  })

  #> download

  download <- reactiveValues(data=NULL)

  observeEvent(
    input$adjust1,
    {
      download$width1 = as.numeric(input$width1)
      download$height1 = as.numeric(input$width1)
    }
  )

  observeEvent(
    input$adjust2,
    {
      download$width2 = as.numeric(input$width2)
      download$height2 = as.numeric(input$width2)
    }
  )

  output$down_tbl_00 = downloadHandler(

    filename = function() {
      if(is.null(ht_obj$phenotable2)){return()}
      "01.trait-haplotype-all.xlsx"
    },
    content = function(file) {
      writexl::write_xlsx(x = ht_obj$phenotable2,path = file)
    }
  )

  output$down_tbl_01 = downloadHandler(

    filename = function() {
      if(is.null(DA_obj$hap_res)){return()}
      paste0("02.",DA_obj$trait,"-haplotype-",DA_obj$da_test,".xlsx")
    },
    content = function(file) {
      writexl::write_xlsx(x = DA_obj$hap_res$DA_result,path = file)
    }
  )

  output$down_tbl_02 = downloadHandler(

    filename = function() {
      if(is.null(DA_obj$hap_res)){return()}
      paste0("04.",DA_obj$trait,"-plotData-",DA_obj$da_test,".xlsx")
    },
    content = function(file) {
      writexl::write_xlsx(x = DA_obj$hap_res$tbl,path = file)
    }
  )

  output$down_tbl_03 = downloadHandler(

    filename = function() {
      if(is.null(snp_obj$snps_tbl)){return()}
      paste0("05.",DA_obj$trait,"-genotype-SelectSNP",".xlsx")
    },
    content = function(file) {
      writexl::write_xlsx(x = snp_obj$snps_tbl$geno_tbl2,path = file)
    }
  )

  output$down_tbl_04 = downloadHandler(

    filename = function() {
      if(is.null(snp_obj$hap_res)){return()}
      paste0("06.",DA_obj$trait,"-haplotype-SelectSNP-DA",".xlsx")
    },
    content = function(file) {
      writexl::write_xlsx(x = DA_obj$hap_res$tbl,path = file)
    }
  )


  output$down_plt_01 = downloadHandler(
    filename = function() {
      paste0("03.",DA_obj$trait,"-haplotype-all-regions.pdf")
    },
    content = function(file) {
      ggsave(filename = file,plot = DA_obj$plt_s,width = download$width1,height = download$height1)
    }
  )

  output$down_plt_02 = downloadHandler(
    filename = function() {
      paste0("07.",DA_obj$trait,"-haplotype-selectSNPs.pdf")
    },
    content = function(file) {
      ggsave(filename = file,plot = snp_obj$da_res$plt,width = download$width2,height = download$height2)
    }
  )

  ##> heatmap plot
  expmatFile <- reactive({
    file5 <- input$ExpMat
    if(is.null(file6)){return()}
    read.delim(file = file6$datapath,header = T)
  })
  snpAnnoFile <- reactive({
    file5 <- input$SnpAnno
    if(is.null(file7)){return()}
    read.delim(file = file7$datapath,header = T)
  })
  p3_obj <- reactiveValues(data=NULL)
  ##>
  observeEvent(
    input$status_check,
    {
      if(is.null(expmatFile())) {return()}
      if(is.null(snpAnnoFile())) {return()} else {ps_obj$snpAnno = snpAnnoFile()}
      p3_obj$expmat = expmatFile()
      ps_obj$regions = input$regions
      ps_obj$snpAnno_extract = list()
      output$ExpmatCheck = renderUI({
        progress_status = c(
          "Step1. Summary of Heatmap database.",
          "Step2. Summary of SNP annotation file",
          "Step3. Extract target region from SNP annotation",
          "Step4. Extract genes expression data by target region"
        )
        withProgress(
          message = "Check Status",value = 0,
          expr = {
            for (i in 1:4) {
              incProgress(1/4,progress_status[i])
              if(i == 1) {

              }
            }
          }
        )
      })
    }
  )



}

shinyApp(ui,server,
         options = list(launch.browser = TRUE)
)
