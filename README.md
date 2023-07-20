Emmax-pipeline

<center>

<font color=salmon size=2pt> 🏡Zhang Lab, Jun 04, 2023 🏡</font>

<font color=green size=2pt> 🔬Multi-omics research center🔬</font>

<font color=green size=2pt> 🔬State Key Laboratory of Crop Stress Adaption and Improvement, HENU, Kaifeng, Henan 🔬</font>

</center>

![Overview](Emmax-pipeline.png)

```r
suppressMessages(library(bruceR))
suppressMessages(library(PCAtools))
suppressMessages(library(tidymass))
suppressMessages(library(tidyverse))
suppressMessages(library(plotly))
```

# 环境配置和软件安装

本流程需要在linux/Unix(MacOS)或者windows WSL环境下配置。建议通过`conda`完成环境搭建及软件安装。下面将讲解详细步骤。

## conda下搭建emmax-gwas分析环境

### 安装conda

首先到清华镜像站下载miniconda安装包：[index of /annconda/miniconda](https://mirrors.tuna.tsinghua.edu.cn/anaconda/miniconda/?C=M&O=D)，根据自己操作系统类型选择**最新日期**对应的安装包。

1.1 通过wget将Miniconda安装包下载到服务器； 1.2 添加执行权限，运行安装脚本，结束的时候注意不要选择随系统启动，以防conda污染环境变量 1.3 创建一个alias快捷启动conda 1.4 替换软件镜像

代码如下：

```bash
#| code-fold: true
#| code-summary: "Show the code"
#> 下载到电脑
wget https://mirrors.tuna.tsinghua.edu.cn/anaconda/miniconda/Miniconda3-py39_23.3.1-0-Linux-x86_64.sh

#> --2023-06-14 09:14:37--  https://mirrors.tuna.tsinghua.edu.cn/anaconda/miniconda/Miniconda3-py39_23.3.1-0-Linux-x86_64.sh
#> Resolving mirrors.tuna.tsinghua.edu.cn (mirrors.tuna.tsinghua.edu.cn)... 101.6.15.130
#> Connecting to mirrors.tuna.tsinghua.edu.cn (mirrors.tuna.tsinghua.edu.cn)|101.6.15.130|:443... connected.
#> HTTP request sent, awaiting response... 200 OK
#> Length: 63969016 (61M) [application/octet-stream]
#> Saving to: ‘Miniconda3-py39_23.3.1-0-Linux-x86_64.sh’

#> Miniconda3-py39_23.3.1-0-Linux-aa 100%[==========================================================>]  61.00M  53.7MB/s    in 1.1s 2023-06-14 09:14:38 (53.7 MB/s) - ‘
#> Miniconda3-py39_23.3.1-0-Linux-x86_64.sh’ saved [63969016/63969016]
#> 添加可执行权限
chmod +x Miniconda3-py39_23.3.1-0-Linux-x86_64.sh
#> 执行安装脚本
bash Miniconda3-py39_23.3.1-0-Linux-x86_64.sh
#> 下面一路yes+回车就可以了 就是在最后看到是否开机启动的时候选择no
#> Do you wish the installer to initialize Miniconda3
#> by running conda init? [yes|no]
#> [no] >>> no
#> 新建一个放置alias的配置文件
touch ~/.bash_alias
#> 将快捷启动conda的conact写入刚创建的~/.bash_alias
echo 'alias conact="source ~/miniconda3/bin/activate"' >> ~/.bash_alias
#> 将source ~/.bash_alias写到用户配置文件下，每次启动shell bash_alias中的快捷命令自动生效
echo 'source ~/.bash_alias' >> ~/.bashrc\
#> 重新加载下配置文件
source ~/.bash_alias
#> 利用conact启动conda
conact
#> 更换镜像源
vim ~/.condarc
#> 全程英文输入模式！！！
#> 键盘输入i进入编辑模式，左下角能看到状态变为insert
#> 将下列代码粘贴进去
channels:
  - defaults
show_channel_urls: true
default_channels:
  - https://mirrors.bfsu.edu.cn/anaconda/pkgs/main
  - https://mirrors.bfsu.edu.cn/anaconda/pkgs/r
  - https://mirrors.bfsu.edu.cn/anaconda/pkgs/msys2
custom_channels:
  conda-forge: https://mirrors.bfsu.edu.cn/anaconda/cloud
  msys2: https://mirrors.bfsu.edu.cn/anaconda/cloud
  bioconda: https://mirrors.bfsu.edu.cn/anaconda/cloud
  menpo: https://mirrors.bfsu.edu.cn/anaconda/cloud
  pytorch: https://mirrors.bfsu.edu.cn/anaconda/cloud
  pytorch-lts: https://mirrors.bfsu.edu.cn/anaconda/cloud
  simpleitk: https://mirrors.bfsu.edu.cn/anaconda/cloud
#> 按ESC退出编辑模式,然后按：号进入命令执行模式，做下角光标闪烁，输入wq回车
conda clean -i #清理index cache
```

### 新建一个emmax分析环境

做分析需要下面的软件：

| 软件      | 作用                                                                    |
|------------------------------------|------------------------------------|
| bcftools  | 操作vcf文件，包括修改染色体prefix，按区间提取，构建索引，过滤snp等      |
| vcftools  | 和bcftools类似                                                          |
| plink     | 将vcf转换为plink格式，然后进行vcf文件操作，emmax要求输入文件格式为plink |
| tassel    | 基因型文件转换格式，pca，生成phylo文件                                  |
| fasttree  | 目标区间快速建立NJ树                                                    |
| admixture | 群体结构分析                                                            |
| beagle    | 缺失基因型填补                                                          |
| tabix     | vcf建立索引                                                             |
| gatk      | snp calling                                                             |
| samtools  | 操作sam，bam文件 对可以结构变异切割bam文件通过IGV验证                   |
| vep       | SNP文件注释                                                             |

代码如下：

```bash
#| code-fold: true
#| code-summary: "Show the code"
conda create -c bioconda -n gwas bcftools vcftools plink tassel admixture beagle fasttree tabix samtools gatk -y
#> 安装成功后可以看到以下提醒
# To activate this environment, use
#
#     $ conda activate gwas
#
# To deactivate an active environment, use
#
#     $ conda deactivate
#> 单独给vep创建一个环境
conda create -c bioconda -n ensembl-vep vep -y
```

### 安装emmax

emmax有可以直接在unbuntu20.04编译好的软件 直接下载，解压使用，不过为了方便，我们还是需要些两个alias方便使用emmax和emmax-kin，或者直接把他们加入环境变量。

代码：

```bash
#| code-fold: true
#| code-summary: "Show the code"
wget http://csg.sph.umich.edu//kang/emmax/download/emmax-beta-07Mar2010.tar.gz
# --2023-06-14 11:05:10--  http://csg.sph.umich.edu//kang/emmax/download/emmax-beta-07Mar2010.tar.gz
# Resolving csg.sph.umich.edu (csg.sph.umich.edu)... 141.211.55.24
# Connecting to csg.sph.umich.edu (csg.sph.umich.edu)|141.211.55.24|:80... connected.
# HTTP request sent, awaiting response... 200 OK
# Length: 519928 (508K) [application/x-gzip]
# Saving to: ‘emmax-beta-07Mar2010.tar.gz’

# emmax-beta-07Mar2010.tar.gz          100%[=====================================================================>] 507.74K   416KB/s    in 1.2s

# 2023-06-14 11:05:13 (416 KB/s) - ‘emmax-beta-07Mar2010.tar.gz’ saved [519928/519928]
tar -xvf emmax-beta-07Mar2010.tar.gz
#> 将emmax文件夹放到家目录的software路径下
mkdir ~/software
mv emmax-beta-07Mar2010 ~/software/
#> 创建alias
echo 'alias emmax="~/software/emmax-beta-07Mar2010/emmax"' >> ~/.bash_alias
echo 'alias emmax-kin="~/software/emmax-beta-07Mar2010/emmax-kin"' >> ~/.bash_alias
source ~/.bash_alias
#> 测试
emmax
#> 出现下面提示说明安装完毕
# Usage: emmax [options]
# Required parameters
# 	-t [tpedf_prefix] : prefix for tped/tfam files
# 	-o [out_prefix]  : output file name prefix
# Likely essential parameters
# 	-p [phenof] : 3-column phenotype file with FAMID, INDID at the first two colmns, in the same order of .tfam file. Not required only with -K option	-k [kinf] : n * n matrix containing kinship values in the individual order consistent to [tpedf].tfam file. [tpedf].kinf will be used if not specified
# 	-c [covf] : multi-column covariate file with FAMID, INDID at the first two colmns, in the same order of .tfam fileOptional parameters
# Optional parameters
# 	-i [in_prefix] : input file name prefix including eigenvectors
# 	-d [# digits]  : precision of the output values (default : 5)
# 	-s [start index of SNP] : start index of SNP (default : 0)
# 	-e [end index of SNP] : end index of SNP (default : #snps)
# 	-w : flag for writing eigenvalues/eigenvector files
# 	-D [delimiters] : delimter string in quotation marks
# 	-P [# heaer cols in tped] : # of column headers in tped file
# 	-F [# heaer cols in tfam] : # of column headers in tfam file
# Segmentation fault (core dumped)
```

### 下载脚本文件

可以通过下载压缩包和通过git等方式将仓库代码下载到服务器，然后给脚本添加alias..

| 步骤                       | 脚本名称              | alias              | 作用                                       |
|------------------|------------------|------------------|------------------|
| 协变量文件格式化           | MakeEmmaxCov.R        | emmax-cov          | 修改协变量文件格式                         |
| 表型数据文件格式化         | split_phenotype.R     | emmax-pheno-split  | 表型分割                                   |
| 开始关联分析及下游文件生成 | run_emmax_pipeline.sh | run_emmax_pipeline | 数据准备好后执行该脚本会将下游分析一键完成 |
| 开始关联分析及下游文件生成 | EMMAX_DownStream.R    | emmax-down         | 已经被集成到run_emmax-pipeline             |
| 区域单倍型划分             | emmax-haplotype.sh    | emmax-haplotype    | 根据候选区间进行单倍型划分                 |
| 区域单倍型热图             | emmax-haploheatmap.R  | emmax-haploheatmap | 倍集成到emmax-haplotype                    |

代码：

```bash
#| code-fold: true
#| code-summary: "Show the code"
cd # 返回家目录
mkdir src && cd src # 创建一个src路径 然后进入该文件夹
git clone https://github.com/ShawnWx2019/emmax-pipeline.git # 通过git将emmax-pipeline下载到本地
echo "alias emmax-pca='Rscript ~/src/emmax-pipeline/code/PCA_var_proportion.R'" >> ~/.bash_alias
echo "alias emmax-cov='Rscript ~/src/emmax-pipeline/code/MakeEmmaxCov.R'" >> ~/.bash_alias
echo "alias emmax-down='Rscript ~/src/emmax-pipeline/code/EMMAX_DownStream.R'" >> ~/.bash_alias
echo "alias emmax-pheno-split='Rscript ~/src/emmax-pipeline/code/split_phenotype.R'" >> ~/.bash_alias
echo "alias run_emmax_pipeline='bash ~/src/emmax-pipeline/code/run_emmax_pipeline.sh'" >> ~/.bash_alias
echo "alias emmax-haplotype='bash ~/src/emmax-pipeline/code/emmax-haplotype.sh'" >> ~/.bash_alias
echo "alias emmax-haploheatmap='Rscript ~/src/emmax-pipeline/code/emmax-haploheatmap.R'" >> ~/.bash_alias
source ~/.bash_alias
```

至此，上游分析环境搭建工作完毕。

# Emmax上游分析

## 基因型文件过滤

首先最原始的基因型文件应该是经过GATK call snp后进行质控后的`.vcf`文件, 可以直接使用`vcftools`进行进一步过滤，对于过滤参数的标准可以参考文献中对应物种的过滤方式。 例如在\@duResequencing243Diploid2018 文章中，陆地棉的过滤条件为：

> only high-quality SNPs (coverage depth ≥3 and ≤50, mapping quality ≥20, missing ratio of samples within population ≤10% (1,980,539 SNPs) or 20% (3,665,030 SNPs), and minor allele frequency (MAF) ≥0.05) were retained for subsequent analyses. SNPs with a miss ratio ≤10% were used in PCA/phylogenetictree/structure analyses, whereas SNPs with a miss ratio ≤20% were used in the rest of the analyses.

<font color=green>**第一步:**</font> 首先emmax要求染色体为纯数字，如果你发现GATK生成的vcf文件中染色体为"chr01"等并非纯数字时需要通过bcftools修改一下染色体。

```bash
#| code-fold: true
#| code-summary: "Show the code"
##> GATK 一般不会给snp加ID，这里用bcftools对vcf文件中snp命名，方式为chr:pos
bcftools annotate --set-id '%CHROM:%POS' raw.vcf.gz
#> 手动写一个配置文件，左边为原始的
vim chr_cn.txt
i ## 进入编辑模式
chr1    1
chr2    2
chr3    3
chr4    4
chr5    5

ESC #退出编辑模式
wq # 保存
## 使用bcftools 修改染色体名称
bcftools annotate --rename-chrs chr_cn.txt raw.vcf.gz | bgzip -c >raw_rename.vcf.gz
## 给新生成的vcf文件建立索引
bcftools index -t raw_rename.vcf.gz
```

<font color=green>**第二步:**</font> 用plink对vcf文件过滤，生成的`Gh_383.maf0.05.int0.8.vcf`

```bash
#| code-fold: true
#| code-summary: "Show the code"
plink --vcf raw_rename.vcf.gz --maf 0.05 --geno 0.2 --recode vcf-fid --out Gh_383.maf0.05.int0.8
#> 生成emmax分析所需要的plink文件
plink --vcf Gh_383.maf0.05.int0.8.vcf --recode 12 transpose --output-missing-genotype 0 --out Gh_383.maf0.05.int0.9
```

<font color=green>**第三步:**</font> 根据LD对SNP筛选，构建运行群体结构分析的基因型文件

```bash
#| code-fold: true
#| code-summary: "Show the code"
#> 过滤LD
plink --vcf Gh_383.maf0.05.int0.8.vcf --indep-pairwise 50 10 0.2 --out Gh_383.maf0.05.int0.8
#> 提取筛选SNP id
plink --vcf Gh_383.maf0.05.int0.8.vcf --make-bed --extract Gh_383.maf0.05.int0.8.prune.in --out Gh_383.maf0.05.int0.8.prune.in
#> 将筛选后的in结果转换为vcf
plink --bfile Gh_383.maf0.05.int0.8.prune.in --recode vcf-fid --out Gh_383.maf0.05.int0.8.prune.in
#> 转换为admixture需要的输入文件格式
plink -bfile Gh_383.maf0.05.int0.8.prune.in --recode 12 --out Gh_383.maf0.05.int0.8.prune.in
```

<font color=green>**第四步:**</font> 利用admixture进行群体结构分析.

通过下面代码进行admixture群体结构分析，根据结果选择合适的K值，然后选择对应的群体结构协变量文件

```bash
#| code-fold: true
#| code-summary: "Show the code"
for i in {1..20}
do
admixture --cv Gh_383.maf0.05.int0.8.prune.in.ped ${i} -j48 >> log.txt
done
##> 运行结束后通过下面代码查看最佳分群
$ grep CV log.txt
##> emmax-cov进行格式转换

emmax-cov -n Gh_383.maf0.05.int0.8 \
          -q admixture.Q3
          -o Cov_Q3_emmax.cov
          -k 3
```

<font color=green>**第五步:**</font> 利用emmax-kin生成亲缘关系矩阵

```bash
#| code-fold: true
#| code-summary: "Show the code"
emmax-kin Gh_383.maf0.05.int0.8 -v -d 10
#> 运行完毕Gh_383.maf0.05.int0.8.BN.kinf 文件
```

至此运行emmax的基因型文件，协变量文件已经具备了。

```yaml
(vep) shawn @ bio-Super-Server: 01.raw_g $ tree
.
├── Cov_P3_emmax.cov # 协变量文件 PCA
├── Cov_Q3_emmax.cov # 协变量文件 群体结构
├── Gh_383.maf0.05.int0.8.bed # 基因型文件plink格式 二进制
├── Gh_383.maf0.05.int0.8.bim # 基因型文件plink格式 二进制
├── Gh_383.maf0.05.int0.8.fam # 基因型文件plink格式 二进制
├── Gh_383.maf0.05.int0.8.hBN.kinf #亲缘关系矩阵
├── Gh_383.maf0.05.int0.8.log # 基因型文件plink格式 二进制
├── Gh_383.maf0.05.int0.8.nosex # 基因型文件plink格式
├── Gh_383.maf0.05.int0.8.tfam # 基因型文件plink格式
├── Gh_383.maf0.05.int0.8.tped # 基因型文件plink格式
├── Gh_383.maf0.05.int0.8.vcf # 基因型文件 vcf格式
```

## SNP注释文件准备

通过vep软件对snp文件进行注释

代码：

```bash

# 格式化
grep -v "#" CRI_Gh_v2.gff3 | sort -k1,1 -k4,4n -k5,5n -t$'\t' | bgzip -c > data.gff.gz
# tabix索引
tabix -p gff data.gff.gz
# 注释
vep -i  raw.vcf.gz  --gff data.gff.gz --fasta CRI_Gh_v2.fa --fork 4 --tab -o Gh_383_anno.xls Gh_383_anno.html
```

## 表型数据准备

首先按照下表准备好输入表型文件`traits.txt`

| sample | trait1 | trait2 | trait3 |
|--------|--------|--------|--------|
| s_001  | 212    | 123    | NA     |
| s_002  | 221    | 1223   | 2      |
| s_004  | NA     | 123    | 3      |

然后运行

```bash
emmax-pheon-split -t Gh_383.maf0.05.int0.8.tfam -p traits.txt
```

然后我们就得到了`trait1.txt trait2.txt trait3.txt`可以用于emmax表型分析的文件

至此，所有运行下游分析的文件都准备完毕

## Emmax一步流程

首先我们看下帮助文档

```bash
-------------------------------

Emmax pipeline . From raw to Manhattan plot

-------------------------------

Usage:

run_emmax_pipline \
 -t [tpedf_prefix] \
 -o [out_prefix] \
 -p [phenof] \
 -k [kinf] \
 -a [anno] \
 -c [covf] (option) \
 -i [img_type] (option) \
 -s [point_size] (option) \
 -w [point_color] (option) \
 -r [rerun] (option) \


-------------------------------

Author Shawn Wang (shawnwang2016@126.com)
Date Web Apr 26, 2023
Version V.0.0.0.99 beta

-------------------------------

Required parameters:

[-t]:tpedf_prefix,   prefix for tped/tfam files
[-o]:out_prefix,  output file name prefix
[-p]:phenof,  3-column phenotype file with FAMID, INDID at the first two colmns, in the same order of .tfam file.
[-k]:kinf,  n * n matrix containing kinship values in the individual order consistent to [tpedf].tfam file. [tpedf].kinf will be used if not specified.
[-a]:anno,  SNP annotation file by vep.

Optional parameters:

[-c]:covf,  multi-column covariate file with FAMID, INDID at the first two colmns, in the same order of .tfam fileOptional parameters.
[-i]:img_type,  output img file type, jpg or pdf.
[-s]:point_size,  point size of manhattan plot, default: cut_off 0_4_6 ==> 0.3_0.5_0.5.
[-w]:point_color,  point color of manhattan plot, default: cut_off 4_6 ==> green_red.
[-r]:rerun,  If you need to rerun, please use this parameter. 1, 2 or 12

-------------------------------
```

-   <font color=blue>**必选参数**</font> -t plink 生成的转置的emmax文件 -o 输出文件的前缀 -p 表型文件 -k 亲缘关系矩阵 -a 注释文件
-   <font color=orange>**可选参数**</font> -c 协变量文件 -i 输出图片文件格式 -s 曼哈顿图点的大小 -w 显著位点点的颜色（-logp 4，6） -r 断点重新运行

```bash
run_emmax_pipeline \
  -t Gh_383.maf0.05.int0.8 \
  -o Morin_raw \
  -p Morin.txt \
  -k Gh_383.maf0.05.int0.8.hBN.kinf \
  -a Gh_383_anno.xls \
  -c Cov_Q3_emmax.cov
```

大概等待10分钟左右就可以出结果，基因型文件较小的话更快。

运行完毕得到的文件如下：

```bash
── LDFilebyChr => #我们根据显染色体上显著位点的个数提取了前5个chr的基因型文件数据用来做LDblock
│   ├── A01_genotype.txt
│   ├── A11_genotype.txt
│   └── D07_genotype.txt
├── Morin_Flower_miss.log
├── Morin_Flower_miss.ps => # emmax 输出的包含位点p值的结果文件
├── Morin_Flower_miss.reml
├── Morin_Flower_miss_result_filter.gwas => #格式化后的结果文件，用与使用IGV看结果
├── Morin_Flower_miss_result_filter.xlsx => #Excel文件，包含显著位点的snp信息，包括了p值和snp注释
├── Morin_Flower_miss_result_genotype_all.txt
├── QQplot.Morin_Flower_miss.jpg => #QQ 图
└── Rectangular-Manhattan.Morin_Flower_miss.jpg => #曼哈顿图


```

# Emmax 下游分析

## 确定LDblock画图范围

首先可以通过Morin_Flower_miss_result_filter.xlsx 文件中连续高信号出现染色体位置划分范围。 其次可以通过用IGV看Morin_Flower_miss_result_filter.gwas 划定范围

## 进行LDblock可视化

通过LDblock进行可视化

经过测试，报错的原因主要是染色体名字不符，许多软件在做gwas分析的时候都会改变chr的前缀，导致基因型文件中染色体标志和gff文件中的不符合，此时需要用bcftools重新修改基因型文件。

```bash
LDBlockShow -InGFF  ~/GeekCloud/01.Database/Cotton/AD1/Gh.gff -Cutline 6 -OutPdf -SeleVar 2 -TopSite -InGenotype A11_Gh_383_miss.geno --OutPut Morin_A11_BGLU --InGWAS A11_genotype.txt -Region A11:91108756:91221656 &
```
