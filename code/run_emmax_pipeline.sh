#!/bin/bash

############################################################
#       Prj: Emmax pipeline
#       Assignment: main pipeline
#       Author: Shawn Wang
#       Date: Apr 26 2023
############################################################

## getopts
set -e ## 报错打断，防止一直错下去

shopt -s expand_aliases
source ~/.bash_alias

start_time=$(date +%s)

## 帮助内容
func(){
    echo -e "\033[32m\n-------------------------------\n\033[0m"
    echo -e "\033[32;1mEmmax pipeline . From raw to Manhattan plot\033[0m"
    echo -e "\033[32m\n-------------------------------\n\033[0m"
    echo -e "\033[32;1mUsage: \n\033[0m"
    echo -e "\033[35mrun_emmax_pipline \ \n \033[31m-t [tpedf_prefix] \ \n \033[31m-o [out_prefix] \ \n \033[31m-p [phenof] \ \n \033[31m-k [kinf] \ \n \033[31m-a [anno] \ \n \033[33m-c [covf] (option) \ \n \033[33m-i [img_type] (option) \ \n \033[33m-s [point_size] (option) \ \n \033[33m-w [point_color] (option) \ \n \033[33m-r [rerun] (option) \ \n\033[0m"
    echo -e "\033[32m\n-------------------------------\n\033[0m"
    echo -e "\033[32mAuthor\033[0m Shawn Wang (shawnwang2016@126.com)"
    echo -e "\033[32mDate\033[0m Web Apr 26, 2023"
    echo -e "\033[32mVersion\033[0m V.0.0.0.99 beta"
    echo -e "\033[32m\n-------------------------------\n\033[0m"
    echo -e "\033[32;1mRequired parameters:\n\033[0m"
    echo -e "\033[31;3m[-t]:tpedf_prefix\033[0m,   prefix for tped/tfam files\033[0m"
    echo -e "\033[31;3m[-o]:out_prefix\033[0m,  output file name prefix\033[0m"
    echo -e "\033[31;3m[-p]:phenof\033[0m,  3-column phenotype file with FAMID, INDID at the first two colmns, in the same order of .tfam file.\033[0m"
    echo -e "\033[31;3m[-k]:kinf\033[0m,  n * n matrix containing kinship values in the individual order consistent to [tpedf].tfam file. [tpedf].kinf will be used if not specified.\033[0m"
    echo -e "\033[31;3m[-a]:anno\033[0m,  SNP annotation file by vep.\n\033[0m"
    echo -e "\033[32;1mOptional parameters:\n\033[0m"
    echo -e "\033[33;3m[-c]:covf\033[0m,  multi-column covariate file with FAMID, INDID at the first two colmns, in the same order of .tfam fileOptional parameters.\033[0m"
    echo -e "\033[33;3m[-i]:img_type\033[0m,  output img file type, jpg or pdf.\033[0m"
    echo -e "\033[33;3m[-s]:point_size\033[0m,  point size of manhattan plot, default: cut_off 0_4_6 ==> 0.3_0.5_0.5.\033[0m"
    echo -e "\033[33;3m[-w]:point_color\033[0m,  point color of manhattan plot, default: cut_off 4_6 ==> green_red.\033[0m"
    echo -e "\033[33;3m[-r]:rerun\033[0m,  If you need to rerun, please use this parameter. 1, 2 or 12\033[0m"
    echo -e "\033[32m\n-------------------------------\n\033[0m"
    exit -1
}


## 默认值
while getopts 'ht:o:p:k:a:c:i::s::w::r::' OPT;
do
    case $OPT in
        t) tpedf_prefix=`echo "$OPTARG"`;;
        o) out_prefix=`echo "$OPTARG"`;;
        p) phenof=`echo "$OPTARG"`;;
        k) kinf=`echo "$OPTARG"`;;
        a) anno=`echo "$OPTARG"`;;
        c) covf=`echo "$OPTARG"`;;
        i) img_type=`echo "$OPTARG"`;;
        s) point_size=`echo "$OPTARG"`;;
        w) point_color=`echo "$OPTARG"`;;
        r) rerun=`echo "$OPTARG"`;;
        h) func
           exit 1
           ;;
        ?) func
           exit 1
           ;;
    esac
done

## 设置默认参数
if [ -z "$tpedf_prefix" ]; then echo -e "\033[31m\nERROR:\033[0m need [-t]:tpedf_prefix parameter,   prefix for tped/tfam files\033[0m"; exit 1; fi
if [ -z "$out_prefix" ]; then echo -e "\033[31m\nERROR:\033[0m need [-o]:out_prefix parameter,  output file name prefix\033[0m"; exit 1; fi
if [ -z "$phenof" ]; then echo -e "\033[31m\nERROR:\033[0m need [-p]:phenof parameter, need 3-column phenotype file with FAMID, INDID at the first two colmns, in the same order of .tfam file.\033[0m"; exit 1; fi
if [ -z "$kinf" ]; then echo -e "\033[31m\nERROR:\033[0m need [-k]:kinf parameter,  n * n matrix containing kinship values in the individual order consistent to [tpedf].tfam file. [tpedf].kinf will be used if not specified.\033[0m"; exit 1; fi
if [ -z "$anno" ]; then echo -e "\033[31m\nERROR:\033[0m need [-a]:anno parameter,  SNP annotation file by vep.\033[0m"; exit 1; fi
if [ -z "$covf" ]; then echo -e "\033[36m\nWarning:\033[0m Not detected [-c]:covf parameter,  Without covariate files, we will not use covariates when performing EMMAX.\033[0m"; fi
if [ -z "$img_type" ]; then img_type='jpg'; fi
if [ -z "$point_size" ]; then point_size='0.3_0.5_0.5'; fi
if [ -z "$point_color" ]; then point_color='green_red'; fi
if [ -z "$rerun" ]; then rerun='3'; fi
## 提示你的参数设置
echo -e "\033[32m\n=====================\033[1m \033[34mYour setting \033[0m \033[32m=====================\033[0m"
echo -e "\033[32mThe prefix of input transposed plink file: \033[0m ${tpedf_prefix}" 
echo -e "\033[32mThe prefix of output file :\033[0m  ${out_prefix}"
echo -e "\033[32mThe phenotype file:\033[0m  ${phenof}"
echo -e "\033[32mThe kinship file:\033[0m  ${kinf}"
echo -e "\033[32mThe SNP annotation file:\033[0m  ${anno}"
echo -e "\033[32mThe covariate file:\033[0m  ${covf:-not set}"
echo -e "\033[32mThe output img type:\033[0m  ${img_type}"
echo -e "\033[32mThe point size of 0-4-6:\033[0m  ${point_size}"
echo -e "\033[32mThe point color of 4-6:\033[0m  ${point_color}"
echo -e "\033[32m=========================================================\033[0m"

## step1 emmax main

echo -e "\033[32m\n==================================================== \033[0m\033[34m\nStep1 run emmax .... \033[0m\n\033[32m====================================================\n\033[0m"

mkdir -p work_dir

if [ -f work_dir/${out_prefix}.ps ] || [ -f work_dir/${out_prefix}_result/${out_prefix}.ps ]
then
    if [ $rerun == "1" ] || [ $rerun == "12" ]
    then
        echo -e "\033[34mRe-run Emmax analysis! \033[0m"
        if [ -z ${covf} ]
        then
            echo -e "\033[34mStep1. Run emmax without covariate .... \n other parameters: -v -d 10 \033[0m"
            emmax -v -d 10 -t ${tpedf_prefix} -o work_dir/${out_prefix} -p ${phenof} -k ${kinf} >/dev/null 2>&1
        else
            echo -e "\033[34mStep1. Run emmax with covariate .... \n other parameters: -v -d 10 \033[0m"
            emmax -v -d 10 -t ${tpedf_prefix} -o work_dir/${out_prefix} -p ${phenof} -k ${kinf} -c ${covf} >/dev/null 2>&1
        fi
        echo -e "\033[34mEmmax analysis finish! Please check your result at: work_dir/${out_prefix}.ps \033[0m"
    else
        echo -e "\033[34mSkip step1. Using previous emmax result .... \033[0m"
    fi
else
    if [ -z ${covf} ]
    then
        echo -e "\033[34mStep1. Run emmax without covariate .... \n other parameters: -v -d 10 \033[0m"
        emmax -v -d 10 -t ${tpedf_prefix} -o work_dir/${out_prefix} -p ${phenof} -k ${kinf}  >/dev/null 2>&1
    else
        echo -e "\033[34mStep1. Run emmax with covariate .... \nother parameters: -v -d 10 \033[0m"
        emmax -v -d 10 -t ${tpedf_prefix} -o work_dir/${out_prefix} -p ${phenof} -k ${kinf} -c ${covf} >/dev/null 2>&1
    fi
    echo -e "\033[34mEmmax analysis finish! Please check your result at: work_dir/${out_prefix}.ps \033[0m"
fi


## 画图
echo -e "\033[32m\n==================================================== \n\033[0m\033[34mStep2 generate QQ and Manhattan plot .... \033[0m\n\033[32m====================================================\n\033[0m"

cd work_dir

if [ -d ${out_prefix}_result ]
then
    if [ $rerun == "2" ] || [ $rerun == "12" ]
    then
        if [ -f ${out_prefix}_result/${out_prefix}.ps ]; then mv ${out_prefix}_result/${out_prefix}.ps ./; fi
        emmax-down -e ${out_prefix}.ps -t ${out_prefix} -a ../${anno} -o ${out_prefix}_result -i ${img_type} -s ${point_size} -w ${point_color}
        echo -e "\033[34mEmmax analysis finish! Please check your result at: work_dir/${out_prefix}_result/ \033[0m"  
        mv ${out_prefix}.* ${out_prefix}_result/
    else
        echo -e "\033[34mSkip step2. Keep previous plots .... \033[0m"
    fi
else
    if [ -f ${out_prefix}_result/${out_prefix}.ps ]; then mv ${out_prefix}_result/${out_prefix}.ps ./; fi
    emmax-down -e ${out_prefix}.ps -t ${out_prefix} -a ../${anno} -o ${out_prefix}_result -i ${img_type} -s ${point_size} -w ${point_color}
    echo -e "\033[34mEmmax analysis finish! Please check your result at: work_dir/${out_prefix}_result/ \033[0m"
    mv ${out_prefix}.* ${out_prefix}_result/
fi

##> check running time
end_time=$(date +%s)
run_time=$((end_time - start_time))


hours=$((run_time / 3600))
minutes=$(((run_time / 60) % 60))
seconds=$((run_time % 60))

printf "\033[34mThe script took %d hours, %d minutes, and %d seconds to run.\033[0m\n" $hours $minutes $seconds