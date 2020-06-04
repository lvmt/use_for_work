#!/bin/bash


#根据目的样本，截取数据量
#需要截取数据的
#sample_list,
#目的样本的名称，
#数据目录


declare -A  info_array  # 定义数组

function get_info_array() {

    while read laneid patientid sampleid libid novoid index path;do
        echo ${laneid}
        if [[ ${laneid} ==  "#"* ]]
        then 
            continue
            
        else 
            info_array[${patientid}]=${path}/${libid}  
        fi
    done < sample_list_B43

}

get_info_array

echo ${info_array[@]}

target_sam="T234"

echo "目的样本是: " ${info_array[${target_sam}]}/*1.fq.gz

function get_target_num_read() {

    target_num_read=$(seqkit seq -j8 -n ${info_array[${target_sam}]}/*1.fq.gz | wc  -l)
    # hello="hello world"

}

# get_target_num_read
# echo ${hello}


function cut_data() {

    seqkit head -j8 -n ${target_num_read}  ${infile} -o  {outdir}/${infile}

}
