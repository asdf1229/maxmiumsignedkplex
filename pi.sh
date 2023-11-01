#!/bin/bash

# 文件列表
files=("dataset/1.edges" "dataset/2.edges" "dataset/3.edges" \
"dataset/4.edges" "dataset/5.edges" "dataset/6.edges")

# 循环遍历文件列表
for file in "${files[@]}"
do
    # 创建唯一的文件名
    output_file="${file}_output.txt"
    
    # 执行可执行文件，并将文件作为参数传递，并将输出写入文件
    5/skce_ok "$file" 3 > "$output_file"
done