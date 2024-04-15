#!/bin/bash

# 创建输出文件夹（如果不存在）
mkdir -p data_out

# 数据文件夹
data_folder="data_in"

# 遍历数据文件夹
for txt_file in "$data_folder"/*.edges; do
    # 提取文件名（不含路径和扩展名）
    filename=$(basename -- "$txt_file")
    filename_noext="${filename%.*}"

    # 遍历不同的 k 值
    for k_value in {1..3}; do
        # 输出文件路径
        output_file="data_out/${filename_noext}_k${k_value}_output.txt"

        # 运行程序
        ./main "$txt_file" "$k_value" > "$output_file"

        echo "Program executed with $txt_file and k=$k_value. Output saved to $output_file"
    done
done
