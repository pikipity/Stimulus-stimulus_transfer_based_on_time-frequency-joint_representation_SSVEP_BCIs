#!/bin/sh
cd ./data

for ((i=1;i<36;i++)); do
    wget http://bci.med.tsinghua.edu.cn/upload/yijun/S$i.mat.7z
    7z e S$i.mat.7z
done