#!/bin/bash
input_file4="particleID400.txt"
input_file7="particleID700.txt"
input_fileList="listruns.txt"
path_file="25/"
listrun=""
#echo "Reading particleID from ${input_file4}"
echo "Reading list from ${input_fileList}"
while IFS=' ' read -r filename; do
        if [[ "$filename" == //* ]]; then
        continue
    fi
    listrun+=\"$path_file$filename\",
done < "$input_fileList"
    echo "test $listrun"
root -l -q 'macros/ragidityVsTransverseMment400.cc({'$listrun'},"'$input_file4'")'
#root -l -q 'macros/ragidityVsTransverseMment700.cc({'$listrun'},"'$input_file7'")'
echo "Done"
