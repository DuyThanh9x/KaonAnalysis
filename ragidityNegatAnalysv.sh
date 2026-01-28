#!/bin/bash
input_file4="particleID400negat.txt"
input_file7="particleID700negat.txt"
input_fileList="listruns.txt"
path_file="25/"
listrun=""
#echo "Reading particleID from ${input_file7}"
echo "Reading list from ${input_fileList}"
while IFS=' ' read -r filename; do
        if [[ "$filename" == //* ]]; then
        continue
    fi
    listrun+=\"$path_file$filename\",
done < "$input_fileList"
    echo "test $listrun"
root -l -q 'macros/ragidityVsTransverseMment400negat.cc({'$listrun'},"'$input_file4'")'
#root -l -q 'macros/ragidityVsTransverseMment700negat.cc({'$listrun'},"'$input_file7'")'
echo "Done"
