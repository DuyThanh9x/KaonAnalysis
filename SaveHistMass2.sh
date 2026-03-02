#!/bin/bash
input_fileList="listruns.txt"
path_file="25/"
listrun=""
echo "Reading list from ${input_fileList}"
while IFS=' ' read -r filename; do
        if [[ "$filename" == //* ]]; then
        continue
    fi
    listrun+=\"$path_file$filename\",
done < "$input_fileList"
    echo "test $listrun"
    #root -l -q 'macros/SaveHistMass2ID400.cc({'$listrun'})'
    #root -l -q 'macros/SaveHistMass2ID400negat.cc({'$listrun'})'
    #root -l -q 'macros/SaveHistMass2ID700.cc({'$listrun'})'
    root -l -q 'macros/SaveHistMass2ID700negat.cc({'$listrun'})'
echo "Done"
