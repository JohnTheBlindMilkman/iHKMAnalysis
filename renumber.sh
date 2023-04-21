#!/bin/bash

# this macro changes the number of given iHKM output file in order to "get rid of" missing files
# usage: specify $namebegin if running from different path than "./" and execute 

i=0
j=0
currFile=""
newFile=""
files=$(ls -1q *.root | wc -l)
namebegin="./"
nameend=".root"
while [ $j -lt $[files-1] ] # while $j < $files-1
do  
    currFile="${namebegin}${i}${nameend}"
    newFile="${namebegin}${j}${nameend}"
    if [ -f $currFile ] # if i-th file exists
    then
        mv $currFile $newFile
        echo "Moving" $currFile "to" $newFile
        j=$[j+1]        
    fi
    i=$[i+1]
done