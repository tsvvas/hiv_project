#!/bin/bash

echo "Hello! I am downloading data. Please wait."
PREFIX=$1
BEGIN=$2
END=$3

for (( VAR=$BEGIN; VAR<=$END; VAR++ ))
do
	touch $PREFIX$VAR.gb
	efetch -db nucleotide -format gb -id $PREFIX$VAR > $PREFIX$VAR.gb
done
