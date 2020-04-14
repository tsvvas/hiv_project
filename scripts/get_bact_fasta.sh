#!/bin/bash

wget -O "CMV_AD169.fasta" "https://www.uniprot.org/uniprot/?include=false&format=fasta&force=true&query=proteome:UP000008991"
wget -O "EBV_AG876.fasta" "https://www.uniprot.org/uniprot/?include=false&format=fasta&force=true&query=proteome:UP000007639"
wget -O "HIV1_HXB2.fasta" "https://www.uniprot.org/uniprot/?include=false&format=fasta&force=true&query=proteome:UP000002241"
wget -O "YFV_17D.fasta" "https://www.uniprot.org/uniprot/?include=false&format=fasta&force=true&query=proteome:UP000000360"
wget -O "HCV_IsoH.fasta" "https://www.uniprot.org/uniprot/?include=false&format=fasta&force=true&query=proteome:UP000000518"
wget -O "CLMD_trachomatis.fasta" "https://www.uniprot.org/uniprot/?include=false&format=fasta&force=true&query=proteome:UP000000431"
wget -O "YERS_pestis.fasta" "https://www.uniprot.org/uniprot/?include=false&format=fasta&force=true&query=proteome:UP000000815"
wget -O "SLML_typhimurium.fasta" "https://www.uniprot.org/uniprot/?include=false&format=fasta&force=true&query=proteome:UP000001014"
wget -O "SHGL_flexneri.fasta" "https://www.uniprot.org/uniprot/?include=false&format=fasta&force=true&query=proteome:UP000001006"
wget -O "CPBT_jejuni.fasta" "https://www.uniprot.org/uniprot/?include=false&format=fasta&force=true&query=proteome:UP000000799"
wget -O "CLOS_difficile.fasta" "https://www.uniprot.org/uniprot/?include=false&format=fasta&force=true&query=proteome:UP000001978"
wget -O "MYPL_pneumoniae.fasta" "https://www.uniprot.org/uniprot/?include=false&format=fasta&force=true&query=proteome:UP000000808"
wget -O "MYBT_smegmatis.fasta" "https://www.uniprot.org/uniprot/?include=false&format=fasta&force=true&query=proteome:UP000000757"
wget -O "KLEB_pneumoniae.fasta" "https://www.uniprot.org/uniprot/?include=false&format=fasta&force=true&query=proteome:UP000000265"
wget -O "MYPL_synoviae.fasta" "https://www.uniprot.org/uniprot/?include=false&format=fasta&force=true&query=proteome:UP000000549"
wget -O "human_proteome.fasta" "https://www.uniprot.org/uniprot/?include=false&format=fasta&force=true&query=proteome:UP000005640"

mkdir -p data/fasta/
mv *fasta data/fasta/
