#!/bin/sh

#SBATCH -J grna
#SBATCH  -p cu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

#python meta_gRNAs.py meta --genomedir /public/home/wangycgroup/wuj/MyDb/Genome/pan_bacteria/genomes/ncbi_dataset/data
python meta_gRNAs.py offtarget --genome /public/home/wangycgroup/wuj/MyDb/Genome/Hsa/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa --bowtie_index /public/home/wangycgroup/wuj/MyDb/Genome/Hsa/GRCh38/bowtie_idx/Homo_sapiens.GRCh38.dna.primary_assembly.fa --gtf /public/home/wangycgroup/public/00_Genome_ref/Homo_sapiens/Homo_sapiens.GRCh38.105.gtf --grnas bm.grna.fa
