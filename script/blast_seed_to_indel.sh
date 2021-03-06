#run in 05_ZhuWC directory
mkdir check_seed
cd check_seed
##blast seed to kinamaze genome
blat /home/lpeng/Gmatic5/genome/Rice_Kinamaze/Kinamaze_genome.fa seed_genes.fasta seed_gene_to_kinamaze.tsv
seqkit fx2tab -l seed_genes.fasta |cut -f 1,4
less -S seed_gene_to_kinamaze.tsv |cut -f 14,16,17,9,10|sed  '1,5d'|awk -v OFS="\t" '{print $3,$4,$5,$1,$2}'>OsRacGEF_gene_chr.bed
bedtools slop -i OsRacGEF_gene_chr.bed -g /home/lpeng/Gmatic5/genome/Rice_Kinamaze/Kinamaze_genome.size -l 3000 -r 3000 >OsRacGEF_gene_chr_3000.bed
less ../vcf/104-1-90-2_combined.flt.vcf|grep -v '#'|awk -v OFS="\t" '{print $1,$2-1,$2,$5}' >../vcf/104-1-90-2_combined.flt.bed
less ../vcf/104-3-1-90-2_combined.flt.vcf|grep -v '#'|awk -v OFS="\t" '{print $1,$2-1,$2,$5}' >../vcf/104-3-1-90-2_combined.flt.bed
bedtools intersect -a OsRacGEF_gene_chr_3000.bed -b ../vcf/104-1-90-2_combined.flt.bed -wa -wb -nonamecheck >104-1-90-2_combined.flt.indel_and_snp_in_OsRacGEF_gene.txt
bedtools intersect -a OsRacGEF_gene_chr_3000.bed -b ../vcf/104-3-1-90-2_combined.flt.bed -wa -wb -nonamecheck >104-3-1-90-2_combined.flt.indel_and_snp_in_OsRacGEF_gene.txt

###
seqkit sliding -s 20 -W 20 outer_seq.fasta |seqkit seq -u>outer_seq_seeds.fasta
