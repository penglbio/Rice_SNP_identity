###check_outer_seq
seqkit sliding -s 30 -W 30 outer_seq.fasta |seqkit seq -u>outer_seq_seeds.fasta
ls ../clean/*_paired.fastq.gz|rush 'seqkit fq2fa {} >./{%:.}.fasta' &
makeblastdb -in 104-1-90-2_combined_R1_paired.fasta -dbtype nucl -out 104-1-90-2_combined_R1 &
makeblastdb -in 104-3-1-90-2_combined_R1_paired.fasta -dbtype nucl -out 104-3-1-90-2_combined_R1 &
makeblastdb -in 104-1-90-2_combined_R2_paired.fasta -dbtype nucl -out 104-1-90-2_combined_R2 &
makeblastdb -in 104-3-1-90-2_combined_R2_paired.fasta -dbtype nucl -out 104-3-1-90-2_combined_R2 &
ls 104*fasta|rush -k 'blastn -query outer_seq_seeds.fasta -db {@(.+)_.*} -outfmt 6 -out outer_to_{@(.+)_pair.*}.tsv' &

#negative control
python random_seq.py outer_seq.fasta outer_random_seq.fasta
seqkit sliding -s 30 -W 30 outer_random_seq.fasta |seqkit seq -u> outer_random_seq_seed.fasta
ls 104*fasta|rush -k 'blastn -query outer_random_seq_seed.fasta -db {@(.+)_.*} -outfmt 6 -out outer_to_{@(.+)_pair.*}.tsv' &
ls *_seed.tsv|rush -k "cat {}|cut -f 1|sort|uniq -c|awk '{print \$2\"\\t\"\$1}' >{@(.+)_seed.tsv}_seed_id_to_read_num.txt" --dry-run
ls *_read_num.txt|rush -k 'less outer_seq_seeds.fasta |seqkit fx2tab|csvtk join -t -H -f 1 - {} >{@outer_to_(.+)_seed_id_to_read_num.txt}_outer_30.tsv'
