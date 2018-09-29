configfile: "config.yaml"

rule all:
	input:
		expand('clean/{sample}_R1_paired.fastq.gz', sample=config['samples']),
		expand('clean/{sample}_R2_paired.fastq.gz', sample=config['samples']),
#		expand('fastqc/raw/{sample}_R1_fastqc.html', sample=config['samples']),
#		expand('fastqc/raw/{sample}_R2_fastqc.html', sample=config['samples']),
		expand('fastqc/clean/{sample}_R1_paired_fastqc.html', sample=config['samples']),
		expand('fastqc/clean/{sample}_R2_paired_fastqc.html', sample=config['samples']),
		expand('stat/fastqc_stat.tsv'),
		expand('bam/{sample}.sort.bam', sample=config['samples']),
		expand('bam/{sample}.rmdup.bam', sample=config['samples']),
		expand('bam/{sample}.bamqc', sample=config['samples']),
#		expand('bam/{sample}.infer_strand_spec', sample=config['samples']),
		expand('stat/bamqc_stat.tsv'),
		expand('vcf/{sample}.raw.vcf',sample=config['samples']),
		expand('vcf/{sample}.flt.vcf',sample=config['samples']),
		expand('vcf/{sample}.flt.gz',sample=config['samples']),
		expand('vcf/{sample}.flt.gz.csi',sample=config['samples']),
		expand('result/{group1}_{group2}',group1=config['group1'],group2=config['group2']),
		['result/{group1}_vs_{group2}'.format(group1=x[0],group2=x[1]) for x in zip(config['group1'],config['group2'])],

rule fastqc_raw_PE:
	input:
		config['path']+'/{sample}_R1.fastq.gz',
		config['path']+'/{sample}_R2.fastq.gz'
	output:
		'fastqc/raw/{sample}_R1_fastqc.html',
		'fastqc/raw/{sample}_R2_fastqc.html'
	shell:
		'fastqc -t 2 -o fastqc/raw {input}'

rule trimmomatic_PE:
	input:
		r1 = config['path']+'/{sample}_R1.fastq.gz',
		r2 = config['path']+'/{sample}_R2.fastq.gz'
	output:
		r1_paired = 'clean/{sample}_R1_paired.fastq.gz',
		r2_paired = 'clean/{sample}_R2_paired.fastq.gz',
		r1_unpaired = 'clean/{sample}_R1_unpaired.fastq.gz',
		r2_unpaired = 'clean/{sample}_R2_unpaired.fastq.gz'
	params:
		adapter = config['adapter']
	shell:
		'trimmomatic PE -threads 3 -phred33 {input.r1} {input.r2} {output.r1_paired} {output.r1_unpaired} {output.r2_paired} {output.r2_unpaired} ILLUMINACLIP:{params.adapter}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36'

rule fastqc_clean_PE:
	input:
		'clean/{sample}_R1_paired.fastq.gz',
		'clean/{sample}_R2_paired.fastq.gz'
	output:
		'fastqc/clean/{sample}_R1_paired_fastqc.html',
		'fastqc/clean/{sample}_R2_paired_fastqc.html'
	shell:
		'fastqc -t 2 -o fastqc/clean {input}'

rule fastqc_stat_PE:
	input:
		['fastqc/raw/{sample}_R1_fastqc.html'.format(sample=x) for x in config['samples']],
		['fastqc/raw/{sample}_R2_fastqc.html'.format(sample=x) for x in config['samples']],
		['fastqc/clean/{sample}_R1_paired_fastqc.html'.format(sample=x) for x in config['samples']],
		['fastqc/clean/{sample}_R2_paired_fastqc.html'.format(sample=x) for x in config['samples']]
	output:
		'stat/fastqc_stat.tsv'
	params:
		Rscript = config['Rscript_path']
	shell:
		'{params.Rscript} script/reads_stat_by_fastqcr.R'

rule map_to_genome_PE:
	input:
		r1='clean/{sample}_R1_paired.fastq.gz',
		r2='clean/{sample}_R2_paired.fastq.gz'
	output:
		bam='bam/{sample}.sort.bam'
	params:
		prefix = 'bam/{sample}',
		index=config['index'],
		cpu=config['cpu'],
	shell:
		'bowtie2 -p {params.cpu} -x {params.index} -1 {input.r1} -2 {input.r2} |samtools view -@ {params.cpu} -Shub|samtools  sort -@ {params.cpu} --output-fmt BAM - -T {params.prefix} -o {output.bam}'

rule bam_idx:
	input:
		bam = 'bam/{sample}.sort.bam'
	output:
		bai = 'bam/{sample}.sort.bam.bai'
	shell:
		'samtools index {input.bam} {output.bai}'

rule bam_rmdup:
	input:
		bam = 'bam/{sample}.sort.bam'
	output:
		rmdup = 'bam/{sample}.rmdup.bam'
	shell:
		'samtools rmdup {input.bam} {output.rmdup}'

rule bam_mpileup:
	input:
		bam='bam/{sample}.rmdup.bam'
	output:
		vcf='vcf/{sample}.raw.vcf'
	params:
		genome=config['genome']
	shell:
		'samtools mpileup -t DP,SP -Eg -f {params.genome} {input.bam}|bcftools call -mv > {output.vcf}'

rule vcf_filter:
	input:
		'vcf/{sample}.raw.vcf'
	output:
		'vcf/{sample}.flt.vcf'
	shell:
		'bcftools view {input} | vcfutils.pl varFilter -d 20 -Q 20 - >{output}'

rule vcf_compress:
	input:
		'vcf/{sample}.flt.vcf'
	output:
		'vcf/{sample}.flt.gz'
	shell:
		'bcftools view -Oz {input} > {output}'

rule vcf_index:
	input:
		'vcf/{sample}.flt.gz'
	output:
		'vcf/{sample}.flt.gz.csi'
	shell:
		'bcftools index {input} -o {output}'

rule vcf_isec:
	input:
		i1='vcf/{group1}-90-2_combined.flt.gz',
		i2='vcf/{group2}_combined.flt.gz'
	output:
		'result/{group1}_{group2}'
	shell:
		'bcftools isec -p {output} {input.i1} {input.i2}'

rule bam_qc:
	input:
		bam = 'bam/{sample}.sort.bam'
	output:
		bamqc_dir = 'bam/{sample}.bamqc',
		bamqc_html = 'bam/{sample}.bamqc/qualimapReport.html'
	params:
		cpu = config['cpu']
	shell:
		"qualimap bamqc --java-mem-size=10G -nt {params.cpu} -bam {input.bam} -outdir {output.bamqc_dir}"

rule bam_qc_stat:
	input:
		['bam/{sample}.bamqc/qualimapReport.html'.format(sample=x) for x in config['samples']]
	output:
		'stat/bamqc_stat.tsv'
	params:
		Rscript = config['Rscript_path']		
	shell:
		"{params.Rscript} script/mapping_stat_by_bamqc.R"


