configfile: 'config.yaml'

INPUT_DIR = config['dirs']['input']
OUTPUT_DIR = config['dirs']['output']
REFERENCE_DIR = config['dirs']['reference'] 

rule all:
    input: OUTPUT_DIR + '/04-multiqc/multiqc.html'

rule umi_tools_extract:
    input: 
        read1 = INPUT_DIR + '/{sample}.R1.fastq.gz',
        read2 = INPUT_DIR + '/{sample}.R2.fastq.gz'
    output: 
        read1 = OUTPUT_DIR + '/02-umitools_markDuplicates/processed.{sample}.R1.fastq.gz',
        read2 = OUTPUT_DIR + '/02-umitools_markDuplicates/processed.{sample}.R2.fastq.gz'
    log: OUTPUT_DIR + '/02-umitools_markDuplicates/.log/{sample}.umitools_markDuplicates.log'
    shell:'''(
umi_tools extract \
--stdin={input.read1} \
--extract-method=string \
--bc-pattern=NNNNNNNN \
--read2-in={input.read2} \
--stdout={output.read1} \
--read2-out={output.read2}
) 2>&1 | tee {log}'''

rule fastqc_seq:
    input: INPUT_DIR + '/{sample}.{read}.fastq.gz'
    output: OUTPUT_DIR + '/01-fastqc_seq/{sample}.{read}_fastqc.html'
    params:
        fastqc_dir= OUTPUT_DIR+ '/01-fastqc_seq/'
    shell: 'fastqc {input} -o {params.fastqc_dir}'

rule star_align:
    input: OUTPUT_DIR + '/02-umitools_markDuplicates/processed.{sample}.R1.fastq.gz',
           OUTPUT_DIR + '/02-umitools_markDuplicates/processed.{sample}.R2.fastq.gz',
    output: OUTPUT_DIR + '/02-star_align/{sample}.star_align.bam'
    log: OUTPUT_DIR + '/02-star_align/.log/{sample}.star_align.log'
    threads: 12
    params:
        genomeDir = REFERENCE_DIR + '/' + config['genome_reference']['star_genome_dir'],
        sjdbGTFfile = REFERENCE_DIR + '/' + config['genome_reference']['gtf'],
        outFileNamePrefix = OUTPUT_DIR + '/02-star_align/{sample}.',
        star_bam_file = '{sample}.Aligned.sortedByCoord.out.bam',
        output_bam_file = '{sample}.star_align.bam'
    shell: '''(
STAR \
--genomeDir {params.genomeDir} \
--runThreadN {threads} \
--readFilesIn {input} \
--readFilesCommand 'gunzip -c' \
--outFileNamePrefix {params.outFileNamePrefix} \
--sjdbGTFfile {params.sjdbGTFfile} \
--runMode alignReads \
--outSAMtype BAM SortedByCoordinate \
--limitBAMsortRAM 64000000000 \
--outFilterType BySJout \
--outFilterMultimapNmax 20 \
--alignSJoverhangMin 8 \
--alignSJDBoverhangMin 1 \
--alignIntronMin 20 \
--outFilterMismatchNmax 999 \
--outFilterScoreMinOverLread 0 \
--outFilterMatchNminOverLread 0 \
--alignMatesGapMax 1000000 \
--alignIntronMax 1000000 \
--chimSegmentMin 20

pushd {OUTPUT_DIR}/02-star_align/
ln --symbolic {params.star_bam_file} {params.output_bam_file}
popd
) 2>&1 | tee {log}
    '''

rule umi_tools_dedup:
    input: OUTPUT_DIR + '/02-star_align/{sample}.star_align.bam'
    output: OUTPUT_DIR + '/03-umi_tools_dedup/{sample}_dedup.bam'
    shell: 'umi_tools dedup -I {input} --paired -S {output}'

rule fastqc_align:
    input: OUTPUT_DIR + '/02-star_align/{sample}.star_align.bam'
    output: OUTPUT_DIR + '/03-fastqc_align/{sample}.fastqc_align.html'
    shell: 'cat {input} >  {output}'

rule multiqc:
    input: expand(OUTPUT_DIR + '/01-fastqc_seq/{sample}.{read}_fastqc.html', \
                  sample=config['samples'], \
                  read=['R1', 'R2']),
           expand(OUTPUT_DIR + '/03-fastqc_align/{sample}.fastqc_align.html', sample=config['samples']),
    output: OUTPUT_DIR + '/04-multiqc/multiqc.html'
    shell: 'cat {input} > {output}'
