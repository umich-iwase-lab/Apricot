configfile: 'config.yaml'

INPUT_DIR = config['dirs']['input']
OUTPUT_DIR = config['dirs']['output']
REFERENCE_DIR = config['dirs']['reference'] 

rule all:
    input: OUTPUT_DIR + '/03-umi_tools_dedup/SampleB_dedup.bam',
           expand(OUTPUT_DIR + '/02-fastqc_seq/processed.{sample}.{read}_fastqc.html', sample=config['samples'], read=['R1', 'R2']),
           expand(OUTPUT_DIR + '/04-fastqc_align/{sample}_dedup_fastqc.html', sample=config['samples']),
           OUTPUT_DIR + '/05-multiqc/multiqc_report.html',

rule umi_tools_extract:
    input: 
        read1 = INPUT_DIR + '/{sample}.R1.fastq.gz',
        read2 = INPUT_DIR + '/{sample}.R2.fastq.gz'
    output: 
        read1 = OUTPUT_DIR + '/01-umitools_markDuplicates/processed.{sample}.R1.fastq.gz',
        read2 = OUTPUT_DIR + '/01-umitools_markDuplicates/processed.{sample}.R2.fastq.gz'
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
    input: OUTPUT_DIR + '/01-umitools_markDuplicates/processed.{sample}.{read}.fastq.gz',
    output: OUTPUT_DIR + '/02-fastqc_seq/processed.{sample}.{read}_fastqc.html', 
    params:
        fastqc_dir= OUTPUT_DIR+ '/02-fastqc_seq/'
    shell: 'fastqc {input} -o {params.fastqc_dir}'

rule star_align:
    input: OUTPUT_DIR + '/01-umitools_markDuplicates/processed.{sample}.R1.fastq.gz',
           OUTPUT_DIR + '/01-umitools_markDuplicates/processed.{sample}.R2.fastq.gz',
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
rule indexing:
    input: OUTPUT_DIR + '/02-star_align/{sample}.star_align.bam'
    output: OUTPUT_DIR + '/02-star_align/{sample}.star_align.bam.bai',
    shell: 'samtools index {input}'
rule umi_tools_dedup:
    input: 
        bam = OUTPUT_DIR + '/02-star_align/{sample}.star_align.bam',
        bai = OUTPUT_DIR + '/02-star_align/{sample}.star_align.bam.bai'
    log: OUTPUT_DIR + '/03-umi_tools_dedup/.log/{sample}_dedup.log'
    output: OUTPUT_DIR + '/03-umi_tools_dedup/{sample}_dedup.bam'
    shell: '(umi_tools dedup -I {input.bam} --paired -S {output}) 2>&1 | tee {log}'

rule fastqc_align:
    input: OUTPUT_DIR + '/03-umi_tools_dedup/{sample}_dedup.bam'
    output: OUTPUT_DIR + '/04-fastqc_align/{sample}_dedup_fastqc.html'
    params:
        fastqc_dir= OUTPUT_DIR+ '/04-fastqc_align/'
    shell: 'fastqc {input} -o {params.fastqc_dir}'

rule multiqc:
    input: OUTPUT_DIR + '/04-fastqc_align/{sample}_dedup_fastqc.html',
           OUTPUT_DIR + '/02-fastqc_seq/processed.{sample}.{read}_fastqc.html',
    output: OUTPUT_DIR + '/05-multiqc/multiqc_report.html'
    params:
        multiqc_dir= OUTPUT_DIR + '/05-multiqc/'
    shell: 'multiqc {OUTPUT_DIR} -o {params.multiqc_dir}'

rule align_stringtie:
    input:
#        reference_checksum = CONFIG_CHECKSUMS_DIR + 'config-references.watermelon.md5',
        bams = OUTPUT_DIR + '/02-star_align/{sample}.star_align.bam',
        gtf = REFERENCE_DIR + '/' + config['genome_reference']['gtf'],
    output:
        gtf = OUTPUT_DIR + '06-stringtie/{sample}.gtf',
        ballgown = expand(OUTPUT_DIR + '06-stringtie/ballgown/{sample}/{bgown_prefixes}.ctab',#
                          sample='{sample}',
                          bgown_prefixes=['e2t','e_data','i2t','i_data','t_data']),
        prep_input = OUTPUT_DIR + '06-stringtie/{sample}_prepDE_input.txt'
    log:
        OUTPUT_DIR + '06-stringtie/.log/{sample}.align_stringtie.log'
    benchmark:
        'benchmarks/{sample}.align_stringtie.benchmark.txt'
    conda:
        'envs/align_hisat2_stringtie.yaml'
    params:
        sample = '{sample}',
#        strand_flag = rnaseq_snakefile_helper.strand_option_stringtie(config),
        ballgown_path = OUTPUT_DIR + '06-stringtie/ballgown/{sample}/'
    threads: 8
    shell:
        '''(
        echo {params.sample} {output.gtf} > {output.prep_input}

        stringtie {input.bams} \
            -p {threads} \
            -G {input.gtf} \
            -b {params.ballgown_path} \
            -e \
            {params.strand_flag} \
            -o {output.gtf}
        ) 2>&1 | tee {log}'''
