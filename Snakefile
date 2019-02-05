configfile: 'config.yaml'

INPUT_DIR = config['dirs']['input']
OUTPUT_DIR = config['dirs']['output']
REFERENCE_DIR = config['dirs']['reference'] 

rule all:
    input: OUTPUT_DIR + '/03-umi_tools_dedup/SampleB_dedup.bam',
           expand(OUTPUT_DIR + '/02-fastqc_seq/processed.{sample}.{read}_fastqc.html', sample=config['samples'], read=['R1', 'R2']),
           expand(OUTPUT_DIR + '/04-fastqc_align/{sample}_dedup_fastqc.html', sample=config['samples']),
           OUTPUT_DIR + '/05-multiqc/multiqc_report.html',
           expand(OUTPUT_DIR + '/06-stringtie/{sample}.gtf',sample=config['samples']),
           gene_counts = OUTPUT_DIR + '/06-stringtie/gene_count_matrix.csv',
           transcript_counts = OUTPUT_DIR + '/06-stringtie/transcript_count_matrix.csv',

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
    input: expand(OUTPUT_DIR + '/04-fastqc_align/{sample}_dedup_fastqc.html',sample=config['samples']),
           expand(OUTPUT_DIR + '/02-fastqc_seq/processed.{sample}.{read}_fastqc.html', sample=config['samples'], read=['R1', 'R2']),
    output: OUTPUT_DIR + '/05-multiqc/multiqc_report.html'
    params:
        multiqc_dir= OUTPUT_DIR + '/05-multiqc/'
    shell: 'multiqc {OUTPUT_DIR} -o {params.multiqc_dir}'

rule stringtie:
    input:
        bams = OUTPUT_DIR + '/02-star_align/{sample}.star_align.bam',
        gtf = REFERENCE_DIR + '/' + config['genome_reference']['gtf'],
    output:
        gtf = OUTPUT_DIR + '/06-stringtie/{sample}.gtf',
        prep_input = OUTPUT_DIR + '/06-stringtie/{sample}_prepDE_input.txt'
    log:
        OUTPUT_DIR + '/06-stringtie/.log/{sample}.stringtie.log'
    benchmark:
        OUTPUT_DIR + '/benchmarks/{sample}.stringtie.benchmark.txt'
    params:
        sample = '{sample}',
        strand_flag = config['stringtie_strand_flag'],
    threads: 8
    shell:
        '''(
        echo {params.sample} {output.gtf} > {output.prep_input}

        stringtie {input.bams} \
            -p {threads} \
            -G {input.gtf} \
            -e \
            {params.strand_flag} \
            -o {output.gtf}
        ) 2>&1 | tee {log}'''

rule stringtie_prepDE:
    input:
        gtfs = expand(OUTPUT_DIR + '/06-stringtie/{sample}.gtf',
                      sample=config['samples']),
        prep_input = expand(OUTPUT_DIR + '/06-stringtie/{sample}_prepDE_input.txt',
                            sample=config['samples']),
    output:
        prep_output = OUTPUT_DIR + '/06-stringtie/prepDE_input.txt',
        gene_counts = OUTPUT_DIR + '/06-stringtie/gene_count_matrix.csv',
        transcript_counts = OUTPUT_DIR + '/06-stringtie/transcript_count_matrix.csv',
    log:
        OUTPUT_DIR + '/06-stringtie/.log/stringtie_prepDE.log'
    benchmark:
        OUTPUT_DIR + '/benchmarks/stringtie_prepDE.benchmark.txt'
    params:
        length = config['stringtie_read_length'],
        prep_input = ' '.join(expand(OUTPUT_DIR + '/06-stringtie/{sample}_prepDE_input.txt',
                                     sample=config['samples']))
    shell:
        '''(
        cat {params.prep_input} > {output.prep_output}

        prepDE.py \
            -i {output.prep_output} \
            -g {output.gene_counts} \
            -t {output.transcript_counts} \
            -l {params.length}
        ) 2>&1 | tee {log}'''

