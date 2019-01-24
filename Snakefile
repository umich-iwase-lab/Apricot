configfile: 'config.yaml'

INPUT_DIR = config['dirs']['input']
OUTPUT_DIR = config['dirs']['output']

rule all:
    input: OUTPUT_DIR + '/04-multiqc/multiqc.html'

rule fastqc_seq:
    input: INPUT_DIR + '/{sample}.{read}.fastq.gz'
    output: OUTPUT_DIR + '/01-fastqc_seq/{sample}.{read}_fastqc.html'
    params: 
        fastqc_dir= OUTPUT_DIR+ '/01-fastqc_seq/'
    shell: 'fastqc {input} -o {params.fastqc_dir}'

rule star_align:
    input: INPUT_DIR + '/{sample}.R1.fastq.gz', INPUT_DIR + '/{sample}.R2.fastq.gz',
    output: OUTPUT_DIR + '/02-star_align/{sample}.star_align.bam'
    shell: 'cat {input} > {output}'

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
