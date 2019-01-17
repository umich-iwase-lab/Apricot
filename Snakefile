# samples = ['SampleA','SampleB']
configfile: 'config.yaml'
rule all:
    input: 'outputs/04-multiqc/multiqc.html'

rule fastqc_seq:
    input: 'inputs/{Sample}.R1.fastq.gz','inputs/{Sample}.R2.fastq.gz'
    output:'outputs/01-fastqc_seq/{Sample}.fastqc_seq.html'
    shell: 'cat {input} > {output}'
rule align:
    input: 'inputs/{sample_name}.R1.fastq.gz','inputs/{sample_name}.R2.fastq.gz',
    output:'outputs/02-align/{sample_name}.align.bam'
    shell: 'cat {input} > {output}'
rule fastqc_align:
    input:'outputs/02-align/{sample}.align.bam'
    output: 'outputs/03-fastqc_align/{sample}.fastqc_align.html'
    shell: 'cat {input} >  {output}'
rule multiqc:
    input: expand('outputs/01-fastqc_seq/{sample}.fastqc_seq.html', sample=config['samples']),
           expand('outputs/03-fastqc_align/{sample}.fastqc_align.html', sample=config['samples']),
    output: 'outputs/04-multiqc/multiqc.html'
    shell: 'cat {input} > {output}'
