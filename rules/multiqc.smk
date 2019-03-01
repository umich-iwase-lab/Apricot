ALL.append(OUTPUT_DIR + '/90-multiqc/multiqc_report.html')

rule multiqc:
    input: expand(OUTPUT_DIR + '/07-fastqc_align/{sample}_fastqc.html',sample=config['samples']),
           expand(OUTPUT_DIR + '/02-fastqc_seq/processed.{sample}.{read}_fastqc.html', sample=config['samples'], read=['R1', 'R2']),
           OUTPUT_DIR + '/10-ballgown/gene_fpkms.txt',
           OUTPUT_DIR + '/10-ballgown/iso_fpkms.txt',
           OUTPUT_DIR + '/09-stringtie_prepDE/gene_count_matrix.stringtie_prepDE.csv',
           OUTPUT_DIR + '/09-stringtie_prepDE/transcript_count_matrix.stringtie_prepDE.csv',
    output: OUTPUT_DIR + '/90-multiqc/multiqc_report.html'
    benchmark: OUTPUT_DIR + '/benchmarks/multiqc.benchmark.txt'
    params:
        multiqc_dir= OUTPUT_DIR + '/90-multiqc/'
    shell: 'multiqc {OUTPUT_DIR} -o {params.multiqc_dir}'