if not config:
    raise ValueError('Please specify --configfile when launching snakemake')

INPUT_DIR = config['dirs']['input']
OUTPUT_DIR = config['dirs']['output']
REFERENCE_DIR = config['dirs']['reference']

ALL = []

include: 'rules/normalize_read_names.smk'
include: 'rules/umi_tools_extract.smk'
include: 'rules/fastqc_seq.smk'
include: 'rules/star_genome_generate.smk'
include: 'rules/star_align.smk'
include: 'rules/star_align_index.smk'
include: 'rules/star_align_flagstat.smk'
include: 'rules/umi_tools_dedup.smk'
include: 'rules/umi_tools_dedup_flagstat.smk'
include: 'rules/fastqc_align.smk'
include: 'rules/stringtie.smk'
include: 'rules/stringtie_prepDE.smk'
include: 'rules/ballgown_build_fpkm_matrices.smk'

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
ALL.append(OUTPUT_DIR + '/90-multiqc/multiqc_report.html')

rule all:
    input:
        *ALL
