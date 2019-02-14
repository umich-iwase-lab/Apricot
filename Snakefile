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

rule all:
    input:
        *ALL












# rule multiqc:
#     input: expand(OUTPUT_DIR + '/05-fastqc_align/{sample}_dedup_fastqc.html',sample=config['samples']),
#            expand(OUTPUT_DIR + '/02-fastqc_seq/processed.{sample}.{read}_fastqc.html', sample=config['samples'], read=['R1', 'R2']),
#     output: OUTPUT_DIR + '/05-multiqc/multiqc_report.html'
#     benchmark: OUTPUT_DIR + '/benchmarks/multiqc.benchmark.txt'
#     params:
#         multiqc_dir= OUTPUT_DIR + '/05-multiqc/'
#     shell: 'multiqc {OUTPUT_DIR} -o {params.multiqc_dir}'
