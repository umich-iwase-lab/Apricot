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
include: 'rules/multiqc.smk'
include: 'rules/deeptools.smk'

rule all:
    input:
        *ALL
