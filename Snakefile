from snakemake.logging import logger
if not config:
    raise ValueError('Please specify --configfile when launching snakemake')

INPUT_DIR = config['dirs']['input']
OUTPUT_DIR = config['dirs']['output']
REFERENCE_DIR = config['dirs']['reference']

_read_counts = list(set([len(reads) for sample, reads in config['samples'].items()]))
if len(_read_counts) != 1:
    raise ValueError('All samples must be all single-end or paired-end')
PAIRED_END = _read_counts[0] == 2
print('paired = {}'.format(PAIRED_END))
READS = ['R1', 'R2'] if PAIRED_END else ['R1']

ALL = []

include: 'rules/normalize_read_names.smk'
include: 'rules/umi_tools_extract.smk'
include: 'rules/fastqc_seq.smk'
include: 'rules/rsem_star_genome_generate.smk'
include: 'rules/rsem_star_align.smk'
include: 'rules/rsem_star_align_index.smk'
include: 'rules/rsem_star_align_flagstat.smk'
include: 'rules/samtools_exclude_non_primary.smk'
include: 'rules/umi_tools_dedup.smk'
include: 'rules/umi_tools_dedup_flagstat.smk'
include: 'rules/fastqc_align.smk'
include: 'rules/rsem_calculate_expression_using_dedup.smk'
include: 'rules/combine_counts_after_dedup.smk'
include: 'rules/combine_counts.smk'
include: 'rules/convert_dedup_for_rsem.smk'
include: 'rules/multiqc.smk'
include: 'rules/samtools_sort_index_after_dedup.smk'
include: 'rules/deeptools_bigwig.smk'

rule all:
    input:
        *ALL

def email(subject_prefix):
    msg = 'config file:\n{}\nlog file:\n{}'.format(workflow.overwrite_configfile, logger.get_logfile())
    email_config = config.get('email')
    if not email_config:
        print(subject_prefix, msg)
    else:
        command = "echo '{msg}' | mutt -s '{subject_prefix}{subject}' {to}".format(
                to=email_config['to'],
                subject_prefix=subject_prefix,
                subject=email_config['subject'],
                msg=msg,
                )
        shell(command)


onstart:
    email('Apricot started: ')

onsuccess:
    email('Apricot completed ok: ')

onerror:
    email('Apricot ERROR: ')
