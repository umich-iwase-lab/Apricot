from os.path import join
_defaultSjdbOverhang = config['sequencing_parameters']['read_length'] - 1
_star_config = config['genome_reference']['star']

#ALL.append(expand(OUTPUT_DIR + '/04-star_align_flagstat/{sample}.star_align_flagstat.flagstat',
#                  sample=config['samples']))
rule rsem_star_genome_generate:
    input:
        fasta = REFERENCE_DIR + '/' + config['genome_reference']['fasta'],
        gtf = REFERENCE_DIR + '/' + config['genome_reference']['gtf'],
    output:
        parametersFiles = join(\
                REFERENCE_DIR,
                _star_config['genome_dir'],
                'genomeParameters.txt'
                ),
    log:
        REFERENCE_DIR + '/' + _star_config['genome_dir'] + '/.log/rsem_star_genome_generate.log'
    benchmark:
        OUTPUT_DIR + '/benchmarks/rsem_star_genome_generate.benchmark.txt'
    threads: 12
    params:
        genomeDir = REFERENCE_DIR + '/' + _star_config['genome_dir'],
        sjdbOverhang = _star_config.get('sjdbOverhang', _defaultSjdbOverhang),
    shell:'''
(rsem-prepare-reference --gtf {input.gtf} \
                    --star \
                    --star-path ~/miniconda3/envs/apricot/bin \
                    --star-sjdboverhang {params.sjdbOverhang} \
                    -p {threads} \
                    {input.fasta} \
                    {params.genomeDir}/RSEM_ref
) 2>&1 | tee {log}
    '''
#     shell: '''(
# mkdir -p {params.genomeDir}.tmp
# STAR \
#   --runMode genomeGenerate \
#   --runThreadN {threads} \
#   --genomeFastaFiles {input.fasta} \
#   --genomeDir {params.genomeDir}.tmp \
#   --sjdbGTFfile {input.gtf} \
#   --sjdbOverhang {params.sjdbOverhang}
# mv {params.genomeDir}.tmp/* {params.genomeDir}
# rmdir {params.genomeDir}.tmp
# ) 2>&1 | tee {log}'''
