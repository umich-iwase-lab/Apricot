# ALL.append(expand(OUTPUT_DIR + '/01-umi_tools_extract/processed.{sample}.{read}.fastq.gz',
#               sample=config['samples'],
#               read=READS))

def _build_inputs():
    return expand(OUTPUT_DIR + '/normalize_read_names/{{sample}}.{read}.fastq.gz',
                  read=READS)

def _build_outputs():
    return expand(OUTPUT_DIR + '/umi_tools_extract/processed.{{sample}}.{read}.fastq.gz',
                  read=READS)

def _build_read_flags(sample):
    flag = {'r1i': '--stdin',
            'r1o': '--stdout',
            'r2i': '--read2-in',
            'r2o': '--read2-out',
            }
    reads = [('r{}i'.format(i + 1), name.format(sample=sample)) for i, name in enumerate(_build_inputs())] \
            + [('r{}o'.format(i + 1), name.format(sample=sample)) for i, name in enumerate(_build_outputs())]

    return ' '.join(['{}={}'.format(flag[read], file) for read, file in reads])

rule umi_tools_extract:
    input:
        _build_inputs()
    output:
        _build_outputs()
    log: OUTPUT_DIR + '/umi_tools_extract/.log/{sample}.umi_tools_extract.log'
    benchmark: OUTPUT_DIR + '/benchmarks/umi_tools_extract.{sample}.benchmark.txt'
    params:
        read_flags = lambda wildcards: _build_read_flags(wildcards.sample)
    shell:'''(
umi_tools extract \
    --extract-method=string \
    --bc-pattern=NNNNNNNN \
    {params.read_flags}
) 2>&1 | tee {log}'''
