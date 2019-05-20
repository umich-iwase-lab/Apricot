ALL.extend(expand(OUTPUT_DIR + '/deeptools_bigwig/{sample}.bw',
                  sample=config['samples']))

rule deeptools_bigwig:
    input: OUTPUT_DIR + '/samtools_sort_index_after_dedup/{sample}.genome.umi_tools_dedup.sorted.bam'
    output: OUTPUT_DIR + '/deeptools_bigwig/{sample}.bw'
    benchmark: OUTPUT_DIR + '/benchmarks/deeptools_bigwig.{sample}.benchmark.txt'
    log:
        OUTPUT_DIR + '/deeptools_bigwig/.log/{sample}.deeptools_bigwig.log'
    shell: '(bamCoverage -b {input} -o {output}) 2>&1 | tee {log}'
