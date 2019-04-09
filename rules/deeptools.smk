ALL.append(expand(OUTPUT_DIR + '/10-deeptools/{sample}.bw',sample=config['samples']))

rule deeptools:
    input: OUTPUT_DIR + '/10-sort_index_after_dedup/{sample}.genome.umi_tools_dedup.sorted.bam'
    output: OUTPUT_DIR + '/10-deeptools/{sample}.bw'
    benchmark: OUTPUT_DIR + '/benchmarks/deeptools.{sample}.benchmark.txt'
    log:
        OUTPUT_DIR + '/10-deeptools/.log/{sample}.deeptools.log'
    shell: '(bamCoverage -b {input} -o {output}) 2>&1 | tee {log}'
