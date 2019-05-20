rule convert_dedup_for_rsem:
    input:
        OUTPUT_DIR + '/umi_tools_dedup/{sample}.{type}.umi_tools_dedup.bam'
    output:
        OUTPUT_DIR + '/convert_dedup_for_rsem/{sample}.{type}.convert_dedup_for_rsem.bam'
    benchmark:
        OUTPUT_DIR + '/benchmarks/{sample}.{type}.convert_dedup_for_rsem.benchmark.txt'
    log:
        OUTPUT_DIR + '/convert_dedup_for_rsem/.log/{sample}.{type}.convert_dedup_for_rsem.log'
    params:
        output_name = OUTPUT_DIR + '/convert_dedup_for_rsem/{sample}.{type}.convert_dedup_for_rsem'
    threads:
        12
    shell:'''(
    convert-sam-for-rsem \
    -p {threads} \
    {input} \
    {params.output_name}
    )2>&1 | tee {log}
    '''
