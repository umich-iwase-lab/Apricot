rule convert_dedup_for_rsem:
    input:
        OUTPUT_DIR + '/06-umi_tools_dedup/{sample}.umi_tools_dedup.bam'
    output:
        OUTPUT_DIR + '/07-convert_dedup_for_rsem/{sample}.umi_tools_dedup.converted.bam'
    benchmark:
        OUTPUT_DIR + '/benchmarks/{sample}.convert_dedup_for_rsem.benchmark.txt'
    log:
        OUTPUT_DIR + '/07-convert_dedup_for_rsem/.log/{sample}.convert_dedup_for_rsem.log'
    params:
        output_name = OUTPUT_DIR + '/07-convert_dedup_for_rsem/{sample}.umi_tools_dedup.converted'
    threads:
        12
    shell:'''(
    convert-sam-for-rsem \
    -p {threads} \
    {input} \
    {params.output_name}
    )2>&1 | tee {log}
    '''
