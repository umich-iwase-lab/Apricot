rule stringtie:
    input:
        bam = OUTPUT_DIR + '/05-umi_tools_dedup/{sample}.umi_tools_dedup.bam',
        gtf = REFERENCE_DIR + '/' + config['genome_reference']['gtf'],
    output:
        gtf = OUTPUT_DIR + '/08-stringtie/{sample}.gtf',
        prep_input = OUTPUT_DIR + '/08-stringtie/{sample}_prepDE_input.txt'
    log:
        OUTPUT_DIR + '/08-stringtie/.log/{sample}.stringtie.log'
    benchmark:
        OUTPUT_DIR + '/benchmarks/{sample}.stringtie.benchmark.txt'
    params:
        sample = '{sample}',
        strand_flag = config['stringtie']['strand_flag'],
    threads: 8
    shell:'''(
echo {params.sample} {output.gtf} > {output.prep_input}
stringtie {input.bam} \
    -p {threads} \
    -G {input.gtf} \
    -e \
    {params.strand_flag} \
    -o {output.gtf}
) 2>&1 | tee {log}'''
