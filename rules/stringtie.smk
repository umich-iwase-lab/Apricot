rule stringtie:
    input:
        bam = OUTPUT_DIR + '/05-umi_tools_dedup/{sample}.umi_tools_dedup.bam',
        gtf = REFERENCE_DIR + '/' + config['genome_reference']['gtf'],
    output:
        gtf = OUTPUT_DIR + '/08-stringtie/{sample}.gtf',
        prep_input = OUTPUT_DIR + '/08-stringtie/{sample}_prepDE_input.txt',
        ballgown_path = directory(OUTPUT_DIR + '/08-stringtie/ballgown/{sample}'),
    log:
        OUTPUT_DIR + '/08-stringtie/.log/{sample}.stringtie.log'
    benchmark:
        OUTPUT_DIR + '/benchmarks/stringtie.{sample}.benchmark.txt'
    params:
        sample = '{sample}',
        strand_flag = config['sequencing_parameters']['strand_flag'],
    threads: 8
    shell:'''(
echo {params.sample} {output.gtf} > {output.prep_input}
stringtie {input.bam} \
    -p {threads} \
    -G {input.gtf} \
    -e \
    -b {output.ballgown_path} \
    {params.strand_flag} \
    -o {output.gtf}
) 2>&1 | tee {log}'''
