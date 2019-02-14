rule umi_tools_dedup:
    input:
        bam = OUTPUT_DIR + '/03-star_align/{sample}.star_align.bam',
        bai = OUTPUT_DIR + '/03-star_align/{sample}.star_align.bam.bai'
    output:
        OUTPUT_DIR + '/05-umi_tools_dedup/{sample}.umi_tools_dedup.bam'
    log:
        OUTPUT_DIR + '/05-umi_tools_dedup/.log/{sample}.umi_tools_dedup.log'
    benchmark:
        OUTPUT_DIR + '/benchmarks/umi_tools_dedup.{sample}.benchmark.txt'
    shell: '''(
umi_tools dedup \
    -I {input.bam} \
    --paired \
    -S {output}\
) 2>&1 | tee {log}'''
