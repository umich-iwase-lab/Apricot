# The downstream dedup script has a weak notion of read pairs and under
# certain condistions it will pick the "best" read pair by combining a primary
# alignment with a secondary alignment. This is unfortunate and causes problems
# with RSEM. To avoid this problem, we exclude all but the primary alignment.
rule samtools_exclude_non_primary:
    input:
        bam = OUTPUT_DIR + '/03-rsem_star_align/{sample}.{type}.bam',
    output:
        bam = OUTPUT_DIR + '/05-samtools_exclude_non_primary/{sample}.{type}.samtools_exclude_non_primary.bam',
        bai = OUTPUT_DIR + '/05-samtools_exclude_non_primary/{sample}.{type}.samtools_exclude_non_primary.bam.bai',
        flagstat = OUTPUT_DIR + '/05-samtools_exclude_non_primary/{sample}.{type}.samtools_exclude_non_primary.flagstat',
    log:
        OUTPUT_DIR + '/05-samtools_exclude_non_primary/.log/{sample}.{type}.samtools_exclude_non_primary.log'
    benchmark:
        OUTPUT_DIR + '/benchmarks/samtools_exclude_non_primary.{sample}.{type}.txt'
    params:
        filtered = OUTPUT_DIR + '/05-samtools_exclude_non_primary/{sample}.{type}.filtered.bam',
    threads:
        12
    shell:'''(
   samtools view \
       -b \
       -F0x900 \
       {input.bam} \
       > {params.filtered}
    samtools sort \
        -@ {threads} \
        -m 1G \
        -o {output.bam} \
        {params.filtered}
    rm {params.filtered}
    samtools flagstat {output.bam} > {output.flagstat}
    samtools index {output.bam}
    ) 2>&1 | tee {log}'''
