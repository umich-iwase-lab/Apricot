rule umi_tools_extract:
    input:
        read1 = OUTPUT_DIR + '/00-normalize_read_names/{sample}.R1.fastq.gz',
        read2 = OUTPUT_DIR + '/00-normalize_read_names/{sample}.R2.fastq.gz',
    output:
        read1 = OUTPUT_DIR + '/01-umi_tools_extract/processed.{sample}.R1.fastq.gz',
        read2 = OUTPUT_DIR + '/01-umi_tools_extract/processed.{sample}.R2.fastq.gz',
    log: OUTPUT_DIR + '/01-umi_tools_extract/.log/{sample}.umi_tools_extract.log'
    benchmark: OUTPUT_DIR + '/benchmarks/umi_tools_extract.{sample}.benchmark.txt'
    shell:'''(
umi_tools extract \
--stdin={input.read1} \
--extract-method=string \
--bc-pattern=NNNNNNNN \
--read2-in={input.read2} \
--stdout={output.read1} \
--read2-out={output.read2}
) 2>&1 | tee {log}'''
