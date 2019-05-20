rule fastqc_seq:
    input:
        OUTPUT_DIR + '/umi_tools_extract/processed.{sample}.{read}.fastq.gz',
    output:
        OUTPUT_DIR + '/fastqc_seq/processed.{sample}.{read}_fastqc.html',
    log:
        OUTPUT_DIR + '/fastqc_seq/.log/processed.{sample}.{read}_fastqc.log',
    benchmark:
        OUTPUT_DIR + '/benchmarks/fastqc_seq.{sample}.{read}.benchmark.txt'
    params:
        fastqc_dir= OUTPUT_DIR+ '/fastqc_seq/'
    shell: '(fastqc {input} -o {params.fastqc_dir}) 2>&1 | tee {log}'
