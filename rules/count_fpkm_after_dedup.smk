ALL.append([OUTPUT_DIR + '/08-rsem_calculate_expression/gene_count.txt',
            OUTPUT_DIR + '/08-rsem_calculate_expression/gene_fpkm.txt',
            OUTPUT_DIR + '/08-rsem_calculate_expression/isoform_count.txt',
            OUTPUT_DIR + '/08-rsem_calculate_expression/isoform_fpkm.txt',
            ])
rule count_fpkm_after_dedup:
    input:
        genes = expand(OUTPUT_DIR + '/08-rsem_calculate_expression/{sample}.genes.results',
                      sample=config['samples']),
        isoforms = expand(OUTPUT_DIR + '/08-rsem_calculate_expression/{sample}.isoforms.results',
                      sample=config['samples']),
    output:
        OUTPUT_DIR + '/08-rsem_calculate_expression/gene_count.txt',
        OUTPUT_DIR + '/08-rsem_calculate_expression/gene_fpkm.txt',
        OUTPUT_DIR + '/08-rsem_calculate_expression/isoform_count.txt',
        OUTPUT_DIR + '/08-rsem_calculate_expression/isoform_fpkm.txt',
    benchmark:
        OUTPUT_DIR + '/benchmarks/count_fpkm_after_dedup.benchmark.txt'
    log:
        OUTPUT_DIR + '/count_fpkm/.log/count_fpkm_after_dedup.log'
    params:
        gene_count_script = srcdir('../scripts/combine_gene_count.py'),
        gene_fpkm_script = srcdir('../scripts/combine_gene_FPKM.py'),
        isoform_count_script = srcdir('../scripts/combine_isoform_count.py'),
        isoform_fpkm_script = srcdir('../scripts/combine_isoform_FPKM.py'),
    shell:'''(
pushd {OUTPUT_DIR}/08-rsem_calculate_expression/
python {params.gene_count_script} `ls *.genes.results`
python {params.gene_fpkm_script} `ls *.genes.results`
python {params.isoform_count_script} `ls *.isoforms.results`
python {params.isoform_fpkm_script} `ls *.isoforms.results`
popd)2>&1 | tee {log}
    '''
