ALL.append([OUTPUT_DIR + '/03-rsem_star_align/gene_count.txt',
            OUTPUT_DIR + '/03-rsem_star_align/gene_fpkm.txt',
            OUTPUT_DIR + '/03-rsem_star_align/isoform_count.txt',
            OUTPUT_DIR + '/03-rsem_star_align/isoform_fpkm.txt',
            ])
rule count_fpkm:
    input:
        genes = expand(OUTPUT_DIR + '/03-rsem_star_align/{sample}.rsem.genes.results',
                      sample=config['samples']),
        isoforms = expand(OUTPUT_DIR + '/03-rsem_star_align/{sample}.rsem.isoforms.results',
                      sample=config['samples']),
    output:
        OUTPUT_DIR + '/03-rsem_star_align/gene_count.txt',
        OUTPUT_DIR + '/03-rsem_star_align/gene_fpkm.txt',
        OUTPUT_DIR + '/03-rsem_star_align/isoform_count.txt',
        OUTPUT_DIR + '/03-rsem_star_align/isoform_fpkm.txt',
    benchmark:
        OUTPUT_DIR + '/benchmarks/count_fpkm.benchmark.txt'
    log:
        OUTPUT_DIR + '/count_fpkm/.log/count_fpkm.log'
    params:
        gene_count_script = srcdir('../scripts/combine_gene_count.py'),
        gene_fpkm_script = srcdir('../scripts/combine_gene_FPKM.py'),
        isoform_count_script = srcdir('../scripts/combine_isoform_count.py'),
        isoform_fpkm_script = srcdir('../scripts/combine_isoform_FPKM.py'),
    shell:'''(
pushd {OUTPUT_DIR}/03-rsem_star_align/
python {params.gene_count_script} `ls *.rsem.genes.results`
python {params.gene_fpkm_script} `ls *.rsem.genes.results`
python {params.isoform_count_script} `ls *.rsem.isoforms.results`
python {params.isoform_fpkm_script} `ls *.rsem.isoforms.results`
popd)2>&1 | tee {log}
    '''
