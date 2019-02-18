# Apricot
Iwase lab RNA-seq pipeline

(1 paragraph overview)

Quickstart
----------

- Setup new machine or new user
  - Clone from github
  - Build conda environment
  - Review versions

- Add a new genome reference
  - Download
  - Add ERCC spike-ins
    - https://www.thermofisher.com/order/catalog/product/4456739

- Setup project
  - inputs

- Setup configfile
  - reference, input, output
  - Adding samples/reads
    - Paired end vs Single end
  - Sequencing params

- Run workflow
  - cd to project dir
  - Setup screen
    - Launch screen session
    - To detach
    - To reattach
  - Show rule graph
    ```
    snakemake \
        --snakefile ../Apricot/Snakefile \
        --configfile config.project_vallianatos.mm10_ercc.yaml \
        --rulegraph | dot -Tpdf > rulegraph.pdf
    ```
  - Run workflow
    - Dry-run to check rules and jobs that will be executed:
    ```
    snakemake \
        --snakefile ../Apricot/Snakefile \
        --configfile config.project_vallianatos.mm10_ercc.yaml \
        --dryrun --quiet
    ```

    - Run
      - Set memory, cores
      ```
      cd {project}
      snakemake \
          --snakefile ../Apricot/Snakefile \
          --configfile config.project_vallianatos.mm10_ercc.yaml \
          --resources mem_gb=105 \
          --cores 36 \
          -p
      ```
    - Run on Flux

- Key outputs
  - Deliverables
  - Logs
    - Snakemake run logs: all details as printed to screen. Note if running jobs
      in parallel, job details will sometimes be interleaved.
    ```
    cd {your project dir}
    cd .snakemake/log
    ls #list of logs ordered by run timestamp
    less `ls -tr | tail -1` #log details of last run
    ```
    - Rule logs
    - Benchmarks
      - All jobs:
      ```
      awk 'BEGIN {OFS="\t"} NR==1 {print "filename", $0} FNR>1 {print FILENAME, $0}' benchmarks/*.txt
      ```
      - Slowest 20 jobs:
        ```
        awk 'BEGIN {OFS="\t"} NR==1 {print "filename", $0} FNR>1 {print FILENAME, $0}' benchmarks/*.txt | sort -k2,2nr | head -20
        ```

Citations
---------

- Genome references
  - Ensembl

- Tools
  - Ballgown: Fu J, Frazee AC, Collado-Torres L, Jaffe AE, Leek JT (2019). ballgown: Flexible, isoform-level differential expression analysis. R package version 2.14.1.
  - fastqc: Andrews S. (2010). FastQC: a quality control tool for high throughput  sequence data. Available online at: http://www.bioinformatics.babraham.ac.uk/projects/fastqc
  - multiqc: MultiQC: Summarize analysis results for multiple tools and samples in a single report. Philip Ewels, Måns Magnusson, Sverker Lundin and Max Käller. Bioinformatics (2016). doi: 10.1093/bioinformatics/btw354. PMID: 27312411
  - samtools: Li H., Handsaker B., Wysoker A., Fennell T., Ruan J., Homer N., et al. (2009). The sequence alignment / map format and SAMtools. Bioinformatics 25 2078–2079.
  - snakemake: Köster, Johannes and Rahmann, Sven. “Snakemake - A scalable bioinformatics workflow engine”. Bioinformatics 2012.
  - STAR: Dobin, A et al. STAR: ultrafast universal RNA-seq aligner. Bioinformatics29, 15–21 (2013).
  - Stringtie: Pertea M, Pertea GM, Antonescu CM, Chang TC, Mendell JT	& Salzberg SL. StringTie enables improved reconstruction of a transcriptome from RNA-seq reads Nature Biotechnology 2015, doi:10.1038/nbt.3122
  - UMITools: Smith TS, Heger A, Sudbery I (2017) UMI-tools: modelling sequencing errors in Unique Molecular Identifiers to improve quantification accuracy. Genome Res 27: 491 – 499.

Basic orientation
-----------------
- Snakemake
  - https://snakemake.readthedocs.io/en/stable/index.html
  - https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html

- Conda
  - https://conda.io/docs/
