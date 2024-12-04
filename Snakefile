Samples=["SRX7080612","SRX7080613","SRX7080614","SRX7080615","SRX7080616","SRX7080617"]
reference_suffixes = ["1", "2", "3", "4", "rev.1", "rev.2"]

rule all: 
    input:
        expand("results/01_raw_data/{SRA_id}.fastq.gz", SRA_id=Samples),
        expand("results/02_Trimming_results/{SRA_id}_trimmed.fastq.gz", SRA_id=Samples),
        "results/03_Reference_Genome/reference.fasta", 
        "results/04_Genome_Annotation/reference.gff",
        expand("results/03_Reference_Genome/reference.{suffixe}.ebwt", suffixe=reference_suffixes),
        expand("results/05_Mapping_results/{SRA_id}.bam", SRA_id=Samples),
        expand("results/05_Mapping_results/{SRA_id}.bai", SRA_id=Samples),
        "results/06_Counting_results/counts.txt",
        "assets/tRNA_synthetases_genes.txt",
        "assets/translation_genes.txt",
        "results/07_Final_results/DESeq2_Results/MA_plot.png",
        "results/07_Final_results/DESeq2_Results/MA_plot_translation1.png",
        "results/07_Final_results/DESeq2_Results/MA_plot_translation2.png",
        "results/07_Final_results/DESeq2_Results/Volcano_plot.png"
             
# Rule that downloads the FASTQ files.
rule fasterq_dump:
    output:
        "results/01_raw_data/{SRA_id}.fastq.gz"
    container:
        "https://zenodo.org/records/13994690/files/fasterq-dump.img?download=1"
    log:
        "logs/fasterq-dump/{SRA_id}.log"
    threads: 20
    shell:
        """
        fasterq-dump --threads {threads} --progress {wildcards.SRA_id} -O results/01_raw_data/ &> {log}
        gzip -f results/01_raw_data/{wildcards.SRA_id}.fastq
        """

# Rule that performs the trimming of FASTQ files.
rule cutadapt:
    input:
        "results/01_raw_data/{SRA_id}.fastq.gz"
    output: 
        trimmed_fastq="results/02_Trimming_results/{SRA_id}_trimmed.fastq.gz",
    container:
        "https://zenodo.org/records/13994994/files/cutadapt_v1.11.img?download=1"
    log:
        "logs/cutadapt/{SRA_id}_trimmed.log"
    threads: 20
    shell:
        """
        cutadapt -o results/02_Trimming_results/{wildcards.SRA_id}_trimmed.fastq.gz {input} -a AGATCGGAAGAGC -q 20 -m 25 > results/02_Trimming_results/{wildcards.SRA_id}_trimming_report.txt
        """
# Rule that downloads the reference genome in FASTA format.
rule reference_genome:
    output:
        "results/03_Reference_Genome/reference.fasta"
    shell:
        """
        wget -q -O  results/03_Reference_Genome/reference.fasta "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=CP000253.1&rettype=fasta"
        """
# Rule that downloads the genome annotation.
rule genome_annotation :
    output:
        "results/04_Genome_Annotation/reference.gff"
    log: 
        "logs/Genome_Annotation"
    shell:
        """
        wget -O {output} "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id=CP000253.1"  &> {log}       
        """


# Rule that creates the genome index
rule genome_index:
    input:
        "results/03_Reference_Genome/reference.fasta"
    output:
        "results/03_Reference_Genome/reference.{suffixe}.ebwt"
    container:
        "https://zenodo.org/records/13994994/files/bowtie_v0.12.7_samtools.img?download=1"
    shell: 
        """
        bowtie-build {input} results/03_Reference_Genome/reference
        """

# Rule that maps the genome and creates an index
rule mapping:
    input:
        index =  expand("results/03_Reference_Genome/reference.{suffixe}.ebwt", suffixe=reference_suffixes),
        fastq_files = "results/02_Trimming_results/{SRA_id}_trimmed.fastq.gz"
    output:
        bam = "results/05_Mapping_results/{SRA_id}.bam",
        bai = "results/05_Mapping_results/{SRA_id}.bai"
    container:
        "https://zenodo.org/records/13994994/files/bowtie_v0.12.7_samtools.img?download=1"
    threads: 8
    shell: 
        """
        bowtie -p {threads} -S results/03_Reference_Genome/reference <(gunzip -c {input.fastq_files}) | samtools sort -@ {threads} -o {output.bam}
        samtools index -@ {threads} -o {output.bai} {output.bam}
        """

# Rule that counts the number of reads mapped on each gene.
rule counting:
    input:
        bamf = expand("results/05_Mapping_results/{SRA_id}.bam", SRA_id=Samples),
        genome_annotation = "results/04_Genome_Annotation/reference.gff"
    output:
        "results/06_Counting_results/counts.txt"
    container:
        "https://zenodo.org/records/13994994/files/featureCounts_v1.4.6-p3.img?download=1"
    log: 
        "logs/counting.log"
    threads: 20
    shell:
        """
        featureCounts -t gene -g ID -s 1 -T {threads} -a {input.genome_annotation} -o {output} {input.bamf}
        """
rule Download_KEGG_Genes:
    output:
        "assets/tRNA_synthetases_genes.txt",
        "assets/translation_genes.txt"
    container:
        "images/R_new_version.img"
    shell:
        """
        Rscript bin/Download_KEGG_Genes.R
        """
rule DESeq2_analysis:
    input:
        counts="results/06_Counting_results/counts.txt",
        translation_genes="assets/translation_genes.txt",
        tRNA_genes="assets/tRNA_synthetases_genes.txt"
    output:
        ma_plot="results/07_Final_results/DESeq2_Results/MA_plot.png",
        ma_plot_translation1="results/07_Final_results/DESeq2_Results/MA_plot_translation1.png",
        ma_plot_translation2="results/07_Final_results/DESeq2_Results/MA_plot_translation2.png",
        volcano_plot="results/07_Final_results/DESeq2_Results/Volcano_plot.png"
    container:
        "images/DESeq2_v1.16.1.img"
    shell:
        """
        Rscript bin/DESeq2_analysis.R {input.counts} {input.translation_genes} {input.tRNA_genes} results/07_Final_results/DESeq2_Results/
        """
