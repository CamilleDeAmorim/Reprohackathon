Samples=["SRX7080612","SRX7080613","SRX7080614","SRX7080615","SRX7080616","SRX7080617"]
reference_suffixes = ["1", "2", "3", "4", "rev.1", "rev.2"]

rule all: 
    input:
        expand("results/01_raw_data/{SRA_id}.fastq.gz", SRA_id=Samples),
        expand("results/02_Trimming_results/{SRA_id}_trimmed.fq.gz", SRA_id=Samples),
        "results/03_Reference_Genome/reference.fasta", 
        "results/04_Genome_Annotation/reference.gff",
        expand("results/03_Reference_Genome/reference.{suffixe}.ebwt", suffixe=reference_suffixes),
        expand("results/05_Mapping_results/{SRA_id}.bam", SRA_id=Samples),
        expand("results/05_Mapping_results/{SRA_id}.bai", SRA_id=Samples),
        "results/06_Counting_results/counts.txt"
             
# Rule that downloads the FASTQ files.
rule fasterq_dump:
    output:
        "results/01_raw_data/{SRA_id}.fastq.gz"
    container:
        "images/fasterq-dump.img"
    log:
        "logs/fasterq-dump/{SRA_id}.log"
    threads: 60
    shell:
        """
        fasterq-dump --threads {threads} --progress {wildcards.SRA_id} -O results/01_raw_data/ &> {log}
        gzip -f results/01_raw_data/{wildcards.SRA_id}.fastq
        """

# Rule that performs the trimming of FASTQ files.
rule trim_galore:
    input:
        "results/01_raw_data/{SRA_id}.fastq.gz"
    output: 
        trimmed_fastq="results/02_Trimming_results/{SRA_id}_trimmed.fq.gz",
        reports="results/02_Trimming_results/{SRA_id}_trimming_report.txt"
    container:
        "images/TrimGalor.img"
    log:
        "logs/TrimGalor/{SRA_id}_trimmed.log"
    threads: 60
    shell:
        """
        trim_galore -q 20 --phred33 --length 25 {input} -o results/02_Trimming_results/ &> {log}
        mv results/02_Trimming_results/{wildcards.SRA_id}.fastq.gz_trimming_report.txt {output.reports}

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
        "images/bowtie.img"
    shell: 
        """
        bowtie-build {input} results/03_Reference_Genome/reference
        """

# Rule that maps the genome and creates an index
rule mapping:
    input:
        index =  expand("results/03_Reference_Genome/reference.{suffixe}.ebwt", suffixe=reference_suffixes),
        fastq_files = "results/02_Trimming_results/{SRA_id}_trimmed.fq.gz"
    output:
        bam = "results/05_Mapping_results/{SRA_id}.bam",
        bai = "results/05_Mapping_results/{SRA_id}.bai"
    container:
        "images/bowtie_samtools.img"
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
        "images/featureCounts.img"
    log: 
        "logs/counting.log"
    threads: 40
    shell:
        """
        featureCounts -t gene -g ID -s 1 -T {threads} -a {input.genome_annotation} -o {output} {input.bamf}
        """

