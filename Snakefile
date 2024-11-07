Samples = ["SRX7080612", "SRX708061", "SRX7080614", "SRX7080615", "SRX7080616", "SRX7080617"]
reference_suffixes = ["1", "2", "3", "4", "rev.1", "rev.2"]

# Rule that downloads all necessary files.
rule all: 
    input:
        expand("results/01_raw_data/{SRA_id}.fastq.gz", SRA_id=Samples),
        expand("results/02_Trimming_results/{SRA_id}_trimmed.fq.gz", SRA_id=Samples),
        "results/03_Reference_Genome/reference.fasta", 
        "results/04_Genome_Annotation/reference.gff",
        expand("results/03_Reference_Genome/reference.{suffixe}.ebwt", suffixe=reference_suffixes),
        expand("results/05_mapping/{SRA_id}.sam", SRA_id=Samples)
        

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
    threads: 40
    shell:
        """
        trim_galore -q 20 --phred33 --length 25 {input} -o results/02_Trimming_results/ &> {log}
        
        """
#mv results/02_Trimming_results/{wildcards.SRA_id}.fastq.gz_trimming_report.txt {output.reports}
# Rule that downloads the reference genome in FASTA format.
rule reference_genome:
    output:
        "results/03_Reference_Genome/reference.fasta"
    shell:
        """
        wget -q -O {output} "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=CP000253.1&rettype=fasta"
        """

# Rule that downloads the genome annotation."
rule genome_annotation:
    output:
        "results/04_Genome_Annotation/reference.gff"
    log: 
        "logs/Genome_Annotation"
    shell:
        """
        wget -O results/04_Genome_Annotation/reference.gff "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id=CP000253.1"  &> {log}       
        """

# Rule that creates the genome index"
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

# Rule that map the genome
# attention, à finir avec les commandes samtools
rule mapping:
    input:
        index =  expand("results/03_Reference_Genome/reference.{suffixe}.ebwt", suffixe=reference_suffixes),
        fastq_files = "results/02_Trimming_results/{SRA_id}_trimmed.fq.gz"
    output:
        "results/05_mapping/{SRA_id}.sam"
    container:
        "images/bowtie_samtools.img"
    threads: 40
    shell: 
        """
        bowtie -p 8 -S results/03_Reference_Genome/reference <(gunzip -c {input.fastq_files}) > {output}
        """

#attention, non vérifiée ! 
rule counting:
    input:
        bamf = expand("results/O5_mapping/{SRA_id}.bam", SRA_id=Samples),
        genome_annotation = "results/04_Genome_Annotation/reference.gff"
    output:
        "results/counts.txt"
    container:
        "images/featureCounts.img"
    log: 
        "logs/counting.log"
    threads: 40
    shell:
        """
        featureCounts --extraAttributes Name -t gene -g ID -F GTF -T {threads} -a {input.genome_annotation} -o {output} {input.bamf}
        """
