
# rule to map reads to the reference genome using STAR
rule pass_one:
    input:
        fq1="work/sortmerna/{sample}.{run}.1.fastq.gz",
        fq2="work/sortmerna/{sample}.{run}.2.fastq.gz",
        genome="work/star_index/Genome"
    output:
        "work/alignment/pass_one/{sample}.{run}.Aligned.out.bam",
        "work/alignment/pass_one/{sample}.{run}.SJ.out.tab"
    threads:
        config["threads"]
    resources:
        tmpdir=config["tmpdir"]
    log:
        "logs/alignment/pass_one/{sample}.{run}.log"
    benchmark:
        "benchmarks/alignment/pass_one/{sample}.{run}.txt"
    conda:
        "../envs/alignment.yaml"
    shell:
        ("STAR \
            --runThreadN {threads} \
            --genomeDir work/star_index/ \
            --readFilesIn {input.fq1} {input.fq2} \
            --outFileNamePrefix work/alignment/pass_one/{wildcards.sample}.{wildcards.run}. \
            --outSAMtype BAM Unsorted \
            --readFilesCommand gunzip -c \
            --outSAMunmapped Within > {log} 2>&1")

rule filter_junctions:
    input:
        get_sjdbs(data)
    output:
        "work/alignment/junctions.txt"
    shell:
        ("cat {input} | awk '{{ if ($7 >= 3) print $0}}' | sort -k1,1 -k2,2n | uniq > {output}")

rule pass_two:
    input:
        fq1="work/sortmerna/{sample}.{run}.1.fastq.gz",
        fq2="work/sortmerna/{sample}.{run}.2.fastq.gz",
        genome="work/star_index_2/Genome"
    output:
        "work/alignment/pass_two/{sample}.{run}.Aligned.out.bam"
    threads:
        config["threads"]
    resources:
        tmpdir=config["tmpdir"]
    log:
        "logs/alignment/pass_two/{sample}.{run}.log"
    benchmark:
        "benchmarks/alignment/pass_two/{sample}.{run}.txt"
    conda:
        "../envs/alignment.yaml"
    shell:
        ("STAR \
            --runThreadN {threads} \
            --genomeDir work/star_index_2/ \
            --readFilesIn {input.fq1} {input.fq2} \
            --outFileNamePrefix work/alignment/pass_two/{wildcards.sample}.{wildcards.run}. \
            --outSAMtype BAM Unsorted \
            --readFilesCommand gunzip -c \
            --outSAMunmapped Within > {log} 2>&1")

rule merge_and_sort_bams:
    input:
        lambda wildcards: expand("work/alignment/pass_two/{sample}.{run}.Aligned.out.bam",
                                 sample=wildcards.sample,
                                 run=data[wildcards.sample].keys())
    output:
        "work/merged_bams/{sample}.bam"
    log:
        "logs/merge_and_sort_bams/{sample}.log"
    threads:
        config["threads"]
    resources:
        tmpdir=config["tmpdir"]
    benchmark:
        "benchmarks/merge_and_sort_bams/{sample}.txt"
    conda:
        "../envs/alignment.yaml"
    shell:
        """
        samtools merge -@ {threads} -u - {input} | \
        samtools sort -@ {threads} -o {output} > {log} 2>&1
        """

rule index_bam:
    input:
        "work/merged_bams/{sample}.bam"
    output:
        "work/merged_bams/{sample}.bam.bai"
    log:
        "logs/index_bam/{sample}.log"
    threads:
        config["threads"]
    benchmark:
        "benchmarks/index_bam/{sample}.txt"
    conda:
        "../envs/alignment.yaml"
    shell:
        "samtools index {input} > {log} 2>&1"
