rule feature_counts:
    input:
        bam=expand("work/merged_bams/{sample}.bam", sample=data.keys()),
        gtf=config["reference"]["gtf"]
    output:
        "work/count_table/feature_counts.txt"
    threads:
        config["threads"]
    log:
        "logs/feature_counts.log"
    benchmark:
        "benchmarks/feature_counts.txt"
    conda:
        "../envs/quantify.yaml"
    shell:
        ("featureCounts -T {threads} \
         -a {input.gtf} \
         -o {output} \
         -s 0 -p -t exon -g gene_id \
         {input.bam} > {log} 2>&1")

rule parse_feature_counts:
    input:
        "work/count_table/feature_counts.txt"
    output:
        "work/count_table/feature_counts_parsed.txt"
    log:
        "logs/parse_feature_counts.log"
    shell:
        ("cut -f 1,7- {input} > {output}")

rule htseq_counts:
    input:
        bam=expand("work/merged_bams/{sample}.bam", sample=data.keys()),
        gtf=config["reference"]["gtf"]
    output:
        "work/count_table/htseq_counts.txt"
    threads:
        config["threads"]
    log:
        "logs/htseq_counts.log"
    benchmark:
        "benchmarks/htseq_counts.txt"
    conda:
        "../envs/quantify.yaml"
    shell:
        ("htseq-count \
            --format bam \
            --order pos \
            --stranded no \
            --type exon \
            --idattr gene_id \
            --mode intersection-nonempty \
            --nprocesses {threads} \
            --with-header \
         {input.bam} {input.gtf} > {output} 2> {log}")
