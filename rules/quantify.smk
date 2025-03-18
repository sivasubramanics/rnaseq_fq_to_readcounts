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
        ("cut -f 1,7- {input} | \
            tail -n +2 | \
            sed 's/work\/merged_bams\///g' | \
            sed 's/\.bam//g' > {output} 2> {log}")

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

rule parse_htseq_counts:
    input:
        "work/count_table/htseq_counts.txt"
    output:
        "work/count_table/htseq_counts_parsed.txt"
    log:
        "logs/parse_htseq_counts.log"
    shell:
        ("""sed 's/work\/merged_bams\///g' {input} | \
            sed 's/\.bam//g' | \
            awk 'BEGIN{{FS=OFS=\"\\t\"}}{{if(NR==1){{print \"Geneid\",$0}}else{{print}}}}' > {output} 2> {log}""")
