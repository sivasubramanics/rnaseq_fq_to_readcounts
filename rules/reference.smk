# Reference genome preparation
rule star_index:
    input:
        fasta=config["reference"]["fasta"],
        gtf=config["reference"]["gtf"]
    output:
        "work/star_index/Genome"
    params:
        sjdbOverhang=config["sjdbOverhang"],
        genomeSAindexNbases=config["genomeSAindexNbases"]
    threads:
        config["threads"]
    log:
        "logs/star_index.log"
    benchmark:
        "benchmarks/star_index.txt"
    conda:
        "../envs/reference.yaml"
    shell:
        ("STAR --runThreadN {threads} \
         --runMode genomeGenerate \
         --genomeDir work/star_index/ \
         --genomeFastaFiles {input.fasta} \
         --sjdbGTFfile {input.gtf} \
         --sjdbOverhang {params.sjdbOverhang} \
         --genomeSAindexNbases {params.genomeSAindexNbases} \
         > {log} 2>&1")

rule samtools_faidx:
    input:
        config["reference"]["fasta"]
    output:
        config["reference"]["fasta"] + ".fai"
    log:
        "logs/samtools_faidx.log"
    benchmark:
        "benchmarks/samtools_faidx.txt"
    conda:
        "../envs/reference.yaml"
    shell:
        "samtools faidx {input.ref_fasta} > {log} 2>&1"

rule star_index_two:
    input:
        fasta=config["reference"]["fasta"],
        gtf=config["reference"]["gtf"],
        junctions="work/alignment/junctions.txt"
    output:
        "work/star_index_2/Genome"
    params:
        sjdbOverhang=config["sjdbOverhang"],
        genomeSAindexNbases=config["genomeSAindexNbases"]
    threads:
        config["threads"]
    log:
        "logs/star_index_2.log"
    benchmark:
        "benchmarks/star_index_2.txt"
    conda:
        "../envs/reference.yaml"
    shell:
        ("STAR --runThreadN {threads} \
         --runMode genomeGenerate \
         --genomeDir work/star_index_2/ \
         --genomeFastaFiles {input.fasta} \
         --sjdbGTFfile {input.gtf} \
         --sjdbFileChrStartEnd {input.junctions} \
         --sjdbOverhang {params.sjdbOverhang} \
         --genomeSAindexNbases {params.genomeSAindexNbases} \
         > {log} 2>&1")
