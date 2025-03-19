# Snakemake rule to run SortMeRNA on the filtered fastq files
def get_rRNA_references(ref_db_dir):
    list_files = []
    for file in os.listdir(ref_db_dir):
        if file.endswith(".fasta") or file.endswith(".fa") or file.endswith(".fna"):
            list_files.append(os.path.join(ref_db_dir, file))
    rRNAs_ref = ""
    for _file in list_files:
        if rRNAs_ref == "":
            rRNAs_ref = "--ref " + _file
        else:
            rRNAs_ref += " --ref " + _file
    print(rRNAs_ref)
    return rRNAs_ref


rule sortmerna:
    input:
        "work/fastp/{sample}.{run}.1.fastq.gz",
        "work/fastp/{sample}.{run}.2.fastq.gz"
    output:
        "work/sortmerna/{sample}.{run}.1.fastq.gz",
        "work/sortmerna/{sample}.{run}.2.fastq.gz",
        "work/sortmerna/{sample}.{run}.rrna.1.fastq.gz",
        "work/sortmerna/{sample}.{run}.rrna.2.fastq.gz"
    params:
        workdir = "work/sortmerna/{sample}.{run}",
        ref = get_rRNA_references("data/rRNA_databases/")
    threads:
        config["threads"]
    log:
        "logs/sortmerna/{sample}.{run}.log"
    benchmark:
        "benchmarks/sortmerna/{sample}.{run}.txt"
    conda:
        "../envs/sortmerna.yaml"
    shell:
        ("mkdir -p {params.workdir} && "
         "sortmerna {params.ref} --reads {input[0]} --reads {input[1]} "
         "--workdir {params.workdir} --threads {threads} --paired_in --out2 --other "
         "--fastx > {log} 2>&1 && "
            "mv {params.workdir}/out/other_fwd.fq.gz {output[0]} && "
            "mv {params.workdir}/out/other_rev.fq.gz {output[1]} && "
            "mv {params.workdir}/out/aligned_fwd.fq.gz {output[2]} && "
            "mv {params.workdir}/out/aligned_rev.fq.gz {output[3]} && "
            "rm -rf {params.workdir}")
