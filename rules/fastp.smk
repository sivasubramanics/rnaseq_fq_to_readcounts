rule fastp:
    input:
        lambda wildcards: data[wildcards.sample][wildcards.run]
    output:
        "work/fastp/{sample}.{run}.1.fastq.gz",
        "work/fastp/{sample}.{run}.2.fastq.gz",
        "work/fastp/{sample}.{run}.html",
        "work/fastp/{sample}.{run}.json",
    params:
        min_len=config["min_len"],
        min_qual=config["min_qual"],
    threads:
        config["threads"]
    log:
        "logs/fastp/{sample}.{run}.log"
    benchmark:
        "benchmarks/fastp/{sample}.{run}.txt"
    conda:
        "../envs/fastp.yaml"
    shell:
        ("fastp --in1 {input[0]} --in2 {input[1]} \
        --out1 {output[0]} --out2 {output[1]} \
        --html {output[2]} --json {output[3]} \
        --length_required {params.min_len} \
        --qualified_quality_phred {params.min_qual} \
        --unqualified_percent_limit 40 \
        --thread {threads} \
        --detect_adapter_for_pe \
        --cut_tail \
        --cut_tail_window_size 1 \
        --cut_tail_mean_quality 20 \
        --cut_front \
        --cut_front_window_size 1 \
        --cut_front_mean_quality 20 \
        --correction > {log} 2>&1")
