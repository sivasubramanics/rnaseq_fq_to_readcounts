# Snakemake pipeline to process the RNAseq data and do read count per gene and isoforms based on STAR alignment
from snakemake.utils import min_version
min_version('8.0')

from collections import defaultdict

configfile: "config/config.yaml"

include: "rules/common.smk"
include: "rules/fastp.smk"
include: "rules/reference.smk"
include: "rules/alignment.smk"
include: "rules/quantify.smk"
include: "rules/sortmerna.smk"
# include: "rules/merge_counts.smk"


rule all:
    input:
        "work/count_table/feature_counts_parsed.txt",
        "work/count_table/htseq_counts_parsed.txt"
