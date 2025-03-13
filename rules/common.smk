# Common functions for the RNA-seq pipeline
from collections import defaultdict

def parse_metadata(metadata_file):
    """
    Parse the metadata file and return a dictionary with the sample name as key and a dictionary with the sample information as value
    """
    data = defaultdict(dict)
    with open(metadata_file) as f:
        header = f.readline().strip().split()
        print(header)
        if header != ["run", "sample", "fq1", "fq2"]:
            raise ValueError("Invalid metadata file. The header should be 'run', 'sample', 'fq1', 'fq2'")
        for line in f:
            fields = line.strip().split()
            if len(fields) != 4:
                raise ValueError("Invalid metadata file. Each line should have 4 fields")
            data[fields[1]][fields[0]] = fields[2:]
        print(data)
    return data

def get_sjdbs(data):
    """
    Return the SJ.out.tab files for the sample and run
    """
    sjlist = []
    for sample in data:
        sjdbfiles=expand("work/alignment/pass_one/{sample}.{run}.SJ.out.tab", sample=sample, run=[x for x in data[sample]])
        sjlist.extend(sjdbfiles)
    return sjlist


data = parse_metadata(config["metadata"])
