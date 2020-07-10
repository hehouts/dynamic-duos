

input_data = "vir_SAMPLEIDs.txt"

samples= [x.strip().split(".tar")[0] for x in open(input_data, 'r')]

#  #takes first 5 samples
#  samples=samples[:5]

#print(samples)

rule all:
    input: 
        # uncompress_genome
        #expand("raw_data/zipped_data/{sample}_{ext}.bz2", sample =samples, ext = [1,2])
        #expand("raw_data/zipped_data/{sample}_1_sequence.txt.bz2", sample =samples),
        #expand("raw_data/zipped_data/{sample}_2_sequence.txt.bz2", sample =samples),
        #expand("fastp/{sample}_1.trimmed.gz", sample=samples),
        #expand("fastp/{sample}_2.trimmed.gz", sample=samples),
        #expand("sourmash/sig/{sample}.sig", sample=samples),
        expand("sourmash/gather/{sample}_gather.csv",sample=samples)

rule download_data:
    output: "raw_data/zipped_data/{sample}.tar"
    params: URL= lambda wildcards: "https://ibdmdb.org/tunnel/static/HMP2/Viromics/1732/" + wildcards.sample + ".tar"
    shell: "wget -O {output} {params.URL}"


rule uncompress_genome:
    input: "raw_data/zipped_data/{sample}.tar"
    output: 
        "raw_data/zipped_data/{sample}_1_sequence.txt.bz2",
        "raw_data/zipped_data/{sample}_2_sequence.txt.bz2" 
    shell:
        "tar xvf {input} -C raw_data/zipped_data/"
rule bzip_to_gzip:
    input: 
        F1="raw_data/zipped_data/{sample}_1_sequence.txt.bz2",
        F2="raw_data/zipped_data/{sample}_2_sequence.txt.bz2"
    output:
        F1="raw_data/zipped_data/{sample}_1.gz", 
        F2="raw_data/zipped_data/{sample}_2.gz"
    log:
        "logs/gzip/{sample}_gzip.log"  
    benchmark:
        "logs/gzip/{sample}_gzip.benchmark"  
    resources:
        # minutes to allocate for sbatch in slurm
        runtime = 20,
        # memory to allocate in mb
        mem_mb = 1000
    threads: 1
    shell:
        """
        bzcat {input.F1} | gzip -c -9 > {output.F1} 2> {log}
        bzcat {input.F2} | gzip -c -9 > {output.F2} 2>> {log}
        """

rule fastp:
    input:
        R1="raw_data/zipped_data/{sample}_1.gz", 
        R2="raw_data/zipped_data/{sample}_2.gz"
    output:
        R1="fastp/{sample}_1.trimmed.gz", 
        R2="fastp/{sample}_2.trimmed.gz",
        html="fastp/{sample}_fastp.html",
        json="fastp/{sample}_fastp.json"
    log:
        "logs/fastp/{sample}_fastp.log"  
    benchmark:
        "logs/fastp/{sample}_fastp.benchmark"  
    conda:
        "virHMP_env.yml"
    shell:
        """
        fastp -i {input.R1} -I {input.R2}  \
        -o {output.R1} -O {output.R2} \
        -h {output.html} -j {output.json} > {log}
        """

rule sourmash_compute:
    input:
        R1="fastp/{sample}_1.trimmed.gz",
        R2="fastp/{sample}_2.trimmed.gz",
    output:
        "sourmash/sig/{sample}.sig"
    log:
        "logs/sourmash/sig/{sample}_sig.log"
    benchmark:
        "logs/sourmash/sig/{sample}_sig.benchmark"
    conda:
        "virHMP_env.yml"
    shell:"""
        sourmash compute -k 21,31,51 --scaled 500 --merge {wildcards.sample} --track-abundance -o {output} {input.R1} {input.R2}
    """

rule download_gather_genbank:
    output: "gather_databases/genbank-d2-k31.tar.gz"
    params:
        dl_link="https://s3-us-west-2.amazonaws.com/sourmash-databases/2018-03-29/genbank-d2-k31.tar.gz"
    shell:
        """
        wget -O {output} {params.dl_link}
        """
rule untar_genbank:
    input:  "gather_databases/genbank-d2-k31.tar.gz"    
    output: "gather_databases/genbank-d2-k31.sbt.json"
    params: 
        outdir = "gather_databases"
    shell: 
        """
        tar xf {input} -C {params.outdir}
        """


rule sourmash_gather:
    input:
        query="sourmash/sig/{sample}.sig"
        database="gather_databases/genbank-d2-k31.sbt.json"
    output:
        csv="sourmash/gather/{sample}_gather.csv",
    # add these, include in shell using `sourmash --help`
    # (also add other databases)
    #specify scaled and k size
        #matches="",
        #unassigned=""
    log:
        "logs/sourmash/gather/{sample}_gather.log"
    benchmark:
        "logs/sourmash/gather/{sample}_gather.benchmark"
    conda:
        "virHMP_env.yml"
    shell:
        """
        sourmash gather {input.query} {input.database} -o {output.csv}
        """




