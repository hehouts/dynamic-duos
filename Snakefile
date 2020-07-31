

input_data = "vir_SAMPLEIDs.txt"

samples= [x.strip().split(".tar")[0] for x in open(input_data, 'r')]

k_sizes= ["21","31","51"]

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
        #expand("sourmash/compare/virHMP_compare_k{ksize}_unfiltered.csv", ksize=k_sizes),
        #expand("sourmash/gather/{sample}_gather.csv",sample=samples)
        expand("sourmash/plots/virHMP_compare_k{ksize}_unfiltered.numpy.dendro.pdf", ksize=k_sizes),
        expand("sourmash/plots/virHMP_compare_k{ksize}_unfiltered.numpy.hist.pdf", ksize=k_sizes),
        expand("sourmash/plots/virHMP_compare_k{ksize}_unfiltered.numpy.matrix.pdf", ksize=k_sizes)


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

#### Sourmashy Stuff ####

rule sourmash_compute:
    input:
        R1="fastp/{sample}_1.trimmed.gz",
        R2="fastp/{sample}_2.trimmed.gz",
    output:
        "sourmash/sig/{sample}.sig"
    params:
        ksizes=",".join(k_sizes)
    log:
        "logs/sourmash/sig/{sample}_sig.log"
    benchmark:
        "logs/sourmash/sig/{sample}_sig.benchmark"
    conda:
        "virHMP_env.yml"
    shell:"""
        sourmash compute -k {params.ksizes} --scaled 500 --merge {wildcards.sample} --track-abundance -o {output} {input.R1} {input.R2}
    """

rule sourmash_compare:
    input: expand("sourmash/sig/{sample}.sig", sample=samples)
    output:
        csv="sourmash/compare/virHMP_compare_k{ksize}_unfiltered.csv",
        numpy="sourmash/compare/virHMP_compare_k{ksize}_unfiltered.numpy"
        
    shell:"""
        sourmash compare {input} -k {wildcards.ksize} -o {output.numpy} --csv {output.csv}
    """

### add a sm plot rule, takes numpy output of ^ and makes a heatmap (pdf) ( sourmash plot --help)

rule sourmash_plot:
#    input: expand("sourmash/compare/virHMP_compare_k{ksize}_unfiltered.numpy", ksize=k_sizes)
    input: "sourmash/compare/virHMP_compare_k{ksize}_unfiltered.numpy"
    output:
        "sourmash/plots/virHMP_compare_k{ksize}_unfiltered.numpy.dendro.pdf",
        "sourmash/plots/virHMP_compare_k{ksize}_unfiltered.numpy.hist.pdf",
        "sourmash/plots/virHMP_compare_k{ksize}_unfiltered.numpy.matrix.pdf"

    shell:"""
        sourmash plot --pdf --labels {input} --output sourmash/plots/
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


rule download_gather_almeida:
    output: "gather_databases/almeida-mags-k31.tar.gz"
    params:
        dl_link="https://osf.io/5jyzr/download"
    shell:
        """
        wget -O {output} {params.dl_link} 
        """
rule untar_almeida:
    input:"gather_databases/almeida-mags-k31.tar.gz"
    output:"gather_databases/almeida-mags-k31.sbt.json"
    params:
        outdir = "gather_databases"
    shell:
        """
        tar xf {input} -C {params.outdir}
        """

rule download_gather_pasolli:
    output: "gather_databases/pasolli-mags-k31.tar.gz"
    params: dl_link="https://osf.io/3vebw/download"
    shell:"""
    wget -O {output} {params.dl_link}
    """
rule untar_pasolli:
    output: "gather_databases/pasolli-mags-k31.sbt.json"
    input: "gather_databases/pasolli-mags-k31.tar.gz"
    params: outdir="gather_databases"
    shell:"""
    tar xf {input} -C {params.outdir}
    """


rule download_gather_nayfach:
    output: "gather_databases/nayfach-k31.tar.gz"
    params: dl_link="https://osf.io/y3vwb/download"
    shell:
        """
        wget -O {output} {params.dl_link}
        """
rule untar_nayfach:
    output:"gather_databases/nayfach-k31.sbt.json"
    input: "gather_databases/nayfach-k31.tar.gz"
    params: outdir="gather_databases"
    shell:"""
    tar xf {input} -C {params.outdir}
    """

rule download_gather_refseq:
    output: "gather_databases/refseq-d2-k31.tar.gz"
    params: dl_link="https://s3-us-west-2.amazonaws.com/sourmash-databases/2018-03-29/refseq-d2-k31.tar.gz"
    shell:"""
    wget -O {output} {params.dl_link}
    """
rule untar_refseq:
    output: "gather_databases/refseq-d2-k31.sbt.json"
    input:  "gather_databases/refseq-d2-k31.tar.gz"
    params: outdir = "gather_databases"
    shell: """
    tar xf {input} -C {params.outdir}
    """



rule sourmash_gather:
    input:
        query="sourmash/sig/{sample}.sig",
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




