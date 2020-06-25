input_data = "vir_SAMPLEIDs.txt"
samples= [x.strip().split(".tar")[0] for x in open(input_data, 'r')]
print(samples)

rule all:
    input: expand("raw_data/zipped_data/{sample}_{ext}.bz2", sample =samples, ext = [1,2])


rule uncompress_genome:
    input: "raw_data/zipped_data/{sample}.tar"
    output: 
        "raw_data/zipped_data/{sample}_1.bz2",
        "raw_data/zipped_data/{sample}_2.bz2" 
    shell:
        "tar xvf {input}"
