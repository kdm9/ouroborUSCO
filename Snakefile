configfile: "config.yml"
wildcard_constraints:
    run="[^/]+",
    lib="[^/]+",
    aligner="[^/]+",
    sample="[^/]+",
    ref="[^/]+",
    type="[^/]+",


rule all:
    input:
        expand("data/{iter}/{ref}/corrected_sequences/{sample}.fasta", iter=config["iterations"], ref=config["origrefs"], sample=config["samples"])

def lastiter(wc, path):
    return expand(path, iter=(int(wc.iter) - 1))

def getref(wc):
    ref = None
    if wc.iter == "1":
        ref = config['origrefs'][wc.ref]
    else:
        ref = lastiter(wc, "data/{iter}/" + wc.ref + "/corrected_reference/" + wc.sample + ".fasta")[0]
    return ref

def busco2regions(resulttab):
    from collections import defaultdict
    buscolocs = defaultdict(list)
    with open(resulttab) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("#"):
                continue
            fields = line.split('\t')
            if fields[1] != "Complete":
                continue
            busco, status, region = fields[:3]
            buscolocs[busco].append(region)
    goodregions = []
    for busco, regions in buscolocs.items():
        if len(regions) == 1:
            goodregions.extend(regions)
    return goodregions


rule fasta_index:
    input:
        "{path}.fasta"
    output:
        "{path}.fasta.fai"
    shell:
        "samtools faidx {input}"

rule bam_index:
    input:
        "{path}.bam"
    output:
        "{path}.bam.bai"
    shell:
        "samtools index {input}"


checkpoint busco_firstrun:
    input:
        ref=lambda wc: config['origrefs'][wc.ref]
    output:
        fulltable="data/1/{ref}/firstbusco/{ref}/run_eudicots_odb10/full_table.tsv",
    log:
        "data/1/{ref}/log/initialbusco.log"
    threads: 16
    params:
        busco_lineage=config["busco"]["lineage"]
    shell:
        "(cd data/1/{wildcards.ref}/firstbusco &&"
        " busco"
        "   --in {input.ref}"
        "   --out {wildcards.ref}"
        "   --lineage {params.busco_lineage}"
        "   --offline"
        "   --mode genome"
        "   --cpu {threads}"
        "   --force"
        ") >{log} 2>&1"


rule ngmidx:
    input:
        "{path}.fasta"
    output:
        "{path}.fasta-enc.2.ngm"

rule ngmap:
    input:
        reads=lambda wc: "rawdata/reads/{sample}.fastq.gz" if wc.iter == "1" else lastiter(wc, "data/{iter}/{{ref}}/readsout/ngm/{{sample}}.fastq.gz"),
        ref=getref,
        idx=lambda wc: getref(wc) + "-enc.2.ngm"
    output:
        bam=temp("data/{iter}/{ref}/rawbam/ngm/{sample}.bam"),
    log:
        "data/{iter}/{ref}/log/ngm/{sample}.log"
    threads:
        8
    params:
        sensitivity=config["mapping"]["ngm"]["sensitivity"],
    shell:
        "( ngm"
        "   -q {input.reads}"
        "   --paired --broken-pairs"
        "   -r {input.ref}"
        "   -t {threads}"
        "   --rg-id {wildcards.iter}_{wildcards.sample}"
        "   --rg-sm {wildcards.sample}"
        "   --sensitivity {params.sensitivity}" # this is the mean from a bunch of different runs
        "| samtools view -Suh - >{output.bam}"
        " ) >{log} 2>&1"


rule bam_markdups_sort:
    input:
        bam="data/{iter}/{ref}/rawbam/{aligner}/{sample}.bam",
        ref=getref,
    output:
        bam=temp("data/{iter}/{ref}/bam/{aligner}/{sample}.bam"),
    threads: 4
    log: "data/{iter}/{ref}/log/markdup/{aligner}/{sample}.log"
    shell:
        "( samtools fixmate "
        "   -m"
        "   -@ {threads}"
        "   --output-fmt bam,level=0"
        "   {input.bam}"
        "   -"
        " | samtools sort"
        "   -T ${{TMPDIR:-/tmp}}/{wildcards.sample}_sort_$RANDOM"
        "   --output-fmt bam,level=0"
        "   -@ {threads}"
        "   -m 1g"
        "   -"
        " | samtools markdup"
        "   -T ${{TMPDIR:-/tmp}}/{wildcards.sample}_markdup_$RANDOM"
        "   -s" # report stats
        "   -@ {threads}"
        "   --output-fmt bam,level=3"
        "   -"
        "   {output.bam}"
        " ) >{log} 2>&1"


def getbuscoregions(wc):
    if wc.iter == "1":
        buscofile = checkpoints.busco_firstrun.get(ref=wc.ref).output[0]
    else:
        last = str((int(wc.iter) - 1))
        buscofile = checkpoints.busco_persamp.get(ref=wc.ref, sample=wc.sample, iter=last).output[0]
    regions = busco2regions(buscofile)
    return " ".join(regions)

rule getreads:
    input:
        bam="data/{iter}/{ref}/bam/{aligner}/{sample}.bam",
        bai="data/{iter}/{ref}/bam/{aligner}/{sample}.bam.bai",
    output:
        "data/{iter}/{ref}/readsout/{aligner}/{sample}.fastq.gz"
    params:
        regions=getbuscoregions,
    shell:
        "cat"
        "   <(samtools view -u -f 4 {input.bam} | samtools collate -O -u | samtools fastq)" 
        "   <(samtools view -u {params.regions} | samtools collate -O -u | samtools fastq)" 
        "| gzip > {output}"

rule getbuscoseqs:
    input:
        ref=getref,
    output:
        "data/{iter}/{ref}/corrected_sequences/{sample}.fasta"
    params:
        regions=getbuscoregions,
    shell:
        "samtools faidx {input.ref} {params.regions} >{output}"


rule mpileup:
    input:
        bam="data/{iter}/{ref}/bam/{aligner}/{sample}.bam",
        bai="data/{iter}/{ref}/bam/{aligner}/{sample}.bam.bai",
        ref=getref,
        fai=lambda wc: getref(wc) + ".fai"
    output:
        bcf="data/{iter}/{ref}/variants/raw_split/mpileup~{aligner}~{ref}~{sample}/{region}.bcf",
    log:
        "data/{iter}/{ref}/log/variants/raw_split/mpileup~{aligner}~{ref}~{sample}/{region}.log",
    params:
        theta=config["varcall"].get("theta_prior", 0.01),
        minmq=lambda wc: config["varcall"]["minmapq"].get(wc.aligner, 5),
        minbq=config["varcall"]["minbq"],
    priority: 1  # get them done earlier, normalisation is super quick
    shell:
        "( bcftools mpileup"
        "   --adjust-MQ 50"
        "   --redo-BAQ"
        "   --max-depth 20000" # the default per file max (250x) is insane, i.e. <1x for most sets. new limit of 20000x  equates to a max. of 20x across all samples.
        "   --min-MQ {params.minmq}"
        "   --min-BQ {params.minbq}"
        "   --fasta-ref {input.ref}"
        "   --annotate FORMAT/DP,FORMAT/AD,FORMAT/SP,INFO/AD" #output extra tags
        "   --region '{wildcards.region}'"
        "   --output-type u" # uncompressed bam
        "   {input.bam}"
        " | bcftools call"
        "   --targets '{wildcards.region}'" # might not be needed
        "   --multiallelic-caller"
        "   --prior {params.theta}"
        "   -O b" # compressed bam
        "   -o {output.bcf}"
        " ) >{log} 2>&1"


rule bcfnorm:
    input:
        bcf="data/{iter}/{ref}/variants/raw_split/mpileup~{aligner}~{ref}~{sample}/{region}.bcf",
        ref=getref,
        fai=lambda wc: getref(wc) + ".fai"
    output:
        bcf=temp("data/{iter}/{ref}/variants/norm/mpileup~{aligner}~{ref}~{sample}/{region}.bcf"),
    log:
        "data/{iter}/{ref}/log/variants/norm/mpileup~{aligner}~{ref}~{sample}/{region}.log",
    shell:
        "( bcftools norm"
        "   --fasta-ref {input.ref}"
        "   -O u"
        "   {input.bcf}" # SKIP VT FOR NOW " | vt decompose_blocksub + -o -" # decompose MNP to multipe SNPs
        " | bcftools norm" # Split multi-alleics
        "   --fasta-ref {input.ref}"
        "   --do-not-normalize"
        "   --multiallelics -snps"
        "   -O u  -o {output.bcf}"
        " ) >{log} 2>&1"

rule bcffilter:
    input:
        bcf="data/{iter}/{ref}/variants/norm/{caller}~{aligner}~{ref}~{sample}/{region}.bcf",
        ref=getref,
        fai=lambda wc: getref(wc) + ".fai"
    output:
        bcf=temp("data/{iter}/{ref}/variants/filtered/{caller}~{aligner}~{ref}~{sample}_filtered~{filter}/{region}.bcf"),
    log:
        "data/{iter}/{ref}/log/variants/filtered/{caller}~{aligner}~{ref}~{sample}_filtered~{filter}/{region}.log",
    params:
        filtarg=lambda wc: config["varcall"]["filters"][wc.filter].replace('\n', ' ')
    shell:
        "( bcftools view"
        "   {params.filtarg}"
        "   -O u"
        "   {input.bcf}"
        " | bcftools norm" # We normalise here to re-join multi-allelic sites, after filtering with multi-allelics split
        "   --fasta-ref {input.ref}"
        "   --do-not-normalize"
        "   --multiallelics +snps" # Split multi-alleic sites
        "   -O b  -o {output.bcf}"
        " ) >{log} 2>&1"


def getbuscobcfs(wc):
    if wc.iter == "1":
        buscofile = checkpoints.busco_firstrun.get(ref=wc.ref).output[0]
    else:
        last = str((int(wc.iter) - 1))
        buscofile = checkpoints.busco_persamp.get(ref=wc.ref, sample=wc.sample, iter=last).output[0]
    regions = busco2regions(buscofile)
    return expand("data/{iter}/{ref}/variants/filtered/{caller}~{aligner}~{ref}~{sample}_filtered~{filter}/{region}.bcf",
                  iter=wc.iter, caller=wc.caller, aligner=wc.aligner, ref=wc.ref, sample=wc.sample, filter=wc.filter, region=regions)

localrules: bcfmerge_fofn
rule bcfmerge_fofn:
    input:
        bcf=getbuscobcfs,
    output:
        fofn=temp("data/{iter}/{ref}/variants/filtered/{caller}~{aligner}~{ref}~{sample}_filtered~{filter}.INPUT_FOFN"),
    run:
        with open(output[0], "w") as fh:
            for s in sorted(input):
                print(s, file=fh)

rule bcfmerge:
    input:
        fofn="data/{iter}/{ref}/variants/filtered/{caller}~{aligner}~{ref}~{sample}_filtered~{filter}.INPUT_FOFN",
    output:
        bcf="data/{iter}/{ref}/variants/final/{caller}~{aligner}~{ref}~{sample}_filtered~{filter}.bcf",
    log:
        "data/{iter}/{ref}/log/variants/final/{caller}~{aligner}~{ref}~{sample}_filtered~{filter}.log",
    threads: 4
    shell:
        "( bcftools concat"
        "   --threads {threads}"
        "   -O b"
        "   -o {output.bcf}"
        "   --file-list {input.fofn}"
        " ) >{log} 2>&1"

rule refupdate:
    input:
        ref=getref,
        fai=lambda wc: getref(wc) + ".fai",
        bcf="data/{iter}/{ref}/variants/final/mpileup~ngm~{ref}~{sample}_filtered~default.bcf",
    output:
        newref="data/{iter}/{ref}/corrected_reference/{sample}.fasta",
    shell:
        "bcftools consensus"
        "   --fasta-ref {input.ref}"
        "   {input.bcf}"
        ">{output.newref}"
        "2>{log}"


checkpoint busco_persamp:
    input:
        ref="data/{iter}/{ref}/corrected_reference/{sample}.fasta",
    output:
        fulltable="data/{iter}/{ref}/busco/{sample}/run_eudicots_odb10/full_table.tsv"
    log:
        "data/{iter}/{ref}/log/busco/{sample}.log"
    threads: 8
    params:
        busco_lineage=config["busco"]["lineage"]
    shell:
        "(cd data/{wildcards.iter}/{wildcards.ref}/busco/ &&"
        " busco"
        "   --in {input.ref}"
        "   --out {wildcards.sample}"
        "   --lineage {params.busco_lineage}"
        "   --mode genome"
        "   --cpu {threads}"
        "   --force"
        ") >{log} 2>&1"

