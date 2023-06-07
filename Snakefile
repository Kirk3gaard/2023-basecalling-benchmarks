configfile: "config.yaml"

rule all:
    input:
        expand("temp/{npID}_NanoPlot-data.tsv.gz", npID=config["inputNPfiles"]),
        expand("data/{npID}.fastq.gz", npID=config["inputNPfiles"]),
        expand("data/{npID}-{nreads}.fastq.gz", npID=config["inputNPfiles"], nreads=config["npreads"]),
        expand("results/{npID}-{nreads}.{asmtype}.fa.gz", npID=config["inputNPfiles"], nreads=config["npreads"], asmtype=config["asmtypes"]),
        expand("temp/{npID}-{nreads}.contigIDs.{asmtype}.txt", npID=config["inputNPfiles"], nreads=config["npreads"], asmtype=config["asmtypes"]),
        "results/fastani.tsv",
        "README.md",

rule get_data:
    input:
        "/user_data/rhk/zymo_hmw_mock/{npID}.fastq.gz"
    output:
        "data/{npID}.fastq.gz"
    threads: 1
    resources:
        mem_mb=10000,
        node_type="basic",
    shell:
        """
        ln -s {input} {output}
        """
        
rule getNPsubsets:
    input:
        NPreads="data/{npID}.fastq.gz",
    output:
        "data/{npID}-{nreads}.fastq.gz"
    threads: 2
    resources:
        mem_mb=100000,
        node_type="basic",
    shell:
        """
        seqtk sample -s100 {input} {wildcards.nreads} | gzip > {output}
        """

rule get_new_refs:
    output:
        expand("data/hmw_individual/{refID}.fasta", refID=config["reffiles"]),
    threads: 1
    resources:
        mem_mb=10000,
        node_type="basic",
    shell:
        """
        mkdir -p data/hmw_individual/
        for f in data/hifirefs/*.fasta ; do
        ## There is an additional \ to avoid python/snakemake from causing problems https://stackoverflow.com/questions/70324411/snakemake-truncating-shell-codes
            NAME=$(echo "$f" | sed -E 's/.*\\/(.*).fasta/\\1/');
            seqname=$(cat $f | bioawk -c fastx '{{ print length($seq), $name }}' | sort -k1,1rn | head -1 | cut -f2)
            samtools faidx $f $seqname > data/hmw_individual/$NAME.fasta
        done
        """

rule NPQC_ref:
    input:
        NPreads="data/{npID}.fastq.gz",
        ref="data/ref.fasta",
    output:
        "temp/{npID}_NanoPlot-data.tsv.gz"
    threads: 30
    resources:
        mem_mb=100000,
        node_type="himem",
    shell:
        """
        # Run mapping
        minimap2 -ax map-ont --secondary=no -t {threads} {input.ref} {input.NPreads} |\
          samtools view --threads {threads} -Sb -F 0x104 - |\
          samtools sort -T temp/ --threads {threads} - > temp/{wildcards.npID}.cov.bam
        samtools index temp/{wildcards.npID}.cov.bam
        bamtools split -in temp/{wildcards.npID}.cov.bam -reference
        for f in temp/{wildcards.npID}*REF*bam;
            do 
                NAME=$(echo $f | sed 's/.*REF_//' | sed 's/.bam//')
                samtools index $f
                NanoPlot --threads {threads} -o temp/NP_QC_ref_{wildcards.npID}/ --no_static --raw --tsv_stats --bam  $f
                cp temp/NP_QC_ref_{wildcards.npID}/NanoPlot-data.tsv.gz temp/{wildcards.npID}_ref_$NAME.NanoPlot-data.tsv.gz
                rm -rf temp/NP_QC_ref_{wildcards.npID}/
                rm $f
            done
        cat temp/{wildcards.npID}_ref_*.NanoPlot-data.tsv.gz > {output}
        """

rule flye:
    input:
        NPreads="data/{npID}.fastq.gz"
    output:
        asm="results/{npID}.flye.fa.gz",
        asminfo="results/{npID}.assembly_info.txt",
    threads: config["assembly_threads"]
    resources:
        mem_mb=config["assembly_mb"],
        node_type="basic",
    shell:
        """
        basecall_mode=$(echo {wildcards.npID} | sed -E 's/.*_([a-z]+)\@.*/\1/')
        if [[ $basecall_mode == "fast" ]]
        then
        flye --nano-raw {input.NPreads} --threads {threads} --meta --out-dir temp/flye.{wildcards.npID}
        else
        flye --nano-hq {input.NPreads} --threads {threads} --meta --out-dir temp/flye.{wildcards.npID}
        fi
        cat temp/flye.{wildcards.npID}/assembly.fasta | gzip > {output.asm}
        cp temp/flye.{wildcards.npID}/assembly_info.txt {output.asminfo}
        """

rule medaka1x:
    input:
        asm="results/{npID}.flye.fa.gz",
        NPreads="data/{npID}.fastq.gz"
    output:
        "results/{npID}.flye.medaka1x.fa.gz"
    threads: config["polish_threads"]
    resources:
        mem_mb=config["polish_mb"],
        node_type="himem",
    params:
        medaka_model=config["medaka_model"],
    shell:
        """
        ## There is an additional \ to avoid python/snakemake from causing problems https://stackoverflow.com/questions/70324411/snakemake-truncating-shell-codes
        basecall_mode=$(echo {wildcards.npID} | sed -E 's/.*_([a-z]+)\@.*/\\1/')
        basecall_version=$(echo {wildcards.npID} | sed -E 's/.*_[a-z]+\@v(.*)-.*/\\1/')
        if [ "$basecall_mode" = "fast" ]; then
            medaka_model="r1041_e82_400bps_fast_g632" # There does not appear to be a 4.0.0, 4.1.0 or 4.2.0 medaka model for fast mode
        else 
            medaka_model="r1041_e82_400bps_"$basecall_mode"_v$basecall_version"
        fi
        zcat {input.asm} > temp/{wildcards.npID}-unpolished.fa
        medaka_consensus -i {input.NPreads} -d temp/{wildcards.npID}-unpolished.fa -o temp/{wildcards.npID}-medaka -t {threads} -m $medaka_model
        cat temp/{wildcards.npID}-medaka/consensus.fasta | gzip > {output}
        """

rule extract_contigs:
    input:
        asm="results/{npID}.{asmtype}.fa.gz",
        asminfo="results/{npID}.assembly_info.txt"
    output:
        "temp/{npID}.contigIDs.{asmtype}.txt"
    threads: 2
    resources:
        mem_mb=5000,
        node_type="basic",
    shell:
        """
        module load seqtk
        mkdir -p results/bins/
        cat {input.asminfo} | cut -f 1 > {output}
        # Separate contigs as individual "bins" (https://twitter.com/lh3lh3/status/1453374559084765185)
        while read -r line;
        do
        echo $line > temp/{wildcards.npID}.{wildcards.asmtype}thisID
        seqtk subseq {input.asm} temp/{wildcards.npID}.{wildcards.asmtype}thisID | gzip > results/bins/{wildcards.npID}.{wildcards.asmtype}.$line.fa.gz
        done < {output}
        module unload seqtk
        """

rule split_quast_QC:
    input:
        IDfiles=expand("temp/{npID}-{nreads}.contigIDs.{asmtype}.txt", npID=config["inputNPfiles"], nreads=config["npreads"], asmtype=config["asmtypes"]),
        genomes=expand("results/{npID}-{nreads}.{asmtype}.fa.gz", npID=config["inputNPfiles"], nreads=config["npreads"], asmtype=config["asmtypes"]),
        ref="data/hmw_individual/{refID}.fasta",
    output:
        "results/quast_{refID}.tsv"
    threads: 16
    resources:
        mem_mb=10000,
        node_type="basic",
    shell:
        """
        quast.py --threads {threads} --output-dir temp/quast_ref_{wildcards.refID} -r {input.ref} --no-plots --no-html results/bins/*.fa.gz {input.ref}
        cat temp/quast_ref_{wildcards.refID}/transposed_report.tsv > {output}
        """

rule split_ANI:
    input:
        expand("temp/{npID}-{nreads}.contigIDs.{asmtype}.txt", npID=config["inputNPfiles"],nreads=config["npreads"], asmtype=config["asmtypes"]),
        expand("results/{npID}-{nreads}.{asmtype}.fa.gz", npID=config["inputNPfiles"], nreads=config["npreads"], asmtype=config["asmtypes"]),
        expand("data/hmw_individual/{refID}.fasta", refID=config["reffiles"]),
    output:
        "results/fastani.tsv"
    threads: 16
    resources:
        mem_mb=10000,
        node_type="basic",
    shell:
        """
        ls results/bins/*.fa.gz > temp/draft_genome_list.txt
        ls data/hmw_individual/*.fasta > temp/ref_genome_list.txt
        fastANI --threads {threads} --queryList temp/draft_genome_list.txt --refList temp/ref_genome_list.txt -o {output}
        """


rule knitRMD:
    input:
        "README.Rmd",
        "results/fastani.tsv",
        expand("results/quast_{refID}.tsv", refID=config["reffiles"]),
        expand("results/{npID}-{nreads}.{asmtype}.fa.gz", npID=config["inputNPfiles"], nreads=config["npreads"], asmtype=config["asmtypes"]),
    output:
        "README.md"
    threads: 4
    resources:
        mem_mb=10000,
        node_type="basic",
    shell:
        """
        R -e "rmarkdown::render('README.Rmd')"
        """