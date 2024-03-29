from urllib.parse import urlparse
import os
from os.path import join, basename

configfile: "pipeline_config.yml"

STEPS = {
    0: "step00.raw_reads",
    1: "step01.fetch_genome",
    2: "step02.remove_duplicates",
    3: "step03.trimmomatic",
    4: "step04.cutadapt",
    5: "step05.STAR_rRNA",
    6: "step06.STAR_genome",
    7: "step07.collect_tags",
    8: "step08.genes_intersect",
    9: "step09.pairwise_interactions"
}

ADAPTERS = join(STEPS[3], 'adapters.fa')
RRNA_FASTA = join(STEPS[5], "rRNA.fa")
GENOME = {}
for filetype, path in config['genome'].items():
    GENOME[filetype] = join(STEPS[6], basename(path)).removesuffix('.gz')
VOTING_BED = join(STEPS[8], basename(GENOME['gff3']).removesuffix('.gff3') + '.genes.voting.bed')

rule all:
    input:
        raw_multiqc=join(STEPS[0], "multiqc", "multiqc_report.html"),
        contacts=join(STEPS[9], 'merge.intermolecular.pairwise.tsv')


rule fetch_reads:
    output:
        join(STEPS[0], "{sample}_{mate}.fastq.gz")
    log:
        log=join(STEPS[0], "{sample}_{mate}.log")
    resources:
        cpus=1,
        mem_gb=2
    run:
        path = config['samples'][wildcards.sample][wildcards.mate]
        parsed_path = urlparse(path)
        with open(log.log, 'wt') as log_file:
            if parsed_path.scheme and parsed_path.netloc:
                print(f'Downloading file from {path}...', file=log_file)
                shell(f'wget {path} --no-verbose -a {log.log} --output-document={output}')
            elif os.path.exists(path):
                print(f'Processing local file: {path}', file=log_file)
                shell(f'ln -s {path} {output}')
            else:
                print(f'File not found: {path}', file=log_file)


def fastqc_input(wildcards):
    input_dir = wildcards.directory
    if input_dir == STEPS[4]:
        postfix = "fq"
    elif input_dir == STEPS[0]:
        postfix = "fastq.gz"
    else:
        raise ValueError(f'Incorrect directory for fastqc: {input_dir}')
    return join(f"{input_dir}", f"{wildcards.prefix}.{postfix}")


rule fastqc:
    output:
        join("{directory}", "fastqc", "{prefix}_fastqc.html")
    input:
        fastq=fastqc_input
    params:
        outdir=join("{directory}", "fastqc")
    conda:
        "preprocessing.yml"
    log:
        join("{directory}", "fastqc", "{prefix}_fastqc.log")
    threads: config["fastqc_threads"]
    resources:
        cpus=config["fastqc_threads"],
        mem_gb=4
    shell:
        """
        fastqc -o {params.outdir} --extract -t {threads} \
            -f fastq {input.fastq} &> {log}
        """


rule raw_multiqc:
    input:
        expand(
            join(STEPS[0], "fastqc", "{sample}_{mate}_fastqc.html"),
            sample=config["samples"],
            mate=config["mates"],
        )
    output:
        join(STEPS[0], "multiqc", "multiqc_report.html"),
        directory(join(STEPS[0], "multiqc", "multiqc_data"))
    log:
        join(STEPS[0], "multiqc", "multiqc.log")
    resources:
        cpus=1,
        mem_gb=1
    params:
        outdir=join(STEPS[0], "multiqc"),
        target=join(STEPS[0], "fastqc")
    conda:
        "preprocessing.yml"
    shell:
        """
        multiqc -s --outdir {params.outdir} {params.target} &> {log}
        """


rule deduplicate:
    input:
        mate1=join(STEPS[0], "{sample}_1.fastq.gz"),
        mate2=join(STEPS[0], "{sample}_2.fastq.gz")
    output:
        mate1=join(STEPS[2], "{sample}_1.dedup.fastq.gz"),
        mate2=join(STEPS[2], "{sample}_2.dedup.fastq.gz")
    conda:
        "preprocessing.yml"
    log:
        join(STEPS[2], "{sample}.dedup.log")
    resources:
        cpus=1,
        mem_gb=2
    shell:
        """
        clumpify.sh\
            in={input.mate1} in2={input.mate2}\
            out={output.mate1} out2={output.mate2}\
            dedupe &> {log}
        """


def get_file(path: str, log_path: str, output_path: str):
    parsed_path = urlparse(path)
    with open(log_path, 'wt') as log_file:
        if parsed_path.scheme and parsed_path.netloc:
            print(f'Downloading file from {path}...', file=log_file)
            shell(f'wget {path} --no-verbose -a {log_path} -O {output_path}')
        elif os.path.exists(path):
            print(f'Processing local file: {path}', file=log_file)
            shell(f'ln -s {path} {output_path}')
        else:
            print(f'File not found: {path}', file=log_file)


def get_file2(path: str, log_path: str, output_path: str):
    parsed_path = urlparse(path)
    with open(log_path, 'wt') as log_file:
        if parsed_path.scheme and parsed_path.netloc:
            print(f'Downloading file from {path}...', file=log_file)
            shell(f'wget {path} --no-verbose -a {log_path} -O {output_path}')
        elif os.path.exists(path):
            print(f'Processing local file: {path}', file=log_file)
            shell(f'ln -s {path} {output_path}')
        else:
            print(f'File not found: {path}', file=log_file)
    


rule fetch_adapters:
    output:
        ADAPTERS
    params:
        path=config['adapters']
    log:
        path=join(STEPS[3], "fetch_adapters.log")
    resources:
        cpus=1,
        mem_gb=1
    run:
        get_file(params.path, log.path, output)


rule trimmomatic:
    input:
        mate1=join(STEPS[2], "{sample}_1.dedup.fastq.gz"),
        mate2=join(STEPS[2], "{sample}_2.dedup.fastq.gz"),
        adapters=ADAPTERS
    output:
        mate1_p=join(STEPS[3], "{sample}_1.dedup.clean.pair.fq"),
        mate1_up=join(STEPS[3], "{sample}_1.dedup.clean.unpair.fq"),
        mate2_p=join(STEPS[3], "{sample}_2.dedup.clean.pair.fq"),
        mate2_up=join(STEPS[3], "{sample}_2.dedup.clean.unpair.fq")
    conda:
        "preprocessing.yml"
    threads: config['trimmomatic']['threads']
    params:
        encoding=config['trimmomatic']['encoding'],
        clip=config['trimmomatic']['clip'].format(ADAPTERS)
    resources:
        cpus=config['trimmomatic']['threads'],
        mem_gb=4
    log:
        join(STEPS[3], "{sample}.dedup.clean.log")
    shell:
        """
        trimmomatic PE {params.encoding} -threads {threads}\
            {input.mate1} {input.mate2}\
            {output.mate1_p} {output.mate1_up}\
            {output.mate2_p} {output.mate2_up}\
            {params.clip} &> {log}
        """


rule cutadapt:
    input:
        mate1_p=join(STEPS[3], "{sample}_1.dedup.clean.pair.fq"),
        mate1_up=join(STEPS[3], "{sample}_1.dedup.clean.unpair.fq"),
        mate2_p=join(STEPS[3], "{sample}_2.dedup.clean.pair.fq"),
        mate2_up=join(STEPS[3], "{sample}_2.dedup.clean.unpair.fq")
    output:
        mate1_p=join(STEPS[4], "{sample}_1.dedup.clean.pair.cut.fq"),
        mate1_up=join(STEPS[4], "{sample}_1.dedup.clean.unpair.cut.fq"),
        mate2_p=join(STEPS[4], "{sample}_2.dedup.clean.pair.cut.fq"),
        mate2_up=join(STEPS[4], "{sample}_2.dedup.clean.unpair.cut.fq")
    log:
        pair=join(STEPS[4], "{sample}.cut.log"),
        unpair1=join(STEPS[4],"{sample}_1.unpair.cut.log"),
        unpair2=join(STEPS[4],"{sample}_2.unpair.cut.log")
    conda:
        "preprocessing.yml"
    threads: config["cutadapt"]['threads']
    resources:
        cpus=config["cutadapt"]['threads'],
        mem_gb=4
    params:
        cut1=lambda wc: config["cutadapt"]['cut1'],
        cut2=lambda wc: config["cutadapt"]['cut2'],
        params=config['cutadapt']['params']
    shell:
        """
        cutadapt -j {threads} {params.cut1} {params.cut2} {params.params}\
            -o {output.mate1_p} -p {output.mate2_p}\
            {input.mate1_p} {input.mate2_p} &> {log.pair}

        cutadapt -j {threads} {params.cut1} {params.params}\
            -o {output.mate1_up} {input.mate1_up} &> {log.unpair1}

        cutadapt -j {threads} {params.cut1} {params.params}\
            -o {output.mate2_up} {input.mate2_up} &> {log.unpair2}
        """


rule fetch_rRNA:
    output:
        RRNA_FASTA
    params:
        rrna_ac=config["rRNA_AC"]
    resources:
        cpus=1,
        mem_gb=1
    conda:
        "preprocessing.yml"
    shell:
        """
        efetch -db nuccore -id {params.rrna_ac} -format fasta > {output}
        """


rule STAR_index:
    input:
        fa=lambda wc: RRNA_FASTA if wc.type == 'rRNA' else GENOME['fasta'],
        gff3=GENOME['gff3']
    output:
        directory(join("{directory}", "{type}_index"))
    wildcard_constraints:
        directory=f'({STEPS[5]}|{STEPS[6]})',
        type='(rRNA|genome)'
    threads: config["STAR"]["indexing_threads"]
    resources:
        cpus=config["STAR"]["indexing_threads"],
        mem_gb=60
    conda:
        "preprocessing.yml"
    log:
        join('{directory}', "{type}_index.log")
    params:
        size=lambda wc: config['STAR'][f'{wc.type}_genomeSAindexNbases'],
        sjdb=lambda wc, input: '' if wc.type == 'rRNA' else f'--sjdbGTFfile {input.gff3} {config["STAR"]["sjdb"]}'
    shell:
        """
        STAR --runMode genomeGenerate\
            --runThreadN {threads}\
            --genomeSAindexNbases {params.size}\
            --genomeDir {output}\
            --genomeFastaFiles {input.fa} {params.sjdb} &> {log}
        """


def STAR_input(wc):
    if wc.target == 'rRNA':
        result = dict(
            index=join(STEPS[5], 'rRNA_index'),
            mate1_p=join(STEPS[4], f"{wc.sample}_1.dedup.clean.pair.cut.fq"),
            mate1_up=join(STEPS[4], f"{wc.sample}_1.dedup.clean.unpair.cut.fq"),
            mate2_p=join(STEPS[4], f"{wc.sample}_2.dedup.clean.pair.cut.fq"),
            mate2_up=join(STEPS[4], f"{wc.sample}_2.dedup.clean.unpair.cut.fq")
        )
    else:
        assert wc.target == 'genome'
        result = dict(
            index=join(STEPS[6], "genome_index"),
            mate1_p=join(STEPS[5], f"{wc.sample}_pair_rRNA_Unmapped.out.mate1"),
            mate1_up=join(STEPS[5], f"{wc.sample}_1_unpair_rRNA_Unmapped.out.mate1"),
            mate2_p=join(STEPS[5], f"{wc.sample}_pair_rRNA_Unmapped.out.mate2"),
            mate2_up=join(STEPS[5], f"{wc.sample}_2_unpair_rRNA_Unmapped.out.mate1")
        )
    return result


rule STAR_alignment:
    output:
        join('{directory}','{sample}_pair_{target}_Aligned.out.bam'),
        expand(
            join('{directory}','{sample}_{mate}_unpair_{target}_{suffix}'),
            mate=config['mates'],
            suffix=('Unmapped.out.mate1', 'Aligned.out.bam'),
            allow_missing=True
        ),
        expand(
            join('{directory}', '{sample}_pair_{target}_Unmapped.out.mate{mate}'),
            mate=config['mates'],
            allow_missing=True
        )
    wildcard_constraints:
        directory=f'({STEPS[5]}|{STEPS[6]})',
        target=r'(rRNA|genome)'
    input:
        unpack(STAR_input)
    params:
        prefix_p=join('{directory}', "{sample}_pair_{target}_"),
        prefix1_up=join('{directory}', "{sample}_1_unpair_{target}_"),
        prefix2_up=join('{directory}', "{sample}_2_unpair_{target}_"),
        args=config["STAR"]['alignment_args'],
        overlap=config["STAR"]['peOverlapNbasesMin']
    log: join('{directory}', "{sample}_{target}.log"),
    threads: config["STAR"]["alignment_threads"]
    resources:
        cpus=config["STAR"]["alignment_threads"],
        mem_gb=60
    conda:
        "preprocessing.yml"
    shell:
        """
        STAR --runMode alignReads\
            --genomeDir {input.index}\
            --readFilesIn {input.mate1_p} {input.mate2_p}\
            --outFileNamePrefix {params.prefix_p}\
            --peOverlapNbasesMin {params.overlap}\
            --runThreadN {threads} {params.args} &> {log}

        STAR --runMode alignReads\
            --genomeDir {input.index}\
            --readFilesIn {input.mate1_up}\
            --outFileNamePrefix {params.prefix1_up}\
            --outSAMflagOR 64\
            --runThreadN {threads} {params.args} &>> {log}

        STAR --runMode alignReads\
            --genomeDir {input.index}\
            --readFilesIn {input.mate2_up}\
            --outFileNamePrefix {params.prefix2_up}\
            --outSAMflagOR 128\
            --runThreadN {threads} {params.args} &>> {log}
        """


def get_gz_file(path: str, log_path: str, outdir: str):
    output_path = join(outdir, basename(path))
    parsed_path = urlparse(path)
    with open(log_path, 'wt') as log_file:
        if parsed_path.scheme and parsed_path.netloc:
            print(f'Downloading file from {path}...', file=log_file)
            shell(f'wget {path} --no-verbose -a {log_path} -O {output_path}')
        elif os.path.exists(path):
            print(f'Processing local file: {path}',file=log_file)
            shell(f'ln -s {path} {output_path}')
        else:
            print(f'File not found: {path}',file=log_file)
        if output_path.endswith('.gz'):
            print(f'Unpacking file: {output_path}', file=log_file)
            shell(f'gunzip -k {output_path}')


for filetype in 'fasta', 'gff3':
    rule:
        name: f'fetch_genome_{filetype}'
        output: GENOME[filetype]
        params:
            path = config['genome'][filetype]
        resources:
            cpus=1,
            mem_gb=2
        log:
            path=join(STEPS[6], f"fetch_{filetype}.log")
        run:
            get_gz_file(params.path, log.path, STEPS[6])



rule collect_tags:
    output:
        ambiguous = join(STEPS[7], '{sample}_tags.ambiguous.bed'),
        bed = join(STEPS[7], '{sample}_tags.bed')
    input:
        script = 'scripts/collect_tags.py',
        p=join(STEPS[6],"{sample}_pair_genome_Aligned.out.bam"),
        up1=join(STEPS[6],"{sample}_1_unpair_genome_Aligned.out.bam"),
        up2=join(STEPS[6],"{sample}_2_unpair_genome_Aligned.out.bam")
    conda:
        "postprocessing.yml"
    params:
        script_threads = config['collect_tags_threads'],
        quality = config['alignment_quality']
    threads:
        config['collect_tags_threads'] - 1
    resources:
        cpus=config['collect_tags_threads'],
        mem_gb=10
    log:
        join(STEPS[7],'{sample}_tags.log')
    shell:
        """
        samtools cat -@ {threads} {input.p} {input.up1} {input.up2} 2>> {log}\
            | samtools view -@ {threads} --keep-tag jM {params.quality} - 2>> {log}\
            | python3 {input.script} --threads {params.script_threads}\
                -i stdin -o {output.bed} -a {output.ambiguous} &>> {log}
        """

rule voting:
    input:
        gff3 = GENOME['gff3'],
        tags = expand(
            join(STEPS[7], '{sample}_tags.bed'),
            sample=config['samples']
        )
    output: VOTING_BED
    conda:
        "postprocessing.yml"
    resources:
        cpus=1,
        mem_gb=60
    threads: 1
    shell:
        """
        awk -F'\t' '$3 == "gene"' {input.gff3}\
            | gff2bed -d | cut -f1-6\
            | bedtools intersect -c -S -F 0.9 -nonamecheck -a stdin -b {input.tags}\
            | awk '{{ print $0"\t"($NF / ($3 - $2)) }}' > {output}
        """

rule genes_intersect:
    output:
        intra=join(STEPS[8], '{sample}_tags.intramolecular.bed'),
        inter=join(STEPS[8], '{sample}_tags.intermolecular.bed'),
        without_gene=join(STEPS[8], '{sample}_tags.without_gene.bed')
    input:
        annotation=VOTING_BED,
        tags=join(STEPS[7], '{sample}_tags.bed'),
        script=join('scripts', 'separate_inter.py')
    params:
        prefix=join(STEPS[8], '{sample}_tags')
    threads: 1
    resources:
        cpus=1,
        mem_gb=10
    conda:
        'postprocessing.yml'
    log:
        join(STEPS[8],'{sample}_tags.genes_intersect.log')
    shell:
        """
        bedtools intersect -f 0.90 -S -wao -nonamecheck\
            -a {input.tags}\
            -b {input.annotation} 2>> {log}\
            | python3 {input.script}\
            -i stdin -p {params.prefix} -n 3 -g 11 &>> {log}
        """

rule choose_alignment:
    output:
        pairwise=join(STEPS[9], '{sample}_tags.intermolecular.pairwise.bed'),
        choice=join(STEPS[9], '{sample}_tags.intermolecular.choice.bed')
    input:
        tags=join(STEPS[8], '{sample}_tags.intermolecular.bed'),
        script=join('scripts', 'choose_annotation.py')
    threads: 1
    resources:
        cpus=1,
        mem_gb=60
    conda:
        'postprocessing.yml'
    log:
        join(STEPS[9],'{sample}_tags.choice.log')
    shell:
        """
        python3 {input.script} -i {input.tags} -o {output.choice} -p {output.pairwise} &> {log}
        """

rule significant_interactions:
    input:
        contacts=expand(
            join(STEPS[9], '{sample}_tags.intermolecular.pairwise.bed'),
            sample=config['samples']
        ),
        script=join('scripts', 'significant.py'),
        annotation=GENOME['gff3']
    output:
        join(STEPS[9], 'merge.intermolecular.pairwise.tsv')
    threads: 1
    resources:
        cpus=1,
        mem_gb=60
    conda:
        'postprocessing.yml'
    log:
        join(STEPS[9],'merge.significant.log')
    shell:
        """
        tail -n +2 {input.contacts}\
            | python3 {input.script} -i stdin -o {output} -a {input.annotation} &> {log}
        """
