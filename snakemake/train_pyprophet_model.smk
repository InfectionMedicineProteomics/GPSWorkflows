#!/usr/bin/env python3
SAMPLES, = glob_wildcards(f"{config['base_file_paths']['osw']}/{{sample}}.osw")

rule all:
    input:
        model = f"{config['base_file_paths']['results']}/pyprophet/models/scoring_model.osw"
        

rule pyprophet_subsample:
    input:
        osw = f"{config['base_file_paths']['osw']}/{{sample}}.osw",
    output:
        subsampled = f"{config['base_file_paths']['results']}/pyprophet/osw/{{sample}}.osws"
    container:
        config['containers']['pyprophet']
    threads: 5
    shell:
        (
            "pyprophet subsample "
            "--in {input.osw} "
            "--out {output.subsampled}"
        )


rule pyprophet_merge:
    input:
        template = config['libraries']['spectral_library'],
        osws = expand(f"{config['base_file_paths']['results']}/pyprophet/osw/{{sample}}.osws", sample=SAMPLES)
    output:
        merged = f"{config['base_file_paths']['results']}/pyprophet/osw/merged.osw"
    wildcard_constraints:
        loop="[0-9]+"
    container:
        config['containers']['pyprophet']
    threads: 1
    shell:
        (
            "pyprophet merge "
            "--template {input.template} "
            "--out {output.merged} "
            "-- {input.osws}"
        )


rule pyprophet_learn:
    input:
        merged = f"{config['base_file_paths']['results']}/pyprophet/osw/merged.osw"
    output:
        model = f"{config['base_file_paths']['results']}/pyprophet/models/scoring_model.osw"
    container:
        config['containers']['pyprophet']
    threads: config['pyprophet']['final_threads']
    params:
        level = 'ms2',
        classifier = "LDA"
    shell:
        (
            "pyprophet score "
            "--in {input} "
            "--out {output.model} "
            "--level {params.level} "
            "--classifier {params.classifier} "
            "--threads {threads}"
        )