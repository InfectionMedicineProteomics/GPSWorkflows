#!/usr/bin/env python3
SAMPLES, = glob_wildcards(f"{config['base_file_paths']['osw']}/{{sample}}.osw")

rule all:
    input:
        combined_output = f"{config['base_file_paths']['results']}/models/combined.pin.tsv",
        model_weights = f"{config['base_file_paths']['results']}/models/{config['percolator']['model_name']}"
        
    
rule gscore_export_training_data:
    input:
        osw_file = f"{config['base_file_paths']['osw']}/{{sample}}.osw",
    output:
        pin = f"{config['base_file_paths']['results']}/training_data/{{sample}}.pin.tsv"
    params:
        percolator_exe = config['percolator']['exe']
    threads:
        config["percolator"]["threads"]
    shell:
        (
            "gscore score "
            "--percolator-output "
            "--input {input.osw_file} "
            "--output {output.pin} "
            "--threads {threads} "
            "--percolator-exe {params.percolator_exe}"
        )

rule gscore_train_model:
    input:
        pins = expand(f"{config['base_file_paths']['results']}/training_data/{{sample}}.pin.tsv", sample=SAMPLES)
    output:
        combined_output = f"{config['base_file_paths']['results']}/models/combined.pin.tsv",
        model_weights = f"{config['base_file_paths']['results']}/models/{config['percolator']['model_name']}"
    threads:
        config["threads"]
    params:
        percolator_exe = config['percolator']['exe']
    shell:
        (
            "gscore train "
            "--train-percolator-model "
            "--percolator-exe {params.percolator_exe} "
            "--input {input.pins} "
            "--percolator-output {output.combined_output} "
            "--model-output {output.model_weights} "
            "--threads {threads}"
        )