#!/usr/bin/env python3
SAMPLES, = glob_wildcards(f"{config['base_file_paths']['mzml']}/{{sample}}.mzML")

print(SAMPLES)

rule all:
    input:
        f"{config['base_file_paths']['results']}/gps_75_pit/quantification/{config['gps']['quantification_file_name']}"

rule openswath_init:
    input:
        spectral_library = config['libraries']['spectral_library'],
        rt_library = config['libraries']['rt_library'],
        mzml = f"{config['base_file_paths']['mzml']}/{{sample}}.mzML"
    output:
        osw = f"{config['base_file_paths']['results']}/osw/{{sample}}.osw"
    container:
        config['containers']['openms']
    params:
        rt_normalization_alignment_method = 'lowess',
        scoring_stop_report_after_feature = config['openswath']['stop_after_features'],
        swath_windows_file = config['openswath']['swath_windows_file'],
        min_rsq = config['openswath']['min_rsq'],
        batch_size = config['openswath']['batch_size'],
    threads:
        config["openswath"]['threads']
    shell:
        (
            "OpenSwathWorkflow "
            "-in {input.mzml} "
            "-tr {input.spectral_library} "
            "-out_osw {output.osw} "
            "-out_chrom {output.chrom} "
            "-swath_windows_file {params.swath_windows_file} "
            "-tr_irt {input.rt_library} "
            "-threads {threads} "
            "-enable_ms1 true "
            "-enable_ipf true "
            "-Scoring:Scores:use_uis_scores "
            "-Scoring:stop_report_after_feature {params.scoring_stop_report_after_feature} "
            "-Scoring:TransitionGroupPicker:compute_peak_quality "
            "-RTNormalization:alignmentMethod {params.rt_normalization_alignment_method} "
            "-batchSize {params.batch_size}"
        )

rule gps_score_run_final:
    input:
        osw = f"{config['base_file_paths']['results']}/osw/{{sample}}.osw"
    output:
        score_file = f"{config['base_file_paths']['results']}/gps_75_pit/scored/{{sample}}.scored.tsv"
    params:
        scoring_model = config['models']['scoring_model'],
        scaler = config['models']['scaler'],
    threads: config["gps"]["score"]["threads"]
    shell:
        (
            "gps score "
            "--input {input.osw} "
            "--output {output.score_file} "
            "--scoring-model {params.scoring_model} "
            "--scaler {params.scaler} "
            "--threads {threads} "
            "--vote-percentage 0.75 "
            "--estimate-pit"
        )


rule gps_build_peptide_model:
    input:
        score_files = expand(f"{config['base_file_paths']['results']}/gps_75_pit/scored/{{sample}}.scored.tsv", sample=SAMPLES)
    output:
        peptide_model = f"{config['base_file_paths']['results']}/gps_75_pit/models/peptide.model"
    shell:
        (
            "gps build "
            "--level peptide "
            "--input {input.score_files} "
            "--output {output.peptide_model} "
            "--estimate-pit"
        )

rule gps_build_protein_model:
    input:
        score_files = expand(f"{config['base_file_paths']['results']}/gps_75_pit/scored/{{sample}}.scored.tsv", sample=SAMPLES)
    output:
        protein_model = f"{config['base_file_paths']['results']}/gps_75_pit/models/protein.model"
    shell:
        (
            "gps build "
            "--level protein "
            "-i {input.score_files} "
            "--output {output.protein_model} "
            "--estimate-pit"
        )

rule gps_combine:
    input:
        score_files = expand(f"{config['base_file_paths']['results']}/gps_75_pit/scored/{{sample}}.scored.tsv", sample=SAMPLES),
        peptide_model = f"{config['base_file_paths']['results']}/gps_75_pit/models/peptide.model",
        protein_model = f"{config['base_file_paths']['results']}/gps_75_pit/models/protein.model"
    output:
        quant_matrix = f"{config['base_file_paths']['results']}/gps_75_pit/quantification/{config['gps']['quantification_file_name']}"
    params:
        max_peakgroup_q_value = config['gps']['combine']['max_peakgroup_q_value']
    shell:
        (
            "gps combine "
            "--input-files {input.score_files} "
            "--peptide-model {input.peptide_model} "
            "--protein-model {input.protein_model} "
            "--output {output.quant_matrix} "
            "--max-peakgroup-q-value {params.max_peakgroup_q_value}"
        )

