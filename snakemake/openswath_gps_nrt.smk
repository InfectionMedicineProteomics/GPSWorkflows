#!/usr/bin/env python3
SAMPLES, = glob_wildcards(f"{config['base_file_paths']['mzml']}/{{sample}}.mzml")

print(SAMPLES)

rule all:
    input:
        f"{config['base_file_paths']['results']}/gps_nrt/final/quantification/{config['gps']['quantification_file_name']}"


rule openswath_init:
    input:
        spectral_library = config['libraries']['spectral_library'],
        rt_library = config['libraries']['rt_library'],
        mzml = f"{config['base_file_paths']['mzml']}/{{sample}}.mzml"
    output:
        osw = f"{config['base_file_paths']['results']}/gps_nrt/init/osw/{{sample}}.osw",
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
            "-swath_windows_file {params.swath_windows_file} "
            "-tr_irt {input.rt_library} "
            "-threads {threads} "
            "-enable_ms1 true "
            "-enable_ipf true "
            "-min_rsq {params.min_rsq} "
            "-Scoring:Scores:use_uis_scores "
            "-Scoring:stop_report_after_feature {params.scoring_stop_report_after_feature} "
            "-Scoring:TransitionGroupPicker:compute_peak_quality "
            "-RTNormalization:alignmentMethod {params.rt_normalization_alignment_method} "
            "-min_rsq {params.min_rsq} "
            "-batchSize {params.batch_size}"
        )

rule gps_score_init:
    input:
        osw = f"{config['base_file_paths']['results']}/gps_nrt/init/osw/{{sample}}.osw",
    output:
        score_file = f"{config['base_file_paths']['results']}/gps_nrt/init/scored/{{sample}}.scored.tsv"
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
            "--threads {threads}"
        )

rule gps_createlib_rt_middle:
    input:
        predicted_files = f"{config['base_file_paths']['results']}/gps_nrt/init/scored/{{sample}}.scored.tsv",
        spectral_library = config['libraries']['spectral_library'],
    output:
        new_spectral_library = f"{config['base_file_paths']['results']}/gps_nrt/middle/libraries/{{sample}}.rt_library.tsv"
    params:
        num_rt_bins = 35,
        num_peptides_per_bin = 5
    threads: 10
    shell:
        (
            "gps createlib "
            "--input {input.predicted_files} "
            "--output {output.new_spectral_library} "
            "--pqp {input.spectral_library} "
            "--library-type rt "
            "--num-rt-bins {params.num_rt_bins} "
            "--num-peptides-per-bin {params.num_peptides_per_bin}"
        )


rule convert_library_middle:
    input:
        new_spectral_library = f"{config['base_file_paths']['results']}/gps_nrt/middle/libraries/{{sample}}.rt_library.tsv"
    output:
        pqp = f"{config['base_file_paths']['results']}/gps_nrt/middle/libraries/{{sample}}.rt_library.pqp"
    container:
        config['containers']['openms']
    shell:
        (

            "TargetedFileConverter -in {input.new_spectral_library} -out {output.pqp}"
        )

rule openswath_middle:
    input:
        spectral_library = config['libraries']['spectral_library'],
        rt_library = f"{config['base_file_paths']['results']}/gps_nrt/middle/libraries/{{sample}}.rt_library.pqp",
        mzml = f"{config['base_file_paths']['mzml']}/{{sample}}.mzml"
    output:
        osw = f"{config['base_file_paths']['results']}/gps_nrt/middle/osw/{{sample}}.osw",
    container:
        config['containers']['openms']
    params:
        rt_normalization_alignment_method = 'lowess',
        scoring_stop_report_after_feature = config['openswath']['stop_after_features'],
        swath_windows_file = config['openswath']['swath_windows_file'],
        min_rsq = config['openswath']['min_rsq'],
        batch_size = config['openswath']['batch_size'],
        num_rt_bins = 25,
        min_bins_filled = 20
    threads:
        config["openswath"]['threads']
    shell:
        (
            "OpenSwathWorkflow "
            "-in {input.mzml} "
            "-tr {input.spectral_library} "
            "-out_osw {output.osw} "
            "-swath_windows_file {params.swath_windows_file} "
            "-tr_irt {input.rt_library} "
            "-threads {threads} "
            "-enable_ms1 true "
            "-enable_ipf true "
            "-min_rsq {params.min_rsq} "
            "-Scoring:Scores:use_uis_scores "
            "-Scoring:stop_report_after_feature {params.scoring_stop_report_after_feature} "
            "-RTNormalization:NrRTBins {params.num_rt_bins} "
            "-RTNormalization:MinBinsFilled {params.min_bins_filled} "
            "-Scoring:TransitionGroupPicker:compute_peak_quality "
            "-RTNormalization:alignmentMethod {params.rt_normalization_alignment_method} "
            "-batchSize {params.batch_size}"
        )

rule gps_score_middle:
    input:
        osw = f"{config['base_file_paths']['results']}/gps_nrt/middle/osw/{{sample}}.osw",
    output:
        score_file = f"{config['base_file_paths']['results']}/gps_nrt/middle/scored/{{sample}}.scored.tsv"
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
            "--threads {threads}"
        )

rule gps_createlib_rt_final:
    input:
        predicted_files = f"{config['base_file_paths']['results']}/gps_nrt/middle/scored/{{sample}}.scored.tsv",
        spectral_library = config['libraries']['spectral_library'],
    output:
        new_spectral_library = f"{config['base_file_paths']['results']}/gps_nrt/final/libraries/{{sample}}.rt_library.tsv"
    params:
        num_rt_bins = 50,
        num_peptides_per_bin = 5
    threads: 10
    shell:
        (
            "gps createlib "
            "--input {input.predicted_files} "
            "--output {output.new_spectral_library} "
            "--pqp {input.spectral_library} "
            "--library-type rt "
            "--num-rt-bins {params.num_rt_bins} "
            "--num-peptides-per-bin {params.num_peptides_per_bin}"
        )

rule convert_library_final:
    input:
        new_spectral_library = f"{config['base_file_paths']['results']}/gps_nrt/final/libraries/{{sample}}.rt_library.tsv"
    output:
        pqp = f"{config['base_file_paths']['results']}/gps_nrt/final/libraries/{{sample}}.rt_library.pqp"
    container:
        config['containers']['openms']
    shell:
        (

            "TargetedFileConverter -in {input.new_spectral_library} -out {output.pqp}"
        )

rule openswath_final:
    input:
        spectral_library = config['libraries']['spectral_library'],
        rt_library = f"{config['base_file_paths']['results']}/gps_nrt/final/libraries/{{sample}}.rt_library.pqp",
        mzml = f"{config['base_file_paths']['mzml']}/{{sample}}.mzml"
    output:
        osw = f"{config['base_file_paths']['results']}/gps_nrt/final/osw/{{sample}}.osw",
    container:
        config['containers']['openms']
    params:
        rt_normalization_alignment_method = 'lowess',
        scoring_stop_report_after_feature = config['openswath']['stop_after_features'],
        swath_windows_file = config['openswath']['swath_windows_file'],
        min_rsq = config['openswath']['min_rsq'],
        batch_size = config['openswath']['batch_size'],
        num_rt_bins = 50,
        min_bins_filled = 15
    threads:
        config["openswath"]['threads']
    shell:
        (
            "OpenSwathWorkflow "
            "-in {input.mzml} "
            "-tr {input.spectral_library} "
            "-out_osw {output.osw} "
            "-swath_windows_file {params.swath_windows_file} "
            "-tr_irt {input.rt_library} "
            "-threads {threads} "
            "-enable_ms1 true "
            "-enable_ipf true "
            "-min_rsq {params.min_rsq} "
            "-Scoring:Scores:use_uis_scores "
            "-Scoring:stop_report_after_feature {params.scoring_stop_report_after_feature} "
            "-RTNormalization:NrRTBins {params.num_rt_bins} "
            "-RTNormalization:MinBinsFilled {params.min_bins_filled} "
            "-Scoring:TransitionGroupPicker:compute_peak_quality "
            "-mz_correction_function weighted_quadratic_regression_delta_ppm "
            "-RTNormalization:alignmentMethod {params.rt_normalization_alignment_method} "
            "-batchSize {params.batch_size}"
        )

rule gps_score_final:
    input:
        osw = f"{config['base_file_paths']['results']}/gps_nrt/final/osw/{{sample}}.osw",
    output:
        score_file = f"{config['base_file_paths']['results']}/gps_nrt/final/scored/{{sample}}.scored.tsv"
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
        score_files = expand(f"{config['base_file_paths']['results']}/gps_nrt/final/scored/{{sample}}.scored.tsv", sample=SAMPLES)
    output:
        peptide_model = f"{config['base_file_paths']['results']}/gps_nrt/final/models/peptide.model"
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
        score_files = expand(f"{config['base_file_paths']['results']}/gps_nrt/final/scored/{{sample}}.scored.tsv", sample=SAMPLES)
    output:
        protein_model = f"{config['base_file_paths']['results']}/gps_nrt/final/models/protein.model"
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
        score_files = expand(f"{config['base_file_paths']['results']}/gps_nrt/final/scored/{{sample}}.scored.tsv", sample=SAMPLES),
        peptide_model = f"{config['base_file_paths']['results']}/gps_nrt/final/models/peptide.model",
        protein_model = f"{config['base_file_paths']['results']}/gps_nrt/final/models/protein.model"
    output:
        quant_matrix = f"{config['base_file_paths']['results']}/gps_nrt/final/quantification/{config['gps']['quantification_file_name']}"
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

