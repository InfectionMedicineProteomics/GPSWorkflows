#!/usr/bin/env python3
SAMPLES, = glob_wildcards(f"{config['base_file_paths']['mzml']}/{{sample}}.mzML")

print(SAMPLES)

rule all:
    input:
        quant_matrix = f"{config['base_file_paths']['results']}/pyprophet_gps_export/quantification/{config['gps']['quantification_file_name']}"

rule openswath_init:
    input:
        spectral_library = config['libraries']['spectral_library'],
        rt_library = config['libraries']['rt_library'],
        mzml = f"{config['base_file_paths']['mzml']}/{{sample}}.mzML"
    output:
        osw = f"{config['base_file_paths']['results']}/osw/{{sample}}.osw",
        chrom = f"{config['base_file_paths']['results']}/chromatograms/{{sample}}.sqMass"
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


rule pyprophet_subsample:
    input:
        osw = f"{config['base_file_paths']['results']}/osw/{{sample}}.osw",
    output:
        subsampled = f"{config['base_file_paths']['results']}/pyprophet_no_nrt/osw/{{sample}}.osws"
    container:
        config['containers']['pyprophet']
    threads: config['openswath']['threads']
    params:
        subsample_ratio = 0.1
    shell:
        (
            "pyprophet subsample "
            "--subsample_ratio {params.subsample_ratio} "
            "--in {input.osw} "
            "--out {output.subsampled}"
        )


rule pyprophet_merge:
    input:
        template = config['libraries']['spectral_library'],
        osws = expand(f"{config['base_file_paths']['results']}/pyprophet_no_nrt/osw/{{sample}}.osws", sample=SAMPLES)
    output:
        merged = f"{config['base_file_paths']['results']}/pyprophet_no_nrt/osw/merged.osw"
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
        merged = f"{config['base_file_paths']['results']}/pyprophet_no_nrt/osw/merged.osw"
    output:
        model = f"{config['base_file_paths']['results']}/pyprophet_no_nrt/models/scoring_model.osw"
    container:
        config['containers']['pyprophet']
    threads: config['pyprophet']['final_threads']
    params:
        level = 'ms2',
        ss_initial_fdr=0.15,
        ss_iteration_fdr=0.05,
        classifier = "LDA"
    shell:
        (
            "pyprophet score "
            "--in {input} "
            "--out {output.model} "
            "--level {params.level} "
            "--classifier {params.classifier} "
            "--ss_initial_fdr {params.ss_initial_fdr} "
            "--ss_iteration_fdr {params.ss_iteration_fdr} "
            "--threads {threads}"
        )

rule pyprophet_apply:
    input:
        scoring_model = f"{config['base_file_paths']['results']}/pyprophet_no_nrt/models/scoring_model.osw",
        osw = f"{config['base_file_paths']['results']}/osw/{{sample}}.osw",
    output:
        score_placeholder = f"{config['base_file_paths']['results']}/pyprophet_no_nrt/scored/{{sample}}"
    container:
        config['containers']['pyprophet']
    params:
        group_id = 'feature_id',
        level = 'ms2'
    threads: config['openswath']['threads']
    shell:
        (
            "pyprophet score "
            "--in {input.osw} "
            "--group_id {params.group_id} "
            "--apply_weights {input.scoring_model} "
            "--level {params.level} && "

            "touch {output.score_placeholder}"
        )

rule gps_export_pyprophet:
    input:
        osw = f"{config['base_file_paths']['results']}/osw/{{sample}}.osw",
        score_placeholder = f"{config['base_file_paths']['results']}/pyprophet_no_nrt/scored/{{sample}}"
    output:
        scored_pyprophet = f"{config['base_file_paths']['results']}/pyprophet_gps_export/scored/{{sample}}.scored.tsv"
    threads: config['gps']['score']['threads']
    shell:
        (
            "gps export "
            "--output-format pyprophet "
            "--input {input.osw} "
            "--output {output.scored_pyprophet}"
        )

rule gps_build_peptide_model:
    input:
        score_files = expand(f"{config['base_file_paths']['results']}/pyprophet_gps_export/scored/{{sample}}.scored.tsv", sample=SAMPLES)
    output:
        peptide_model = f"{config['base_file_paths']['results']}/pyprophet_gps_export/models/peptide.model"
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
        score_files = expand(f"{config['base_file_paths']['results']}/pyprophet_gps_export/scored/{{sample}}.scored.tsv", sample=SAMPLES)
    output:
        protein_model = f"{config['base_file_paths']['results']}/pyprophet_gps_export/models/protein.model"
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
        score_files = expand(f"{config['base_file_paths']['results']}/pyprophet_gps_export/scored/{{sample}}.scored.tsv", sample=SAMPLES),
        peptide_model = f"{config['base_file_paths']['results']}/pyprophet_gps_export/models/peptide.model",
        protein_model = f"{config['base_file_paths']['results']}/pyprophet_gps_export/models/protein.model"
    output:
        quant_matrix = f"{config['base_file_paths']['results']}/pyprophet_gps_export/quantification/{config['gps']['quantification_file_name']}"
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

