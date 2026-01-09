#!/bin/bash nextflow
nextflow.enable.dsl=2

process MUNGE_SUMSTATS {
    module 'any/ldsc/1.0.1'

    input:
    tuple val(study), file(file_to_be_munged), val(metabolite)

    output:
    tuple val(metabolite), val(study), file("${study}_${metabolite}_munged.sumstats.gz")

    script:
    """
    munge_sumstats.py --out ${study}_${metabolite}_munged --sumstats ${file_to_be_munged}.tsv --merge-alleles ${params.snplist}
    """
}

//calculates genetic correlation and parses the log file for the result matrix
process CALCULATE_CORRELATION_BETWEEN_STUDIES {
    label 'correlation'
    module 'any/ldsc/1.0.1'

    input:
    tuple val(metabolite), file(submission_files), val(correlation_input), val(n_comparisons), path(ld_folder)

    output:
    file("${metabolite}_${n_comparisons}_corr_matrix.txt")

    script:
    """
    ldsc.py --rg $correlation_input --ref-ld-chr ${ld_folder}/ --w-ld-chr ${ld_folder}/ --out $metabolite
    grep p1 -A ${n_comparisons} ${metabolite}.log | tail -n+2 > ${metabolite}_${n_comparisons}_corr_matrix.txt
    """
}

process CALCULATE_CORRELATION_BETWEEN_METABOLITES {
    label 'correlation'
    module 'any/ldsc/1.0.1'

    input:
    tuple val(study), file(submission_files), val(correlation_input), val(n_comparisons), path(ld_folder)

    output:
    file("${study}_${n_comparisons}_corr_matrix.txt")

    script:
    """
    ldsc.py --rg $correlation_input --ref-ld-chr ${ld_folder}/ --w-ld-chr ${ld_folder}/ --out $study
    grep p1 -A ${n_comparisons} ${study}.log | tail -n+2 > ${study}_${n_comparisons}_corr_matrix.txt
    """
}

workflow between_studies {
    study_ch = Channel.fromPath(params.studies)
                    .splitCsv(header: true, sep: '\t', strip: true)
                    .map{row -> [row.study, row.regenie_folder]}
    metabolite_ch = Channel.fromPath(params.metabolites)
                        .splitCsv(header: false, sep: ",", strip: true)
                        .flatten()
    combined_ch = study_ch
                    .combine(metabolite_ch)
                    .map{ row ->
                        def metabolite_file = row[1]+row[2]
                        tuple(row[0], metabolite_file, row[2])
                    }
    
    MUNGE_SUMSTATS(combined_ch)

    ld_folder = Channel.fromPath(params.ld_folder)
                .collect()

    correlation_ch = MUNGE_SUMSTATS.out
        .groupTuple(by: 0)
        .map{ row ->
            def study_files = row[2].sort()
            def submission_files = []
            def correlation_inputs = []
            def n_comparisons = []
            def n_studies = study_files.size()
            for (i in 2..n_studies){
                def correlation_input = ""
                submission_files.add(study_files)
                n_comparisons.add(n_studies-i+1)
                for (study_file in study_files){
                    if (correlation_input == ""){
                        correlation_input = correlation_input + study_file
                    } else {
                        correlation_input = correlation_input + "," + study_file
                    }
                }
                correlation_inputs.add(correlation_input)

                study_files = study_files.tail()
            }
            tuple(row[0], submission_files, correlation_inputs, n_comparisons)
        }
        .transpose()
        .combine(ld_folder)

    CALCULATE_CORRELATION_BETWEEN_STUDIES(correlation_ch)

    CALCULATE_CORRELATION_BETWEEN_STUDIES.out.collectFile(name: "$params.output_file", storeDir: "$baseDir")
}

workflow between_metabolites{
    study_ch = Channel.fromPath(params.studies)
                    .splitCsv(header: true, sep: '\t', strip: true)
                    .map{row -> [row.study, row.regenie_folder]}
    metabolite_ch = Channel.fromPath(params.metabolites)
                        .splitCsv(header: false, sep: ",", strip: true)
                        .flatten()
    combined_ch = study_ch
                    .combine(metabolite_ch)
                    .map{ row ->
                        def metabolite_file = row[1]+row[2]
                        tuple(row[0], metabolite_file, row[2])
                    }
    
    MUNGE_SUMSTATS(combined_ch)

    ld_folder = Channel.fromPath(params.ld_folder)
            .collect()

    correlation_ch = MUNGE_SUMSTATS.out
        .groupTuple(by: 1)
        .map{ row ->
            def study_files = row[2].sort()
            def submission_files = []
            def correlation_inputs = []
            def n_comparisons = []
            def n_studies = study_files.size()
            for (i in 2..n_studies){
                def correlation_input = ""
                submission_files.add(study_files)
                n_comparisons.add(n_studies-i+1)
                for (study_file in study_files){
                    if (correlation_input == ""){
                        correlation_input = correlation_input + study_file
                    } else {
                        correlation_input = correlation_input + "," + study_file
                    }
                }
                correlation_inputs.add(correlation_input)

                study_files = study_files.tail()
            }
            tuple(row[1], submission_files, correlation_inputs, n_comparisons)
        }
        .transpose()
        .combine(ld_folder)

    CALCULATE_CORRELATION_BETWEEN_METABOLITES(correlation_ch)

    CALCULATE_CORRELATION_BETWEEN_METABOLITES.out.collectFile(name: "$params.output_file", storeDir: "$baseDir")
}