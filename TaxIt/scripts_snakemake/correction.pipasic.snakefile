### Count Corrections ###

# Using Pipasic relies on Python 2 interpreter and libraries. Those are not
# installed via conda yet since conda doesn’t support different versions of python
# in the same environment. However, they can be automatically installed in a
# separated environment on pipeline execution using the "--use-conda" parameter.
# Otherwise, the pipeline will expect a python2 executable available via your PATH
# environment (see [environments/python2.yaml](environments/python2.yaml) for
# additional library dependencies).


## Pipasic Abundance Correction ##

# digestion
rule step5_protein_digestion:
    input:
        tagT = "iteration{i}/taxonomy/tag.01.taxonomy",
        tagF = "iteration{i}/abundance/tag.01.filter",
        taxids = "iteration{i}/abundance/taxids.txt",
    output:
        tag = temp("iteration{i}/taxonomy/tag.02.digested")
    params:
        dir = "iteration{i}/taxonomy",
        wfpath = WORKFLOW_PATH
    conda:
        srcdir("../environments/python2.yaml")
    script:
        "../scripts_python/pipasic.digestion.py"

# merge non-unique psm, i.e. with same score and sequence
rule step5_merge_non_uniques:
    input:
        tag = "iteration{i}/taxonomy/tag.01.taxonomy",
        taxids = "iteration{i}/abundance/taxids.txt"
    output:
        tag = temp("iteration{i}/taxonomy/tag.03.merged")
    params:
        engine = SEARCH_ENGINE,
        wfpath = WORKFLOW_PATH
    run:
        taxids = [taxid.strip("\n") for taxid in open(input.taxids).readlines()]
        for taxid in taxids:
            tsv_in = "iteration" + wildcards.i + "/taxonomy/" + taxid + ".tsv"
            tsv_out = "iteration" + wildcards.i + "/taxonomy/" + taxid + ".merged.tsv"
            R(r'''
            merge.separaters <- list(msgfplus = ";", xtandem = "&&&")
            PARAMETERS <- list(search.engine = "{params.engine}")
            PARAMETERS$merge.separater <- merge.separaters[[PARAMETERS$search.engine]]
            source(paste0("{params.wfpath}", "/scripts_r/fdr.R"))

            psm.in <- psm.read("{tsv_in}")
            psm.out <- psm.merge(psm.in)
            psm.write(psm.out, "{tsv_out}")
            ''')
        shell("touch " + output.tag)

# similarity matrix
rule step5_pipasic_matrix:
    input:
        tag2 = "iteration{i}/taxonomy/tag.02.digested",
        tag3 = "iteration{i}/taxonomy/tag.03.merged",
        taxids = "iteration{i}/abundance/taxids.txt"
    output:
        matrix = "iteration{i}/abundance/pipasic.matrix.txt"
    params:
        wfpath = WORKFLOW_PATH,
        log = "logs/iteration{i}.pipasic.matrix.log"
    benchmark:
        "logs/iteration{i}.pipasic.matrix.txt"
    log:
        "logs/iteration{i}.pipasic.matrix.log"
    conda:
        srcdir("../environments/python2.yaml")
    script:
        "../scripts_python/pipasic.matrix.py"

# correction
rule step5_pipasic_correction:
    input:
        sample  = iteration_sample,
        matrix = "iteration{i}/abundance/pipasic.matrix.txt",
        taxids = "iteration{i}/abundance/taxids.txt"
    output:
        rel = "iteration{i}/abundance/results.pipasic.relative.txt",
        cor = "iteration{i}/abundance/results.pipasic.corrected.txt"
    params:
        wfpath = WORKFLOW_PATH,
        log = "logs/iteration{i}.pipasic.correction.log"
    benchmark:
        "logs/iteration{i}.pipasic.correction.txt"
    log:
        "logs/iteration{i}.pipasic.correction.log"
    conda:
        srcdir("../environments/python2.yaml")
    script:
        "../scripts_python/pipasic.correction.py"
