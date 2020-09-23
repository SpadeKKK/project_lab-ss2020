### Iteration Control ###

# sample per iteration
# use full sample each iteration so far
# or host-filtered sample, if a host db is provided
def iteration_sample(wc):
    # host filter iteration
    if wc.i == "F":
        return SAMPLE[:-3] + "0free.mgf"
    # normal iterations
    else:
        return "iterationF/sample_new.mgf" if DATABASE_HOST else SAMPLE[:-3] + "0free.mgf"

# target database per iteration
# either first database (e.g. NCBI RefSeq) or new database created by last iteration
def iteration_database(wc):
    # host filter iteration
    if wc.i == "F":
        return DATABASE_HOST + ".clean"
    # normal iterations
    else:
        return "iteration" + str(int(wc.i)+1) + "/output/db_targets.fasta.clean" \
        if int(wc.i) < ITERATIONS else DATABASE_REF + ".clean"

# corresponding decoy database
def iteration_database_decoys(wc):
    return iteration_database(wc) + ".decoys"

# STEP 1:
# execute spectral search
# (see search.*.snakefile)

# STEP 2:
# calculate FDR and clean up PSMs
# (see processing.snakefile)

# STEP 3:
# split up PSMs according to taxonomy
# (see taxonomy.snakefile)

# STEP 4:
# filter taxa with low abundance
rule step4_abundance_filter:
    input:
        taxids = "iteration{i}/taxonomy/taxids.txt",
        counts = "iteration{i}/abundance/results.original.txt"
    output:
        taxids = "iteration{i}/abundance/taxids.txt",
        counts = "iteration{i}/abundance/results.filtered.txt"
    params:
        type = FILTER
    script:
        "../scripts_r/filter.R"

# check whether any PSMs are left after filtering
rule step4_abundance_check:
    input:
        counts = "iteration{i}/abundance/results.filtered.txt"
    output:
        tag = "iteration{i}/abundance/tag.01.filter"
    run:
        if os.stat(input.counts).st_size == 0:
            print("ERROR:\tno PSMs left after filter " + FILTER + ", no more analysis possible!")
            os.remove(input.counts)
            sys.exit(1)
        shell("touch {output.tag}")

# STEP 5:
# pipasic abundance correction
# (see abundance.snakefile)

# STEP 6:
# select present taxa/organism(s)
RELATIVE = {"weighted" : "iteration{i}/abundance/results.weighted.relative.txt",
            "pipasic"  : "iteration{i}/abundance/results.pipasic.relative.txt",
            "uniques"  : "iteration{i}/abundance/results.uniques.relative.txt"}
rule step6_organism_selection:
    input:
        taxids = "iteration{i}/abundance/taxids.txt",
        org = "iteration{i}/abundance/results.filtered.txt",
        rel = RELATIVE[CORRECTION],
        names = "iteration{i}/abundance/taxnames.txt",
        clean = "iteration{i}/taxonomy/tag.03.cleanup"
    output:
        taxids = "iteration{i}/results/taxids.txt",
        names = "iteration{i}/results/taxnames.txt",
        plot = "iteration{i}/results/abundance.png",
        plot10org = "iteration{i}/results/abundance_top10org.png",
        plot10rel = "iteration{i}/results/abundance_top10rel.png",
    params:
        mode = MODE
    script:
        "../scripts_r/selection.R"

# STEP 7:
# download and create databases for next iteration
# (see download.snakefile)

# collect final results
rule final_results:
    input:
        "iteration1/results/taxids.txt",
        "iteration1/results/taxnames.txt",
        "iteration1/results/abundance.png",
        "iteration1/results/abundance_top10org.png",
        "iteration1/results/abundance_top10rel.png",
        "iteration1/results/abundance_stacked.pdf",
        "iteration1/results/abundance_stacked.png"
    output:
        "results/identified.taxid.txt",
        "results/identified.names.txt",
        "results/abundance.png",
        "results/abundance_stacked.pdf",
        "results/abundance_stacked.png",
    shell:
        "cp -pr iteration1/results/* results && " \
        "mv results/taxids.txt results/identified.taxid.txt && " \
        "mv results/taxnames.txt results/identified.names.txt"
