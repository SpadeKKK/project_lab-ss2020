### FDR and PSM Control ###

# additional decoy per iteration
# either no decoys or decoys from previous iterations (including host filtering)
def iteration_decoys(wc):
    if (wc.i == "F" or not PASS_DECOYS):
        decoys = []
    elif (int(wc.i) == ITERATIONS):
        decoys = "iterationF/psms/decoys.fasta" if DATABASE_HOST else []
    else:
        decoys = "iteration" + str(int(wc.i)+1) + "/psms/decoys.fasta"
    return decoys

# STEP 2:
# calculate FDR and clean up PSMs
rule step2_psm_processing:
    input:
        psms_in = lambda wc: "iteration" + wc.i + "/search/" + SEARCH_RESULTS[SEARCH_ENGINE]
    output:
        psms_out = "iteration{i}/psms/psms.txt",
        psms_dec = "iteration{i}/psms/tmp/psms.decoys.tsv",
        psms_tmp = "iteration{i}/psms/tmp/"
    params:
        temp = "iteration{i}/psms/tmp/",
        engine = SEARCH_ENGINE,
        wfpath = WORKFLOW_PATH,
        fdr = FDR_CUT,
    run:
        R('''
        DEBUGGING <- FALSE

        merge.separaters <- list(msgfplus = ";", xtandem = "&&&")
        PARAMETERS <- list(search.engine = "{params.engine}")
        PARAMETERS$merge.separater <- merge.separaters[[PARAMETERS$search.engine]]

        WFPATH <- "{params.wfpath}"
        print(WFPATH)
        source(paste0(WFPATH, "/scripts_r/fdr.R"))
        source(paste0(WFPATH, "/scripts_r/decoys.R"))

        in.file <- "{input.psms_in}"
        out.file <- "{output.psms_out}"
        print(in.file)
        print(out.file)

        psms.org <- psm.read(in.file)

        wd.org <- getwd()
        wd.tmp <- "{params.temp}"
        print(wd.tmp)
        dir.create(wd.tmp, showWarnings = TRUE, recursive = TRUE)
        setwd(wd.tmp)
        psms.new <- psm.parse(psms.org, {params.fdr})

        setwd(wd.org)
        if (nrow(psms.new) > 0 ){{
            psm.write(psms.new, out.file)
        }}
        ''')
        # check whether any PSMs were exported
        if not os.path.isfile(output.psms_out):
            print("ERROR:\tno PSMs left after TDA FDR adjustment, no more analysis possible!")
            sys.exit(1)

# extract identified decoy title
rule step2_psm_decoy_titles:
    input:
        dec_psms = "iteration{i}/psms/tmp/psms.decoys.tsv"
    output:
        dec_title = "iteration{i}/psms/decoys.title"
    run:
        with open(input.dec_psms, 'r') as f_in:
            title_idx = f_in.readline().rstrip("\n").split("\t").index("Protein")
            with open(output.dec_title, 'w') as f_out:
                for line in f_in:
                    f_out.write(line.rstrip("\n").split("\t")[title_idx] + "\n")

# extract identified decoy sequences
rule step2_psm_decoy_fasta:
    input:
        dec_title = "iteration{i}/psms/decoys.title",
        dec_db    = iteration_database_decoys, # current decoy sequences
        dec_prev  = iteration_decoys # previous decoy sequences
    output:
        dec_fasta = "iteration{i}/psms/decoys.fasta"
    run:
        # get decoy title
        decoys = set()
        with open(input.dec_title, 'r') as f_in:
            for line in f_in:
                decoys.add(line.rstrip("\n"))

        # iterate proteins
        with open(output.dec_fasta, 'w') as f_out:
            with open(input.dec_db, 'r') as f_in:
                is_decoy = False
                for line in f_in:
                    if line.startswith(">"):
                        title = line[1:].strip("\n")
                        is_decoy = title in decoys
                    if is_decoy:
                        f_out.write(line)
            if input.dec_prev:
                with open(input.dec_prev, 'r') as f_in:
                    is_decoy = False
                    for line in f_in:
                        if line.startswith(">"):
                            title = line[1:].strip("\n")
                            is_decoy = title in decoys
                        if is_decoy:
                            f_out.write(line)
