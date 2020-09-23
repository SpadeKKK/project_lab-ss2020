### MS-GF+ Spectral Search ###

# execute MS-GF+
rule step1_msgfplus_execute:
    input:
        spectra = iteration_sample,
        database = "iteration{i}/search/msgfplus/db_merged.fasta",
        mods = lambda wc: MSGFPLUS_MODS if MSGFPLUS_MODS else []
    output:
        mzid = "iteration{i}/search/msgfplus/msgfplus.output.mzid",
        temp1 = temp("iteration{i}/search/msgfplus/db_merged.canno"),
        temp2 = temp("iteration{i}/search/msgfplus/db_merged.cnlcp"),
        temp3 = temp("iteration{i}/search/msgfplus/db_merged.csarr"),
        temp4 = temp("iteration{i}/search/msgfplus/db_merged.cseq")
    params:
        pmt = MSGFPLUS_PMT,
        tda = 0, # decoys are created separately
        inst = MSGFPLUS_INST,
        minc = MSGFPLUS_MINC,
        maxc = MSGFPLUS_MAXC,
        mods = lambda wc: "-mod " + MSGFPLUS_MODS if MSGFPLUS_MODS else [],
        exec = MSGFPLUS_EXEC,
        add = MSGFPLUS_PARAS
    threads:
        THREADS_MAX
    benchmark:
        "logs/iteration{i}.msgfplus.txt"
    log:
        "logs/iteration{i}.msgfplus.log"
    shell:
        "{params.exec} " \
        "-s {input.spectra} -o {output.mzid} " \
        "-d {input.database} " \
        "-t {params.pmt} " \
        "-thread {threads} " \
        "-tda {params.tda} " \
        "-inst {params.inst} " \
        "-minCharge {params.minc} " \
        "-maxCharge {params.maxc} " \
        "{params.mods} " \
        "{params.add} " \
        "2>&1 | tee {log}"

# convert MS-GF+ results to tsv
rule step1_msgfplus_mzid_to_tsv:
    input:
        "iteration{i}/search/msgfplus/msgfplus.output.mzid"
    output:
        "iteration{i}/search/msgfplus/msgfplus.convert.tsv"
    params:
        exec = "java -Xmx" + MEMORY_MAX + " -cp " + MSGFPLUS_JAR if MSGFPLUS_JAR else "msgf_plus -Xmx" + MEMORY_MAX
    shell:
        "{params.exec} edu.ucsd.msjava.ui.MzIDToTsv -i {input} -o {output} -showDecoy 1"

# remove digestion infos (such as "(pre=K,post=L)") from Protein column
rule step1_msgfplus_clean_id:
    input:
        tsv = "iteration{i}/search/msgfplus/msgfplus.convert.tsv"
    output:
        tsv = "iteration{i}/search/msgfplus/msgfplus.output.tsv"
    run:
        import re
        protein_index = 10 # column 11 in MS-GF+ tsv output
        cleavage_regex = re.compile("\(pre=[A-Z\-]{1},post=[A-Z\-]{1}\)")
        with open(input.tsv, 'r') as f_in:
            with open(output.tsv, 'w') as f_out:
                for line in f_in:
                    splits = line.split("\t")
                    proteins = splits[protein_index]
                    proteins = cleavage_regex.sub("", proteins)
                    splits[protein_index] = proteins
                    f_out.write("\t".join(splits))

# merge search database for usage with MS-GF+
rule step1_msgfplus_merge_dbs:
    input:
        db_main_targets = iteration_database,
        db_main_decoys = iteration_database_decoys,
        db_conts_targets = DATABASE_CONTS + ".clean",
        db_conts_decoys = DATABASE_CONTS + ".clean.decoys",
        prev_decoys = iteration_decoys
    output:
        db_merged = temp("iteration{i}/search/msgfplus/db_merged.fasta")
    shell:
        "cat {input} > {output.db_merged}"
