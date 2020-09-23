### Host Spectra Filtering ###

# 1. search with xtandem (via iteration rules)

# 2. FDR and PSMs clean up (via iteration rules)

# 3. extract spectra title
rule hostfilter_spectra_title:
    input:
        psms = "iterationF/psms/psms.txt"
    output:
        titles = "iterationF/host_spectra_title.txt"
    run:
        titles = set()
        with open(input.psms, 'r') as f:
            header = f.readline().strip().split("\t")
            idx_title = header.index("Title")
            out = open(output.titles, "w")
            for line in f:
                splits = line.split("\t")
                titles.add(splits[idx_title])
        with open(output.titles, "w") as out:
            out.writelines([title + "\n" for title in titles])

# 4. remove host spectra
rule hostfilter_spectra_selection:
    input:
        sample = SAMPLE[:-3] + "0free.mgf",
        titles = "iterationF/host_spectra_title.txt"
    output:
        new = temp("iterationF/sample_new.mgf"),
        rej = temp("iterationF/sample_rej.mgf")
    params:
        jar = WORKFLOW_PATH + "/modules_java/spectra-filter/target/spectra-filter-1.0-SNAPSHOT.jar",
        memory  = "-Xmx" + MEMORY_MAX
    shell:
        "ln -sfr {input.sample} iterationF/sample.mgf && " \
        "java {params.memory} -jar {params.jar} -i iterationF/sample.mgf -f {input.titles} && " \
        "rm iterationF/sample.mgf"
