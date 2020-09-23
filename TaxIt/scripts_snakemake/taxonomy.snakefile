### Taxonomic Rules ###

# STEP 3:
# identify organism (taxid) per PSM
rule step3_taxonomic_identification:
    input:
        psms = lambda wc: "iteration" + wc.i + "/psms/psms.txt",
        tax_mapping = lambda wc: TAX_MAPPING_NCBI if int(wc.i) == 1 and ITERATIONS > 1 else TAX_MAPPING,
    output:
        tax = "iteration{i}/taxonomy/protein2taxid.txt"
    params:
        jar = WORKFLOW_PATH + "/modules_java/taxonomic-identification/target/taxonomic-identification-1.0-SNAPSHOT-jar-with-dependencies.jar",
        outdir = "iteration{i}/taxonomy",
        memory  = "-Xmx" + MEMORY_MAX
    benchmark:
        "logs/iteration{i}.taxonomic.identification.txt"
    log:
        "logs/iteration{i}.taxonomic.identification.log"
    shell:
        "java {params.memory} -jar {params.jar} " \
        "{input.psms} {input.tax_mapping} {output.tax} 2>&1 | tee {log}"

# split up PSMs according to taxonomy
rule step3_taxonomic_splitting:
    input:
        psms = lambda wc: "iteration" + wc.i + "/psms/psms.txt",
        fasta = iteration_database,
        tax_mapping_ref = lambda wc: TAX_MAPPING_NCBI if int(wc.i) == 1 and ITERATIONS > 1 else TAX_MAPPING,
        tax_mapping = "iteration{i}/taxonomy/protein2taxid.txt",
        tax_nodes = TAX_NODES,
        tax_names = TAX_NAMES
    output:
        tag = temp("iteration{i}/taxonomy/tag.01.taxonomy")
    params:
        jar = WORKFLOW_PATH + "/modules_java/taxonomic-separation/target/taxonomic-separation-1.0-SNAPSHOT-jar-with-dependencies.jar",
        iter = "{i}",
        outdir = "iteration{i}/taxonomy",
        memory  = "-Xmx" + MEMORY_MAX
    benchmark:
        "logs/iteration{i}.taxonomic.splitting.txt"
    log:
        "logs/iteration{i}.taxonomic.splitting.log"
    shell:
        "java {params.memory} -jar {params.jar} {params.iter} {input.psms} {input.fasta} " \
        "{input.tax_mapping_ref} {input.tax_mapping} {input.tax_nodes} {input.tax_names} {params.outdir} 2>&1 | tee {log} && touch {output.tag}"

# taxid collection
rule taxid_collection:
    input:
        tag = "iteration{i}/taxonomy/tag.01.taxonomy"
    output:
        taxids = "iteration{i}/taxonomy/taxids.txt"
    params:
        dir = "iteration{i}/taxonomy"
    run:
        taxids, = glob_wildcards(params.dir + "/{taxid}.fasta")
        taxids = [str(taxid) + "\n" for taxid in sorted([int(i) for i in taxids])]
        f_out = open(output.taxids, 'w')
        f_out.writelines(taxids)
        f_out.close()

# count hits per taxid
rule count_hits:
    input:
        tag = "iteration{i}/taxonomy/tag.01.taxonomy",
        taxids = "iteration{i}/taxonomy/taxids.txt"
    output:
        counts = "iteration{i}/abundance/results.original.txt"
    params:
        dir = "iteration{i}/taxonomy"
    run:
        taxids = [taxid.strip("\n") for taxid in open(input.taxids).readlines()]
        tsvs =  [params.dir + "/" + taxid + ".tsv" for taxid in taxids]

        with open(output.counts, 'w') as f_out:
            for tsv in tsvs:
                hit_count = -1 # subtract header
                with open(tsv, 'r') as f_in:
                    for line in f_in: hit_count += 1
                f_out.write(str(hit_count) + "\n")

# get scientific names
rule taxname_collection:
    input:
        taxids = "iteration{i}/abundance/taxids.txt",
        names = TAX_NAMES
    output:
        names = "iteration{i}/abundance/taxnames.txt"
    run:
        taxids = set()
        with open(input.taxids, 'r') as f_in:
            [taxids.add(line.strip("\n")) for line in f_in]
        with open(input.names, 'r') as f_in:
            with open(output.names, 'w') as f_out:
                for line in f_in:
                    splits = line.strip("\n").split("\t|\t")
                    taxid = splits[0]
                    name_class = splits[3]
                    if taxid in taxids and name_class.startswith("scientific name"):
                        name = splits[1].replace('/', ' ')
                        f_out.write(taxid + "\t" + name + "\n")


### Cleanup Taxa Files ###

# data for counting and correcting hits takes a lot of disc space,
# thus it is removed after calculations
rule step5_cleanup_taxa_data:
    input:
        tag1 = "iteration{i}/taxonomy/tag.01.taxonomy",
        tag2 = "iteration{i}/taxonomy/tag.02.digested" if CORRECTION == "pipasic" else [],
        tag3 = "iteration{i}/taxonomy/tag.03.merged" if CORRECTION == "pipasic" else [],
        rel = RELATIVE[CORRECTION] # see iterations.snakefile STEP 6
    output:
        tag = "iteration{i}/taxonomy/tag.03.cleanup"
    params:
        rem = "iteration{i}/taxonomy/*.tsv iteration{i}/taxonomy/*.fasta iteration{i}/taxonomy/*.merged " \
              "iteration{i}/taxonomy/*.digest" if CORRECTION == "pipasic" else
              "iteration{i}/taxonomy/*.tsv iteration{i}/taxonomy/*.fasta"
    shell:
        "rm {params.rem} && touch {output.tag}"
