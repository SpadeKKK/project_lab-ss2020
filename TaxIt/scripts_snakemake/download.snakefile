### Download Rules ###

# STEP 7:
# download and create databases for next iteration
rule step7_database_downloads:
    input:
        taxids = "iteration{i}/results/taxids.txt",
        visual = "iteration{i}/results/abundance_stacked.pdf",
        tax_nodes = TAX_NODES
    output:
        fasta = "iteration{i}/output/all.fasta"
    params:
        jar = WORKFLOW_PATH + "/modules_java/resource-downloads/target/resource-downloads-1.0-SNAPSHOT-jar-with-dependencies.jar",
        outdir = "iteration{i}/output/",
        memory  = "-Xmx" + MEMORY_MAX,
        proxy = lambda wc: "-Dhttps.proxyHost=" + HTTPS_PROXY_URL +
                           " -Dhttps.proxyPort=" + HTTPS_PROXY_PORT if HTTPS_PROXY_URL else ""
    benchmark:
        "logs/iteration{i}.downloads.txt"
    log:
        "logs/iteration{i}.downloads.log"
    run:
        # read in identified taxids (after correction)
        taxids = [taxid.strip("\n") for taxid in open(input.taxids).readlines()]

        taxid2children = dict()
        for line in open(input.tax_nodes, "r"):
            splits = line.split("\t|\t")
            taxid_child = splits[0]
            taxid_parent = splits[1]
            if taxid_parent in taxid2children:
                taxid2children[taxid_parent].append(taxid_child)
            else:
                taxid2children[taxid_parent] = [taxid_child]

        def getChildren(parent, descendants):
            if parent in taxid2children:
                for child in taxid2children[parent]:
                    descendants.append(child)
                    getChildren(child, descendants)

        print("looking for descendants of " + str(taxids) + "...")
        descendants = []
        for taxid in taxids:
            getChildren(taxid, descendants)

        # check whether any strains were available
        if descendants:
            print("found descendants: " + str(descendants))
        else:
            print("ERROR:\tfound no descendants for " + str(taxids) + ", no more iterations possible!\n" +
                  "\tintermediate results (species) can be found in iteration2/results.")
            sys.exit(1)

        taxid_string = ','.join(descendants)
        shell("java {params.memory} {params.proxy} -jar {params.jar} " + taxid_string + " {params.outdir} all 2>&1 | tee {log}")

        fastas =  ' '.join(["iteration" + wildcards.i + "/output/" + descendant + "_exp_all.fasta" for descendant in descendants])
        shell("cat " + fastas + " > {output.fasta}")

        # check whether any strain proteins were available
        if os.stat(output.fasta).st_size == 0:
            print("ERROR:\tfound no strain proteins for " + str(taxids) + ", no more iterations possible!\n" +
                  "\tintermediate results (species) can be found in iteration2/results.")
            sys.exit(1)

# remove duplicated sequences (shared proteins)
# based on SeqKit: http://bioinf.shenwei.me/seqkit/usage/#rmdup
rule step7_database_duplicate_remove:
    input:
        all = "iteration{i}/output/all.fasta"
    output:
        uni = "iteration{i}/output/db_targets.fasta"
    shell:
        "seqkit rmdup {input.all} > {output.uni}"
