### Count Corrections ###

## Uniques Only ##

# count hits per taxid, if they are unique for this taxid
rule step5_unique_count_hits:
    input:
        tagT = "iteration{i}/taxonomy/tag.01.taxonomy",
        tagF = "iteration{i}/abundance/tag.01.filter",
        taxids = "iteration{i}/abundance/taxids.txt"
    output:
        counts = "iteration{i}/abundance/results.uniques.txt"
    params:
        dir = "iteration{i}/taxonomy"
    run:
        taxids = [taxid.strip("\n") for taxid in open(input.taxids).readlines()]
        tsvs =  [params.dir + "/" + taxid + ".tsv" for taxid in taxids]

        # count spectra occurrence, but only once per taxid
        # (thus taxa will not be filtered based on own ambiguous hits)
        counts = dict()
        for tsv in tsvs:
            seen = set() # to remember if spectra got counted already
            with open(tsv, 'r') as f_in:
                idx_title = f_in.readline().rstrip("\n").split("\t").index("Title")
                for line in f_in:
                    spectra = line.rstrip("\n").split("\t")[idx_title]
                    # if spectra not counted yet in this taxa, add to counts
                    if spectra not in seen:
                        seen.add(spectra)
                        if spectra not in counts:
                            counts[spectra] = 0
                        counts[spectra] += 1

        # check whether any unique PSMs are available
        if not any(c == 1 for c in counts.values()):
            print("WARNING:\tno unique PSMs available, no more analysis possible!")

        # count taxa, i.e. 1 per unique hit
        with open(output.counts, 'w') as f_out:
            for tsv in tsvs:
                taxa_count = float(0)
                with open(tsv, 'r') as f_in:
                    idx_title = f_in.readline().rstrip("\n").split("\t").index("Title")
                    for line in f_in:
                        spectra = line.rstrip("\n").split("\t")[idx_title]
                        if counts[spectra] == 1:
                            taxa_count += 1
                f_out.write(str(taxa_count) + "\n")

# calculate relative unique counts
rule step5_unique_relative_hits:
    input:
        tag = "iteration{i}/abundance/tag.01.filter",
        org = "iteration{i}/abundance/results.uniques.txt"
    output:
        rel = "iteration{i}/abundance/results.uniques.relative.txt"
    run:
        R(r'''
        counts <- as.numeric(readLines("{input.org}"))
        if (sum(counts) > 0) counts <- counts / sum(counts)
        writeLines(as.character(counts), "{output.rel}")
        ''')
