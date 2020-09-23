### Count Corrections ###

## Spectrum Frequency Weighting ##

# count hits per taxid, weighted by spectra occurrence/frequency
rule step5_weighted_count_hits:
    input:
        tagT = "iteration{i}/taxonomy/tag.01.taxonomy",
        tagF = "iteration{i}/abundance/tag.01.filter",
        taxids = "iteration{i}/abundance/taxids.txt"
    output:
        counts = "iteration{i}/abundance/results.weighted.txt"
    params:
        dir = "iteration{i}/taxonomy"
    run:
        taxids = [taxid.strip("\n") for taxid in open(input.taxids).readlines()]
        tsvs =  [params.dir + "/" + taxid + ".tsv" for taxid in taxids]

        # count spectra occurrence
        counts = dict()
        for tsv in tsvs:
            with open(tsv, 'r') as f_in:
                idx_title = f_in.readline().rstrip("\n").split("\t").index("Title")
                for line in f_in:
                    spectra = line.rstrip("\n").split("\t")[idx_title]
                    if spectra not in counts:
                        counts[spectra] = 0
                    counts[spectra] += 1

        # count taxa, i.e. 1 per hit / # spectrum occurrence
        with open(output.counts, 'w') as f_out:
            for tsv in tsvs:
                taxa_count = float(0)
                with open(tsv, 'r') as f_in:
                    idx_title = f_in.readline().rstrip("\n").split("\t").index("Title")
                    for line in f_in:
                        spectra = line.rstrip("\n").split("\t")[idx_title]
                        taxa_count += 1 / counts[spectra]
                f_out.write(str(taxa_count) + "\n")

# calculate relative weighted abundance
rule step5_weighted_relative_hits:
    input:
        org = "iteration{i}/abundance/results.weighted.txt"
    output:
        rel = "iteration{i}/abundance/results.weighted.relative.txt"
    run:
        R(r'''
        counts.org <- as.numeric(readLines("{input.org}"))
        counts.rel <- counts.org / sum(counts.org)
        writeLines(as.character(counts.rel), "{output.rel}")
        ''')
