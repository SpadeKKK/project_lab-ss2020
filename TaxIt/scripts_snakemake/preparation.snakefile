### Input Data Preparations ###

# remove descriptions from fasta header, keep IDs only
rule fasta_header_clean:
    input:
        prots = "{proteins}.fasta"
    output:
        prots = "{proteins}.fasta.clean"
    run:
        f = open(output.prots,'w')
        for line in open(input.prots, 'r'):
            if line.startswith(">"):
                prefix = line.strip("\n").split(" ")[0]
                f.write(prefix + "\n")
            else:
                f.write(line)
        f.close()

# remove descriptions from fasta header, keep IDs only
rule faa_header_clean:
    input:
        prots = "{proteins}.faa"
    output:
        prots = "{proteins}.faa.clean"
    run:
        f = open(output.prots,'w')
        for line in open(input.prots, 'r'):
            if line.startswith(">"):
                prefix = line.strip("\n").split(" ")[0]
                f.write(prefix + "\n")
            else:
                f.write(line)
        f.close()

# remove 0 intensities from mgf file
# as well as comment lines
rule mgf_intensity_filter:
    input:
        mgf = "{spectra}.mgf"
    output:
        mgf = "{spectra}.0free.mgf"
    shell:
        r"grep -Pv '(^\d+(\.\d+)?[ \t]0+(\.0+)?$)|(^#)' {input.mgf} > {output.mgf}"

# create decoys
# based on: http://search.cpan.org/~alexmass/InSilicoSpectro-Databanks-0.0.43/scripts/fasta-decoy.pl
rule fasta_decoys:
    input:
        db_target = "{fasta}.clean"
    output:
        db_decoys = "{fasta}.clean.decoys",
    shell:
        WORKFLOW_PATH + "/tools/decoys/fasta-decoy.pl --ac-prefix=XXX_ --in={input.db_target} --method=reverse --out={output.db_decoys}"

# remove descriptions from fasta header, keep IDs only + add contaminant label
rule fasta_header_clean_conts:
    input:
        unlabeled = DATABASE_CONTS
    output:
        labeled = DATABASE_CONTS + ".clean"
    run:
        with open(input.unlabeled, 'r') as f_in:
            with open(output.labeled, 'w') as f_out:
                for line in f_in:
                    if line.startswith(">"):
                        prefix = line[1:].strip("\n").split(" ")[0]
                        f_out.write(">contaminant|" + prefix + "\n")
                    else:
                        f_out.write(line)
