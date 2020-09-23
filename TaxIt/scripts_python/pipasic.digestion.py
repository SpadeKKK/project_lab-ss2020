# code for pipasic rule step5_protein_digestion in correction.pipasic.snakefile
import subprocess

taxids = [taxid.strip("\n") for taxid in open(snakemake.input.taxids).readlines()]
fastas = ["iteration" + snakemake.wildcards.i + "/taxonomy/" + taxid + ".fasta" for taxid in taxids]
# fastas = sorted(glob(params.dir + "/*.fasta"))
print(fastas)
for fasta in fastas:
    # digestion based on: https://github.com/yafeng/trypsin
    exec_string = "python2 " + snakemake.params.wfpath + "/tools/trypsin/trypsin.py --input \"" + fasta + "\" --output \"" + fasta + ".tmp\" --miss 1"
    subprocess.check_call(exec_string, shell=True)

    # convert digested peptides to proper fasta
    f_out = open(fasta + ".digest", 'w')
    for line in open(fasta + ".tmp", 'r'):
        splits = line.rstrip("\n").split("\t")
        # filter peptides by a minimum size of 6 amino acids
        if (len(splits[1]) >= 6):
            f_out.write(">" + splits[0] + "\n")
            f_out.write(splits[1] + "\n")
    f_out.close()
subprocess.check_call("rm " + snakemake.params.dir + "/*.fasta.tmp && touch " + snakemake.output.tag, shell=True)
