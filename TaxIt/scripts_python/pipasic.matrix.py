# code for pipasic rule step5_pipasic_matrix in correction.pipasic.snakefile
import subprocess

taxids = [taxid.strip("\n") for taxid in open(snakemake.input.taxids).readlines()]
digests = ','.join(["iteration" + snakemake.wildcards.i + "/taxonomy/" + taxid + ".fasta.digest" for taxid in taxids])
fastas =  ','.join(["iteration" + snakemake.wildcards.i + "/taxonomy/" + taxid + ".fasta" for taxid in taxids])
tsvs =  ','.join(["iteration" + snakemake.wildcards.i + "/taxonomy/" + taxid + ".merged.tsv" for taxid in taxids])
input_string = "\"" + digests + "\" \"" + tsvs + "\" \"" + fastas + "\""
exec_string = "python2 " + snakemake.params.wfpath + "/modules_python/pipasic/WeightedSimilarityMatrix.py " + \
              input_string + " -o " + snakemake.output.matrix + " 2>&1 | tee " + snakemake.params.log
subprocess.check_call(exec_string, shell=True)
subprocess.check_call("rm weightedPeptides*.txt", shell=True)