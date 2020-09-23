# code for pipasic rule step5_pipasic_correction in correction.pipasic.snakefile
import subprocess

taxids = [taxid.strip("\n") for taxid in open(snakemake.input.taxids).readlines()]
tsvs =  ','.join(["iteration" + snakemake.wildcards.i + "/taxonomy/" + taxid + ".tsv" for taxid in taxids])
input_string = "\"" + tsvs + "\""
exec_string = "python2 " + snakemake.params.wfpath + "/modules_python/pipasic/AbundanceCorrection.py " + \
              snakemake.input.sample + " " + input_string + " " + snakemake.input.matrix + " -r " + \
              snakemake.output.rel + " -c " + snakemake.output.cor + " 2>&1 | tee " + snakemake.params.log
subprocess.check_call(exec_string, shell=True)