# code for Python 2 version rule in versions.snakefile
import subprocess

def software_version(call):
    # execute software version call on shell
    if not call:
        return "not available"
    try:
        return str(subprocess.check_output(call, shell=True, stderr=subprocess.STDOUT, universal_newlines=True)).rstrip()
    except subprocess.CalledProcessError as e:
        return "not available"

calls = {m : 'python2 -c "import ' + m + '; print ' + m + '.__version__"' for m in snakemake.params.libs}
calls["Python2"] = "which python2 >/dev/null && python2 -c \"import sys; print(sys.version.split(' ')[0])\""
with open(snakemake.output.ver, 'w') as f_out:
    [f_out.write(s + ":\t" + software_version(calls[s]) + "\n") for s in sorted(calls.keys())]