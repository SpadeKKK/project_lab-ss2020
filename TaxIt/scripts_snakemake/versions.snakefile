### Collect Pipeline, Software & Package Version ###

# merge pipeline, software and R package version
rule collect_versions_all:
    input:
        p = "reports/version.pipeline.txt",
        s = "reports/versions.software.txt",
        r = "reports/versions.R-packages.txt",
        p2 = "reports/versions.python2.txt" if CORRECTION == "pipasic" else []
    output:
        ver = "reports/versions.txt"
    shell:
        '''
        echo "Pipeline:" > {output} &&
        cat {input.p} >> {output} &&
        echo "\n\nSoftware:" >> {output} &&
        cat {input.s} >> {output} &&
        echo "\nR packages:" >> {output} &&
        cat {input.r} >> {output} &&
        if [ -n "{input.p2}" ]; then
            echo "\nPython 2:" >> {output} &&
            cat {input.p2} >> {output};
        fi
        '''


### Pipeline Versioning ###

# get version of snakemake workflow (aka git tag + commit hash)
# NOTE: make sure to commit latest workflow to get an updated version number!
rule collect_version_pipeline_git:
    output:
        temp("reports/version.pipeline.txt")
    run:
        try:
            version = str(subprocess.check_output("cd " + WORKFLOW_PATH +
                  " && git describe --tags --long", shell=True,
                  stderr=subprocess.STDOUT, universal_newlines=True)).rstrip()
        except subprocess.CalledProcessError as e:
            version = "not available"
        with open(output[0], "w") as f_out:
            f_out.write(version)

### Software Versions ###

# dictionary of software and their version parameter call
# (incl. cutting, trimming, cleaning etc. if necessary)
dict_software_version = {
    "GNU R" : "which R >/dev/null && R --version | head -n1 | cut -d ' ' -f3",
    "MSGFPlus" : MSGFPLUS_EXEC + " | head -n1 | cut -d \" \" -f 3 | grep -Po \"(?<=\().*(?=\))\"",
    "Python3" : "which python3 >/dev/null && python3 -c \"import sys; print(sys.version.split(' ')[0])\"",
    "seqkit" : "which seqkit  >/dev/null && seqkit version | head -1 | cut -d \" \" -f2",
    "Snakemake" : "which snakemake >/dev/null && snakemake --version",
    "X!Tandem" : "which tandem.exe >/dev/null && tandem.exe -v | grep \"TANDEM\" | grep -Po \"(?<=\().*(?=\))\""
}

def software_version(call):
    # execute software version call on shell
    if not call:
        return "not available"
    try:
        return str(subprocess.check_output(call, shell=True, stderr=subprocess.STDOUT, universal_newlines=True)).rstrip()
    except subprocess.CalledProcessError as e:
        return "not available"

# get version of software
rule collect_versions_software:
    output:
        ver = temp("reports/versions.software.txt")
    run:
        with open(output.ver, 'w') as f_out:
            [f_out.write(s + ":\t" + software_version(dict_software_version[s]) + "\n")
             for s in sorted(dict_software_version.keys())]


### Python Libraries ###
list_python2_libraries = [
    "Bio",
    "numpy",
    "scipy"
]

# get versions of Python 2 libraries
rule collect_versions_python2:
    output:
        ver = temp("reports/versions.python2.txt")
    conda:
        srcdir("../environments/python2.yaml")
    params:
        libs = list_python2_libraries
    script:
        "../scripts_python/python2.versions.py"


### R-Package Versions ###

# list of all R packages which version should be reported
list_R_packages = [
    "ggplot2"
]

# get versions of R packages
rule collect_versions_R_packages:
    output:
        ver = temp("reports/versions.R-packages.txt")
    run:
        R(r'''
        packages <- strsplit("{list_R_packages}", " ")[[1]]
        lines <- vector(mode = "character", length = length(packages))
        for (i in 1:length(packages)){{
            if (require(packages[i], character.only = TRUE)){{
              version <- as.character(packageVersion(packages[i]))
            }} else {{
              version <- "not available"
            }}
            lines[i] <- paste0(packages[i], ":\t", version)
        }}
        writeLines(lines, "{output.ver}")
        ''')
