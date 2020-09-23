### X!Tandem ###

XTANDEM_OUTPUT = "iteration{i}/search/xtandem/tandem.output.xml"

# locate default X!Tandem executable
# look for xtandem (conda) first, then tandem.exe (manual installation)
def xtandem_executable():
    import shutil
    if shutil.which("xtandem"):
        return("xtandem")
    elif shutil.which("tandem.exe"):
        return ("tandem.exe")
    else:
        print("neither xtandem or tandem.exe is available", file=sys.stderr)
        sys.exit(1)

# execute X!Tandem
rule step1_xtandem_execute:
    input:
        spectra = iteration_sample,
        paras = "iteration{i}/search/xtandem/tandem.input.xml",
        taxon = "iteration{i}/search/xtandem/tandem.taxonomy.xml"
    output:
        XTANDEM_OUTPUT
    threads:
        THREADS_MAX
    benchmark:
        "logs/iteration{i}.xtandem.txt"
    log:
        "logs/iteration{i}.xtandem.log"
    params:
        exec = xtandem_executable()
    shell:
        "{params.exec} {input.paras} 2>&1 | tee {log}"

# convert X!Tandem results
rule step1_xtandem_convert:
    input:
        XTANDEM_OUTPUT
    output:
        "iteration{i}/search/xtandem/tandem.output.tsv"
    params:
        memory  = "-Xmx" + MEMORY_MAX,
        jar     = WORKFLOW_PATH + "/modules_java/xtandem-conversion/target/xtandem-conversion-1.0-SNAPSHOT-jar-with-dependencies.jar"
    shell:
        "java {params.memory} -jar {params.jar} {input} {output}"

# create taxonomy file for X!Tandem
rule step1_xtandom_input_taxonomy:
    input:
        db_main_targets = iteration_database,
        db_main_decoys = iteration_database_decoys,
        db_conts_targets = DATABASE_CONTS + ".clean",
        db_conts_decoys = DATABASE_CONTS + ".clean.decoys",
        prev_decoys = iteration_decoys
    output:
        "iteration{i}/search/xtandem/tandem.taxonomy.xml"
    run:
        lines = ["<?xml version=\"1.0\"?>",
        "<bioml label=\"x! taxon-to-file matching list\">",
        "<taxon label=\"iteration" + wildcards.i + "\">",
        "<file format=\"peptide\" URL=\"" + input.db_main_targets + "\" />",
        "<file format=\"peptide\" URL=\"" + input.db_main_decoys + "\" />",
        "<file format=\"peptide\" URL=\"" + input.db_conts_targets + "\" />",
        "<file format=\"peptide\" URL=\"" + input.db_conts_decoys + "\" />"]
        if input.prev_decoys:
            lines.append("<file format=\"peptide\" URL=\"" + input.prev_decoys + "\" />")
        lines.append("</taxon>")
        lines.append("</bioml>")
        f_out = open(output[0], "w")
        f_out.writelines([line + "\n" for line in lines])
        f_out.close()

# create input file for X!Tandem
rule step1_xtandem_input_parameters:
    input:
        spectra = iteration_sample,
        taxonomy = "iteration{i}/search/xtandem/tandem.taxonomy.xml"
    output:
        "iteration{i}/search/xtandem/tandem.input.xml"
    params:
        xtandem_output = XTANDEM_OUTPUT
    threads:
        THREADS_MAX
    run:
        f_out = open(output[0], "w")

        lines = ["<?xml version=\"1.0\"?>",
                "<bioml>"]

        lines.append(header_note("list path parameters"))
        lines.append(input_note("list path, default parameters", XTANDEM_DEFAULT))
        lines.append(input_note("list path, taxonomy information", input.taxonomy))

        lines.append(header_note("spectrum parameters"))
        lines.append(input_note("spectrum, path", input.spectra))

        lines.append(input_note("spectrum, fragment monoisotopic mass error", XTANDEM_FMME))
        lines.append(input_note("spectrum, fragment monoisotopic mass error units", XTANDEM_FMMEU))

        lines.append(input_note("spectrum, parent monoisotopic mass error plus", XTANDEM_PMMEP))
        lines.append(input_note("spectrum, parent monoisotopic mass error minus", XTANDEM_PMMEM))
        lines.append(input_note("spectrum, parent monoisotopic mass error units", XTANDEM_PMMEU))

        lines.append(header_note("spectrum conditioning parameter"))
        lines.append(input_note("spectrum, threads", threads))

        if XTANDEM_MODS_FIX or XTANDEM_MODS_VAR:
            lines.append(header_note("residue modification parameters"))
        if XTANDEM_MODS_FIX:
            lines.append(input_note("residue, modification mass", XTANDEM_MODS_FIX))
        if XTANDEM_MODS_VAR:
            lines.append(input_note("residue, potential modification mass", XTANDEM_MODS_VAR) )
        if XTANDEM_MODS_VAR_NTERM:
            lines.append(header_note("model refinement parameters"))
            lines.append(input_note("refine", "yes"))
            lines.append(input_note("refine, potential N-terminus modifications", XTANDEM_MODS_VAR_NTERM))

        lines.append(header_note("protein parameters"))
        lines.append(input_note("protein, taxon", "iteration" + wildcards.i))
        # enzyme is trypsin by default

        lines.append(header_note("scoring parameters"))
        lines.append(input_note("scoring, include reverse", "no"))

        lines.append(header_note("output parameters"))
        lines.append(input_note("output, path", params.xtandem_output))
        lines.append(input_note("output, path hashing", "no"))  # no date/time tag in output file name
        lines.append(input_note("output, message", ""))

        if XTANDEM_PARAS:
            lines.append(header_note("additional parameters"))
            lines.extend([input_note(param, value) for param, value in XTANDEM_PARAS.items()])

        lines.append("</bioml>")

        f_out.writelines([line + "\n" for line in lines])
        f_out.close()

# create input note string
def input_note(label, value):
    return "<note type=\"input\" label=\"" + str(label) + "\">" + str(value) + "</note>"

# create heading note string
def header_note(value):
    return "\n<note type=\"heading\">" + str(value) + "</note>"
