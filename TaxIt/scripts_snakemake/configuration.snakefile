### Configuration ###

# Global thread config for rules
THREADS_MAX = config.get("threads", 4)

# Global Java memory maximum
MEMORY_MAX = config.get("memory", "4G")

# Search iterations
ITERATIONS = config.get("iterations", 2)

# Sample
SAMPLE = config["sample"]
MODE = config.get("mode", "single")

# Databases (initial)
DATABASE_REF = config["database_reference"]
DATABASE_CONTS = config.get("database_contaminants", WORKFLOW_PATH + "/resources/crap.fasta")
DATABASE_HOST = config.get("database_host", "")

# Taxonomy
TAX_MAPPING_NCBI = config["taxonomy_mapping_ncbi"]
TAX_MAPPING = config["taxonomy_mapping"]
TAX_NODES = config["taxonomy_nodes"]
TAX_NAMES = config["taxonomy_names"]

# Decoy passing
PASS_DECOYS = config.get("pass_decoys", True)

# Search engine
SEARCH_ENGINE = config.get("search_engine", "xtandem")
SEARCH_RESULTS = {
    "msgfplus" : "msgfplus/msgfplus.output.tsv",
    "xtandem"  : "xtandem/tandem.output.tsv"
}

FDR_CUT= config.get("fdr_cutoff", 0.01)

# Filter type: none (default), ratioR (e.g. ratio0.005), topN (e.g. top10), minN (e.g. min2), maxR (e.g. max0.005)
FILTER = config.get("filter", ["none"])

# Correction type: weighted (default), pipasic, uniques
CORRECTION = config.get("correction", "weighted")

# Visualization
# max number of taxa labels in stacked bar plot legend
LEGEND_MAX = config.get("stacked_legend_max", 15)

# MS-GF+ parameter
# indicate jar file only if MS-GF+ is not installed via conda (e.g. "~/software/MSGFPlus/MSGFPlus.jar")
MSGFPLUS_JAR = config.get("msgfplus_jar", "")
MSGFPLUS_EXEC = "java -Xmx" + MEMORY_MAX + " -jar " + MSGFPLUS_JAR if MSGFPLUS_JAR else "msgf_plus -Xmx" + MEMORY_MAX
MSGFPLUS_MODS = config.get("msgfplus_mods", "")

MSGFPLUS_PMT  = config.get("msgfplus_PrecursorMassTolerance", "10ppm")
MSGFPLUS_INST = config.get("msgfplus_MS2DetectorID", 1)
MSGFPLUS_MINC = config.get("msgfplus_MinCharge", 2)
MSGFPLUS_MAXC = config.get("msgfplus_MaxCharge", 4)
# additional parameters (string)
MSGFPLUS_PARAS = config.get("msgfplus_add_params", "")

# X!Tandem parameter
XTANDEM_DEFAULT = config.get("xtandem_default", WORKFLOW_PATH + "/resources/xtandem_default_input.xml")
XTANDEM_FMME = config.get("xtandem_fmme", 0.4)
XTANDEM_FMMEU = config.get("xtandem_fmmeu", "DA")
XTANDEM_PMMEP = config.get("xtandem_pmmep", 100)
XTANDEM_PMMEM = config.get("xtandem_pmmem", 100)
XTANDEM_PMMEU = config.get("xtandem_pmmeu", "ppm")
# X!Tandem PTMs (comma separated if more than one)
XTANDEM_MODS_FIX = config.get("xtandem_mods_fixed", "57@C") # default:
XTANDEM_MODS_VAR = config.get("xtandem_mods_variable", "") # e.g. "16@M"
XTANDEM_MODS_VAR_NTERM = config.get("xtandem_mods_variable_nterm", "") # e.g. "+42.0@[" for N-terminal acetylation
# additional parameters (dict: param name -> value)
XTANDEM_PARAS = config.get("xtandem_add_params", dict())

# optional proxy for protein downloads
HTTPS_PROXY_URL = config.get("https_proxy_url", "")
if HTTPS_PROXY_URL:
    HTTPS_PROXY_PORT = config["https_proxy_port"]
