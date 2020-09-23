### Taxa Count Visualizations ###

# STEP 6:
# plot taxa counts stacked per stage (ie. original/filtered and corrected)
rule step6_taxa_vis_stacked:
    input:
        org_counts = "iteration{i}/abundance/results.original.txt",
        org_taxids = "iteration{i}/taxonomy/taxids.txt",
        red_taxids = "iteration{i}/abundance/taxids.txt",
        red_names  = "iteration{i}/abundance/taxnames.txt",
        red_counts = "iteration{i}/abundance/results.filtered.txt",
        cor_counts = RELATIVE[CORRECTION]
    output:
        plot_pdf = "iteration{i}/results/abundance_stacked.pdf",
        plot_png = "iteration{i}/results/abundance_stacked.png"
    params:
        wfpath = WORKFLOW_PATH,
        cor = CORRECTION,
        fil = FILTER,
        topL = LEGEND_MAX
    script:
        "../scripts_r/vis.stacked.R"
