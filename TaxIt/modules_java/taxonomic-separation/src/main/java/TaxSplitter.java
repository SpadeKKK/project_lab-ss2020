/* Copyright (c) 2016,
 * Mathias Kuhring, KuhringM@rki.de, Robert Koch Institute, Germany,
 * All rights reserved. For details, please note the license.txt.
 */

import files.Fasta;
import files.TSV;
import tools.MemoryMonitor;
import tools.Timer;

import java.io.File;
import java.util.HashMap;
import java.util.HashSet;

// parse PSMs (TSV) from a spectral search against several organism,
// calculate taxonomic grouping of protein IDs per organism
//
public class TaxSplitter {

	public static void main(String[] args){

		MemoryMonitor memory = new MemoryMonitor();

		Timer totalTime = new Timer();
		totalTime.start();
		Timer timer = new Timer();

		// parameter
		if (args.length != 8){
			System.out.println("wrong number of parameters!");
			System.out.println("iteration, msgf_results, fasta_file, id_taxid_ref_file, id_taxid_file, tax_node_file, tax_name_file, outputFolder");
			System.exit(1);
		}

		int index = 0;
		// files
		int iteration = Integer.parseInt(args[index++]);
		String msgf_results = args[index++];
		String fasta_file = args[index++];
		String id_taxid_ref_file = args[index++];
		String id_taxid_file = args[index++];
		String tax_node_file = args[index++];
		String tax_name_file = args[index++];
		String outputFolder = args[index++];

		// create output folder (mkdir checks itself whether it already exist)
		new File(outputFolder).mkdir();

		// Import PSM/hit Ids
		TSV psms = null;
		if (id_taxid_ref_file.endsWith("gi_taxid_prot.dmp")){
			psms = new TSV(msgf_results);
		}
		else if (id_taxid_ref_file.endsWith("prot.accession2taxid")){
			psms = new TSV(msgf_results, "accession");
		}
		else if (id_taxid_ref_file.endsWith("idmapping_selected.tab")){
			psms = new TSV(msgf_results, "uniref");
		}
		else{
			System.err.println("error: unexpected mapping file:\n\t" + id_taxid_ref_file);
			System.exit(1);
		}

		psms.parseProteinIDs();
		HashSet<String> matches = psms.getProteinIDs();

		// Import Taxonomy: IDMapping, Nodes, Names
		timer.start();
		Taxonomy tax = new Taxonomy(id_taxid_file, tax_node_file, tax_name_file);
		timer.stop(); System.out.println();

		// Remove ids with no taxid/organism mapping
		// TODO: report removed ids
		timer.start();
		tax.cleanIdMap(matches);
		timer.stop(); System.out.println();

		// Down propagation of redundant proteins and
		// species assignment of strain level hits (1. iteration)
		// TODO: probably not useful for UniRef cluster member hits
		timer.start();
		tax.reassignIds(iteration);
		timer.stop(); System.out.println();

		// Identify Orgs, use identified Ids to create a taxid Set
		HashSet<Integer> taxids = tax.getTaxids(matches);
		System.out.println("\n" + taxids);

		// select all IDs per taxid and create ID to taxid map
		// (multiple taxids per ID are possible, e.g. in case of UniRef hits)
		HashMap<String, HashSet<Integer>> idTaxidMap = tax.getIdTaxidMap(taxids);

		// Iterate over input Fasta
		timer.start();
		System.out.print("iterate input fasta and tsv... ");
		
		Fasta fasta = null;
		if (id_taxid_ref_file.endsWith("gi_taxid_prot.dmp")){
			fasta = new Fasta(fasta_file);
		}
		else if (id_taxid_ref_file.endsWith("prot.accession2taxid")){
			fasta = new Fasta(fasta_file, "accession");
		}
		else if (id_taxid_ref_file.endsWith("idmapping_selected.tab")){
			fasta = new Fasta(fasta_file, "uniref");
		}
		else{
			System.err.println("error: unexpected mapping file:\n\t" + id_taxid_ref_file);
			System.exit(1);
		}

		fasta.divideFile(idTaxidMap, outputFolder);
		psms.divideFile(idTaxidMap, outputFolder);
		timer.stop(); System.out.println();

		System.out.print("total runtime: ");
		totalTime.stop(); System.out.println();
	}
}
