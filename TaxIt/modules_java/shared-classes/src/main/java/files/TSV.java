/* Copyright (c) 2016,
 * Mathias Kuhring, KuhringM@rki.de, Robert Koch Institute, Germany,
 * All rights reserved. For details, please note the license.txt.
 */

package files;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

// handle TSV files, i.e. PSM reports from e.g. MS-GF+ or X!Tandem (after conversion)
// in particular ID parsing and taxonomic partitioning
public class TSV {

	private String filename;
	private FastaHeaderIdExtractor idExtractor;
	private HashSet<String> proteinIDs;

	// indicate input tsv (expecting NCBI GIs as proteins IDs)
	public TSV(String filename){
		this.filename = filename;
		this.proteinIDs = new HashSet<String>();
		this.idExtractor = new FastaHeaderIdExtractor("gi");
	}

	// indicate input tsv and ID type (e.g. NCBI: "gi", "accession", UniProt: "uniprot", "uniref")
	public TSV(String filename, String idType){
		this(filename);
		this.idExtractor = new FastaHeaderIdExtractor(idType);
	}

	public HashSet<String> getProteinIDs(){
		return proteinIDs;
	}

	// read tsv, extract protein hit IDs
	// shared hits must be separated on one line each!
	public void parseProteinIDs(){
		try {
			BufferedReader br = new BufferedReader(new FileReader(new File(filename)));

			try {
				// read header
				String line = br.readLine();

				// get protein position
				int proteinIndex = 0;
				for (String split : line.split("\t")){
					if (split.matches("Protein")) break;
					else proteinIndex++;
				}

				// read hit IDs
				while ((line = br.readLine()) != null){
						proteinIDs.add(idExtractor.extractProteinId(line.split("\t")[proteinIndex]));
				}
			}
			finally {
				br.close();
			}

		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (NumberFormatException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	// split tsv file into several files based on different names (e.g. scientific organism)
	// input is a map with an ID to name mapping, files are numbered in alphabetical order
	public void divideFile(HashMap<String, HashSet<Integer>> idTaxidMap, String outputFolder) {

		try {
			// open input fasta
			BufferedReader br = new BufferedReader(new FileReader(new File(this.filename)));

			// get list of organism names and sort them
			HashSet<Integer> merged = new HashSet<Integer>();
			for (HashSet<Integer> set : idTaxidMap.values()){
				merged.addAll(set);
			}
			List<Integer> taxids = new ArrayList<Integer>(merged);
			Collections.sort(taxids);

			// create name to file map and initiate a buffer for each file
			HashMap<Integer, BufferedWriter> bwPerTaxid = new HashMap<Integer, BufferedWriter>();
			int length = String.valueOf(taxids.size()).length();
			int number = 1;
			for (Integer taxid : taxids){
				// String outfile = outputFolder + "/" + String.format("%0" + length + "d", number++) + " " + name + ".tsv";
				String outfile = outputFolder + "/" + taxid + ".tsv";
				bwPerTaxid.put(taxid, new BufferedWriter(new FileWriter(new File(outfile))));
			}
			
			try {
				String line  = br.readLine(); // read header

				for (Integer taxid : taxids){
					bwPerTaxid.get(taxid).write(line); // write header
					bwPerTaxid.get(taxid).newLine();
				}
				
				// get protein position (column) in tsv file
				int proteinIndex = 0;
				for (String split : line.split("\t")){
					if (split.matches("Protein")) break;
					else proteinIndex++;
				}

				String id = "0";

				// read tsv file
				while ((line = br.readLine()) != null){
					// extract ID
					id = idExtractor.extractProteinId(line.split("\t")[proteinIndex]);

					// if ID is in idNameMap, write to corresponding name file
					if (idTaxidMap.containsKey(id)){
						for (Integer taxid : idTaxidMap.get(id)){
							bwPerTaxid.get(taxid).write(line);
							bwPerTaxid.get(taxid).newLine();
						}
					}
				}
			}
			finally {
				// close all buffer
				br.close();
				for (Integer taxid : taxids){
					bwPerTaxid.get(taxid).close();
				}
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}
