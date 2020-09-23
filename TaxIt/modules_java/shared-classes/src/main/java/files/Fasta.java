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

// handle fasta files, in particular taxonomic partioning
public class Fasta {

	private String filename;
	private FastaHeaderIdExtractor idExtractor;

	// indicate input fasta (expecting NCBI GIs as proteins IDs)
	public Fasta(String filename){
		this.filename = filename;
		this.idExtractor = new FastaHeaderIdExtractor("gi");
	}

	// indicate input fasta and ID type (e.g. NCBI: "gi", "accession", UniProt: "uniprot", "uniref")
	public Fasta(String filename, String idType){
		this(filename);
		this.idExtractor = new FastaHeaderIdExtractor(idType);
	}
	
	// split fasta file into several files based on different names (e.g. scientific organism)
	// input is a map with an ID to name mapping, files are numbered in alphabetical order
	public void divideFile(HashMap<String, HashSet<Integer>> idTaxidMap, String outputFolder) {

		try {
			// open input fasta
			BufferedReader fastaBR = new BufferedReader(new FileReader(new File(this.filename)));

			// get list of organism names and sort them
			HashSet<Integer> merged = new HashSet<Integer>();
			for (HashSet<Integer> set : idTaxidMap.values()){
				merged.addAll(set);
			}
			List<Integer> taxids = new ArrayList<Integer>(merged);
			Collections.sort(taxids);

			// create taxid-to-file map and initiate a buffer for each file
			HashMap<Integer, BufferedWriter> bwPerTaxid = new HashMap<Integer, BufferedWriter>();
			int length = String.valueOf(taxids.size()).length();
			int number = 1;
			for (Integer taxid : taxids){
				// String outfile = outputFolder + "/" + String.format("%0" + length + "d", number++) + " " + name + ".fasta";
				String outfile = outputFolder + "/" + taxid + ".fasta";
				bwPerTaxid.put(taxid, new BufferedWriter(new FileWriter(new File(outfile))));
			}

			try {
				String line;
				String id = "0";

				// read fasta file
				while ((line = fastaBR.readLine()) != null){
					// if fasta header
					if (line.startsWith(">")){
						// extract ID
						id = idExtractor.extractProteinId(line);

						// if ID is in idTaxidMap, write to corresponding taxid files
						if (idTaxidMap.containsKey(id)){
							for (Integer taxid : idTaxidMap.get(id)){
								bwPerTaxid.get(taxid).write(line);
								bwPerTaxid.get(taxid).newLine();
							}
						}
					}
					// else fasta body (i.e. peptide sequence)
					else{
						// if ID is in idNameMap, write to corresponding name file
						if (idTaxidMap.containsKey(id)){
							for (Integer taxid : idTaxidMap.get(id)){
								bwPerTaxid.get(taxid).write(line);
								bwPerTaxid.get(taxid).newLine();
							}
						}
					}
				}
			}
			finally {
				// close all buffer
				fastaBR.close();
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
