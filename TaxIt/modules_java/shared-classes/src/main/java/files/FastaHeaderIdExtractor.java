/* Copyright (c) 2016,
 * Mathias Kuhring, KuhringM@rki.de, Robert Koch Institute, Germany,
 * All rights reserved. For details, please note the license.txt.
 */

package files;

// extractor for different ID/accession types in fasta files
// support so far (header examples, extracted id -->ID<--):
// NCBI GI	>gi|-->302596923<--|ref|YP_003848705.1| glycoprotein [Rift Valley fever virus]
// NCBI Acc	>-->YP_009333364.1<-- hypothetical protein 9 [Wuhan heteroptera virus 3]
// UniProt	>sp|-->P62258<--|1433E_HUMAN 14-3-3 protein epsilon OS=Homo sapiens GN=YWHAE PE=1 SV=1
// UniRef	>-->UniRef50_Q6GZX3<-- Uncharacterized protein 002L n=19 Tax=Ranavirus RepID=002L_FRG3G
public class FastaHeaderIdExtractor {

	private String idType;

	// initiate by indicating the ID type to extract
	public FastaHeaderIdExtractor(String idType){
		this.idType = idType;
	}

	// provide a fasta header and get the ID
	public String extractProteinId(String proteinHeader){

		// remove ">" if still present
		if (proteinHeader.startsWith(">")){
			proteinHeader = proteinHeader.substring(1);
		}

		// extract ID depending on type
		switch (idType){
			case "gi":
				return proteinHeader.split("\\|")[1];
			case "accession":
				// NCBI protein fasta header seem to be rather inconsistent at the moment
				// some will just contain the plain accession as prefix
				// but some are more detailed and include e.g. source database
				// (separated with "|", similar to old GI format: e.g. gb|AOR81611.1|, sp|P03518.2|GP_RVFV)
				String id = proteinHeader.split(" ")[0];
				if (id.contains("|")) return id.split("\\|")[1];
				else return id;
			case "uniprot":
				return proteinHeader.split("\\|")[1];
			case "uniref":
				return proteinHeader.split(" ")[0];
			default:
				return null;
		}
	}
}
