import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

import org.jsoup.Jsoup;
import org.jsoup.nodes.Document;
import org.jsoup.nodes.Element;
import org.jsoup.select.Elements;

// Access to specific NCBI databases via Entrez, includes so far:
// NCBI Protein: non-RefSeq proteins per TAXID, GI to TAXID
public class Entrez {

	// NBCI Entrez eutils base urls
	private static final String URL_ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?";
	private static final String URL_EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?";
	private static final String URL_ELINK =  "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?";

	// maximal number of entries per fetch iteration
	// some data has to be downloaded in chunks due to NCBI server limits/policy
	private static final int ITER_MAX = 200;

	// maximal number of IDs per fetch iteration
	// to many ID requests result in a time out
	// (NCBI error message: "Search Backend failed: read request has timed out")
	private static final int ID_MAX = 10000000;

	// so far the main functions downloads NCBI proteins for a given taxid
	public static void main(String[] args) throws IOException {

		// simple parameter parsing: a list of taxids, a output folder for fasta results,
		// a NCBI database (refseq, non-refseq, all) and expansion control (optional, default=exp)
		String[] taxids = args[0].split(",");
		String folder = args[1];
		String database = args[2];
		boolean expansion = !(args.length >= 4 && args[3].matches("noexp"));

		// download non-refseq proteins
		Entrez entrez = new Entrez();

		for (String taxid : taxids){
			entrez.downloadProteins(Integer.parseInt(taxid), folder, database, expansion);
		}
	}

	// Downloads an Entrez XML via the provided query url and returns it as Jsoup Document
	private Document xmlDownload(String entrezQueryUrl){
		// download file
		String data = Download.downloadFile(entrezQueryUrl);
		// create Jsoup Document from
		Document doc = Jsoup.parse(data);
		return doc;
	}

	// query number of NCBI proteins for a given TAXID and DB (refseq, non-refseq, all)
	private int getNumberOfProteins(int taxid, String dbString, String expString){
		// create query url
		String url	= URL_ESEARCH + "db=protein&" + "retmax=0&"
				+ "term=txid" + taxid + expString + dbString;
		// download search results
		Document doc = xmlDownload(url);
		// return count
		return Integer.parseInt(doc.select("eSearchResult > Count").first().text());
	}

	// query all GIs of NCBI proteins for a given TAXID and DB (refseq, non-refseq, all)
	private ArrayList<Integer> getIDsOfProteins(int taxid, String dbString, String expString){
		// get protein count
		int count = getNumberOfProteins(taxid, dbString, expString);

		// id list
		ArrayList<Integer> ids = new ArrayList<Integer>();

		// iterative query
		int retmax = ID_MAX;
		for (int retstart = 0; retstart < count; retstart += retmax){
			// create query url with retmax=count, thus all IDs will be returned
			String url	= URL_ESEARCH + "db=protein&" + "retmax=" + retmax + "&"
					+ "retstart=" + retstart + "&" + "term=txid" + taxid + expString + dbString;
			System.out.println(url);
			// download search results
			Document doc = xmlDownload(url);

			// check for errors
			Elements errors = doc.select("eSearchResult > ERROR");
			if (errors.size() > 0){
				for (Element error : errors){
					System.err.println("error: " + error.text());
				}
				System.exit(1);
			}

			// extract ids
			Elements elems = doc.select("eSearchResult > IdList > Id");
			for (Element elem : elems){
				ids.add(Integer.parseInt(elem.text()));
			}
		}

		return ids;
	}

	// download all NCBI proteins for a given TAXID and DB (refseq, non-refseq, all) to a fasta file
	public void downloadProteins(int taxid, String folder, String database, boolean expansion){
		// expansion, output file and info
		String expString, filename;
		if (expansion) {
			expString = "[Organism:exp]";
			filename = folder + "/" + taxid + "_exp_" + database + ".fasta";
			System.out.println("Searching " + database + " proteins for txid " + taxid + " (and subtaxa)");
		}
		else {
			expString = "[Organism:noexp]";
			filename = folder + "/" + taxid + "_noexp_" + database + ".fasta";
			System.out.println("Searching " + database + " proteins for txid " + taxid + " (exclusively)");
		}
		System.out.println("\t(" + new Date().toString() + ")");

		// db string
		String dbString = "";
		switch (database) {
			case "refseq":
				dbString = "+AND+refseq[filter]";
				break;
			case "non-refseq":
				dbString = "+NOT+refseq[filter]";
				break;
			case "all":
				break;
			default:
				System.err.println("error: db selection \"" + database + "\" unavailable (should be either \"refseq\", \"non-refseq\" or \"all\")");
				System.exit(1);
		}

		// get ids of all non-refseq proteins for the taxid organism
		ArrayList<Integer> ids = getIDsOfProteins(taxid, dbString, expString);

		// proteins are downloaded in chunks of ITER_MAX due to NCBI server limits/policy
		// iteration parameter
		int iterations = (int) Math.ceil((double) ids.size() / (double) ITER_MAX);
		int rest = ids.size() % ITER_MAX;

		System.out.println("Found " + ids.size() + " proteins. Need " + iterations + " iterations");

		String fasta = null;
		int ret_count, ret_start, ret_end;
		double progress;
		String url;

		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(filename));
			try{
				for (int i=0; i<iterations; i++){
					// determine chunk
					ret_start = i * ITER_MAX;
					ret_end = Math.min(ret_start + ITER_MAX, ids.size());
					System.out.println("Iterarion: " + (i+1) + " - " + "id index: " + ret_start);

					// create query url
					url = URL_EFETCH + "db=protein&" + "rettype=fasta&" + "retmode=text&" //+ "showgi=1&"
							+ "id=" + idString(ids.subList(ret_start, ret_end));
					System.out.println(url);

					// repeat current download if expected protein number is not matched
					ret_count = 0;
					while ((i<iterations-1 && ret_count!=ITER_MAX) || 	// normal iteration
							(i==iterations-1 && ret_count!=rest)){		// last iteration

						// download fasta
						fasta = Download.downloadFile(url);
						// count fasta entries aka proteins
						ret_count = countFastaHeaders(fasta);
					}

					// add proteins to fasta file
					bw.write(fasta);

					progress = Math.round((double)(ret_end)/(double) ids.size() * 10000.) / 100.;
					System.out.println("\tdownloaded " + ret_count + " sequences"
							+ " (" + progress + " %)");
				}
			}
			finally{
				bw.close();
				System.out.println("Done");
			}
		} catch (IOException e) {
			e.printStackTrace();
			System.err.println("error: unable to write file " + filename);
			System.exit(1);
		}
	}

	// count the number of fasta headers in a String list (one String per fasta line)
	private int countFastaHeaders(String fasta){
		int counter = 0;
		for (String line : fasta.split("\n")){
			if (line.startsWith(">")){
				counter++;
			}
		}
		return counter;
	}

	// create proper id string for efetch from a list of ids
	private String idString(List<Integer> list){
		// combine Integers into String via toString(), e.g.: "[1234, 3214, 5325]"
		String idString = list.toString();
		// remove brackets and spaces
		idString = idString.substring(1, idString.length()-1).replace(" ", "");
		return idString;
	}

	// get the corresponding TAXID for a given ACCESSION.VERSION
	public int accToTaxid(String accession){
		// get protein entry as TinySeq XML
		String url = URL_EFETCH + "db=protein&rettype=fasta&retmode=xml&id=" + accession;

		Document doc = xmlDownload(url);

		String taxid = doc.select("TSeqSet > TSeq > TSeq_taxid").first().text();
		return Integer.parseInt(taxid);
	}

	// get Map of corresponding TAXIDs for a Set of Accessions (NCBI Protein)
	public HashMap<String,Integer> accsToTaxids(HashSet<String> accessions){
		HashMap<String,Integer> accToTaxid = new HashMap<String,Integer>();

		ArrayList<String> acc = new ArrayList<String>(accessions);

		// IDs are queried in chunks of ITER_MAX due to NCBI server limits/policy
		// iteration parameter
		int iterations = (int) Math.ceil((double) acc.size() / (double) ITER_MAX);
		int rest = acc.size() % ITER_MAX;

		System.out.println("Got " + acc.size() + " accessions. Need " + iterations + " iterations");

		int ret_count, ret_start, ret_end;
		double progress;
		String url;

		for (int i=0; i<iterations; i++){
			// determine chunk
			ret_start = i * ITER_MAX;
			ret_end = Math.min(ret_start + ITER_MAX, acc.size());
			System.out.println("Iterarion: " + (i+1) + " - " + "acc index: " + ret_start);

			// create query url
			url = URL_EFETCH + "db=protein&" + "rettype=fasta&" + "retmode=xml&"
					+ "id=" + String.join(",", acc.subList(ret_start, ret_end));

			// repeat current download if expected protein number is not matched
			ret_count = 0;
			while ((i<iterations-1 && ret_count!=ITER_MAX) || 	// normal iteration
					(i==iterations-1 && ret_count!=rest)){		// last iteration

				// download search results
				Document doc = xmlDownload(url);

				// get proteins
				Elements elems = doc.select("TSeqSet > TSeq");

				// count proteins
				ret_count = elems.size();

				// iterate proteins, if expected number was downloaded
				if ((i<iterations-1 && ret_count==ITER_MAX) || 	// normal iteration
						(i==iterations-1 && ret_count==rest)) {    // last iteration
					for (Element elem : elems) {
						// extract taxids
						accToTaxid.put(elem.select("TSeq_accver").first().text(),
								Integer.parseInt(elem.select("TSeq_taxid").first().text()));
					}
				}
			}

			progress = Math.round((double)(ret_end)/(double) acc.size() * 10000.) / 100.;
			System.out.println("\tdownloaded " + ret_count + " sequences"
					+ " (" + progress + " %)");
		}

		return accToTaxid;
	}

	// get the corresponding TAXID for a given GI via elink
	public int giToTaxid(int gi){
		String url = URL_ELINK + "dbfrom=protein&db=taxonomy&id=" + gi;
		Document doc = xmlDownload(url);

		// selection string should actually be: "eLinkResult > LinkSet > LinkSetDb > Link > Id"
		// however, the xml file is not parsed properly by JSoup ("Link" element is not closed)
		// (don't know whether a JSoup bug or elink xml file not formatted properly in first place)
		String taxid = doc.select("eLinkResult > LinkSet > LinkSetDb > Id").first().text();
		return Integer.parseInt(taxid);
	}

	// get Map of corresponding TAXIDs for a Set of GIs
	public HashMap<Integer,Integer> gisToTaxids(HashSet<Integer> gis){
		HashMap<Integer,Integer> giToTaxid = new HashMap<Integer,Integer>();
		for (int gi : gis){
			giToTaxid.put(gi, giToTaxid(gi));
		}
		return giToTaxid;
	}
}
