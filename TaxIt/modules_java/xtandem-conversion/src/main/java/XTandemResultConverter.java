import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.net.URISyntaxException;
import java.util.ArrayList;
import javax.xml.parsers.ParserConfigurationException;
import org.xml.sax.SAXException;

import de.proteinms.xtandemparser.xtandem.Domain;
import de.proteinms.xtandemparser.xtandem.Peptide;
import de.proteinms.xtandemparser.xtandem.Spectrum;
import de.proteinms.xtandemparser.xtandem.XTandemFile;

// for comparison, MS-GF+ tsv output:
// #SpecFile	SpecID	ScanNum	Title	FragMethod	Precursor	IsotopeError	PrecursorError(ppm)	Charge	Peptide	Protein	DeNovoScore	MSGFScore	SpecEValue	EValue	QValue	PepQValue
// cowpox_new.mgf	index=16779	26724	File3284 Spectrum20595 scans: 26724	CID	1064.976	0	-3.3240445	2	YVDGSASEDAADDTSLINSAK	gi|20178372|ref|NP_619792.1|(pre=K,post=L);gi|20178585|ref|NP_620006.1|(pre=K,post=L)	272	269	7.259532E-25	6.107548E-17	0.0	0.0

// simple programm to convert XTandem result files (XML) to TSV format
// output is either printed to the screen (standard output) or writen to file
public class XTandemResultConverter {

	public static void main(String[] args) throws SAXException, ParserConfigurationException, URISyntaxException {

		if (args.length < 1) {
			System.out.println("need more arguments: XTandem results file, output file (optional)");
		}
		else if (args.length > 2){
			System.out.println("to many arguments, only need: XTandem results file, output file (optional)");
		}
		else {
			String fileXTandemResults = args[0];

			// fetch necessary spectra data from XTandem file
			ArrayList<String> tsvLines = extractSpectra(fileXTandemResults);

			if (args.length == 1) {
				// print the data to standard output
				for (String line : tsvLines) {
					System.out.println(line);
				}
			}
			else {
				// write the data to a tsv file
				String fileConvertedResults = args[1];
				System.out.println("write spectra XTandem results (" + fileXTandemResults + ") to TSV (" + fileConvertedResults + ")");
				writeTSV(tsvLines, fileConvertedResults);
			}
		}
	}

	// extract spectra from an XTandemFile and return it as a TSV format String with header + one spectrum per line
	// data included so far: Spectrum File, Spectrum Title, Peptide, Protein, Hyperscore, EValue
	public static ArrayList<String> extractSpectra(String fileXTandemResults){

		// open XTandem file via XTandemParser
		XTandemFile xt = null;
		try {
			xt = new XTandemFile(fileXTandemResults);
		} catch (SAXException e) {
			e.printStackTrace();
			System.err.println("error: unable to read file " + fileXTandemResults);
			System.exit(1);
		} catch (ParserConfigurationException e) {
			e.printStackTrace();
			System.err.println("error: unable to read file " + fileXTandemResults);
			System.exit(1);
		}

		ArrayList<String> sb = new ArrayList<String>();

		// need at least: SpecFile, Title, Peptide, Protein, EValue (using titles similar to MS-GF+ TSV)
		sb.add("#SpecFile\tTitle\tPeptide\tProtein\tHyperscore\tEValue");

		if (xt.getSpectraList().size() > 0) {
			// iterate over spectra
			for (Spectrum spec : xt.getSpectraList()) {
				// their peptide matches
				for (Peptide pep : xt.getPeptideMap().getAllPeptides(spec.getSpectrumNumber())) {
					// and their corresponding domains (hits within the peptide?)
					for (Domain dom : pep.getDomains()) {
						sb.add(
								// XTandemFile provides the file path
								xt.getInputParameters().getSpectrumPath()
										// SupportData provide the spectrum title
										+ "\t" + xt.getSupportData(spec.getSpectrumNumber()).getFragIonSpectrumDescription().split("\\n* RTINSECONDS=")[0]
										// Domain provides ...
										// ... the matched peptide sequence
										+ "\t" + dom.getDomainSequence()
										// ... the ID of protein origin,
										// (can't split(" ") to get IDs only because this will remove decoy suffix!)
										+ "\t" + dom.getProteinKey()
										// ... the XTandem hyper score for this match
										+ "\t" + dom.getDomainHyperScore()
										// ... the eValue for this match
										+ "\t" + dom.getDomainExpect()
						);
					}
				}
			}
		}
		else{
			System.err.println("error: no spectra in " + fileXTandemResults);
			System.exit(1);
		}
		return (sb);
	}

	// write TSV lines (header + one spectrum per line) to file
	public static void writeTSV(ArrayList<String> data, String filename){
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(filename));

			try{
				for (String line : data){
					bw.write(line);
					bw.newLine();
				}
			}
			finally{
				bw.close();
			}
		} catch (IOException e) {
			e.printStackTrace();
			System.err.println("error: unable to write file " + filename);
			System.exit(1);
		}
	}
}