import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashSet;
import java.util.Vector;

public class Main {

	private static boolean verbose = false;
	private static boolean countonly = false;
	private static boolean replaceTitle = false;
	private static boolean filterByTitle = false;
	private static boolean keep = false;
	private static String filename, titlesFile;
	private static HashSet<String> titles;

	/**
	 * @param args
	 */
	public static void main(String[] args) {

		printHeader();

		parseParas(args);

		String input = filename;
		String output = new StringBuffer(filename).insert(filename.length()-4, "_new").toString();
		String rejected = new StringBuffer(filename).insert(filename.length()-4, "_rej").toString();

		if (countonly){
			new SpectrumCounter(new File(filename)).countSpectra();
		}
		else{
			if (filterByTitle){
				readTitles();
			}
			readMGF(input, output, rejected);			
		}

		System.out.println("done");
	}

	

	private static void printHeader(){
		String title = "mgfres - Rejection of Empty Spectra (MS/MS) in MGF files\n";
		System.out.println(title);
	}

	private static void printHelp(){
		String help = 	
				"- rejects empty spectra as well spectra with missing TITLE, CHARGE or PEPMASS\n" +
				"- converts integer values to double (concerns ions, PEPMASS, RTINSECONDS)\n" +
				"- checks CHARGE suffixes (+/-) and attaches \"+\" if no suffix is present\n" +
				"- optional: renames TITLES with \"spectrum\" + increasing number\n" +
				"- optional: keep/reject spectra based on file with titles\n" +
				"\n" + "example call:" + "\t" + "mgfref -i myfile.mgf" + "\n" + 
				"example output:" + "\t" + "myfile_new.mgf and myfile_rej.mgf" + "\n\n" + 
				"mgfres accepts following parameters:" + "\n\n" +
				"\t" + "-h" + "\t" + "prints this help message" + "\n" +
				"\t" + "-c" + "\t" + "simply counts the spectra (count mode)" + "\n" +
				"\t" + "-v" + "\t" + "itemizes each spectrum (verbose mode)" + "\n" +
				"\t" + "-t" + "\t" + "replace titles with increasing numbers (e.g. spectrum#)" + "\n" +
				"\t" + "-i" + "\t" + "indicate the mgf file to process" + "\n" + 
				"\t" + "-f" + "\t" + "indicate a file with spectra (titles) to reject (default)" + "\n" +
				"\t" + "-k" + "\t" + "keep indicated spectra, reject the rest";
		System.out.println(help);
		System.exit(0);
	}

	private static void parseParas(String[] args){
		if (args.length==0){
			printHelp();
		}
		else{
			for (int i=0; i<args.length; i++){
				if (args[i].matches("-h")){
					printHelp();
				}
				else if (args[i].matches("-v")){
					verbose = true;
				}
				else if (args[i].matches("-c")){
					countonly = true;
				}
				else if (args[i].matches("-t")){
					replaceTitle = true;
				}
				else if (args[i].matches("-i")){
					filename = args[i+1];
					i++;
				}
				else if (args[i].matches("-f")){
					titlesFile = args[i+1];
					i++;
					filterByTitle = true;
				}
				else if (args[i].matches("-k")){
					keep = true;
				}
				else{
					System.out.println("unknown parameter");
					System.exit(0);
				}
			}
		}
	}

	private static void readTitles() {
		titles = new HashSet<String>();
		
		try {
			BufferedReader tf = new BufferedReader(new FileReader(new File(titlesFile)));	
			try {
				String line;
				while ((line = tf.readLine()) != null){
					titles.add(line);
				}
			}
			finally {
				tf.close();
			}
		}
		catch (FileNotFoundException e) {
			e.printStackTrace();
		} 
		catch (IOException e) {
			e.printStackTrace();
		}
				
	}
	
	private static void readMGF(String input, String output, String rejected){

		try {
			BufferedReader mgf_input = new BufferedReader(new FileReader(new File(input)));				
			BufferedWriter mgf_output = new BufferedWriter(new FileWriter(new File(output)));
			BufferedWriter mgf_rejected = new BufferedWriter(new FileWriter(new File(rejected)));

			try {
				String line;
				int speccount = 0;
				int rejcount = 0;
				int ioncount;
				Vector<String> speclines;
				
				boolean hasTITLE, hasPEPMASS, hasCHARGE, keepSpectrum;
				String title;

				while ((line = mgf_input.readLine()) != null){

					if (line.startsWith("BEGIN IONS")){
						speccount++;

						// read spectrum (complete)
						speclines = new Vector<String>();
						ioncount = 0;

						hasTITLE = false;
						hasPEPMASS = false;
						hasCHARGE = false;
						keepSpectrum = !filterByTitle;
						
						speclines.add(line);	// "BEGIN IONS"

						while (!(line = mgf_input.readLine()).startsWith("END IONS")){
							if (!(line.isEmpty() || line.contains("="))){
								ioncount++;
								line = assureDoubles(line);
							}
							if (line.startsWith("TITLE")){
								hasTITLE = true;
								
								if (filterByTitle){
									title = line.replaceFirst("^TITLE=", "");
									if (keep) keepSpectrum = titles.contains(title);
									else keepSpectrum = !titles.contains(title);
								}
								
								if (replaceTitle){
									line = "TITLE=spectrum" + speccount;
								}
							}
							if (line.startsWith("PEPMASS")){
								hasPEPMASS = true;
								line = correctValues(line);
							}
							if (line.startsWith("CHARGE")){
								hasCHARGE = true;
								line = correctCharge(line);
							}
							if (line.startsWith("RTINSECONDS")){
								line = correctValues(line);
							}
							speclines.add(line);
						}
						speclines.add(line);	// "END IONS"

						// write spectrum
						if (ioncount > 0 && hasTITLE && hasPEPMASS && hasCHARGE && keepSpectrum){
							if (verbose) System.out.println("accepting spectrum " + speccount + " (" + ioncount + " ions; " + 
									hasTITLE + " " + hasPEPMASS + " " + hasCHARGE + ")");

							for (String sl : speclines){
								mgf_output.write(sl);
								mgf_output.newLine();
							}
						}
						else{
							if (verbose) System.out.println("rejecting spectrum " + speccount + " (" + ioncount + " ions; " + 
									hasTITLE + " " + hasPEPMASS + " " + hasCHARGE + ")");
							rejcount++;

							for (String sl : speclines){
								mgf_rejected.write(sl);
								mgf_rejected.newLine();
							}
							mgf_rejected.newLine();
						}
					}
					else{
						mgf_output.write(line);
						mgf_output.newLine();
					}
				}

				System.out.println("rejected " + rejcount + " of " + speccount + " spectra");
			}
			finally {
				mgf_input.close();
				mgf_output.close();
				mgf_rejected.close();
			}
		} 
		catch (FileNotFoundException e) {
			e.printStackTrace();
		} 
		catch (IOException e) {
			e.printStackTrace();
		}
	}

	// make sure the CHARGE value ends with "+" or "-", if not add "+"
	private static String correctCharge(String line) {
		if (!(line.endsWith("+") || line.endsWith("-"))){
			line = line + "+";
		}
		return line;
	}

	// make sure parameter (e.g. PEPMASS) have double values
	private static String correctValues(String line) {
		String[] keyValueSplits = line.split("=");
		return keyValueSplits[0] + "=" + assureDoubles(keyValueSplits[1]);
	}

	// make sure ions have double values
	private static String assureDoubles(String line) {
		String[] splits = line.split(" ");
		StringBuffer newLine = new StringBuffer();
		double value;
		for (String split : splits){
			value = Double.parseDouble(split);
			if (value == (long) value){ 
				newLine.append((long) value + ".0 ");
			}
			else{
				newLine.append(split + " ");				
			}
		}
		return newLine.substring(0, newLine.length()-1);
	}
}
