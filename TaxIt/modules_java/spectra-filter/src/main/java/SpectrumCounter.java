import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

public class SpectrumCounter {

	private File mgfFile;
	
	public SpectrumCounter(File mgfFile){
		this.mgfFile = mgfFile;
	}
	
	public void countSpectra(){

		try {
			BufferedReader mgf_input = new BufferedReader(new FileReader(mgfFile));
			
			try {
				String line;
				int speccount = 0;
				
				while ((line = mgf_input.readLine()) != null){
					
					if (line.startsWith("BEGIN IONS")){
						speccount++;
					}
				}
				
				System.out.println("counted " + speccount + " spectra");
			}
			finally {
				mgf_input.close();
			}
		} 
		catch (FileNotFoundException e) {
			e.printStackTrace();
		} 
		catch (IOException e) {
			e.printStackTrace();
		}
	}
	
}
