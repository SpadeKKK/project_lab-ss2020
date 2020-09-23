import java.io.*;
import java.util.HashSet;

/**
 * Class to import and provide NCBI accession.version to taxid mappings
 */
public class TaxidMappingParserNCBIAcc extends TaxidMappingParser{

    public TaxidMappingParserNCBIAcc(){
        super("NCBI accession");
    }

    // import mapping from text file
    protected void parseMapping(String filename, HashSet<String> filter){

        try {
            BufferedReader br = new BufferedReader(new FileReader(new File(filename)));

            try {
                String line;
                String[] splits;

                // header line
                line = br.readLine();
                splits = line.split("\t");

                // find relevant columns
                int idx_accver = -1;
                int idx_taxid = -1;
                for (int i=0; i<splits.length; i++){
                    if (splits[i].matches("accession.version")){
                        idx_accver = i;
                    }
                    if (splits[i].matches("taxid")){
                        idx_taxid = i;
                    }
                }

                // parse data lines
                while ((line = br.readLine()) != null){
                    splits = line.split("\t");
                    if (filter.contains(splits[idx_accver])){
                        if (!protein2taxid.containsKey(splits[idx_accver])){
                            protein2taxid.put(splits[idx_accver], new HashSet<Integer>());
                        }
                        protein2taxid.get(splits[idx_accver]).add(Integer.parseInt(splits[idx_taxid]));
                    }
                }
            }
            finally {
                br.close();
            }
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
