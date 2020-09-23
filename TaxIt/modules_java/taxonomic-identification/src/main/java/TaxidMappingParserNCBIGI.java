import java.io.*;
import java.util.HashSet;

/**
 * Class to import and provide NCBI accession.version to taxid mappings
 */
public class TaxidMappingParserNCBIGI extends TaxidMappingParser{

    public TaxidMappingParserNCBIGI(){
        super("NCBI GI");
    }

    // import mapping from text file
    protected void parseMapping(String filename, HashSet<String> filter){

        try {
            BufferedReader br = new BufferedReader(new FileReader(new File(filename)));

            try {
                String line;
                String[] splits;

                // relevant columns
                int idx_gi = 0;
                int idx_taxid = 1;

                // parse data lines
                while ((line = br.readLine()) != null){
                    splits = line.split("\t");
                    if (filter.contains(splits[idx_gi])){
                        if (!protein2taxid.containsKey(splits[idx_gi])){
                            protein2taxid.put(splits[idx_gi], new HashSet<Integer>());
                        }
                        protein2taxid.get(splits[idx_gi]).add(Integer.parseInt(splits[idx_taxid]));
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
