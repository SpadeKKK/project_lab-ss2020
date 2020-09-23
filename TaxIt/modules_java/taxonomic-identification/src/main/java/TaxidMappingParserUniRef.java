import java.io.*;
import java.util.HashSet;

/**
 * Class to import and provide UniRef to taxid mappings
 *
 * Based on the UniProt ID Mapping file idmapping_selected.tab
 * (ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/README)
 *
 * Relevant fields: 8. UniRef100, 9. UniRef90, 10. UniRef50, 13. NCBI-taxon
 *
 */
public class TaxidMappingParserUniRef extends TaxidMappingParser{

    public TaxidMappingParserUniRef(){
        super("UniRef");
    }

    // import mapping from text file
    protected void parseMapping(String filename, HashSet<String> filter){

        // identify UniRef version of input IDs
        String uniref = filter.iterator().next().split("_")[0];
        System.out.print("identified UniRef database from PSMs: " + uniref);

        try {
            BufferedReader br = new BufferedReader(new FileReader(new File(filename)));

            try {
                String line;
                String[] splits;

                // find relevant columns
                int idx_uniref = 0;
                switch (uniref){
                    case "UniRef50":
                        idx_uniref = 9; // 10-1
                        break;
                    case "UniRef90":
                        idx_uniref = 8; // 9-1
                        break;
                    case "UniRef100":
                        idx_uniref = 7; // 8-1
                        break;
                    default:
                        System.err.println("error: unrecognized UniRef database: " + uniref);
                        System.exit(1);
                }
                int idx_taxid = 12; // 13-1, NCBI-taxon

                // parse data lines
                while ((line = br.readLine()) != null){
                    splits = line.split("\t");
                    if (filter.contains(splits[idx_uniref])){
                        if (!protein2taxid.containsKey(splits[idx_uniref])){
                            protein2taxid.put(splits[idx_uniref], new HashSet<Integer>());
                        }
                        protein2taxid.get(splits[idx_uniref]).add(Integer.parseInt(splits[idx_taxid]));
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
