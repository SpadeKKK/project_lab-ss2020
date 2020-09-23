
import java.io.*;
import java.util.HashMap;
import java.util.HashSet;

/**
 * Class for importing NCBI taxid to name mappings
 */
public class Names {

    // taxid to names
    private HashMap<Integer,String> taxidToName;

    public Names(){
        taxidToName = new HashMap<Integer,String>();
    }

    /**
     *
     * @param taxid The taxid of the organism of interest
     * @return Scientific name of organism of interest
     */
    public String getName(int taxid){
        return taxidToName.get(taxid);
    }

    /**
     * Import taxid name mappings, keep scientific names for indicated taxids only
     *
     * @param filename Name of the taxid name mapping file
     * @param filterTaxids Set of taxid for which names should be imported
     */
    public void parseNames(String filename, HashSet<Integer> filterTaxids){

        try {
            BufferedReader br = new BufferedReader(new FileReader(new File(filename)));

            try {
                String line;
                String[] splits;
                Integer taxid;
                String name, name_class;
                while ((line = br.readLine()) != null){
                    splits = line.split("\t\\|\t");
                    taxid = Integer.parseInt(splits[0]);
                    name_class = splits[3];
                    if (filterTaxids.contains(taxid) && name_class.startsWith("scientific name")){
                        name = splits[1].replace('/', ' ');
                        taxidToName.put(taxid, name);
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

    /**
     * Import taxid name mappings, keep all scientific names
     *
     * @param filename Name of the taxid name mapping file
     */
    public void parseNames(String filename){

        try {
            BufferedReader br = new BufferedReader(new FileReader(new File(filename)));

            try {
                String line;
                String[] splits;
                Integer taxid;
                String name, name_class;
                while ((line = br.readLine()) != null){
                    splits = line.split("\t\\|\t");
                    taxid = Integer.parseInt(splits[0]);
                    name_class = splits[3];
                    if (name_class.startsWith("scientific name")){
                        name = splits[1].replace('/', ' ');
                        taxidToName.put(taxid, name);
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
