
import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;

/**
 * Class for importing id to NCBI taxid mappings
 * e.g. GIs to TAXID, UniRef Cluster ID to TAXID
 */
public class IDMapping {

    // the relation of protein ids and taxids is represented by two HashMaps:
    private HashMap<String, HashSet<Integer>> idToTaxid;
    private HashMap<Integer, HashSet<String>> taxidToId;

    public IDMapping(){
        idToTaxid = new HashMap<String, HashSet<Integer>>();
        taxidToId = new HashMap<Integer, HashSet<String>>();
    }

    public HashMap<String, HashSet<Integer>> getIdToTaxid() {
        return idToTaxid;
    }

    public HashMap<Integer, HashSet<String>> getTaxidToId() {
        return taxidToId;
    }

    public synchronized HashSet<Integer> getTaxids(HashSet<String> matches){
        HashSet<Integer> taxids = new HashSet<Integer>();
        for (String id : matches){
            taxids.addAll(idToTaxid.get(id));
        }
        return(taxids);
    }

    /**
     * Import protein taxid mappings from file
     * (as created by the taxonomic-identification module)
     *
     * @param filename Name of the protein mapping file
     */
    public void importMapping(String filename){

        try {
            BufferedReader br = new BufferedReader(new FileReader(new File(filename)));

            try {
                String line;
                String[] splits;
                String id;
                Integer taxid;
                while ((line = br.readLine()) != null){
                    splits = line.split("\t");
                    id = splits[0];
                    taxid = Integer.parseInt(splits[1]);

                    if (!idToTaxid.containsKey(id)){
                        idToTaxid.put(id, new HashSet<Integer>());
                    }
                    idToTaxid.get(id).add(taxid);

                    if (!taxidToId.containsKey(taxid)){
                        taxidToId.put(taxid, new HashSet<String>());
                    }
                    taxidToId.get(taxid).add(id);
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

    /*
    * Remove ids from a given set which have no corresponding taxid
    */
    public ArrayList<String> removeIdsWithNoTaxidMapping(HashSet<String> matches){
        ArrayList<String> removed = new ArrayList<String>();

        for (Iterator<String> i = matches.iterator(); i.hasNext();) {
            String id = i.next();
            if (!idToTaxid.containsKey(id)) {
                removed.add(id);
                i.remove();
            }
        }

        return(removed);
    }

    /*
    * Remove taxids from taxidToId/idToTaxid which are not available nodeTaxids (e.g. NCBI tax nodes)
    */
    public void removeTaxidsWithoutNodes(HashSet<Integer> nodeTaxids){

        System.out.println("clean up mapped taxid... ");

        ArrayList<Integer> removed = new ArrayList<Integer>();

        // for each taxid in the mapping
        HashSet<Integer> mappingTaxids = new HashSet<Integer>(taxidToId.keySet());
        for (int taxid : mappingTaxids){
            // if it is not in the filter
            if (!nodeTaxids.contains(taxid)){
                // iterate corresponding ids...
                for (String id : taxidToId.get(taxid)){
                    // ... to remove id2taxid mappings
                    idToTaxid.get(id).remove(taxid);
                    // remove id completely if no mapping is left
                    if (idToTaxid.get(id).size() == 0){
                        idToTaxid.remove(id);
                    }
                }
                // remove taxid
                taxidToId.remove(taxid);
                removed.add(taxid);
            }
        }

        if (removed.size() > 0){
            System.out.println("removed " + removed.size() + " taxids without corresponding nodes... ");
            System.out.println("\tfirst 3 removed taxids: " + removed.subList(0, Math.min(3, removed.size())));
        }
    }

    public HashSet<Integer> getAllTaxids(){
        return new HashSet<Integer>(taxidToId.keySet());
    }

    public HashSet<String> getIds(int taxid){
        return new HashSet<String>(taxidToId.get(taxid));
    }

    public int countTaxids() {
        return taxidToId.size();
    }

    public boolean hasGis(int taxid) {
        return taxidToId.containsKey(taxid);
    }
}
