
import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

/**
 * Class for importing NCBI taxid nodes and representing the tree structure
 */
public class Nodes {

    // the NCBI taxonomy tree is represented by 3 HashMaps.
    // each node (represented by a taxid) in the taxonomy is linked to its
    // parent (taxidToParent) and to its taxonomic rank (taxidToRank).
    private HashMap<Integer, Integer> taxidToParent;
    private HashMap<Integer, String> taxidToRank;
    // an additional HashMap contains all nodes with a certain taxonomic rank
    private HashMap<String, HashSet<Integer>> rankToTaxids;

    public Nodes(){
        taxidToParent = new HashMap<Integer, Integer>();
        taxidToRank = new HashMap<Integer, String>();
        rankToTaxids = new HashMap<String, HashSet<Integer>>();
    }

    public HashMap<Integer, Integer> getTaxidToParent() {
        return taxidToParent;
    }

    public HashMap<Integer, String> getTaxidToRank() {
        return taxidToRank;
    }

    public HashMap<String, HashSet<Integer>> getRankToTaxids() {
        return rankToTaxids;
    }

    public synchronized String getRank(int taxid){
        return(taxidToRank.get(taxid));
    }

    // parse a ncbi taxonomy node file, keep taxid, parent_taxid and rank
    // fills HashMaps with:
    // taxid->parent_taxid
    // taxid->rank
    // example line, note the separater is "\t|\t" and the end of line is "\t|\n"
    // "1	|	1	|	no rank	|		|	8	|	0	|	1	|	0	|	0	|	0	|	0	|	0	|		|"
    public void parseNodes(String filename){

        try {
            BufferedReader br = new BufferedReader(new FileReader(new File(filename)));

            try {
                String line;
                String[] splits;
                Integer taxid;
                Integer parent_taxid;
                String rank;
                while ((line = br.readLine()) != null){
                    splits = line.split("\t\\|\t");
                    taxid = Integer.parseInt(splits[0]);
                    parent_taxid = Integer.parseInt(splits[1]);
                    rank = splits[2];

                    taxidToParent.put(taxid, parent_taxid);
                    taxidToRank.put(taxid, rank);
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

    public ArrayList<String> getRanks() {
        return new ArrayList<String>(rankToTaxids.keySet());
    }
}
