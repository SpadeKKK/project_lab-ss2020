import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;

/**
 * Abstract class for importing/parsing and exporting protein-taxid mappings
 */
public abstract class TaxidMappingParser {

    private final String idType;

    protected HashMap<String, HashSet<Integer>> protein2taxid;

    public TaxidMappingParser(String idType){
        this.idType =  idType;
        protein2taxid = new HashMap<String, HashSet<Integer>>();
    }

    public void parse(String filename, HashSet<String> filter){
        // parse mapping file
        System.out.println("parse " + idType + " to taxid mapping...");
        System.out.println("\tstart: " + new Date().toString());
        parseMapping(filename, filter);

        System.out.println("\tdone: " + new Date().toString());

        System.out.println(protein2taxid);
    }

    abstract protected void parseMapping(String filename, HashSet<String> filter);

    public void export(String filename){

        try {
            BufferedWriter bw = new BufferedWriter(new FileWriter(new File(filename)));
            try {
                for (Map.Entry<String, HashSet<Integer>> entry : protein2taxid.entrySet()) {
                    for (Integer taxid : entry.getValue()) {
                        bw.write(entry.getKey() + "\t" + taxid);
                        bw.newLine();
                    }
                }
            }
            finally {
                bw.close();
            }
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

}
