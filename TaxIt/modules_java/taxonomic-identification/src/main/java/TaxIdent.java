import files.TSV;
import tools.MemoryMonitor;
import tools.Timer;

import java.util.HashSet;

/**
 * Created by kuhringm on 1/5/17.
 */
public class TaxIdent {

    public static void main(String[] args){

        MemoryMonitor memory = new MemoryMonitor();

        Timer totalTime = new Timer();
        totalTime.start();
        Timer timer = new Timer();

        // parameter
        if (args.length != 3){
            System.out.println("wrong number of parameters!");
            System.out.println("need: result_tsv, id_taxid_file, outputFolder");
            System.exit(1);
        }

        int index = 0;
        // files
        String msgf_results = args[index++];
        String id_taxid_file = args[index++];
        String output_file = args[index++];

        // Import PSM/hit IDs
        TSV psms = null;
        if (id_taxid_file.endsWith("gi_taxid_prot.dmp")){
            psms = new TSV(msgf_results);
        }
        else if (id_taxid_file.endsWith("prot.accession2taxid")){
            psms = new TSV(msgf_results, "accession");
        }
        else if (id_taxid_file.endsWith("idmapping_selected.tab")){
            psms = new TSV(msgf_results, "uniref");
        }
        else{
            System.err.println("error: unexpected mapping file:\n\t" + id_taxid_file);
            System.exit(1);
        }

        psms.parseProteinIDs();
        HashSet<String> matches = psms.getProteinIDs();
        memory.printGbUsage();

        System.out.println(matches.size());
        System.out.println(matches);

        // parse taxid mapping, filtered by matches
        TaxidMappingParser parser = TaxidMappingParserFactory.getParser(id_taxid_file);
        parser.parse(id_taxid_file, matches);

        memory.printGbUsage();
        // export mapping
        parser.export(output_file);
    }

}
