/**
 * Class to select the proper TaxidMappingParser
 */
public class TaxidMappingParserFactory {

    public static TaxidMappingParser getParser(String filename){
        if (filename.endsWith("prot.accession2taxid")){
            return (new TaxidMappingParserNCBIAcc());
        }
        else if (filename.endsWith("gi_taxid_prot.dmp")){
            return (new TaxidMappingParserNCBIGI());
        }
        else if (filename.endsWith("idmapping_selected.tab")){
            return (new TaxidMappingParserUniRef());
        }
        else{
            System.err.println("error: unexpected mapping file:\n\t" + filename);
            System.exit(1);
            return(null);
        }
    }
}
