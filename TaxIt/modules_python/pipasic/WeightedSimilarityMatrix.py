import numpy as np
from Bio import SeqIO
from optparse import OptionParser, SUPPRESS_HELP
from re import sub

usage = '''
%prog [options]

Creates a weighted similarity matrix for the pipasic tool
'''

def constructWeightedSimilarityMatrix(digest_files, peptide_files, database_files, output_file, globalN=True):
    '''
    Objective:
    
    construct a weighted similarity matrix
    (return: proteome similarity matrix weighted by dataset in output_file)
    '''
    
    # make sure we have the same amount of datasets from digest and identified peptides
    assert len(digest_files) == len(peptide_files)
    
    similarity_matrix = np.zeros([len(digest_files),len(digest_files)]) # proteome similarity matrix

    globalTrypCountDic = None
    if (globalN): globalTrypCountDic = globalTrypticCount(digest_files)

    for i, files in enumerate(zip(peptide_files, digest_files)):
        peptide_file, digest_file = files
        
        # peptides are already found => identifiedPeptides
        # and the databases are already digested
        weighted_peptide_file = weightPeptides(peptide_file, digest_file,
                                               globalTrypCountDic=globalTrypCountDic,
                                               outfile='weightedPeptides' + str(i) + '.txt')

        for j, digest_file in enumerate(digest_files):
            pMatch, cMatch, cPept = weightedMatching(weighted_peptide_file, digest_file) # may need more arguments
            similarity_matrix[i, j] = pMatch
    
    # write the similarity matrix into a file (human unreadable format)
    np.savetxt(output_file, similarity_matrix, delimiter=',\t')

def globalTrypticCount(digest_files):
    """
    Count overall occurrence of tryptic peptides in all digest files
    (to calculate preliminary weights w~_p with global N_p = number of peptides the spectrum can be identified with)
    :param digest_files:
    :return: a dictionary with the counts of each tryptic peptide
    """

    trypCountDic = dict()  # dictionary of all peptide sequences and corresponding count

    for file in digest_files :
        for pept in SeqIO.parse(open(file, "rU"), "fasta") :

            # peptide sequence
            pept = str(pept.seq)
            pept = sub('Q', 'K', sub('I', 'L', pept))

            # dictionary of peptide sequences indicating unique and shared peptides
            if pept in trypCountDic:
                trypCountDic[pept] += 1
            else:  # add new entry to sequence dictionary & set counts per pept to 0
                trypCountDic[pept] = 1

    return trypCountDic

def weightPeptides(identifiedPeptides, trypticPeptides, globalTrypCountDic=None,
                   outfile='weighthedPeptides.txt', init=0, verbose=True):
    """
    search identified peptides in tryptic peptide list of a proteome in order 
    to weight each peptide
    (return: path to output file with line structure: "weight    sequence")
    
    identifiedPeptides: inspectparser.py output (experimental spectra)
    trypticPeptides:    typticCut result of proteome in question
    init:               initial weight for all peptides
    """
    
    ############### Read the digested Peptides ###############
    
    trypList = [] # list of tryptic peptides
    trypDic = dict() # dictionary of all peptide sequences and corresponding indices in trypList
    protDic = dict() # dictionary of protein names and corresponding indices in trypList

    # dict for peptide match counts
    matchCountDic = dict()

    if verbose: print "reading sequences of tryptic peptides ..."
    
    for idx, pept in enumerate(SeqIO.parse(open(trypticPeptides, "rU"), "fasta")) :
        
        # dictionary of protein names (descriptions)
        p_name = str(pept.description)
        if p_name in protDic: # add peptide (by index) to its protein
            protDic[p_name].append(idx)
        else: # add new protein name (+ first entry) to dictionary
            protDic[p_name] = [idx]
            
        # peptide sequence
        pept = str(pept.seq)
        pept = sub('Q','K', sub('I','L', pept))
        trypList.append(pept)
        
        # dictionary of peptide sequences indicating unique and shared peptides
        if pept in trypDic: # add index of current peptide to already known sequence
            trypDic[pept].append(idx)
        else: # add new entry to sequence dictionary & set counts per pept to 0
            trypDic[pept] = list([idx])
            matchCountDic[pept] = 0

    if verbose: print len(trypDic), "different sequences out of", len(trypList), "peptides from", len(protDic), "proteins"

    # counter for peptides not in tryptic peptide list
    notfound = 0
    
    ############### Comparison to identified peptides ###############
    
    if verbose: print "reading and searching sequences of identified peptides ..."
    
    with open(identifiedPeptides,"r") as inF:
        
        ## idXML parsing ##
        # file processing to extract the sequences
        # text = inF.read()
        # peptide_hits = text.split('<PeptideIdentification ')[1:] # => '... sequence="KJAHGSNNK" ....'
        # peptides = [peptide_hit.split('sequence="')[1].split('"')[0] for peptide_hit in peptide_hits] # => 'KJAHGSNNK'
        
        # the idxml file peptides can contain modifications like 'M(Oxidation)RSDE...' 0> needs to be filtered out 
        # peptides = [sub('\(.*?\)', '', peptide) for peptide in peptides]
        
        ## tsv parsing (MS-GF+ text format)##
        # read from tsv file
        text = inF.readlines()
        # identify peptide column and extract peptides
        peptidePosition = text[0].split('\t').index("Peptide")
        peptides = [peptide_hit.split('\t')[peptidePosition] for peptide_hit in text[1:]]
        # remove modifications
        peptides = [sub('[\+\-]\d+\.\d*', '', peptide) for peptide in peptides]
        # substitute aa with similar/same mol weight
        peptides = [sub('Q','K', sub('I','L', peptide)) for peptide in peptides]

        print peptides, "..."

        for pept in peptides:
            # search each peptide in list of proteome tryptic peptides
            if pept in matchCountDic:
                matchCountDic[pept] += 1
            else:
                notfound += 1

    # array for peptide weights
    weights = [init] * len(trypList)

    # calculate weights from counts
    if verbose: print "not found: %i of %i"%(notfound,sum(matchCountDic.values())+notfound)
    if sum(matchCountDic.values()) > 0:
        if verbose: print "recalculating counts for shared peptides ..."
        for t in trypList: # equally distribute counts for peptides of same sequence
            t_ind = trypDic.get(t)

            # norm counts either by local peptide occurrence (current proteome only)...
            if (globalTrypCountDic is None):
                x = matchCountDic[t] / float(len(t_ind))
            # or global peptide occurrence (all proteomes)
            else:
                x = matchCountDic[t] / float(globalTrypCountDic[t])

            for i in t_ind:
                weights[i] = x
    
        sum_c = sum(weights)
        for p in protDic: # equally distribute counts for peptides of same protein
            # normalized by the sum of all counts (sum_c)
            p_ind = protDic.get(p)
            x = sum([weights[c] for c in p_ind])/float(len(p_ind))/sum_c
            for i in p_ind:
                weights[i] = x
    else:
        print "warning: no peptides identified!"
    if verbose: print "writing resulting weights ..."
    cPos = 0 # number of positively weighted peptides
    outF = open(outfile,"w+")
    # output file line structure: "weight    sequence"
    for idx, seq in enumerate(trypList):
        outF.write(str(weights[idx])+"\t"+seq+"\n")
        if not weights[idx] == 0:
            cPos += 1
    outF.close()
    if verbose: print cPos,"of",len(trypList),"peptides got positive weight."
    return outfile
    
def weightedMatching(weightedPeptFile, digest_file, verbose=True):
    """
    (return: weighted sum of matches, number of matches, number of peptides)
    search all weighted peptide sequences from a given file in a given database (digested)
    
    weightedPeptFile: textfile of weighted peptide sequences 
             (weightPeptides output file line structure: "count    sequence")
    digest_file:  protein database, digested (fasta-format)
    """
    
    DB_proteins = set() # set of protein header
    DB_peptides = set() # set of peptide sequences
    cPept = 0 # number of peptides
    cMatch = 0 # number of matches (weighted peptide <-> DB peptide)
    pMatch = float(0) # weighted sum of matches
    
    if verbose: print "reading DB sequences ..."
    for pept in SeqIO.parse(open(digest_file, "rU"), "fasta") :
        DB_proteins.add(pept.name)
        DB_peptides.add(sub('Q', 'K', sub('I', 'L', str(pept.seq))))

    if verbose: print "reading and mapping peptide sequences ..."
    for line in open(weightedPeptFile, "rU"):
        if float(line.split("\t")[0]) > 0: #this line is NEW (irritates cmatches and cPept)
            cPept += 1
            peptide_seq = line.split("\t")[1].rstrip().upper() # another minor (speed) change
            if peptide_seq in DB_peptides:
                cMatch += 1
                pMatch += float(line.split("\t")[0])
        
    if verbose: print "\r%i of %i peptides found" %(cMatch,cPept)
    return pMatch, cMatch, cPept # weighted sum of matches, #matches, #peptides
    
if __name__ == '__main__':

    op = OptionParser(usage=usage)
    
    op.add_option('-d','--digest', dest='digest_files', action='append', default=None, help='files with the (tryptic) digest of the .fasta protein libaries')
    op.add_option('-p','--peptides', dest='peptide_files', action='append', default=None, help='(.idXML) files of the found peptides for the different species')
    op.add_option('-b','--db', dest='database_files', action='append', default=None, help='(.fasta) database files for the different species')
    op.add_option('-o','--output', dest='output_file', default = None, help='the file the similarity matrix is output to (currently a python pickle file)')
    op.add_option('-g','--global', default=True, help=SUPPRESS_HELP) # secret option to control globalTrypCount aka globalN
    op.add_option('-i','--init', default=0, help = 'initial weight for all peptides')
    
    opts, args = op.parse_args()
    
    ### read file lists (don't use options except -o) ###
    print args
    opts.digest_files = args[0].split(',')
    opts.peptide_files = args[1].split(',')
    opts.database_files = args[2].split(',')
    print opts.digest_files

    constructWeightedSimilarityMatrix(opts.digest_files, opts.peptide_files, opts.database_files, opts.output_file)
