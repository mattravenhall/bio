#!/usr/bin/env python3

import sys
import os

####################################
# Utilities of broader application #
####################################

def namesToFile(fasta, keepStart='N'):
    """Creates a file of names from a given fasta.
    """
    titles = []

    f = open(fasta, "r")

    for line in f:
        if ">" in line:
            if (keepStart.upper() == 'N'):
                titles.append(line.lstrip('>'))
            elif (keepStart.upper() == 'Y'):
                titles.append(line)
            else:
                print("keepStart input unclear, accepts Y/N only")
                return
    f.close()

    n = open(r'names_out.txt', 'w')

    for element in titles:
        n.write(element)
    n.close()

def fastaToDict(fasta):
    """Given a fasta as its harddrive location, returns a dictionary
    where keys are the strain name and values are their associated
    sequences. Note that dictionaries are unordered and that all
    bases will be forced to uppercase.
    """
    strains = {}
    seqcatch = ""
    f = open(fasta, "r")
    
    for line in f: #take name from > to /n, seq from /n to >
        if ">" in line:
            tmp = line[1:-1]
            strains[tmp] = seqcatch.upper()
            seqcatch = ""
        else:
            seqcatch += line.rstrip()
        strains[tmp] = seqcatch.upper()
    f.close()
    return(strains)

def getStats(filename, givenThreshold=500, scaffold=False, returnLens=False):
    """Given a fasta, return the top, bottom and mean lengths.
    Can also optionally return the full list of sequence lengths.
    """

    # NB: Replacing with a 'FASTA' class might be the best future solution

    contigLengths = [0]
    threshold = int(givenThreshold) #bp
    x = -1 # This can probably be replaced with something better suited

    f = open(filename, "r")

    # Get list of sequence lengths
    for line in f:
        if ">" in line:
            # Create new entry in contigLengths ()
            contigLengths.append(0)
            x += 1
        else:
            # Add length of line to value for that contig
            contigLengths[x] += len(line.rstrip('\n'))

    contigLengths = contigLengths[:-1] # Removes excess list entry
    totalMean = sum(contigLengths) / len(contigLengths)

    # Calculate N50 (The smallest contig length that at least half the nucleotides belong to)
    contigLengths.sort(reverse=True)

    tmpTotal = tmpTotalThreshold = 0

    for index, value in enumerate(contigLengths):
        tmpTotal += int(contigLengths[index])
        if tmpTotal > (sum(contigLengths) / 2):
            N50 = contigLengths[index]
            break

    toremove = []
    contigThreshold = []

    # Remove below threshold contigs from contigLengths
    for i in contigLengths:
        if i <= threshold:
            toremove.append(i)
    
    for i in contigLengths:
        if i not in toremove:
            contigThreshold.append(i)

    # This time for the threshold limit only, but only if the threshold removed contigs
    if len(contigThreshold) != len(contigLengths):
        for index, value in enumerate(contigThreshold):
            tmpTotalThreshold += int(contigThreshold[index])
            if tmpTotalThreshold > (sum(contigThreshold) / 2):
                N50threshold = contigThreshold[index]
                break
    else:
        N50threshold = N50

    print(len(contigLengths))
    print(len(contigThreshold))
    print(toremove)

    # Calculate GC
    statsGC = getGC(filename, total=True)

    # Optional: Return full list of contig lengths
    if returnLens:
        print ("Contig lengths: " + str(contigLengths))

    print ("Mean: " + str(totalMean))
    print ("N50: " + str(N50))
    print("N50 (>" + str(threshold) + "bp): " + str(N50threshold))
    print ("GC: " + str(statsGC) + "%")
    print ("Longest: " + str(max(contigLengths)))
    print ("Shortest: " + str(min(contigLengths)))
    print ("No. Contigs: " + str(len(contigLengths)))
    print ("Total length: " + str(sum(contigLengths)))
    ## Include statistics with a bottom threshold (e.g. N50 for contigs > 500bp)
    # Could also add mode, median etc.

#########################
# Biologically Relevant #
#########################

def countNucs(seq):
    """Return the number of A, C, G & Ts within a sequence of DNA given
    as a string.
    """
    # DNA/RNA check
    RNA = DNA = scaffold = False
    seq = seq.upper()

    if "U" in seq:
        RNA = True
    if "T" in seq:
        DNA = True
    if "N" in seq:
        scaffold = True

    if (RNA):
        if (DNA):
            raise Exception("Sequence contains both DNA and RNA.")

    A = C = G = T = U = N = 0
    for x in list(seq):
        if x.upper() == "A":
            A += 1
        elif x.upper() == "C":
            C += 1
        elif x.upper() == "G":
            G += 1
        elif x.upper() == "T":
            T += 1
        elif x.upper() == "U":
            U += 1
        elif x.upper() == "N":
            N += 1
    if (DNA):
        print("A: " + str(A) + "\n" + "C: " + str(C) + "\n" + "G: " + str(G) + "\n" + "T: " + str(T))
    if (RNA):
        print("A: " + str(A) + "\n" + "C: " + str(C) + "\n" + "G: " + str(G) + "\n" + "U: " + str(U))
    if (scaffold):
        print("N: " + str(N))

def transcribe(seq):
    """Converts a sequence of bases, provided as a string, from RNA to
    DNA or DNA to RNA. An automatic check is included to determine if
    the given sequence is DNA or RNA.
    """
    seq.upper()
    newSeq = ""

    #DNA/RNA check
    if "U" in seq:
        switch = "toDNA"
        if "T" in seq:
            print("Invalid sequence, contains both DNA and RNA. " +\
                  "Function aborted.")
            return
    if "T" in seq:
        switch = "toRNA"
        if "U" in seq:
            print("Invalid sequence, contains both DNA and RNA. " +\
                  "Function aborted.")
            return
    
    #Sequence conversion
    if switch == "toDNA":
        for x in list(seq):
            if x == "U":
                newSeq += "T"
            else:
                newSeq += str(x)
    elif switch == "toRNA":
        for x in list(seq):
            if x == "T":
                newSeq += "U"
            else:
                newSeq += str(x)
    print(newSeq)

def getComplement(seq):
    """Returns the reverse complement of a given sequence,
    as a string
    """
    revSeq = seq[::-1]
    newSeq = ""
    for x in list(revSeq):
        if x == "A":
            newSeq += "T"
        elif x == "C":
            newSeq += "G"
        elif x == "G":
            newSeq += "C"
        elif x == "T":
            newSeq += "A"
    print(newSeq)

def getHeteroProb(k, m, n):
        """Returns the probability of gaining a dominate positive
        offspring for a pairing within a population where dominant (k),
        heterozygous (m) and recessive (n) traits are known
        """

        # Probability of parent X being D, H or R
        D = k/(k+m+n)
        H = m/(k+m+n)
        R = n/(k+m+n)

        # Probability of parent Y being D, H or R
        DD = (k-1)/((k+m+n)-1)
        Dx = (k)/((k+m+n)-1)
        HH = (m-1)/((k+m+n)-1)
        Hx = (m)/((k+m+n)-1)
        RR = (n-1)/((k+m+n)-1)
        Rx = (n)/((k+m+n)-1)

        # Probabilities of dominant positive for pairings
        pDD = 1
        pDH = 1
        pDR = 1
        pHH = 0.75
        pHR = 0.5
        pRR = 0

        # Full calculation
        prob = (D*DD*pDD)+(D*Dx*pDH)+(D*Dx*pDR)+(H*HH*pHH)+(H*Hx*pDH)+(H*Hx*pHR)+(R*RR*pRR)+(R*Rx*pDR)+(R*Rx*pHR)
        print(prob)        

def getGC(fasta, total=False):
    """Given the location of a fasta, returns the GC value for each
    strain as a dictionary by default. Alternatively, the total GC will
    be returned if total=True.
    """
    GC = 0
    Ncount = 0
    if (total == False): # Default: return dict of strains with GCs
        strains = fastaToDict(fasta) # Convert fasta into dictionary
        strainsGC = {}
        for key in strains: # Calculate GCs and assign to strains
            for base in strains[key]: # Workhorse
                if base.upper() in ("G", "C"):
                    GC += 1
                if base.upper() in ("N"):
                    Ncount += 1
            GCperc = (GC / (len(strains[key]) - Ncount)) * 100
            strainsGC[key] = GCperc
        return(strainsGC)

    elif (total): # Alternative: return the total GC for all sequences
        totalLength = 0
        f = open(fasta, "r")
        for line in f: # Iterate over sequences, ignore sequences names
            if ">" not in line:
                totalLength += len(line.rstrip('\n'))
                for base in line:
                    if base.upper() in ("G", "C"):
                        GC += 1
        totalGC = (float(GC) / float(totalLength)) * 100
        f.close()
        return(round(totalGC, 2))
        
    
def getHighestGC(fasta):
    """Given a fasta file location, prints the strain with the highest
    GC content and its value.
    """
    strainsGC = getGC(fasta) # Get GCs
    highestGC = 0
    highestGCstrain = ""
    for key in strainsGC: # Find the highest GC, return the strain & value
        if strainsGC[key] > highestGC:
            highestGC = strainsGC[key]
            highestGCstrain = key
    print(str(highestGCstrain))
    print(str(strainsGC[highestGCstrain]))

def calcHamming(A, B):
    """Calculate the Hamming distance between two sequences (as strings)
    of equal length.
    """
    dH = 0
    if len(A) == len(B):
        for index, value in enumerate(A):
            if value != B[index]:
                dH += 1
    if len(A) != len(B):
        print("Error: Sequences must be the same length!")
        return
    print(dH)

def translate(seq):
    """Given an RNA sequence, as a string, returns the protein
    sequence.
    """
    ##Convert RNA to Protein
    #Could include a "check type" function to check if DNA/RNA/Protein
    protein = ""
    end = False
    ticker = 0
    RNAcodons = {"UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L",
                 "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
                 "UAU": "Y", "UAC": "Y", "UAA": "-", "UAG": "-",
                 "UGU": "C", "UGC": "C", "UGA": "-", "UGG": "W",
                 "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
                 "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
                 "CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
                 "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",
                 "AUU": "I", "AUC": "I", "AUA": "I", "AUG": "M",
                 "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
                 "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K",
                 "AGU": "S", "AGC": "S", "AGA": "R", "AGG": "R",
                 "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
                 "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
                 "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
                 "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G"}
    while end == False:
        fc = ticker
        codon = seq[fc:fc+3]
        protein += RNAcodons[codon]
        ticker += 3
        if fc+3 >= len(seq):
            end = True
    print(protein)

def findMotif(motif, seq, vocal=False): 
    # TODO: Optionally supress print output, perhaps behind a debug=False flag
    """When given a motif <string> and a fasta <file location> or 
    sequence <string> returns the locations of that motif as a 
    dictionary if a file is provided or a list if a string was provided.
    """
    # Force motif into uppercase to avoid case issues
    motif = motif.upper()

    # Check is seq is a file or a sequence:
    if (os.path.isfile(seq)): # if seq is a file
        allSeqs = fastaToDict(seq) # convert fasta to dictionary
        locationsDict = {}

        for key in allSeqs: # Iterate over each sequence
            if (vocal):
                print("Searching " + str(key) + "...")
            # Force seq into uppercase to avoid case issues
            seq = allSeqs[key].upper()

            locations = []
            mod = 1
            for index, base in enumerate(seq): # Iterate over each base
                if base == motif[0]:
                    while mod < len(motif) and (index + mod) < len(seq):
                        if seq[index+mod] == motif[mod]:
                            # Base match found, scan for whole motif
                            mod += 1
                        else:
                            mod = 1 # Vital reset of mod on fail
                            break
                    if mod >= len(motif):
                        if (vocal):
                            print("Whole motif found, adding location...")
                        locations.append(index+1)
                        mod = 1
            if locations != []: # Only add to hits to locations if not empty
                locationsDict[key] = locations
        return(locationsDict) # Returns a dictionary of hits according to each strain/contig
    elif (type(seq)) == str: # if seq is a string
        # Force seq into uppercase to avoid case issues
        seq = seq.upper()

        locations = []
        mod = 1
        for index, base in enumerate(seq):
            if base == motif[0]:
                while mod < len(motif) and (index + mod) < len(seq):
                    if seq[index+mod] == motif[mod]:
                        # Base match found, scan for whole motif
                        mod += 1
                    else:
                        mod = 1 # Vital reset of mod on fail
                        break
                if mod >= len(motif):
                    print("Whole motif found, adding location...")
                    locations.append(index+1)
                    mod = 1
        return(locations) # NB: locations is a dictionary for files and a list for strings
    else:
        # If not a file or string, abort with an error.
        # Not that os.path.isfile will probably raise an error before this pops.
        raise TypeError("Inappropriate datatype supplied, findMotif currently only accepts strings and fastas.")    

def countRepeats(motifs, ref):
    """Print the number of times a motif repeats within a larger 
    sequence. Will work when given a string or a list of strings.
    """
    if (type(motifs) == str): # convert a single string into a one element list
        motifs = [motifs]
    for motif in motifs:
        motifDict = findMotif(motif, ref)
        count = sum(len(i) for i in motifDict.itervalues())
        print(motif + ": " + str(count))  

def findConsensus(fasta):
        """Returns the consensus string and profile matrix for a given
        set of strains.
        Strains must be provided as a fasta. NB: Errors will occur when
        there is a base conflict.
        """
        base = ['A', 'C', 'G', 'T']
        baseDict = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
        count = 0
        consensus = ""

        #pull sequences into a dictionary {strain: sequence}
        sequences = fastaToDict(fasta)

        #check sequences are the same length, else abort
        #if check passes, assigns that length to seqLength
        seqLen = []
        for value in sequences.values():
                seqLen.append(len(value))
        if max(seqLen) != min(seqLen): #exit if sequences different lengths
                print("Error: Sequences must be of the same length!")
                print("Sequences current range from " + str(min(seqLen)) +\
                      " to " + str(max(seqLen)) + " bases long.")
                return
        else: #set sequence lengths
                seqLength = max(seqLen)
                
        #create empty matrix
        profile = [[0]*seqLength for i in range(4)]

        #populate matrix
        for value in sequences.values():
                for basePos, baseType in enumerate(value):
                       profile[baseDict[baseType]][basePos] += 1
                       #NB: base position = y, base type = x
        
        for i in range(seqLength):
                colHolder = []
                for value, x in enumerate(profile):
                        colHolder.append(profile[value][i])
                consensus += (base[colHolder.index(max(colHolder))])

        #print consensus string
        print(consensus)
                
        #print profile matrix
        for x in profile:
                row = ' '.join(map(str, x))
                print(base[count] + ": " + row)
                count += 1

def findKmers(length, reference):
    """Given a kmer length and a genome as a fasta, returns a dictionary
    of kmers and their counts.
    """
    refDict = fastaToDict(reference)
    candidates = {} # sequence:count

    for key in refDict:
        for i, base in enumerate(refDict[key]): # Focus on the seq for each contig
            if i > len(refDict[key])-int(length): # Can this be encorporated into the for loop statement?
                break;
            window = refDict[key][i:i + int(length)] # Move across sequence with window
            if window not in candidates:
                candidates[window] = 1
            elif window in candidates:
                candidates[window] += 1
    return candidates

def predictMT(seq):
    """Given a sequence as a string, returns a rough prediction of its 
    melting temperature as an int. Note that an empirically determined 
    melting temperature may be significantly different.
    """
    GC = AT = 0

    # determine GC & AT content of given sequence
    for base in seq.upper():
        if base in ("G", "C"):
            GC += 1
    AT = len(seq) - GC

    # predict melting temperature, differs with seq length
    if len(seq) < 14:
        mt = 4 * GC + 2 * AT
    else:
        mt = 64.9 + 41 * (GC - 16.4)/len(seq)
    return(mt)

####################################################
# For distinguishing command line/import behaviour #
####################################################

# Gateway function
def main(args):
    # args[0] = subfunction to be called, args[1:] = input arguments for args[0].
    # e.g. if 'mt', call predictMT() with arguments
    if len(args) == 1:
        print("Warning: No arguments sent to function.")

    if args[0].lower() == "stats":
        if len(args) == 2:
            getStats(args[1]) # passing additional arguments needs to be optimised, current solution is terrible
        if len(args) == 3:
            getStats(args[1], args[2])
        if len(args) == 4:
            getStats(args[1], args[2], args[3])
    if args[0].lower() == "mt":
        predictMT(args[1])
    if args[0].lower() == "kmers":
        findKmers(args[1], args[2])
    if args[0].lower() == "consensus":
        findConsensus(args[1])
    if args[0].lower() == "translate":
        translate(args[1])
    if args[0].lower() == "transcribe":
        transcribe(arg[1])
    # else:
    #     print("Operation aborted: Function not recognised.")
    #     sys.exit()
    pass

# If being directly executed (ie. not imported)
if __name__ == "__main__":
    if len(sys.argv) <= 1: # ie. if no arguments were passed to biocore
        # Fill this with something useful explaining basic uses of biocore
        print("\nUsage: biocore <command> <arguments>\n\nCommands:\n"
            +"stats\tBasic statistics for fastas\n"
            +"mt\tRough melting temperature prediction\n"
            +"kmers\tFind kmers of given length within a fasta\n"
            +"translate\tTranslate from RNA/DNA to DNA/RNA, auto-detects\n"
            +"transcribe\tTranscribe an RNA sequence to proteins\n"
            +"\nNB: This list is incomplete & will be added to later.\n")
        sys.exit()
    else:
        exit(main(sys.argv[1:])) # Call main() with arguments from the command line