#!/usr/bin/env python3

import sys
import os

# TODO: Build a fastq to fasta converter (base qualities will have to be removed)
# TODO: Update findMotif to support IUPAC ambiguous nucleotides (Y, R, W etc.)
# TODO: Add ambiguous allele support for translate.

####################################
# Utilities of broader application #
####################################

def detectType(filename):
    """Given filename, will attempt to predict and return the filetype.
    Note that this depends entirely on the name of the file given.
    In future this may be replaced with a more dynamic approach.
    """
    if filename[-2:] == "fa" or filename[-5:] == "fasta":
        return("fasta")
    elif filename[-2:] == "fq" or filename[-5:] == "fastq":
        return("fastq")
    elif filename[-2:] == "gz":
        return("gzip")

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

def ToDict(filename):
    """Given a fasta or fastq, returns a dictionary with the format
    {"contig_name": ["sequence"]} if a fasta and
    {"contig_name": ["sequence", "qualities"]} if a fastq.
    Note that dictionaries are unordered and that all bases will be
    forced to uppercase.
    """

    dictionary = {}
    seqcatch = ""
    filetype = detectType(filename)

    f = open(filename, "r")

    if filetype.lower() == "fasta":
        for line in f: #take name from > to /n, seq from /n to >
            if ">" in line:
                tmp = line[1:-1]
                if tmp in dictionary.keys():
                    raise Exception("Warning: Duplicate ToDict() key identified for: "+tmp+" Please change name for this sequence.")
                dictionary[tmp] = seqcatch.upper()
                seqcatch = ""
            else:
                seqcatch += line.rstrip()
            dictionary[tmp] = [seqcatch.upper()]
    elif filetype.lower() == "fastq":
        for line1 in f:
                if line1.rstrip('\n') in dictionary.keys():
                    raise Exception("Warning: Duplicate ToDict() key identified for: "+tmp+" Please change name for this sequence.")
                line2 = next(f)
                line3 = next(f)
                line4 = next(f)
                dictionary[line1.rstrip('\n')] = [line2.upper().rstrip('\n'), line4.rstrip('\n')]
    else:
        raise Exception("Incorrect file format supplied, please indicate fasta or fastq.")

    f.close()
    return(dictionary)

def ToList(filename):
    """Given a fasta file, return a list with the format:
    [contig1, contig2, ..., contigN]
    Note that this function will remove contig labels.
    """

    # TODO: Allow fasta/fastq recognition.

    theList = []
    seqcatch = ''

    with open(filename, "r") as f:
        for line in f: #ignore original > line, catch all else
            if ">" in line:
                pass
            else:
                if 'N' not in line.upper():
                    seqcatch += line.rstrip()
                elif 'N' in line.upper():
                    for base in line.rstrip():
                        if base.upper() != 'N':
                            seqcatch += base
                        elif base.upper() == 'N':
                            if seqcatch != '':
                                theList.append(seqcatch)
                                seqcatch = ''
                            else:
                                pass
    return(theList)

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

def transcribe(seq, asPrint=False):
    """Converts a sequence of bases, provided as a string, from RNA to
    DNA or DNA to RNA. An automatic check is included to determine if
    the given sequence is DNA or RNA.
    """
    seq = seq.upper()
    newSeq = ""

    #DNA/RNA check - this breaks if U or T is not present
    if "U" in seq:
        switch = "toDNA"
        if "T" in seq:
            print("Invalid sequence, contains both DNA and RNA. " +\
                  "Function aborted.")
            return
    elif "T" in seq:
        switch = "toRNA"
        if "U" in seq:
            print("Invalid sequence, contains both DNA and RNA. " +\
                  "Function aborted.")
            return
    else:
        switch = "toAlt"

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
    elif switch == "toAlt":
        newSeq == seq

    #if __name__ == "__main__": # for command line execution
    if asPrint:
        print(switch[2:]+" Sequence: "+newSeq)
    else:
        return(newSeq)

def getComplement(seq, silent=False): # 'silent' is internal use only
    """Returns the reverse complement of a given sequence,
    as a string. DNA-exclusive. Supports ambiguous alleles,
    including unknowns (as 'X' or 'N') and gaps (as '-').
    """
    revSeq = seq[::-1]
    newSeq = ""

    # This can be optimised slightly but functions perfectly
    for x in list(revSeq):
        if x == "A":
            newSeq += "T"
        elif x == "C":
            newSeq += "G"
        elif x == "G":
            newSeq += "C"
        elif x == "T":
            newSeq += "A"
        elif x == "Y":
            newSeq += "R"
        elif x == "R":
            newSeq += "Y"
        elif x == "W":
            newSeq += "W"
        elif x == "S":
            newSeq += "S"
        elif x == "K":
            newSeq += "M"
        elif x == "M":
            newSeq += "K"
        elif x == "D":
            newSeq += "H"
        elif x == "V":
            newSeq += "B"
        elif x == "H":
            newSeq += "D"
        elif x == "B":
            newSeq += "V"
        elif x == 'X':
            newSeq += 'X'
        elif x == 'N':
            newSeq += 'N'
        elif x == '-':
            newSeq += '-'
        else:
            raise Exception('Error: Non-DNA strands cannot be complemented.')
    if __name__ == "__main__" and not silent: # for command line execution
        print("Complement: "+newSeq[::-1])
    return(newSeq[::-1])

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

def calcHamming(seqA, seqB):
    """Calculate the Hamming distance between two sequences (as strings)
    of equal length.
    """
    distH = 0
    if len(seqA) == len(seqB):
        for index, value in enumerate(seqA):
            if value != seqB[index]:
                distH += 1
    elif len(seqA) != len(seqB):
        print("Error: Sequences must be the same length!")
        return
    print(distH)

def translate(seq, silent=False): # 'silent' is internal use only
    """Given an RNA sequence, as a string, returns the protein
    sequence. If given a DNA sequence, will attempt conversion
    to RNA.
    """

    #Detects DNA and converts to RNA
    if "T" in seq.upper():
        if not silent:
            print("Notice: DNA sequence detected, converting to RNA")
    seq = transcribe(seq)

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
        if len(seq[fc:fc+3]) < 3:
            break
        codon = seq[fc:fc+3]
        protein += RNAcodons[codon]
        ticker += 3
        if fc+3 >= len(seq):
            end = True
    if __name__ == "__main__": # for command line execution
        if not silent:
            print("Protein Sequence: "+protein)
    return(protein)

def findMotif(motif, seq, vocal=False, asPrint=False):
    """When given a motif <string> and a fasta <file location> or
    sequence <string> returns the locations of that motif as a
    dictionary if a file is provided or a list if a string was provided.
    """
    # Force motif into uppercase to avoid case issues
    motif = motif.upper()

    # Check is seq is a file or a sequence:
    if (os.path.isfile(seq)): # if seq is a file
        allSeqs = ToDict(seq) # convert fasta to dictionary
        locationsDict = {}

        for key in allSeqs: # Iterate over each sequence
            if (vocal):
                print("Searching " + str(key) + "...")
            # Force seq into uppercase to avoid case issues
            seq = allSeqs[key][0].upper()

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
        if asPrint:
            for keys,values in locationsDict.items():
                print("Motif found starting at (bp):",keys,values)
        else:
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
        if asPrint:
            print("Motif found starting at (bp):",locations)
        else:
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
        sequences = ToDict(fasta)

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

def findKmers(reference, length):
    """Given a kmer length and a genome as a fasta, returns a dictionary
    of kmers and their counts.
    """
    refDict = ToDict(reference)
    candidates = {} # sequence:count

    for key in refDict:
        for i, base in enumerate(refDict[key][0]): # Focus on the seq for each contig
            if i > len(refDict[key][0])-int(length): # Can this be encorporated into the for loop statement?
                break;
            window = refDict[key][0][i:i + int(length)] # Move across sequence with window
            if window not in candidates:
                candidates[window] = 1
            elif window in candidates:
                candidates[window] += 1

    return(candidates)

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
    print(mt)

def simCleaveMulti(genomefile, enzyme, csite):
    """The fasta/fastq multi-sequence gateway to simCleave.
    """
    genomeS = ToDict(genomefile)
    for strain in genomeS:
        print("\nSimulating cleavage of "+strain+" by "+enzyme+"...")
        simCleave(genomeS[strain][0], enzyme, csite)

def simCleave(genome, enzyme, csite):
    """Given a restriction enzyme's recognition site (string), the
    index at which it cleaves (int) and a sequence to cleave (string),
    returns a list of the fragment lengths and a list of the fragment
    sequences.

    E.g. for enzyme 'GCCG' with cleavage site 'GC|CG' csite = 2.
    NB: Non-specific base sites should be indicated with 'N'.
    NB: Both the genome and enzyme sequences should be provided 5'->3'.
    """

    sites = []
    fragLens = []
    csite = int(csite)
    enzyme = enzyme.upper()
    enzymeComp = enzyme[::-1]

    for i, baseG in enumerate(genome):
        if (genome[i].upper() == enzyme[0]) or (enzyme[0] == "N"):
            for j, baseE in enumerate(enzyme):
                if (i+len(enzyme)) > len(genome): # if the enzyme would go beyond the given sequence
                    break
                if (enzyme[j] == genome[i+j].upper()) or (enzyme[j] == "N"):
                    if (j + 1) == len(enzyme): #ie. if there is a complete match
                        sites.append(i + csite)
                else:
                    break
        elif (genome[i].upper() == enzymeComp[0]) or (enzymeComp[0] == "N"):
            # NB: This loop only makes sense in regards to double-strained DNA
            for l, baseC in enumerate(enzymeComp):
                if (i+len(enzymeComp)) > len(genome):
                    break
                if (enzymeComp[l] == genome[i+l].upper()) or enzymeComp[l] == "N":
                    if (l + 1) == len(enzyme):
                        sites.append(len(enzymeComp) - csite)
                else:
                    break
    for k in sites[::-1]: #go over the cleavage sites in reverse
        genome = genome[:k] + '\n' + genome[k:]
    genome = genome.split('\n')

    for fragment in genome:
        fragLens.append(len(fragment))

    print(sorted(fragLens, reverse=True))
    print(sorted(genome, key=len, reverse=True))

def simPCRMulti(genomefile, primer1, primer2, passmark=90):
    """The fasta/fastq multi-sequence gateway to simPCR.
    """
    genomeS = ToDict(genomefile)
    for strain in genomeS:
        print("\nSimulating PCR of "+strain+" by "+primer1+" and "+primer2+"...")
        simPCR(genomeS[strain][0], primer1, primer2, passmark)

def simPCR(sequence, primer1, primer2, passmark=90):
    """Given strings for a base sequence and two primer sequences,
    returns the fragment/s that PCR amplification would produce.
    NB: Both primer sequences should be in the 5' to 3' direction.
    NB: Primer recognition is currently set at >= 80%.
    """

    sequence = sequence.upper()

    primer = primer1.upper()
    complement = getComplement(primer2.upper(),silent=True)[::-1]

    frags = []
    fraglens = []

    csites = []
    fragstmp = []

    for i, base in enumerate(sequence): # for each base in sequence
        count = 0 # reset match counter
        total = len(primer)
        for j, baseP in enumerate(primer): # for each base in primerA
            if i+len(primer) > len(sequence): # prevent overhanging matches
                break
            if primer[j] == sequence[i+j]: # if seq base matches primer base
                count += 1
            if ((j+1)==len(primer)) & (((count/total)*100)>=passmark):
                # if at end of primer & match >= 80%
                csites.append(i) # add cleavage site to list
    for k in csites: # perform cleavage by primer sites
        fragstmp.append(sequence[k:]) # fragments produced from step 1

    for fragment in fragstmp: # apply complement primerB to cleavage on fragments
        csites = []
        for i, base in enumerate(fragment): # for each base in fragment
            count = 0
            total = len(complement)
            for j, baseP in enumerate(complement): # for each base in primerB
                if i+len(complement) > len(fragment): # prevent overhanging matches
                    break
                if complement[j] == fragment[i+j]:
                    count += 1
                if ((j+1)==len(complement)) & (((count/total)*100)>=passmark):
                    # if at end of complement primer & good enough match
                    csites.append(i+j) # add cleavage site to list
        for l in csites:
            frags.append(fragment[:l]) # fragments produced from all cleavage events
    for piece in frags:
        fraglens.append(len(piece))
    print(sorted(fraglens, reverse=True))
    print(sorted(frags, key=len, reverse=True))

def scaffoldToContigs(infile, outfile):
    """Given a scaffolded genome as a fasta, returns a fasta of the
    contigs with a given output name.
    Note that this approach will remove the > identifier.
    """
    wholeGenome = ToList(infile)

    # TODO: Allow fasta and fastq distinction?
    # TODO: Allow input with more than a single > line

    with open(outfile, 'w') as f:
        count = 0
        for contig in wholeGenome:
            #print > line
            f.write('> MKAN_contig_'+str(count)+'\n')
            #print section out
            f.write(contig+'\n')
            count += 1

def AAchange(snp, gene):
    """Given a gene sequence and SNP, return the amino acid change.
    """
    ## As predictMutation('4-A','GATCATGCATGCAGACTAGCATCGA')

    # Set SNP variables
    if '-' not in snp:
        print('Error: SNP must be formatted as <location-allele>.')
        return
    snpLoc = int(snp.split('-')[0]) - 1
    snpAllele = str(snp.split('-')[1])

    # SNP location to AA location
    AAloc = snpLoc // 3

    # Reference amino acid sequence
    refAA = translate(gene, silent=True)

    # Mutation sequence
    newGene = list(gene)
    newGene[snpLoc] = snpAllele
    newGene = "".join(newGene)

    newAA = translate(newGene, silent=True)

    print('AA Change: '+refAA[AAloc]+str(AAloc+1)+newAA[AAloc]) # e.g. AA Change: Y412G

####################################################
# For distinguishing command line/import behaviour #
####################################################

# Gateway function
def main(args):
    # args[0] = subfunction to be called, args[1:] = input arguments for args[0].
    # e.g. if 'mt', call predictMT() with arguments
    if len(args) == 1:
        print("Warning: No arguments sent to function.")
    if args[0].lower() == "predictmt":
        if len(args) >= 2:
            predictMT(args[1])
        else:
            return("Required arguments: <sequence:str>")
    if args[0].lower() == "findkmers":
        if len(args) >= 2:
            findKmers(args[1], args[2])
        else:
            return("Required arguments: <genome:fasta/q> <length:int>")
    if args[0].lower() == "consensus":
        if len(args) >= 2:
            findConsensus(args[1])
        else:
            return("Required arguments: <sequences:fasta/q>")
    if args[0].lower() == "translate":
        if len(args) >= 2:
            translate(args[1])
        else:
            return("Required arguments: <rna:str>")
    if args[0].lower() == "transcribe":
        if len(args) >= 2:
            transcribe(args[1], asPrint=True)
        else:
            return("Required arguments: <sequence:str>")
    if args[0].lower() == "simcleave":
        if len(args) >= 2:
            simCleave(args[1], args[2], args[3])
        else:
            return("Required arguments: <sequence:str> <enzyme:str> <cleavage_site:int>")
    if args[0].lower() == "simcleavemulti":
        if len(args) >= 2:
            simCleaveMulti(args[1], args[2], args[3])
        else:
            return("Required arguments: <sequences:file> <enzyme:str> <cleavage_site:int>")
    if args[0].lower() == "simpcr":
        if len(args) >= 2:
            if len(args) == 4:
                simPCR(args[1], args[2], args[3])
            elif len(args) == 5:
                simPCR(args[1], args[2], args[3], args[4])
        else:
            return("Required arguments: <sequence:str> <primer1:str> <primer2:str>")
    if args[0].lower() == "simpcrmulti":
        if len(args) >= 2:
            if len(args) == 4:
                simPCRMulti(args[1], args[2], args[3])
            elif len(args) == 5:
                simPCRMulti(args[1], args[2], args[3], int(args[4]))
        else:
            return("Required arguments: <sequences:location> <primer1:str> <primer2:str>")
    if args[0].lower() == "complement":
        if len(args) >= 2:
            getComplement(args[1])
        else:
            return("Required arguments: <sequence:str>")
    if args[0].lower() == "calchamming":
        if len(args) >= 3:
            calcHamming(args[1], args[2])
        else:
            return("Required arguments: <sequence1:str> <sequence2:str>")
    if args[0].lower() == "scafftocontigs":
        if len(args) >= 3:
            scaffoldToContigs(args[1], args[2])
        else:
            return("Required arguments: <infile:file_location> <outfile:str>")
    if args[0].lower() == "findmotif":
        if len(args) >= 3:
            findMotif(args[1], args[2], vocal=False, asPrint=True)
        else:
            return("Required arguments: <motif:str> <fasta:file_location>")
    if args[0].lower() == 'aachange':
        if len(args) >= 3:
            AAchange(args[1], args[2])
        else:
            return("Required arguments: <snp:location_int-allele_str> <sequence:str>")
    # else:
    #     print("Operation aborted: Function not recognised.")
    #     sys.exit()
    pass

# If being directly executed (ie. not imported)
if __name__ == "__main__":
    if len(sys.argv) <= 1: # ie. if no arguments were passed to biocore
        # Fill this with something useful explaining basic uses of biocore
        print("\nUsage: biocore <command> <arguments>\n\nCommands:\n"
            +"predictMT\tRough melting temperature prediction\n"
            +"findkmers\tFind kmers of given length within a fasta\n"
            +"translate\tTranslate from DNA/RNA to Protein, auto-detects\n"
            +"transcribe\tTranscribe from RNA/DNA to DNA/RNA, auto-detects\n"
            +"complement\tFind the complement of a DNA sequence\n"
            +"calcHamming\tDetermine the Hamming distance between two sequences\n"
            +"simCleave\tSimulate cleavage of a sequence by a given enzyme\n"
            +"simCleaveMulti\tsimCleave for multiple sequences provided as a fasta/fastq\n"
            +"simPCR\t\tPredict PCR fragments of a given sequence and two primers\n"
            +"simPCRMulti\tsimPCR for multiple sequences provided as a fasta/fastq\n"
            +"scaffToContigs\tConvert single scaffold genome to contigs\n"
            +"findMotif\tGiven a motif, find start positions in fasta file or sequence\n"
            +"AAchange\tPredict AA change from SNP and gene sequence\n"
            +"\nNB: This list is incomplete & will be added to later.\n")
        sys.exit()
    else:
        exit(main(sys.argv[1:])) # Call main() with arguments from the command line
