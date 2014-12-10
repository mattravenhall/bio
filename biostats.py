#!/usr/bin/env python3

import sys
import os

def getStats(filename, givenThreshold=500, scaffold=False, returnLens=False):
    """Given a fasta, return the top, bottom and mean lengths.
    Can also optionally return the full list of sequence lengths.
    """
    
    # Stats is growing big enough to be its own file, ideally with each stat type as its own function

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
    # Consider changing all stats in relation to a threshold
    # Could also add mode, median etc.

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

def main(args):
    # args[0] = subfunction to be called, args[1:] = input arguments for args[0].
    # e.g. if 'mt', call predictMT() with arguments
    if len(args) == 1:
        print("Warning: No arguments sent to function.")

    if args[0].lower() == "full":
        if len(args) == 2:
            getStats(args[1]) # passing additional arguments needs to be optimised, current solution is terrible
        if len(args) == 3:
            getStats(args[1], args[2])
        if len(args) == 4:
            getStats(args[1], args[2], args[3])
    if args[0].lower() == "gc":
    	if len(args) == 2:
    		getGC(args[1])
    	if len(args) == 3:
    		getGC(args[1], args[2])
    if args[0].lower() == "topgc":
    	getHighestGC(args[1])
    # else:
    #     print("Operation aborted: Function not recognised.")
    #     sys.exit()
    pass

# If being directly executed (ie. not imported)
if __name__ == "__main__":
    if len(sys.argv) <= 1: # ie. if no arguments were passed to biocore
        # Fill this with something useful explaining basic uses of biocore
        print("\nUsage: biostats <command> <arguments>\n\nCommands:\n"
            +"full\tOverview statistics for fastas\n"
            +"GC\tPercentage GC for given contigs or whole fasta\n"
            +"topGC\tReturns the contig/strain with the highest GC\n")
        	# Add N50, range, contigs etc. as separate functions
        sys.exit()
    else:
        exit(main(sys.argv[1:])) # Call main() with arguments from the command line