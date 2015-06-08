#!/usr/bin/env python3

import sys
import os
import biocore

######################################
# Functions with a statistical focus #
######################################

def getStats(filename, givenThreshold=500, scaffold=False, returnLens=False):
    """Given a fasta or fastq, returns range of bioinformatic metrics.
    Can also optionally return the full list of sequence lengths.
    """
    
    # getStats can probably eventually become a gateway to call multiple stats functions
    # these additional functions could then be called directly with biostats <fasta> calling getStats by default

    #Auto-detect file name
    filetype = biocore.detectType(filename)

    
    threshold = int(givenThreshold) #bp
    x = -1 # This can probably be replaced with something better suited
    contigLengths = []

    f = open(filename, "r")

    # Get list of sequence lengths, different for fasta/fastq
    if filetype == "fasta":
        contigLengths.append(0)
        for line in f:
            if ">" in line:
                # Create new entry in contigLengths ()
                contigLengths.append(0)
                x += 1
            else:
                # Add length of line to value for that contig
                contigLengths[x] += len(line.rstrip('\n'))
        contigLengths = contigLengths[:-1] # Removes excess list entry
    elif filetype == "fastq":
        for line1 in f:
            line2 = next(f)
            line3 = next(f)
            line4 = next(f)
            contigLengths.append(len(line2.rstrip('\n')))

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
    for i in contigLengths: # Find the lengths under the threshold, add to toremove
        if i <= threshold:
            toremove.append(i)
    
    for i in contigLengths: # Add lengths above the threshold to contigThreshold
        if i not in toremove:
            contigThreshold.append(i)

    N50threshold = 0

    # This time for the threshold limit only, but only if the threshold removed contigs
    if len(contigThreshold) == 0:
        N50threshold = "Error: No contigs longer than " + str(threshold) + "bp!"
    elif len(contigThreshold) != len(contigLengths):
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

def getGC(filename, total=False):
    """Given the location of a fasta or fastq, returns the GC value for
    each strain as a dictionary by default. Alternatively, the total GC 
    will be returned if total=True.
    """

    filetype = biocore.detectType(filename)

    Ncount = 0
    if (total == False): # Default: return dict of strains with GCs
        strains = biocore.ToDict(filename) # Convert file into dictionary
        strainsGC = {}
        for key in strains: # Calculate GCs and assign to strains
            GC = 0
            for base in strains[key][0]: # Workhorse
                if base.upper() in ("G", "C"):
                    GC += 1
                if base.upper() in ("N"):
                    Ncount += 1
            GCperc = (float(GC) / float((len(strains[key][0])) - float(Ncount))) * 100
            strainsGC[key] = round(GCperc, 2)
        return(strainsGC)

    elif (total): # Alternative: return the total GC for all sequences
        GC = 0
        totalLength = 0

        genome = biocore.ToDict(filename)
        for key in genome:
            totalLength += len(genome[key][0])
            for base in genome[key][0]:
                if base.upper() in ("G", "C"):
                    GC += 1
                    
        totalGC = (float(GC) / float(totalLength)) * 100
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
        if len(args) == 1:
            print("\nUsage: biostats full [<threshold:int (default:500)>] [<scaffold:boolean (default:False)>] [<return all lengths:boolean (default:False)>]\n")
        if len(args) == 2:
            getStats(args[1]) # passing additional arguments needs to be optimised, current solution is terrible
        if len(args) == 3:
            getStats(args[1], args[2])
        if len(args) == 4:
            getStats(args[1], args[2], args[3])
    if args[0].lower() == "gc":
        if len(args) == 1:
            print("\nUsage: biostats gc <filename:str> [<Return total:boolean (default:False)>]\n")
        if len(args) == 2:
            getGC(args[1])
        if len(args) == 3:
            getGC(args[1], args[2])
    if args[0].lower() == "topgc":
        if len(args) == 1:
            print("\nUsage: biostats topgc <fasta>\n")
        if len(args) == 2:
            getHighestGC(args[1])
    # else:
    #     print("Operation aborted: Function not recognised.")
    #     sys.exit()
    pass

# If being directly executed (ie. not imported)
if __name__ == "__main__":
    if len(sys.argv) <= 1: # ie. if no arguments were passed to biostats
        # Fill this with something useful explaining basic uses of biostats
        print("\nUsage: biostats <command> <arguments>\n\nCommands:\n"
            +"full\tOverview statistics for fastas\n"
            +"GC\tPercentage GC for given contigs or whole fasta\n"
            +"topGC\tReturns the contig/strain with the highest GC\n")
        	# Add N50, range, contigs etc. as separate functions
        sys.exit()
    else:
        exit(main(sys.argv[1:])) # Call main() with arguments from the command line