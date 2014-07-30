def fastaToDict(fasta):
	"""Given a fasta as its harddrive location, returns a dictionary where keys are the strain name and values are their associated sequences. Note that dictionaries are unordered!"""
	strains = {}
	seqcatch = ""
	f = open(fasta, "r")
	
	for line in f: #take name from > to /n, seq from /n to >
		if ">" in line:
			tmp = line[1:-1]
			strains[tmp] = seqcatch
			seqcatch = ""
		else:
			seqcatch += line.rstrip()
		strains[tmp] = seqcatch
	f.close()
	return(strains)

def findConsensus(fasta):
        
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
                print("Sequences current range from " + str(min(seqLen)) + " to " + str(max(seqLen)) + " bases long.")
                return
        else: #set sequence lengths
                seqLength = max(seqLen)
                
        #create empty matrix
        profile = [[0]*seqLength for i in range(4)]

        #populate matrix
        for value in sequences.values():
                for basePos, baseType in value:
                       profile[baseDict[baseType]][basePos] += 1
                       #NB: base position = y, base type = x

        ##############       
	#make consensus string from profile matrix
        #get max of each column, convert into associated base
        #consensus =
        #use this to iterate over columns within the profile matrix
        #for x in range(len(base)):
        #       print(profile[x][0]) #grabs all first bases of a seq, use [y] to go through a sequence
        #       grab the x for the highest value (0:A, 1:C, 2:G, 3:T)

        #grab a column
        #for value, x in enumerate(profile):
        #       profile[value][COLUMN])
        ##############
        
        for i in range(seqLength):
               colHolder = []
               for value, x in enumerate(profile):
	               colHolder.append(profile[value][i])
	       consensus.append(base[colHolder.index(max(colHolder))])

        #print consensus string
	print(consensus)
                
        #print profile matrix
	for x in profile:
		row = ' '.join(map(str, x))
		print(base[count] + row)
		count += 1
