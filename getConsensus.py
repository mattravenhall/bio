def fastaToMatrix(fasta):
	#given a fasta, returns a matrix as:
	#matrix = [[seqA], [seqB], [seqC]]
	
	strains = []
	strains.append([])
	location = -1
	f = open(fasta, "r")

        #create empty matrix

	#populate matrix
	for line in f: #take name from > to /n, seq from /n to >
	#currently fastaToDict, needs to be ToMatrix
		if ">" in line:
			location += 1
			#add matrix position here? strains[0][location]
		else:
                        #add seq to the previously created position
			strains[location].append(line.rstrip())
		
	f.close()
	return(strains)

        #alt: capture each sequence in a list, then create matrix from that list

def findConsensus(fasta):
	base = ['A', 'C', 'G', 'T']
	count = 0

	matrix = fastaToMatrix(fasta)

	#TEST: print matrix
        #for row in matrix:
                #print(' '.join(row))

	#make profile matrix from matrix
	#make consensus string from profile matrix

	print(consensus)
	for x, y in profile:
		print(base[count], x, y)
		count += 1

