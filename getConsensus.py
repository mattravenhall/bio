def fastaToMatrix(fasta):
	#given a fasta, returns a matrix as:
	#matrix = [[seqA], [seqB], [seqC]]
	
	seqcatch = ""
	location = -1
	f = open(fasta, "r")
	
	for line in f: #take name from > to /n, seq from /n to >
	#currently fastaToDict, needs to be ToMatrix
		if ">" in line:
			location += 1
			seqcatch = ""
		else:
			seqcatch += line.rstrip()
		strains[0][location] = seqcatch
	f.close()
	return(strains)	

def findConsensus(fasta):
	base = ['A', 'C', 'G', 'T']
	count = 0

	matrix = fastaToMatrix(fasta)

	#make profile matrix from matrix
	#make consensus string from profile matrix

	print(consensus)
	for x, y in profile:
		print(base[count], x, y)
		count += 1
