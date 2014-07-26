def countNucs(seq):
	## Returns number of A, C, G & T
	# Add DNA/RNA check
	A = 0
	C = 0
	G = 0
	T = 0
	for x in list(seq):
		if x == "a" or "A":
			A += 1
		elif x == "c" or "C":
			C += 1
		elif x == "g" or "G":
			G += 1
		elif x == "t" or "T":
			T += 1
	print(str(A) + ' ' + str(C) + ' ' + str(G) + ' '+ str(T))

def DNAtoRNA(seq, switch = "toRNA"):
	## Converts DNA string to RNA
	# Add check for RNA
	# Develop this into a DNA to/from RNA two-way convertor with auto-detection and optional switch
	
	seq.upper()
	newSeq = ""

	if "U" in seq:
		switch = "toDNA"
		if "T" in seq:
			print("Invalid sequence, contains both DNA and RNA. Function aborted.")
			return
	
	if switch = "toDNA":
		for x in list(seq):
			if x == "U":
				newSeq += "T"
			else:
				newSeq += str(x)
	elif switch = "toRNA":
		for x in list(seq):
			if x == "T":
				newSeq += "U"
			else:
				newSeq += str(x)
	print(newSeq)

def getComplement(seq):
	## Get reverse complement
	# Reverse given string
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
		# Returns the probability of a dominant positive offspring
		# in a pairing of a population where dominant (k),
		# heterozygous (m) and recessive (n) traits are known

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

def fastaToDict(fasta):
	#given a fasta, returns a dictionary as:
	#strains = {"strainB": "sequenceB", "strainA": "sequenceA", "strainC": "sequenceC"}
	#NB: dictionary entries are not ordered
	
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

def getGC(fasta):
	strains = fastaToDict(fasta) #convert fasta into dictionary
	strainsGC = {}
		for key in strains: #calculate GCs and assign to strains
			GC = 0
			for base in strains[key]:
				if base in ("G", "C"):
					GC += 1
			GCperc = (GC / len(strains[key])) * 100
			strainsGC[key] = GCperc
	return(strainsGC)
	
	
def getHighestGC(fasta):
		#Given a fasta returns the strain with the highest GC content, with that value
		strainsGC = getGC(fasta)

		#get GC
		highestGC = 0
		highestGCstrain = ""
		for key in strainsGC: #find the highest GC, return the strain & value
			if strainsGC[key] > highestGC:
				highestGC = strainsGC[key]
				highestGCstrain = key
		print(str(highestGCstrain))
		print(str(strainsGC[highestGCstrain]))

def calcHamming(A, B):
	#calculate the Hamming distance between two sequences of equal length
	dH = 0
	if len(A) == len(B):
		for index, value in enumerate(A):
			if value != B[index]:
				dH += 1
	print(dH)

def RNAtoPro(seq):
	##Convert RNA to Protein
	#Could include a "check type" function to check if DNA/RNA/Protein
	protein = ""
	end = False
	ticker = 0
	RNAcodons = {"UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L", "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S", "UAU": "Y", "UAC": "Y", "UAA": "-", "UAG": "-", "UGU": "C", "UGC": "C", "UGA": "-", "UGG": "W", "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L", "CCU": "P", "CCC": "P", "CCA": "P", "CCG":"P", "CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q", "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AUU": "I", "AUC": "I", "AUA": "I", "AUG": "M", "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T", "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K", "AGU": "S", "AGC": "S", "AGA": "R", "AGG": "R", "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V", "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A", "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E", "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G"}
	while end == False:
		fc = ticker
		codon = seq[fc:fc+3]
		protein += RNAcodons[codon]
		ticker += 3
		if fc+3 >= len(seq):
			end = True
	print(protein)

def findMotif(motif, seq):
	locations =  []
	mod = 1
	for index, base in enumerate(seq):
		if base == motif[0]:
			while mod < len(motif) and (index + mod) < len(seq):
				if seq[index+mod] == motif[mod]:
					print("Base match found, scanning for whole motif...")
					mod += 1
				else:
					mod = 1 #NB: reset mod on fail
					break
			if mod >= len(motif):
				print("Whole motif found, adding location...")
				print(seq[index:index+mod])
				locations.append(index+1)
				mod = 1
	print("Full sequence scanned, returning locations...")
	print(locations)