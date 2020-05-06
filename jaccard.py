import random
import matplotlib.pyplot as plt


filepath = 'SRR000001.fastq'
my_seqs =[]
sim_idx = [24477, 33654, 99860, 120174, 147319, 284943, 287512, 303812, 313166, 346085]

# Code to read in all sequences from genetic data file
with open(filepath) as fp:
	line = fp.readline()
	cnt = 1
	while line:
   	   	if cnt % 4 == 2:
   	   		my_seqs.append(line.strip())
   	   		#print("Line {}: {}".format(cnt, line.strip()))
   	   	line = fp.readline()
   	   	cnt += 1

# Prints first sequence 
print(my_seqs[0])

"""
Input: two indices and an ngram length
Output: The jaccard similarity using given n gram length for the sequences at given indices
"""
def jaccard(idx0, idx1, ngram):
	set0 = set()
	set1 = set()
	seq0 = my_seqs[idx0]
	seq1 = my_seqs[idx1]
	for i in range(len(seq0)-ngram):
		set0.add(seq0[i:i+ngram])
	for i in range(len(seq1)-ngram):
		set1.add(seq1[i:i+ngram])

	return len(set0.intersection(set1))/len(set0.union(set1))

# Randomly sampling 100 sequences to compare to the 10 test values
random.seed(1998)
controls = random.sample(range(len(my_seqs)), 100)
test_vals = [308597, 319039, 464183, 394249, 439143, 442117, 18615, 153262, 260764, 82983]
control_jaccards = [] # this will have one average jaccard similarity per sequence 
for val in test_vals:
	inner = []
	for control in controls:
		inner.append(jaccard(val, control, 6))
	control_jaccards.append(sum(inner)/len(inner))
print(sum(control_jaccards)/len(control_jaccards))

# Average jaccard similarites between query and results found using LSH (calculated lower down)
lsh_jaccards = [0, 0.1115, 0.1577, 0.0948, 0.0949, 0.1042, 0.0568, 0.1355, 0.1061, 0.0897]
# Average jaccard similarities between query and results found using Levenshtein distance
leven_jaccards = [0.16311541251645245, 0.12537053016486355, 0.16181571365391317, 0.13752382969320054,
				  0.14254201691143042, 0.23696711687257413, 0.07739710534246695, 0.1355598244711999,
				  0.24928708705044508,
				  0.16254122792356118]

# Script to plot 
plt.plot(control_jaccards, label="Randomly Generated")
plt.plot(lsh_jaccards, label="Retrieved Using LSH")
plt.plot(leven_jaccards, label="Retrieved Using Leven Distance")
plt.legend(loc="lower right")
plt.ylabel("Average Jaccard Similarity")
plt.xlabel("Index in List of Test Values")
plt.title("Jaccard Similarity of Retrieved Sequences")
plt.show()

# Index of sequence that was queried
current = 82983
# Indices of returned results
vals = [
455086,
374861,
418330,
465575,
349090,
387229,
402923,
280191,
310101,
408587]

# Calculates average jaccard similarity of queried to returned results
jaccards = []
for val in vals:
	if val is not current: 
		jaccards.append(jaccard(current, val, 6))
print(sum(jaccards)/len(jaccards))
print(len(my_seqs))
