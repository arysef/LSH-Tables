import random


filepath = 'SRR000001.fastq'
my_seqs =[] 
sim_idx = [24477, 33654, 99860, 120174, 147319, 284943, 287512, 303812, 313166, 346085]
with open(filepath) as fp:
	line = fp.readline()
	cnt = 1
	while line:
   	   	if cnt % 4 == 2:
   	   		my_seqs.append(line.strip())
   	   		#print("Line {}: {}".format(cnt, line.strip()))
   	   	line = fp.readline()
   	   	cnt += 1

print(my_seqs[0])

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

print(jaccard(0, 85181, 6))

print(jaccard(0, 444433, 6))

rand = []
for i in range(1, 100):
	rand.append(jaccard(0, i, 6))
	#print(jaccard(0, i, 6))
print(sum(rand)/len(rand))

print("below are calculated ones")
top_k = []
for idx in sim_idx:
	top_k.append(jaccard(0, idx, 6))
	#print(jaccard(0, idx, 6))
print(sum(top_k)/len(top_k))
'''
for idx in sim_idx:
	print(my_seqs[idx])
'''
'''
random.seed(1998)
controls = random.sample(range(len(my_seqs)), 100)
test_vals = [308597, 319039, 464183, 394249, 439143, 442117, 18615, 153262, 260764, 82983]
control_jaccards = []
for val in test_vals:
	for control in controls:
		control_jaccards.append(jaccard(val, control, 6))
print(sum(control_jaccards)/len(control_jaccards))
print(len(control_jaccards))
'''
current = 82983
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
jaccards = []
for val in vals:
	if val is not current: 
		jaccards.append(jaccard(current, val, 6))
print(sum(jaccards)/len(jaccards))
print(len(my_seqs))