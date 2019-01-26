from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem import AllChem
import matplotlib.pyplot as plt
import matplotlib as mpl
import scipy
import scipy.cluster.hierarchy as sch
import scipy.spatial.distance as dist
import numpy
import string
plt.switch_backend('agg')

input1 = '../dt/all_input_compounds.txt'

output1 = '../dt/2d_clustering/input_2D_distance.txt'
output2 = '../dt/2d_clustering/2D_clustering.pdf'
output3 = '../dt/2d_clustering/distance_matrix.txt'
output4 = '../dt/2d_clustering/2D_clusters.txt'

drugs = []
fps = []
d2color = {}
with open(input1,'r') as f:
	next(f)
	for line in f:
		line1 = line.strip().split('\t')
		name = line1[0]
		drugs.append(name)
		smile = line1[4]
		d2color[name] = line1[5]
		m = Chem.MolFromSmiles(smile)
		if m != None:
			fp = AllChem.GetMorganFingerprint(m,2,useFeatures=True)
			fps.append(fp)
		else:
			print name + 'no Fingerprints'
f.close()

n = len(drugs)
x = []
out = open(output1,'w')
for i in range(n-1):
	fp1 = fps[i] 
	for j in range(i+1,n):
		fp2 = fps[j]
		dis = 1- DataStructs.TanimotoSimilarity(fp1,fp2)
		dis = round(dis,4)
		out.write(drugs[i] + '\t' + drugs[j] + '\t' + str(dis)+'\n')
		x.append(dis)
out.close()

D = dist.squareform(x)
numpy.savetxt(output3, D,fmt='%.4f')
Y = sch.linkage(D, method='single', metric='euclidean') ### array-clustering metric - 'average', 'single', 'centroid', 'complete'

fig = plt.figure(figsize = (6,10))
ax = fig.add_axes([0.4, 0.08, 0.55, 0.90],frame_on=True)
Z = sch.dendrogram(Y, labels = drugs, orientation="right",leaf_font_size = 11, distance_sort=True)

# write out clusters
ind = sch.fcluster(Y,0.5*max(Y[:,2]),'distance')
idx = Z['leaves']
ind = ind[idx]
new_drugs = [drugs[i] for i in idx]
out = open(output4,'w')
for i in range(n):
	out.write(str(ind[i]) + '\t'+ new_drugs[i] + '\n')
out.close()

# Create a color palette with 6 color for the 6 drug groups
ax = plt.gca()
xlbls = ax.get_ymajorticklabels()
for lbl in xlbls:
	d = str(lbl).split(',')[2][2:-2]
	val= d2color[d]
	lbl.set_color(val)
print len(xlbls)
plt.savefig(output2)
plt.show()



