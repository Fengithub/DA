from decimal import Decimal
import numpy as np
path_pval_file = 'pathway_pval.txt'

path_adjust_pval_outfile = 'pathway_adjust_pval.txt'

path2pval_known = {}
path2pval_merge = {}
with open(path_pval_file,'r') as f:
	next(f)
	for line in f:
		line = line.strip().split('\t')
		if line[4]!='-':
			path2pval_known[line[1]] = float(line[4])
		path2pval_merge[line[1]] = float(line[5])

def adjust_pval_bh(path2pval_in):

	m = len(path2pval_in)
	keys = []
	vals = []
	path2pval_out = {}
	for k,v in sorted(path2pval_in.items(),key = lambda x:x[1]):
		keys.append(k)
		vals.append(v)
	pvals = np.array(vals)
	adjust_pvals_temp = pvals*m/np.arange(1, m+1)
	adjust_pvals = [np.min([p, adjust_pvals_temp[i+1],1]) if i < m-1 else p for i, p in enumerate(adjust_pvals_temp)]
	for i in range(m):
		path2pval_out[keys[i]] = adjust_pvals[i]
	return path2pval_out

path2adjust_pval_known = adjust_pval_bh(path2pval_known)
path2adjust_pval_merge = adjust_pval_bh(path2pval_merge)

out = open(path_adjust_pval_outfile,'w')
out.write('Pathway index\tPathway name\tNO. of known targets\tNO. of merged targets\tP-value (known)\tAdjusted_P-value (known)\tP-value (merged)\tAdjusted_P-value (merged)\n')
with open(path_pval_file,'r') as f:
	next(f)
	for line in f:
		line = line.strip().split('\t')
		path = line[1]
		if path in path2adjust_pval_known:
			adjust_pval_known = '%.2E' % Decimal(path2adjust_pval_known[path])
		else:
			adjust_pval_known = '-'
		adjust_pval_merge = '%.2E' % Decimal(path2adjust_pval_merge[path])
		out.write('\t'.join(line[0:5])+'\t'+str(adjust_pval_known)+'\t'+line[5]+'\t'+str(adjust_pval_merge)+'\n')
out.close()