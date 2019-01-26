"""
use all selected drug-target interactions (both from DrugBank and STITCH) as input
return pathway enrichment results
"""
from scipy.stats import hypergeom
from decimal import Decimal

path_folder = '/Users/fenpei/Box Sync/BalestraWeb/Module_code/KEGG_pathway/'
kegg_path2gene_file = path_folder+'path2gene.txt' # downloaded from http://rest.kegg.jp/list/pathway/hsa
kegg_path2name_file = path_folder + 'path2name.txt' # downloaded from http://rest.kegg.jp/link/hsa/pathway
work_dir = '/Users/fenpei/Box Sync/DA/dt/result/merge_all/'
known_dt = work_dir+'all_known_interactions1.txt'
merge_dt = work_dir+'all_interactions.txt'

path_enrich_known = work_dir + 'path_enrich_known.txt'
path_enrich_merge = work_dir + 'path_enrich_merge.txt'

path2name = {}
with open(kegg_path2name_file,'r') as f:
	for line in f:
		line = line.strip().split('\t')
		path2name[line[0]] = line[1]

path2gene = {}
gene2path = {}
with open(kegg_path2gene_file,'r') as f:
	for line in f:
		line = line.strip().split('\t')
		name = path2name[line[0]]
		gene = line[1]
		path2gene.setdefault(name,set()).add(gene)
		gene2path.setdefault(gene,set()).add(name)
tot_gene = len(gene2path)
print tot_gene

def get_enrich(in_file,out_file):
	tars = set()
	path2tar = {}
	with open(in_file,'r') as f:
		for line in f:
			line = line.strip().split('\t')
			tar = line[4]
			if tar not in tars and tar in gene2path:
				tars.add(tar)
				paths = gene2path[tar]
				for p in paths:
					path2tar.setdefault(p,set()).add(tar)

	tot_tar = len(tars)
	out = open(out_file,'w')
	for k,v in path2tar.items():
		genes = path2gene[k]
		num_gene = len(genes)
		num_tar = len(v)
		pval = hypergeom.sf(num_tar-1, tot_gene, num_gene, tot_tar)
		pval = '%.2E' % Decimal(pval)
		out.write(k+'\t'+str(num_tar)+'\t'+str(pval)+'\n')
	out.close()
get_enrich(known_dt,path_enrich_known)
get_enrich(merge_dt,path_enrich_merge)



