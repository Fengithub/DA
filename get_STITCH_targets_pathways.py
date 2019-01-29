""" 
	input a list of chemicals
	return the known and predicted chemical-protein interactions, enriched targets, enriched pathways, enriched GO terms
	note: the chemical-protein interactions of human experimental dataset from STITCH are stored in our server database, predictions made by PMF are precalculated and stored on the server.
"""
import time
start_time = time.time()
import os
import psycopg2
import sys
from scipy.stats import hypergeom
import matplotlib.pyplot as plt
import numpy as np
plt.switch_backend('agg')

try:
	conn = psycopg2.connect(database="database", user="username", password="password", host="---", port="---")
except:
	out_log.write("I am unable to connect to the database\n")

cur = conn.cursor()

folder = sys.argv[1]
folder = folder.strip('/')
dataset = sys.argv[2] # stitch_v5
query_type = sys.argv[3] # drugs or proteins
max_pred = int(sys.argv[4])
query_file = sys.argv[5] # input file, including a list of chemical name/DB_id/pubchem id or protein name/Uniprot id/Protein Entry
#similarity cutoff
# table names in database
prefix = 'fen_'
tbl_chem_index = prefix + 'stitch_v5_chem_index'
tbl_pro_index = prefix + dataset + '_pro_index'
tbl_dt_known = prefix + dataset + '_known_d2t_'
tbl_dt_predict = prefix + dataset + '_prediction_d2t_'
tbl_td_known = prefix + dataset + '_known_t2d_'
tbl_td_predict = prefix + dataset + '_prediction_t2d_'
tbl_dp_known = prefix + dataset + '_chem2path'
tbl_path_index = prefix + dataset + '_path_index'

if not os.path.exists(folder):
	os.makedirs(folder)

max_known = 200
cutoff_score_known = 0.4 # cutoff score for known interactions
cutoff_score_pred = 0.7 # cutoff score for predicted interactions

def split_drugs(drug_list):
	table_num2drugs = {}
	for d in drug_list:
		table_num = str(d/60000 + 1)
		if table_num == '6':
			table_num = '5'
		table_num2drugs.setdefault(table_num,[]).append(d)
	return table_num2drugs

def split_tars(tar_list):
	table_num2tars = {}
	for t in tar_list:
		table_num = str(t/5000 + 1)
		table_num2tars.setdefault(table_num,[]).append(t)
	return table_num2tars

def merge_dict(dict1, dict2):
	dict3 = {}
	for k,v in dict1.items():
		dict3[k] = v
	for k,v in dict2.items():
		if k in dict3:
			dict3[k] = dict3[k]|v
		else:
			dict3[k] = v
	return dict3

def enrich_plot(pval_dict, fig_file, item_type):
	""" takes the {item: enrichment score} dictionary, output figure filename, item type (drug, target, pathway, GO terms)
		returns a bar plot of enrichment score of the items
	"""
	sorted_pval_dict = sorted(pval_dict.items(), key = lambda x: x[1], reverse = True)
	items = [x[0] for x in sorted_pval_dict]
	values = [x[1] for x in sorted_pval_dict]  ## enrichment score

	axes1 = [0.60, 0.08, 0.35, 0.90] # for target,pathway and GO terms
	axes2 = [0.35, 0.08, 0.60, 0.90] # for drug
	plt.rcdefaults()
	fig = plt.figure(figsize = (10,8))
	if item_type == 'target' or item_type == 'pathway':
		ax = fig.add_axes(axes1,frame_on=True)
	else:
		ax = fig.add_axes(axes2,frame_on=True)
	y_pos = np.arange(len(items))
	width = 0.4

	rec1 = ax.barh(y_pos, values, width, color = 'gray', alpha = 0.65, edgecolor='k')
	ax.set_ylim([-0.25,19.75])
	ax.set_yticks(y_pos+width/2)
	ax.set_yticklabels(items)
	ax.invert_yaxis()  # labels read top-to-bottom
	ax.set_xlabel('Enrichment Score [-log10(p-value)]')
	plt.savefig(fig_file)

def find_drug(query,ind2drug):
	tbl_chem_index = prefix + 'stitch_v5_chem_index'
	if query.isdigit(): # case1: query pubchem id
		query_id = "\'"+query+"\'"
		cur.execute("select chem_index, chem_id, chem_name, MW, SMILES, num_pros from "+ tbl_chem_index + " where chem_id="+query_id+";")
	else: # case2: query drug_name
		query_like = "\'%"+query+"%\'"
		cur.execute("select chem_index, chem_id, chem_name, MW, SMILES, num_pros from "+ tbl_chem_index + " where chem_name ILIKE"+query_like+";")
	drug = cur.fetchall()
	if drug == []:
		out_log.write('cannot find '+query+' in STITCH_v5\n')
		query_index = None
	else:
		drug = drug[0]
		query_index = drug[0]
		ind2drug[query_index] = drug[1:]
	return ind2drug

def find_target(query, ind2pro):
	query_id = "\'"+query+"\'"
	cur.execute("select pro_index, uniprot_id, protein_name, entry_name, num_chemicals, gene_name, gene_id, PDB, path_ids, path_names, go_function, go_process, go_component,path_inds from "+ tbl_pro_index + " where uniprot_id="+query_id+";")
	pro = cur.fetchall()
	if pro == []:
		cur.execute("select pro_index,uniprot_id, protein_name, entry_name, num_chemicals, gene_name, gene_id, PDB, path_ids, path_names, go_function, go_process, go_component,path_inds from "+ tbl_pro_index + " where gene_name="+query_id+";")
		pro = cur.fetchall()
		if pro == []:
			query_like = "\'%"+query+"%\'"
			cur.execute("select pro_index,uniprot_id, protein_name, entry_name, num_chemicals, gene_name, gene_id, PDB, path_ids, path_names, go_function, go_process, go_component,path_inds from "+ tbl_pro_index + " where protein_name ILIKE "+query_like+";")
			pro = cur.fetchall()
	if pro == []:
		out_log.write('cannot find '+query+' in ' + dataset+'\n')
		query_index = None
	else:
		pro = pro[0]
		query_index = pro[0]
		ind2pro[query_index] = pro[1:]
	return ind2pro

def get_drugs(drug_ind_lst, ind2drug):
	tbl_chem_index = prefix + 'stitch_v5_chem_index'
	cur.execute("select chem_index, chem_id, chem_name, MW, SMILES, num_pros from "+ tbl_chem_index + " where chem_index in "+drug_ind_lst+";")
	rows = cur.fetchall()
	for row in rows:
		#row = ['NA' if x==None else x for x in row]
		ind2drug[row[0]] = row[1:]
	return ind2drug

def get_pros(pro_ind_lst, ind2pro):
	cur.execute("select pro_index, uniprot_id, protein_name, entry_name, num_chemicals, gene_name, gene_id, PDB, path_ids, path_names, go_function, go_process, go_component,path_inds from "+ tbl_pro_index + " where pro_index in "+pro_ind_lst+";")
	rows = cur.fetchall()
	for row in rows:
		ind2pro[row[0]] = row[1:]
	return ind2pro

def get_paths(path_ind_lst):
	ind2path = {}
	cur.execute("select path_index, path_id, path_name, num_genes, class_1st, class_2nd, num_chemicals, num_proteins from " +tbl_path_index+ " where path_index in "+path_ind_lst+";")
	rows = cur.fetchall()    
	for row in rows:
		ind2path[str(row[0])] = row[1:]
	return ind2path

def get_d2t_known(query_index_lst, table_num, dt, td):
	"""
	take chemical indexes as input, search known targets for these chemicals
	return drug2target, target2drug dicts
	"""
	tbl_dt = tbl_dt_known + table_num 
	cur.execute("select chem_index, pro_index_list, score_list from "+ tbl_dt +" where chem_index in "+query_index_lst+";")
	rows = cur.fetchall()
	if rows != []:
		for chem_index, pro_index_list, score_list in rows:
			pro_index_list = pro_index_list.split(';')
			score_list = score_list.split(';')
			for i in range(min(len(pro_index_list), max_known)):
				pro_index = int(pro_index_list[i])
				interaction = float(score_list[i])
				if interaction >= cutoff_score_known:
					dt.setdefault(chem_index,set()).add((pro_index,interaction))
					td.setdefault(pro_index,set()).add((chem_index,interaction))
	return [dt, td]

def get_d2t_predict(query_index_lst, table_num, dt, td):
	"""
	take chemical indexes as input, search predicted targets for these chemicals
	return drug2target, target2drug dicts
	"""
	tbl_dt = tbl_dt_predict + table_num
	cur.execute("select chem_index, pro_index_list, score_list from "+ tbl_dt +" where chem_index in "+query_index_lst+";")
	rows = cur.fetchall()
	if rows != []:
		for chem_index, pro_index_list, score_list in rows:
			pro_index_list = pro_index_list.split(';')
			score_list = score_list.split(';')
			for i in range(min(len(pro_index_list), max_pred)):
				pro_index = int(pro_index_list[i])
				interaction = float(score_list[i])
				if interaction >= cutoff_score_pred:
					dt.setdefault(chem_index,set()).add((pro_index,interaction))
					td.setdefault(pro_index,set()).add((chem_index,interaction))
	return [dt, td]

def get_simi_drug_pairs(query_index_lst, query_index_set, ind2drug, table_num, simi_in, simi_out):
	simi_index_set = set()
	tbl_D2D = prefix + dataset + '_d2d_list_' + table_num
	cur.execute("select chem1_index, chem2_index_list, score_list from "+ tbl_D2D + " where chem1_index in "+query_index_set+";")
	rows = cur.fetchall()
	for chem1_index, chem2_index_list, score_list in rows:
		chem2_index_list = chem2_index_list.split(';')
		score_list = score_list.split(';')
		for i in range(min(max_pred,len(chem2_index_list))):
			chem2_index = int(chem2_index_list[i])
			similarities = score_list[i].split(',')
			similarity = float(similarities[0])
			similarity_2D = similarities[1]
			if chem2_index in query_index_lst:				
				simi_in.setdefault(chem1_index,set()).add((chem2_index,similarity,similarity_2D))
			else:
				simi_index_set.add(chem2_index)
				simi_out.setdefault(chem1_index,set()).add((chem2_index,similarity,similarity_2D))
	if simi_index_set != set():
		simi_index_lst = list(simi_index_set)
		simi_index_set = '(' + str(simi_index_lst)[1:-1] +')'
		ind2drug = get_drugs(simi_index_set, ind2drug)
	return [ind2drug, simi_in, simi_out]

def get_simi_pro_pairs(query_index_lst, query_index_set, ind2pro, table_num, simi_in, simi_out):
	simi_index_set = set()
	tbl_T2T = prefix + dataset + '_t2t_list_' + table_num
	cur.execute("select pro1_index, pro2_index_list, score_list from "+ tbl_T2T + " where pro1_index in "+query_index_set+";")
	rows = cur.fetchall()
	for pro1_index, pro2_index_list, score_list in rows:
		pro2_index_list = pro2_index_list.split(';')
		score_list = score_list.split(';')
		for i in range(min(max_pred,len(pro2_index_list))):
			pro2_index = int(pro2_index_list[i])
			similarity = float(score_list[i])
			if pro2_index in query_index_lst:		
				simi_in.setdefault(pro1_index,set()).add((pro2_index,similarity))
			else:
				simi_index_set.add(pro2_index)
				simi_out.setdefault(pro1_index,set()).add((pro2_index,similarity))
	if simi_index_set != set():
		simi_index_lst = list(simi_index_set)
		simi_index_set = '(' + str(simi_index_lst)[1:-1] +')'
		ind2pro = get_pros(simi_index_set, ind2pro)
	return [ind2pro, simi_in, simi_out, simi_index_lst]

def get_name2go(go_type):
	"""
	take GO term to targets dictionary and GO term type as input
	search the database for the total number of proteins and drugs that are associated with that term
	return a name2go dictionary: {GO term name: total number of proteins, total number of drugs}, used to calculated enrichment p-value
	"""
	#go_lst = '('+str(go2p.keys())[1:-1]+')'
	tbl_go_index = prefix +dataset+ '_go_' + go_type
	name2go = {}
	#cur.execute("select go_name, num_pros, num_chemicals from "+ tbl_go_index+" where go_name in "+go_lst +";")
	cur.execute("select go_name, num_pros, num_chemicals from "+ tbl_go_index+";")	
	rows = cur.fetchall()
	for go_name, num_pros, num_chemicals in rows:
		name2go[go_name] = (num_pros, num_chemicals)
	return name2go

def get_pro2path_go(query_index_lst, ind2pro):
	p2t = {}
	t2p = {}
	go_function2t = {}
	t2go_function = {}
	go_process2t = {}
	t2go_process = {}
	go_component2t = {}
	t2go_component = {}	
	for pro_index in query_index_lst:
		v = ind2pro[pro_index]
		genename = v[4]
		gene_id = v[5]
		path_inds = v[-1]
		go_function = v[-4]
		go_process = v[-3]
		go_component = v[-2]
		if path_inds != 'NA':
			path_inds = path_inds.split(';')
			for path_index in path_inds:
				p2t.setdefault(path_index,set()).add((pro_index,genename,gene_id))
				t2p.setdefault(pro_index,set()).add(path_index)
		if go_function != 'NA':
			go_function = go_function.split('|')
			for go in go_function:
				go_function2t.setdefault(go,set()).add((pro_index,genename,gene_id))
				t2go_function.setdefault(pro_index,set()).add(go)
		if go_process != 'NA':
			go_process = go_process.split('|')
			for go in go_process:
				go_process2t.setdefault(go,set()).add((pro_index,genename,gene_id))
				t2go_process.setdefault(pro_index,set()).add(go)
		if go_component != 'NA':
			go_component = go_component.split('|')
			for go in go_component:
				go_component2t.setdefault(go,set()).add((pro_index,genename,gene_id))
				t2go_component.setdefault(pro_index,set()).add(go)
	return [p2t, t2p, go_function2t,t2go_function,go_process2t,t2go_process,go_component2t,t2go_component]

def write_drug_pairs(out_file, pair_dict, ind2drug):
	out = open(out_file,'w')
	out.write('Chemical1_ID\tChemical1_Name\tChemical1_MW\tChemical1_SMILES\tChemical2_ID\tChemical2_Name\tChemical2_MW\tChemical2_SMILES\tSimilarity_DT\tSimilarity_2D\n')
	for k,v in pair_dict.items():
		a = ind2drug[k]
		for b in v:
			d2 = b[0]
			if d2 in ind2drug:
				c = ind2drug[d2]
				out.write('\t'.join(x for x in a[0:-1]) + '\t' + '\t'.join(x for x in c[0:-1]) +'\t'+ str(b[1])+'\t'+ b[2]+ '\n')
	out.close()

def write_pro_pairs(out_file, pair_dict, ind2pro):
	out = open(out_file,'w')
	out.write('Gene1_Name\tGene2_Name\tSimilarity_DT\n')	
	for k,v in pair_dict.items():
		a = ind2pro[k][4]
		for b in v:
			p2 = b[0]
			if p2 in ind2pro:
				c = ind2pro[b[0]][4]
				out.write(a + '\t' + c +'\t'+ str(b[1]) + '\n')
	out.close()

def write_pros(out_file, pro_lst, ind2pro): 
	out = open(out_file,'w')
	out.write('Uniprot_ID\tProtein_Name\tEntry_Name\tTot_Num_Chemicals\tGene_Name\tGene_ID\tPDB\tPathways_ID\tPathways_Name\tGO_Function\tGO_Process\tGO_Component\n')
	for one in pro_lst:
		pro = ind2pro[one]
		out.write('\t'.join(str(x) for x in pro[0:-1]) + '\n')
	out.close()

def write_dt(out_file, dt_dict,ind2drug,ind2pro):
	out = open(out_file,'w')
	out.write('Pubchem_ID\tChemical_Name\tMolecular_Weight\tSMILES\tUniprot_ID\tProtein_Name\tEntry_Name\tTot_Num_Chemicals\tGene_Name\tGene_ID\tPDB\tPathways_ID\tPathways_Name\tGO_Function\tGO_Process\tGO_Component\tConfidnece_Score\n')
	for k,v in dt_dict.items():
		drug = ind2drug[k]
		v = sorted(list(v), key = lambda x:x[1], reverse= True)
		for one in v:
			pro_ind = one[0]
			pro = ind2pro[pro_ind]
			score = one[1]
			out.write('\t'.join(str(x) for x in drug[0:-1]) + '\t' +'\t'.join(str(x) for x in pro[0:-1])+'\t'+ str(score) + '\n')
	out.close()

def write_td(out_file, td_dict, ind2drug, ind2pro):
	out = open(out_file,'w')
	out.write('Uniprot_ID\tProtein_Name\tEntry_Name\tTot_Num_Chemicals\tGene_Name\tGene_ID\tPDB\tPathways_ID\tPathways_Name\tGO_Function\tGO_Process\tGO_Component\tPubchem_ID\tChemical_Name\tMolecular_Weight\tSMILES\tConfidnece_Score\n')
	for k,v in td_dict.items():
		pro = ind2pro[k]
		v = sorted(list(v), key = lambda x:x[1], reverse= True)
		for one in v:
			drug_ind = one[0]
			drug = ind2drug[drug_ind]
			score = one[1]
			out.write('\t'.join(str(x) for x in pro[0:-1]) + '\t' +'\t'.join(str(x) for x in drug[0:-1])+'\t'+ str(score) + '\n')
	out.close()

def write_tar_rank(tbl_file, fig_file, pro2drug, ind2pro, M, N): # write targets ranked by enrichment p-value
	pro2pval = {}
	for k, v in pro2drug.items():
		pro = ind2pro[k]
		q = pro[3]
		#m = len(v)
		m = sum([x[1] for x in v])
		if m > q: # for predicted case, we keep use the total number of known chemicals interacting with a certain target, if predicted number greater than known, pval is 0.
			pval = hypergeom.sf(m-1, M, m, N)
		else:
			pval = hypergeom.sf(m-1, M, q, N) # Hypergeometric test: the probability of getting more than (m-1) items from N, when the backgroud is q in M.
		pro2pval[k] = pval
	
	n = 1
	out = open(tbl_file,'w')
	out.write('Rank\tUniprot_ID\tProtein_Name\tEntry_Name\tTot_Num_Chemicals\tGene_Name\tGene_ID\tPDB\tPathway_Ids\tPathway_Names\tGO_Function\tGO_Process\tGO_Component\tChemicals\tNum_Chemicals\tP_value\n')
	pro_name2enrich_score = {}
	for k, v in sorted(pro2pval.items(),key=lambda x:x[1]):
		pro = ind2pro[k]
		drugs = [(ind2drug[x[0]][1],x[1]) for x in pro2drug[k]]
		out.write(str(n)+'\t'+'\t'.join(str(x) for x in pro[0:-1])+'\t'+';'.join(str(x) for x in drugs)+'\t'+str(len(drugs))+'\t'+str(v)+'\n')
		if n <= 20:
			pro_name2enrich_score[pro[1]] = -np.log10(v)
		n += 1
	out.close()
	enrich_plot(pro_name2enrich_score, fig_file, 'target')

def write_drug_rank(tbl_file, fig_file, drug2pro, ind2drug, ind2pro, M, N): # write targets ranked by enrichment p-value
	drug2pval = {}
	for k, v in drug2pro.items():
		drug = ind2drug[k]
		q = drug[4]
		#m = len(v)
		m = sum([x[1] for x in v])
		if m > q: # for predicted case, we keep use the total number of known chemicals interacting with a certain target, if predicted number greater than known, pval is 0.
			pval = hypergeom.sf(m-1, M, m, N)
		else:
			pval = hypergeom.sf(m-1, M, q, N) # Hypergeometric test: the probability of getting more than (m-1) items from N, when the backgroud is q in M.
		drug2pval[k] = pval
	
	n = 1
	out = open(tbl_file,'w')
	out.write('Rank\tPubchem_ID\tChemical_Name\tMolecular_Weight\tSMILES\tProteins\tNum_Targets\tP_value\n')
	drug_name2erich_score = {} ## in order to draw the enrichment plot
	for k, v in sorted(drug2pval.items(),key=lambda x:x[1]):
		drug = ind2drug[k]
		pros = [(ind2pro[x[0]][4],x[1]) for x in drug2pro[k]]
		out.write(str(n)+'\t'+'\t'.join(str(x) for x in drug[0:-1])+'\t'+';'.join(str(x) for x in pros)+'\t'+str(len(pros))+'\t'+str(v)+ '\n')
		if n <= 20:
			name = drug[1]
			drug_name2erich_score[name] = -np.log10(v)
		n += 1
	out.close()
	enrich_plot(drug_name2erich_score, fig_file, 'drug')
	
def write_path_rank(query_type, tbl_file, fig_file_t, fig_file_d, path2pro, path2drug, ind2path, N_t, N_d):
	"""
	take the output table file name, output figure file names, pathway-protein dictionary: {pathway_index: [protein index]}, pathway-drug dictionary: {pathway index: [drug]}
	ind2path, the total number of targets (N_t) based on the query, and the total number of input drugs (N_d) as input 
	write the pathways ranked by the enrichment pvalue based on associated targets
	plot the barplots of pathways with the enrichment pvalue
	"""
	path2pval_t = {}
	path2pval_d = {}
	path_name2enrich_score_t = {}
	path_name2enrich_score_d = {}
	for k, v in path2pro.items():
		m_t = len(v)
		q_t = ind2path[k][2]
		if m_t > q_t:
			pval = hypergeom.sf(m_t-1, M_t, m_t, N_t)
		else:
			pval = hypergeom.sf(m_t-1, M_t, q_t, N_t)
		path2pval_t[k] = pval

	if query_type == 'drugs':
		path2pval_d = {}
		for k, v in path2drug.items():
			m_d = len(v)
			q_d = ind2path[k][5]
			if m_d > q_d:
				pval = hypergeom.sf(m_d-1, M_d, m_d, N_d)
			else:
				pval = hypergeom.sf(m_d-1, M_d, q_d, N_d)
			path2pval_d[k] = pval

		n = 1
		out = open(tbl_file,'w')
		out.write('Rank\tPathway_ID\tPathway_Name\tTotal_Gene_Number\tPathway_Class_1st\tPathway_Class_2nd\tGene_Names\tGene_IDs\tNum_Targets\tChemicals\tNum_Chemicals\tP_value_target\tP_value_drug\n')
		for k, v in sorted(path2pval_t.items(),key=lambda x:x[1]):
			path = ind2path[k][0:-2]
			tars = path2pro[k]
			genenames = [x[1] for x in tars]
			genes = [x[2] for x in tars]
			drugs = path2drug[k]
			pval_d = path2pval_d[k]
			out.write(str(n)+'\t'+'\t'.join(str(x) for x in path)+'\t'+';'.join(genenames)+'\t'+';'.join(genes)+'\t'+str(len(genes))+'\t'+';'.join(drugs)+'\t'+str(len(drugs))+'\t'+ str(v)+'\t' + str(pval_d) +'\n')
			if n <= 20:
				path_name2enrich_score_t[path[1]] = -np.log10(v)
			n += 1
		out.close()

		n = 1
		for k,v in sorted(path2pval_d.items(),key=lambda x:x[1]):
			path = ind2path[k]
			if n <= 20:
				path_name2enrich_score_d[path[1]] = -np.log10(v)
			n += 1
		enrich_plot(path_name2enrich_score_d, fig_file_d, 'pathway')

	elif query_type == 'targets':
		n = 1
		out = open(tbl_file,'w')
		out.write('Rank\tPathway_ID\tPathway_Name\tTotal_Gene_Number\tPathway_Class_1st\tPathway_Class_2nd\tGene_Names\tGene_IDs\tNum_Targets\tP_value\n')
		for k, v in sorted(path2pval_t.items(),key=lambda x:x[1]):
			path = ind2path[k][0:-2]
			tars = path2pro[k]
			genenames = [x[1] for x in tars]
			genes = [x[2] for x in tars]
			out.write(str(n)+'\t'+'\t'.join(str(x) for x in path)+'\t'+';'.join(genenames)+'\t'+';'.join(genes)+'\t'+str(len(genes))+'\t'+ str(v)+'\n')
			if n <= 20:
				path_name2enrich_score_t[path[1]] = -np.log10(v)
			n += 1
		out.close()
	enrich_plot(path_name2enrich_score_t, fig_file_t, 'pathway')

def write_go_rank(query_type, tbl_file, fig_file_t, fig_file_d, go2pro, go2drug, name2go, N_t, N_d):
	"""
	take the query GO term type, output table file name, output figure file names, GO-protein dictionary: {GO term name: [protein index]}, GO term-drug dictionary: {GO term index: [drug]}
	name2go, the total number of targets (N_t) based on the query, and the total number of input drugs (N_d) as input 
	write the GO terms ranked by the enrichment pvalue based on associated targets
	plot the barplots of GO terms with the enrichment pvalue
	"""
	go2pval_t = {}
	go2pval_d = {}
	go_name2enrich_score_t = {}
	go_name2enrich_score_d = {}
	for k, v in go2pro.items():
		m_t = len(v)
		q_t = name2go[k][0]
		if m_t > q_t:
			pval = hypergeom.sf(m_t-1, M_t, m_t, N_t)
		else:
			pval = hypergeom.sf(m_t-1, M_t, q_t, N_t)
		go2pval_t[k] = pval

	if query_type == 'drugs':
		go2pval_d = {}
		for k, v in go2drug.items():
			m_d = len(v)
			q_d = name2go[k][1]
			if m_d > q_d:
				pval = hypergeom.sf(m_d-1, M_d, m_d, N_d)
			else:
				pval = hypergeom.sf(m_d-1, M_d, q_d, N_d)
			go2pval_d[k] = pval

		n = 1
		out = open(tbl_file,'w')
		out.write('Rank\tGo_Name\tGene_Names\tGene_IDs\tNum_Targets\tChemicals\tNum_Chemicals\tP_value_target\tP_value_drug\n')
		for k, v in sorted(go2pval_t.items(),key=lambda x:x[1]):
			tars = go2pro[k]
			genenames = [x[1] for x in tars]
			genes = [x[2] for x in tars]
			drugs = go2drug[k]
			pval_d = go2pval_d[k]
			out.write(str(n)+'\t'+k+'\t'+';'.join(genenames)+'\t'+';'.join(genes)+'\t'+str(len(genes))+'\t'+';'.join(x for x in drugs)+'\t'+str(len(drugs))+'\t'+ str(v)+'\t' + str(pval_d) +'\n')
			if n <= 20:
				go_name2enrich_score_t[k] = -np.log10(v)
			n += 1
		out.close()

		n = 1
		for k,v in sorted(go2pval_d.items(),key=lambda x:x[1]):
			if n <= 20:
				go_name2enrich_score_d[k] = -np.log10(v)
			n += 1
		enrich_plot(go_name2enrich_score_d, fig_file_d, 'pathway')

	elif query_type == 'targets':
		n = 1
		out = open(tbl_file,'w')
		out.write('Rank\tGo_Name\tGene_Names\tGene_IDs\tNum_Targets\tP_value\n')
		for k, v in sorted(go2pval_t.items(),key=lambda x:x[1]):
			tars = go2pro[k]
			genenames = [x[1] for x in tars]
			genes = [x[2] for x in tars]
			out.write(str(n)+'\t'+k+'\t'+';'.join(genenames)+'\t'+';'.join(genes)+'\t'+str(len(genes))+'\t'+ str(v)+'\n')
			if n <= 20:
				go_name2enrich_score_t[k] = -np.log10(v)
			n += 1
		out.close()
	enrich_plot(go_name2enrich_score_t, fig_file_t, 'pathway')

def drug_gene_path_go(tar2drug, ind2drug):
	"""
	takes a dictionary of targets to drugs, and a dictionary of drug information
	get the drug-pathway, drug-GO terms relationships
	return 12 dictionaries as follows:
		pro2path: {protein index: pathway indexes}
		path2pro: {pathway index: (protein index, gene name, gene id)}
		p2go_function: {protein index: GO function terms}
		p2go_process: {protein index: GO process terms}
		p2go_component: {protein index: GO component terms}
		go_function2p: {GO function term: (protein index, gene name, gene id)}
		go_process2p: {GO process term: (protein index, gene name, gene id)}
		go_component2p: {GO component term: (protein index, gene name, gene id)}
		path2drug: {pathway index: drug names}
		go_function2d: {GO function term: drug names}
		go_process2d: {GO process term: drug names}
		go_component2d: {GO component term: drug names}
	"""		
	path2drug = {}
	go_function2d = {}
	go_process2d = {}
	go_component2d = {}
	pros_lst = tar2drug.keys()
	[path2pro, pro2path,  go_function2p, p2go_function, go_process2p, p2go_process, go_component2p, p2go_component] = get_pro2path_go(pros_lst, ind2pro)
	for k,v in path2pro.items():
		drugs = []
		for each in v:
			drugs = [ind2drug[x[0]][1] for x in tar2drug[each[0]]] + drugs
		path2drug[k] = set(drugs)
	for k,v in go_function2p.items():
		drugs = []
		for each in v:
			drugs = [ind2drug[x[0]][1] for x in tar2drug[each[0]]] + drugs
		go_function2d[k] = set(drugs)
	for k,v in go_process2p.items():
		drugs = []
		for each in v:
			drugs = [ind2drug[x[0]][1] for x in tar2drug[each[0]]] + drugs
		go_process2d[k] = set(drugs)
	for k,v in go_component2p.items():
		drugs = []
		for each in v:
			drugs = [ind2drug[x[0]][1] for x in tar2drug[each[0]]] + drugs
		go_component2d[k] = set(drugs)		
	return [pro2path, path2pro, p2go_function, p2go_process, p2go_component, go_function2p, go_process2p, go_component2p, path2drug, go_function2d, go_process2d, go_component2d]

M_t = 9457 # total number of chemicals
M_d = 315514 # total number of targets
M_g = 7314 # total number of human genes in KEGG

out_log = open(folder + '/log.txt','w')

if query_type == 'drugs': # write all result files when input drugs
	ind2drug = {}
	with open(query_file,'r') as f:
		for line in f:
			query = line.strip()
			if len(query) > 0:
				ind2drug = find_drug(query, ind2drug)

	if ind2drug == {}:
		out_log.write('cannot find any queried drugs in ' + dataset + '\n')
		exit()
	else:
		query_chemicals_info = folder + '/input_chemicals_details.txt'
		input_simi_outfile = folder + '/input_similarity.txt'  # only similarity within the top 64 is presented
		simi_chemicals_outfile = folder + '/similar_chemical_pairs.txt'
		## drug-target interactions ##
		known_dt_outfile = folder + '/known_interactions.txt'
		pred_dt_outfile = folder + '/predicted_interactions.txt'
		## targets ##
		known_tars_outfile = folder + '/known_targets_rank.txt'
		merge_tars_outfile = folder + '/merge_targets_rank.txt'
		known_tar_enrich_plot = folder + '/known_target_enrichment_plot.png'
		merge_tar_enrich_plot = folder + '/merge_target_enrichment_plot.png'
		## pathways ##
		known_paths_outfile = folder + '/known_pathways_rank.txt'            
		merge_paths_outfile = folder + '/merge_pathways_rank.txt'
		gene_paths_outfile = folder + '/gene_pathway_relationship.txt'
		known_path_enrich_d_plot = folder + '/known_pathway_enrichment_d_plot.png'
		merge_path_enrich_d_plot = folder + '/merge_pathway_enrichment_d_plot.png'
		known_path_enrich_t_plot = folder + '/known_pathway_enrichment_t_plot.png'
		merge_path_enrich_t_plot = folder + '/merge_pathway_enrichment_t_plot.png'
		## GO terms ##
		known_go_function_outfile = folder + '/known_GO_function_rank.txt'            
		merge_go_function_outfile = folder + '/merge_Go_function_rank.txt'
		known_go_function_enrich_d_plot = folder + '/known_Go_function_enrichment_d_plot.png'
		merge_go_function_enrich_d_plot = folder + '/merge_Go_function_enrichment_d_plot.png'
		known_go_function_enrich_t_plot = folder + '/known_Go_function_enrichment_t_plot.png'
		merge_go_function_enrich_t_plot = folder + '/merge_Go_function_enrichment_t_plot.png'
		known_go_process_outfile = folder + '/known_GO_process_rank.txt'            
		merge_go_process_outfile = folder + '/merge_GO_process_rank.txt'
		known_go_process_enrich_d_plot = folder + '/known_GO_process_enrichment_d_plot.png'
		merge_go_process_enrich_d_plot = folder + '/merge_GO_process_enrichment_d_plot.png'
		known_go_process_enrich_t_plot = folder + '/known_GO_process_enrichment_t_plot.png'
		merge_go_process_enrich_t_plot = folder + '/merge_GO_process_enrichment_t_plot.png'
		known_go_component_outfile = folder + '/known_GO_component_rank.txt'            
		merge_go_component_outfile = folder + '/merge_GO_component_rank.txt'
		known_go_component_enrich_d_plot = folder + '/known_GO_component_enrichment_d_plot.png'
		merge_go_component_enrich_d_plot = folder + '/merge_GO_component_enrichment_d_plot.png'
		known_go_component_enrich_t_plot = folder + '/known_GO_component_enrichment_t_plot.png'
		merge_go_component_enrich_t_plot = folder + '/merge_GO_component_enrichment_t_plot.png'

		out = open(query_chemicals_info,'w')
		out.write('Pubchem_ID\tChemical_Name\tMolecular_Weight\tSMILES\n')
		for v in ind2drug.values():
			out.write('\t'.join(x for x in v[0:-1]) + '\n')
		out.close()

		N_d = len(ind2drug)
		query_index_lst = ind2drug.keys()
		table_num2query = split_drugs(query_index_lst)
		simi_in = {}
		simi_out = {}
		ind2pro = {}
		drug2tar_known = {}
		tar2drug_known = {}
		drug2tar_predict = {}
		tar2drug_predict = {}
		for k,v in table_num2query.items():
			query_index_set = '(' + str(v)[1:-1] + ')'
			#similarity
			[ind2drug, simi_in, simi_out] = get_simi_drug_pairs(v, query_index_set, ind2drug, k, simi_in, simi_out)
			[drug2tar_known, tar2drug_known] = get_d2t_known(query_index_set, k, drug2tar_known, tar2drug_known)
			[drug2tar_predict, tar2drug_predict] = get_d2t_predict(query_index_set, k, drug2tar_predict, tar2drug_predict)
		
		if simi_in != {}:
			write_drug_pairs(input_simi_outfile, simi_in, ind2drug)
		else:
			out_log.write('No similarity data stored within the input chemicals')
		write_drug_pairs(simi_chemicals_outfile, simi_out, ind2drug)
	#find targets
		N_t_known = len(tar2drug_known)
		if drug2tar_predict == {}:
			out_log.write('cannot find predicted interactions in '+ dataset+'\n')
			tar2drug_merge = tar2drug_known
		else:
			tar2drug_merge = merge_dict(tar2drug_known, tar2drug_predict)
			N_t_merge = len(tar2drug_merge)
		pros_lst = '(' + str(tar2drug_merge.keys())[1:-1] +')'
		ind2pro = get_pros(pros_lst, ind2pro)

	# pathways
		in_file = known_paths_outfile
		[pro2path_known, path2pro_known, p2go_function_known, p2go_process_known, p2go_component_known, go_function2p_known, go_process2p_known, go_component2p_known, path2drug_known, go_function2d_known, go_process2d_known, go_component2d_known]=drug_gene_path_go(tar2drug_known, ind2drug)
		N_g_known = len(pro2path_known) # number of target genes in KEGG
		if tar2drug_predict!={}:
			[pro2path_predict, path2pro_predict, p2go_function_predict, p2go_process_predict, p2go_component_predict, go_function2p_predict, go_process2p_predict, go_component2p_predict, path2drug_predict, go_function2d_predict, go_process2d_predict, go_component2d_predict]=drug_gene_path_go(tar2drug_predict, ind2drug)
			path2pro_merge = merge_dict(path2pro_known, path2pro_predict)
			pro2path_merge = merge_dict(pro2path_known, pro2path_predict)
			N_g_merge = len(pro2path_merge)
			path2drug_merge = merge_dict(path2drug_known, path2drug_predict)
			p2go_function_merge = merge_dict(p2go_function_known, p2go_function_predict)
			go_function2p_merge = merge_dict(go_function2p_known, go_function2p_predict)
			p2go_process_merge = merge_dict(p2go_process_known, p2go_process_predict)
			go_process2p_merge = merge_dict(go_process2p_known, go_process2p_predict)
			p2go_component_merge = merge_dict(p2go_component_known, p2go_component_predict)
			go_component2p_merge = merge_dict(go_component2p_known, go_component2p_predict)		
			go_function2d_merge = merge_dict(go_function2d_known, go_function2d_predict)
			go_process2d_merge = merge_dict(go_process2d_known, go_process2d_predict)
			go_component2d_merge = merge_dict(go_component2d_known, go_component2d_predict)
			path_ind_lst = '('+str(path2pro_merge.keys())[1:-1]+')'
			ind2path = get_paths(path_ind_lst)
		else:
			path_ind_lst = '('+str(path2pro_known.keys())[1:-1]+')'
			ind2path = get_paths(path_ind_lst)

		name2go_function = get_name2go('function')
		name2go_process = get_name2go('process')
		name2go_component = get_name2go('component')

		# known        
		write_dt(known_dt_outfile, drug2tar_known,ind2drug,ind2pro)
		write_tar_rank(known_tars_outfile, known_tar_enrich_plot, tar2drug_known, ind2pro, M_d, N_d)
		write_path_rank(query_type, known_paths_outfile, known_path_enrich_t_plot, known_path_enrich_d_plot, path2pro_known, path2drug_known, ind2path, N_g_known, N_d)
		write_go_rank(query_type, known_go_function_outfile, known_go_function_enrich_t_plot, known_go_function_enrich_d_plot, go_function2p_known, go_function2d_known, name2go_function, N_t_known, N_d)
		write_go_rank(query_type, known_go_process_outfile, known_go_process_enrich_t_plot, known_go_process_enrich_d_plot, go_process2p_known, go_process2d_known, name2go_process, N_t_known, N_d)
		write_go_rank(query_type, known_go_component_outfile, known_go_component_enrich_t_plot, known_go_component_enrich_d_plot, go_component2p_known, go_component2d_known, name2go_component, N_t_known, N_d)

		# merge        
		if drug2tar_predict != {}:
			in_file = merge_paths_outfile
			write_dt(pred_dt_outfile, drug2tar_predict,ind2drug,ind2pro)
			write_tar_rank(merge_tars_outfile, merge_tar_enrich_plot, tar2drug_merge, ind2pro, M_d, N_d)
			write_path_rank(query_type, merge_paths_outfile, merge_path_enrich_t_plot, merge_path_enrich_d_plot, path2pro_merge, path2drug_merge, ind2path, N_g_merge, N_d)
			write_go_rank(query_type, merge_go_function_outfile, merge_go_function_enrich_t_plot, merge_go_function_enrich_d_plot, go_function2p_merge, go_function2d_merge, name2go_function, N_t_merge, N_d)
			write_go_rank(query_type, merge_go_process_outfile, merge_go_process_enrich_t_plot, merge_go_process_enrich_d_plot, go_process2p_merge, go_process2d_merge, name2go_process, N_t_merge, N_d)
			write_go_rank(query_type, merge_go_component_outfile, merge_go_component_enrich_t_plot, merge_go_component_enrich_d_plot, go_component2p_merge, go_component2d_merge, name2go_component, N_t_merge, N_d)

		# write all entry-pathway link
		out = open(gene_paths_outfile,'w')
		out.write('Gene_Name\tPathway_ID\tPathway_Name\n')
		with open(in_file, 'r') as f:
			next(f)
			for line in f:
				line = line.strip().split('\t')
				path_id = line[1]
				path_name = line[2]
				genenames = line[6].split(';')
				for genename in genenames:
					out.write(genename+'\t'+path_id+'\t'+path_name+'\n')
		out.close()

out_log.close()
conn.close()
print "--- %s seconds ---" % (time.time() - start_time)