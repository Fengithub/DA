""" 
	input a list of drugs
	return the known and predicted drug-target interactions, enriched targets, enriched pathways, enriched GO terms
	note: the drug-target interactions from DrugBank are stored in our server database, predictions made by PMF are precalculated and stored on the server.
"""
import os
import psycopg2
import sys
from scipy.stats import hypergeom
import matplotlib.pyplot as plt
import numpy as np
plt.switch_backend('agg')

try:
	conn = psycopg2.connect(database="balestra", user="fep7", password="Balestra2017", host="127.0.0.1", port="5432")
except:
	out_log.write("I am unable to connect to the database\n")

cur = conn.cursor()

folder = sys.argv[1]
folder = folder.strip('/')
dataset = sys.argv[2] # approved, all
query_type = sys.argv[3] # drugs or proteins
max_pred = int(sys.argv[4])
query_file = sys.argv[5] # input file, including a list of chemical name/DB_id/pubchem id or protein name/Uniprot id/Protein Entry

# table names in database
prefix = 'fen_DB5_1_1_'
tbl_drug_index =  prefix + dataset + '_drug_index'
tbl_pro_index = prefix + dataset + '_pro_index'
tbl_path_index = prefix + 'path_index'
tbl_go_function = prefix + 'go_function'
tbl_go_process = prefix + 'go_process'
tbl_go_component = prefix + 'go_component'

tbl_d2t_predict = prefix + dataset +'_prediction_d2t'
tbl_t2d_predict = prefix + dataset +'_prediction_t2d'
tbl_D2D = prefix + dataset +'_d2d_list'
tbl_T2T = prefix + dataset +'_t2t_list'

# make a results folder #
if not os.path.exists(folder):
	os.makedirs(folder)

dt2pdb = {}
cur.execute("select drug_id, uniprot_id, pdb from fen_DB_dt_pdb;")
lines = cur.fetchall()
for drug_id, uniprot_id, pdb in lines:
	dt2pdb[(drug_id, uniprot_id)] = pdb

def merge_dict(dict1, dict2):
	""" 
		takes two dictinaries of the same kind as input
		return a dictionary that merge the values of each key from two dictionaries
		used to merge known and predicted results
	"""
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

def find_drug(query,ind2drug,dt_known,td_known):
	"""
	takes a query, a dictionary contains drug details (ind2drug), a dictionary of drug to known targets, and a dictionary of target to known drugs
	search the query drug in the database
	return the three input dictionaries updated with the query drug
	"""
	if query.startswith('DB'): # case1: query drugbank_id
		query_id = "\'"+query+"\'"
		cur.execute("select drug_index, drug_id, drug_name, synonyms, drug_type, drug_group, smiles, pubchem_id, chembl_id, conditions, num_pros, pro_inds, num_paths, path_ids, path_names, go_function, go_process, go_component from "+ tbl_drug_index + " where drug_id="+query_id+";")
	else: # case2: query drug_name, case insensative, do not need to write the full name
		query_like = "\'%"+query+"%\'"
		cur.execute("select drug_index, drug_id, drug_name, synonyms, drug_type, drug_group, smiles, pubchem_id, chembl_id, conditions, num_pros,  pro_inds, num_paths, path_ids, path_names, go_function, go_process, go_component  from "+ tbl_drug_index + " where drug_name ILIKE"+query_like+";")
	drug = cur.fetchall()
	if drug == []:
		out_log.write('cannot find '+query+' in DrugBank ' + dataset+' dataset\n')
	else:
		drug = drug[0]
		query_index = drug[0]
		ind2drug[query_index] = drug[1:]
		### find known targets for the query ###
		tars = drug[11].split(';')
		for t in tars:
			t = int(t)
			dt_known.setdefault(query_index,set()).add((t,1))
			td_known.setdefault(t,set()).add((query_index,1))
	return [ind2drug,dt_known,td_known]

def get_drugs(drug_ind_lst, ind2drug):
	"""
	takes a list of drug indexes, a dictionary contains drug details (ind2drug)
	search the list of drugs in the database
	return the input dictionary updated with the input list of drug
	"""
	cur.execute("select drug_index, drug_id, drug_name, synonyms, drug_type, drug_group, smiles, pubchem_id, chembl_id, conditions, num_pros, pro_inds, num_paths, path_ids,path_names, go_function, go_process, go_component from "+ tbl_drug_index + " where drug_index in "+drug_ind_lst+";")
	rows = cur.fetchall()
	for row in rows:
		ind2drug[row[0]] = row[1:]
	return ind2drug

def get_pros(pro_ind_lst, ind2pro):
	"""
	takes a list of protein indexes, a dictionary contains protein details (ind2drug)
	search the list of proteins in the database
	return the input dictionary updated with the input list of proteins
	"""
	cur.execute("select pro_index, uniprot_id, protein_name,  entry_name, num_drug, gene_name, gene_id, PDB, path_ids, path_names, go_function, go_process, go_component from "+ tbl_pro_index + " where pro_index in "+pro_ind_lst+";")
	rows = cur.fetchall()
	for row in rows:
		ind2pro[row[0]] = row[1:]
	return ind2pro

def get_d2t_known(query_index_lst): 
	"""
	takes a list of drug indexes
	search the list of drugs in the database for their known targets
	return the two dictionaries: {drug: known targets} and {target: known drugs} 
	"""	
	dt = {}
	td = {}
	cur.execute("select drug_index, pro_inds from "+ tbl_drug_index +" where drug_index in "+query_index_lst+";")
	rows = cur.fetchall()
	if rows != []:
		for drug_index, pro_inds in rows:
			pro_index_list = pro_inds.split(';')
			for t in pro_index_list:
				pro_index = int(t)
				dt.setdefault(drug_index,set()).add((pro_index,1))
				td.setdefault(pro_index,set()).add((drug_index,1))
	return [dt, td]

def get_d2t_predict(query_index_lst, max_pred):
	"""
	takes a list of drug indexes
	search the list of drugs in the database for their predicted targets
	return the two dictionaries: {drug: predicted targets} and {target: predicted drugs} 
	"""	
	dt = {}
	td = {}
	cur.execute("select drug_index, pro_index_list, score_list from "+ tbl_d2t_predict +" where drug_index in "+query_index_lst+";")
	rows = cur.fetchall()
	if rows != []:
		for drug_index, pro_index_list, score_list in rows:
			pro_index_list = pro_index_list.split(';')
			score_list = score_list.split(';')
			num_pred = len(pro_index_list)
			for i in range(min(num_pred, max_pred)):
				pro_index = int(pro_index_list[i])
				interaction = float(score_list[i])
				dt.setdefault(drug_index,set()).add((pro_index,interaction))
				td.setdefault(pro_index,set()).add((drug_index,interaction))
	return [dt, td]

def get_paths(path2pro):
	"""
	takes a dictionary of {pathway index: proteins}
	search the keys of the dictionary (pathway indexes) in the database
	return a dictionary with the information of pathways
	"""	
	path_ind_lst = '(' + str(path2pro.keys())[1:-1]+')'
	ind2path = {}
	cur.execute("select path_index, path_id, path_name, num_genes, class_1st, class_2nd, num_drugs_all from "+ tbl_path_index + " where path_index in "+path_ind_lst+";")
	rows = cur.fetchall()    
	for row in rows:
		ind2path[str(row[0])] = row[1:]
	return ind2path

def get_simi_drug_pairs(query_index_lst, query_index_set, ind2drug):
	"""
	takes a list of query drug indexes, the set format of these drug indexes, and a dictionary of {drug index: drug information}
	search the similar drugs for the query drugs
	return a dictionary with the updated information of query drugs and similar drugs, a dictionary of the similarities within the query drugs, and a dictionary of the similarities outside the query drugs.
	"""		
	simi_index_set = set()
	simi_in = {}
	simi_out = {}
	cur.execute("select drug1_index, drug2_index_list, score_dt_list, score_2d_list from "+ tbl_D2D + " where drug1_index in "+query_index_set+";")
	rows = cur.fetchall()
	for drug1_index, drug2_index_list, score_dt_list, score_2d_list in rows:
		drug2_index_list = drug2_index_list.split(';')
		score_dt_list = score_dt_list.split(';')
		score_2d_list = score_2d_list.split(';')
		for i in range(min(max_pred,len(drug2_index_list))):
			drug2_index = int(drug2_index_list[i])
			similarity = float(score_dt_list[i])
			similarity_2D = score_2d_list[i]
			if drug2_index in query_index_lst:				
				simi_in.setdefault(drug1_index,set()).add((drug2_index,similarity,similarity_2D))
			else:
				simi_index_set.add(drug2_index)
				simi_out.setdefault(drug1_index,set()).add((drug2_index,similarity,similarity_2D))
	if simi_index_set != set():
		simi_index_lst = list(simi_index_set)
		simi_index_set = '(' + str(simi_index_lst)[1:-1] +')'
		ind2drug = get_drugs(simi_index_set, ind2drug)
	return [ind2drug, simi_in, simi_out]

def get_pros2path_go(query_index_lst):
	"""
	takes a list of query protein indexes
	get the protein-pathway, parotein-GO terms relationships
	return 8 dictionaries as follows:
		pro2path: {protein index: pathway indexes}
		path2pro: {pathway index: (protein index, gene name, gene id)}
		p2go_function: {protein index: GO function terms}
		p2go_process: {protein index: GO process terms}
		p2go_component: {protein index: GO component terms}
		go_function2p: {GO function term: (protein index, gene name, gene id)}
		go_process2p: {GO process term: (protein index, gene name, gene id)}
		go_component2p: {GO component term: (protein index, gene name, gene id)}
	"""		
	pro2path = {}
	path2pro = {}
	p2go_function = {}
	p2go_process = {}
	p2go_component = {}
	go_function2p = {}
	go_process2p = {}
	go_component2p = {}
	cur.execute("select pro_index, gene_name, gene_id, path_inds, go_function, go_process, go_component from "+ tbl_pro_index +" where pro_index in "+query_index_lst+";")
	rows = cur.fetchall()
	if rows != []:
		for pro_index, gene_name, gene_id, path_inds, go_function, go_process, go_component in rows:
			if path_inds != 'NA':
				path_inds = path_inds.split(';')
				for path_index in path_inds:
					pro2path.setdefault(pro_index,set()).add(path_index)
					path2pro.setdefault(path_index,set()).add((pro_index,gene_name,gene_id))			
			if go_function != 'NA':
				go_function = go_function.split('|')
				for one in go_function:
					p2go_function.setdefault(pro_index,set()).add(one)
					go_function2p.setdefault(one,set()).add((pro_index,gene_name,gene_id))
			if go_process != 'NA':
				go_process = go_process.split('|')
				for one in go_process:
					p2go_process.setdefault(pro_index,set()).add(one)
					go_process2p.setdefault(one,set()).add((pro_index,gene_name,gene_id))
			if go_component != 'NA':
				go_component = go_component.split('|')
				for one in go_component:
					p2go_component.setdefault(pro_index,set()).add(one)
					go_component2p.setdefault(one,set()).add((pro_index,gene_name,gene_id))					
	return [pro2path, path2pro, p2go_function, p2go_process, p2go_component, go_function2p, go_process2p, go_component2p]  #pro2path, path2pro

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
	pros_lst = '(' + str(tar2drug.keys())[1:-1] + ')'
	[pro2path, path2pro, p2go_function, p2go_process, p2go_component, go_function2p, go_process2p, go_component2p] = get_pros2path_go(pros_lst)
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

def get_name2go(go_type):
	"""
	take GO term to targets dictionary and GO term type as input
	search the database for the total number of proteins and drugs that are associated with that term
	return a name2go dictionary: {GO term name: total number of proteins, total number of drugs}, used to calculated enrichment p-value
	"""
	tbl_go_index = prefix + 'go_' + go_type
	name2go = {}
	cur.execute("select go_name, num_pros_all, num_drugs_all from "+ tbl_go_index +";")
	rows = cur.fetchall()
	if rows != []:
		for go_name, num_pros_all, num_drugs_all in rows:
			name2go[go_name] = (num_pros_all, num_drugs_all)
	return name2go

def write_drug_pairs(out_file, pair_dict, ind2drug):
	"""
	take the output file name, drug pairs dictionary: {drug_index: [(similar drug index, DT similarity, 2D similarity)]}, ind2drug dictionary as input 
	write the drug-drug pairs with similarities in a file
	"""	
	out = open(out_file,'w')
	out.write('Drug1_ID\tDrug1_Name\tDrug1_Synonyms\tDrug1_Type\tDrug1_Group\tDrug1_SMILES\tDrug1_Pubchem_ID\tDrug1_CHMBEL_ID\tDrug1_Associated_Conditions\tDrug2_ID\tDrug2_Name\tDrug2_Synonyms\tDrug2_Type\tDrug2_Group\tDrug2_SMILES\tDrug2_Pubchem_ID\tDrug2_CHMBEL_ID\tDrug2_Associated_Conditions\tSimilarity_DT\tSimilarity_2D\n')
	for k,v in pair_dict.items():
		a = ind2drug[k]
		for b in v:
			c = ind2drug[b[0]]
			out.write('\t'.join(x for x in a[0:9]) + '\t' + '\t'.join(x for x in c[0:9]) +'\t'+ str(b[1])+'\t'+ b[2]+ '\n')
	out.close()

def write_dt(out_file, dt_dict,ind2drug,ind2pro):
	"""
	take the output file name, drug-target interaction dictionary: {drug_index: [(target index, confidence score)]}, ind2drug, ind2pro as input 
	write the drug-target pairs with detailed drug, target information and confidence score in a file
	"""		
	out = open(out_file,'w')
	out.write('Drug_ID\tDrug_Name\tSynonyms\tDrug_Type\tDrug_Group\tSMILES\tPubchem_ID\tCHMBEL_ID\tAssociated_Conditions\tUniprot_ID\tProtein_Name\tEntry_Name\tTot_Num_Drugs\tGene_Name\tGene_ID\tPDB\tPathways_ID\tPathways_Name\tGO_Function\tGO_Process\tGO_Component\tConfidnece_Score\tPDB_dt\n')
	for k,v in dt_dict.items():
		drug = ind2drug[k]
		v = sorted(list(v), key = lambda x:x[1], reverse= True)
		for one in v:
			pro_ind = one[0]
			pro = ind2pro[pro_ind]
			dt = (drug[0],pro[0])
			if dt in dt2pdb:
				dt_pdb = dt2pdb[dt]
			else:
				dt_pdb = 'NA'
			out.write('\t'.join(str(x) for x in drug[0:9]) + '\t' +'\t'.join(str(x) for x in pro[0:9])+ '\t' + '\t'.join(pro[-3:])+'\t'+ str(one[1]) +'\t' + dt_pdb + '\n')
	out.close()

def write_tar_rank(tbl_file, fig_file, pro2drug, ind2pro, pro2path, ind2path, M, N):
	"""
	take the output table file name, output figure file name, target-drug interaction dictionary: {protein_index: [(drug index, confidence score)]}, ind2drug, ind2pro, the total number of drugs (M_d) in the database, and the total number of input drugs (N_d) as input 
	write the targets ranked by the enrichment pvalue based on interacting drugs
	plot the barplot of targets with the enrichment pvalue based on the interacting drugs
	"""	
	pro2pval = {}
	for k, v in pro2drug.items():
		pro = ind2pro[k]
		q = pro[3]
		m = sum([x[1] for x in v])
		if m > q: # for predicted case, we keep use the total number of known chemicals interacting with a certain target, if predicted number greater than known, pval is 0.
			pval = hypergeom.sf(m-1, M, m, N)
		else:
			pval = hypergeom.sf(m-1, M, q, N) # Hypergeometric test: the probability of getting more than (m-1) items from N, when the backgroud is q in M.
		pro2pval[k] = pval
	
	n = 1
	out = open(tbl_file,'w')
	out.write('Rank\tUniprot_ID\tProtein_Name\tEntry_Name\tTot_Num_Drugs\tGene_Name\tGene_ID\tPDB\tPathways_ID\tPathways_Name\tGO_Function\tGO_Process\tGO_Component\tDrugs\tNum_Drugs\tP_value\n')
	pro_name2enrich_score = {}
	for k, v in sorted(pro2pval.items(),key=lambda x:x[1]):
		pro = ind2pro[k]
		drugs = [(ind2drug[x[0]][1],x[1]) for x in pro2drug[k]]
		out.write(str(n)+'\t'+'\t'.join(str(x) for x in pro[0:9])+ '\t' + '\t'.join(pro[-3:])+'\t'+';'.join(str(x) for x in drugs)+'\t'+str(len(drugs))+'\t'+str(v)+'\n')
		if n <= 20:
			pro_name2enrich_score[pro[1]] = -np.log10(v)
		n += 1
	out.close()
	enrich_plot(pro_name2enrich_score, fig_file, 'target')
	
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
		out.write('Rank\tPathway_ID\tPathway_Name\tTotal_Gene_Number\tPathway_Class_1st\tPathway_Class_2nd\tGene_Names\tGene_IDs\tNum_Targets\tDrugs\tNum_Drugs\tP_value_target\tP_value_drug\n')
		for k, v in sorted(path2pval_t.items(),key=lambda x:x[1]):
			path = ind2path[k][0:-1]
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
			path = ind2path[k][0:-1]
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
		out.write('Rank\tGo_Name\tGene_Names\tGene_IDs\tNum_Targets\tDrugs\tNum_Drugs\tP_value_target\tP_value_drug\n')
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
	
def find_target(query, ind2pro, d2t, t2d):
	"""input: a target query
		search the query in the protein table, get the corresponding protein information
		return query index, a dictionary with protein information, dictionaries of drug2target,target2drug
	"""
	query_id = "\'"+query+"\'"
	cur.execute("select pro_index, uniprot_id, protein_name, entry_name, num_drug, gene_name, gene_id, PDB, path_ids, path_names, path_inds, drug_inds, go_function, go_process, go_component from "+ tbl_pro_index + " where uniprot_id="+query_id+";")
	pro = cur.fetchall()
	if pro == []:
		cur.execute("select pro_index,uniprot_id, protein_name, entry_name, num_drug, gene_name, gene_id, PDB, path_ids, path_names, path_inds, drug_inds, go_function, go_process, go_component from "+ tbl_pro_index + " where gene_names="+query_id+";")
		pro = cur.fetchall()
		if pro == []:
			query_like = "\'%"+query+"%\'"
			cur.execute("select pro_index,uniprot_id, protein_name, entry_name, num_drug, gene_name, gene_id, PDB, path_ids, path_names, path_inds, drug_inds, go_function, go_process, go_component from "+ tbl_pro_index + " where protein_name ILIKE "+query_like+";")
			pro = cur.fetchall()
	if pro == []:
		out_log.write('cannot find '+query+' in ' + DATASET+'\n')
		query_index = None
	else:
		pro = pro[0]
		query_index = pro[0]
		ind2pro[query_index] = pro[1:]
		drugs = pro[11].split(';')
		for d in drugs:
			d = int(d)
			d2t.setdefault(d,set()).add((query_index,1))
			t2d.setdefault(query_index,set()).add((d,1))
	return [query_index, ind2pro, d2t, t2d]

def get_t2d_known(query_index_lst):
	"""input a list of protein indexes, find the known drugs for them, return two dictionaries: drug2target, target2drug """
	dt = {}
	td = {}
	cur.execute("select pro_index, drug_inds from "+ tbl_pro_index +" where pro_index in "+query_index_lst+";")
	rows = cur.fetchall()
	if rows != []:
		for pro_index, drug_inds in rows:
			drug_inds = drug_inds.split(';')
			for d in drug_inds:
				drug_index = int(d)
				dt.setdefault(drug_index,set()).add((pro_index,1))
				td.setdefault(pro_index,set()).add((drug_index,1))
	return [dt, td]

def get_t2d_predict(query_index_lst, max_pred):
	"""input a list of target indexes, find the predicted drugs for them, return two dictionaries: drug2target, target2drug """
	dt = {}
	td = {}
	cur.execute("select pro_index, drug_index_list, score_list from "+ tbl_t2d_predict +" where pro_index in "+query_index_lst+";")
	rows = cur.fetchall()
	if rows != []:
		for pro_index, drug_index_list, score_list in rows:
			drug_index_list = drug_index_list.split(';')
			score_list = score_list.split(';')
			for i in range(min(max_pred,len(drug_index_list))):
				drug_index = int(drug_index_list[i])
				interaction = float(score_list[i])
				dt.setdefault(drug_index,set()).add((pro_index,interaction))
				td.setdefault(pro_index,set()).add((drug_index,interaction))
	return [dt, td]

def get_simi_pro_pairs(query_index_lst, query_index_set, ind2pro):
	"""
	takes a list of query pro indexes, the set format of these pro indexes, and a dictionary of {pro index: pro information}
	search the similar proteins for the query proteins
	return a dictionary with the updated information of query proteins and similar proteins, a dictionary of the similarities within the query proteins, and a dictionary of the similarities outside the query proteins.
	"""		
	simi_index_set = set()
	simi_in = {}
	simi_out = {}
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

def write_pro_pairs(out_file, pair_dict, ind2pro):
	"""
	take the output file name, protein pairs dictionary: {protein_index: [(similar protein index, DT similarity)]}, ind2pro dictionary as input 
	write the protein-protein pairs with similarities in a file
	"""	
	out = open(out_file,'w')
	out.write('Gene1_Name\tGene2_Name\tSimilarity_DT\n')	
	for k,v in pair_dict.items():
		a = ind2pro[k][4]
		for b in v:
			c = ind2pro[b[0]][4]
			out.write(a + '\t' + c +'\t'+ str(b[1]) + '\n')
	out.close()

def write_pros(out_file, pro_lst, ind2pro): 
	"""write the protein information into a txt file"""
	out = open(out_file,'w')
	out.write('Uniprot_ID\tProtein_Name\tEntry_Name\tTot_Num_Drugs\tGene_Name\tGene_ID\tPDB\tPathways_ID\tPathways_Name\tGO_Function\tGO_Process\tGO_Component\n')
	for one in pro_lst:
		pro = ind2pro[one]
		out.write('\t'.join(str(x) for x in pro[0:9])+ '\t' + '\t'.join(pro[-3:]) + '\n')
	out.close()

def write_td(out_file, td_dict, ind2drug, ind2pro):
	"""
	input the output filename, a target-drug dictionary, ind2drug and ind2pro
	write the drug-target interactions into a txt file
	"""
	out = open(out_file,'w')
	out.write('Uniprot_ID\tProtein_Name\tEntry_Name\tTot_Num_Drugs\tGene_Name\tGene_ID\tPDB\tPathways_ID\tPathways_Name\tGO_Function\tGO_Process\tGO_Component\tDrug_ID\tDrug_Name\tSynonyms\tDrug_Type\tDrug_Group\tSMILES\tPubchem_ID\tCHMBEL_ID\tAssociated_Conditions\tConfidnece_Score\tPDB_dt\n')
	for k,v in td_dict.items():
		pro = ind2pro[k]
		v = sorted(list(v), key = lambda x:x[1], reverse= True)
		for one in v:
			drug_ind = one[0]
			drug = ind2drug[drug_ind]
			dt = (drug[0],pro[0])
			if dt in dt2pdb:
				dt_pdb = dt2pdb[dt]
			else:
				dt_pdb = 'NA'
			out.write('\t'.join(str(x) for x in pro[0:9])+ '\t' + '\t'.join(pro[-3:]) + '\t' +'\t'.join(str(x) for x in drug[0:9])+'\t'+ str(one[1]) +'\t'+ dt_pdb +'\n')
	out.close()

def write_drug_rank(tbl_file, fig_file, drug2pro, ind2drug, ind2pro, M, N): 
	"""
		input a drug2pro dictionary
		rank the drugs according to the enrichment pvalue based on targets
		write it into a txt file
	"""
	drug2pval = {}
	for k, v in drug2pro.items():
		drug = ind2drug[k]
		q = drug[9] # where the number of pros
		m = len(v)
		if m > q: # for predicted case, we keep use the total number of known drugs interacting with a certain target, if predicted number greater than known, pval is 0.
			pval = hypergeom.sf(m-1, M, m, N)
		else:
			pval = hypergeom.sf(m-1, M, q, N) # Hypergeometric test: the probability of getting more than (m-1) items from N, when the backgroud is q in M.
		drug2pval[k] = pval
	
	n = 1
	out = open(tbl_file,'w')
	out.write('Rank\tDrug_ID\tDrug_Name\tSynonyms\tDrug_Type\tDrug_Group\tSMILES\tPubchem_ID\tCHMBEL_ID\tAssociated_Conditions\tProteins\tNum_proteins\tP_value\n')
	drug_name2erich_score = {} ## in order to draw the enrichment plot
	for k, v in sorted(drug2pval.items(),key=lambda x:x[1]):
		drug = ind2drug[k]
		pros = [(ind2pro[x[0]][4],x[1]) for x in drug2pro[k]]
		out.write(str(n)+'\t'+'\t'.join(str(x) for x in drug[0:9])+'\t'+';'.join(str(x) for x in pros)+'\t'+str(len(pros))+'\t'+str(v)+ '\n')
		if n <= 20:
			name = drug[1]
			drug_name2erich_score[name] = -np.log10(v)
		n += 1
	out.close()
	enrich_plot(drug_name2erich_score, fig_file, 'drug')

if dataset == 'approved':
	M_t = 2244  # total number of targets
	M_d = 1883  # total number of drugs
elif dataset == 'all':
	M_t = 2807 # total number of targets
	M_d = 5494 # total number of drugs

out_log = open(folder + '/log.txt','w')

if query_type == 'drugs':
	query_id_set = set()
	query_name_set = set()
	ind2drug = {}
	d2t_known = {}
	t2d_known = {}

	with open(query_file,'r') as f:
		for line in f:
			query = line.strip()
			if len(query) >0:
				[ind2drug,d2t_known,t2d_known]=find_drug(query,ind2drug,d2t_known,t2d_known)

	if ind2drug == {}:
		out_log.write('cannot find any queried drugs in DrugBank ' + dataset + ' dataset\n')
		exit()
	else:
		### output filenames ###
		## drugs ##
		query_drugs_info = folder + '/input_drugs_details.txt'
		input_simi_outfile = folder + '/input_similarity.txt' 
		simi_drugs_outfile = folder + '/similar_drug_pairs.txt'
		## drug-target interactions ##
		known_dt_outfile = folder + '/known_interactions.txt'
		pred_dt_outfile = folder + '/predicted_interactions.txt'
		## targets ##
		known_tars_outfile = folder + '/known_targets_rank.txt'
		known_tar_enrich_plot = folder + '/known_target_enrichment_plot.png'
		merge_tars_outfile = folder + '/merge_targets_rank.txt'
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
		merge_go_function_outfile = folder + '/merge_GO_function_rank.txt'
		known_go_function_enrich_d_plot = folder + '/known_GO_function_enrichment_d_plot.png'
		merge_go_function_enrich_d_plot = folder + '/merge_GO_function_enrichment_d_plot.png'
		known_go_function_enrich_t_plot = folder + '/known_GO_function_enrichment_t_plot.png'
		merge_go_function_enrich_t_plot = folder + '/merge_GO_function_enrichment_t_plot.png'
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

		### write out the details of each queried drug ###
		out = open(query_drugs_info,'w')
	
		out.write('Drug_ID\tDrug_Name\tSynonyms\tDrug_Type\tDrug_Group\tSMILES\tPubchem_ID\tCHMBEL_ID\tAssociated_Conditions\tPathway_Name\tGO_Function\tGO_Process\tGO_Component\n')
		for v in ind2drug.values():
			drug_id = v[0]
			drug_name = v[1]
			out.write('\t'.join(v[0:9]) + '\t'+ '\t'.join(v[13:]) + '\n')
		out.close()

		N_d = len(ind2drug)
		query_index_lst = ind2drug.keys()
		query_index_set = '(' + str(query_index_lst)[1:-1] + ')'
	#similarity
		[ind2drug, simi_in, simi_out] = get_simi_drug_pairs(query_index_lst, query_index_set, ind2drug)
		write_drug_pairs(input_simi_outfile, simi_in, ind2drug)
		write_drug_pairs(simi_drugs_outfile, simi_out, ind2drug)

	#find targets
		ind2pro = {}
		N_t_known = len(t2d_known)
		[d2t_predict, t2d_predict] = get_d2t_predict(query_index_set, max_pred)
		if d2t_predict == {}:
			out_log.write('cannot find predicted interactions in DrugBank '+ dataset+' dataset\n')
		else:
			t2d_merge = merge_dict(t2d_known, t2d_predict)
			N_t_merge = len(t2d_merge)
		pros_lst = '(' + str(t2d_merge.keys())[1:-1] +')'
		ind2pro = get_pros(pros_lst, ind2pro)

	# pathways and GO terms
		in_file = known_paths_outfile
		[pro2path_known, path2pro_known, p2go_function_known, p2go_process_known, p2go_component_known, go_function2p_known, go_process2p_known, go_component2p_known, path2drug_known, go_function2d_known, go_process2d_known, go_component2d_known]=drug_gene_path_go(t2d_known, ind2drug)
		if t2d_predict!={}:		
			[pro2path_predict, path2pro_predict, p2go_function_predict, p2go_process_predict, p2go_component_predict, go_function2p_predict, go_process2p_predict, go_component2p_predict, path2drug_predict, go_function2d_predict, go_process2d_predict, go_component2d_predict]=drug_gene_path_go(t2d_predict, ind2drug)
			path2pro_merge = merge_dict(path2pro_known, path2pro_predict)
			pro2path_merge = merge_dict(pro2path_known, pro2path_predict)
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
			ind2path = get_paths(path2pro_merge)
		else:
			ind2path = get_paths(path2pro_known)
		name2go_function = get_name2go('function')
		name2go_process = get_name2go('process')
		name2go_component = get_name2go('component')

		# known        
		write_dt(known_dt_outfile, d2t_known,ind2drug,ind2pro)
		write_tar_rank(known_tars_outfile, known_tar_enrich_plot, t2d_known, ind2pro, pro2path_known, ind2path, M_d, N_d)
		write_path_rank(query_type, known_paths_outfile, known_path_enrich_t_plot, known_path_enrich_d_plot, path2pro_known, path2drug_known, ind2path, N_t_known, N_d)
		write_go_rank(query_type, known_go_function_outfile, known_go_function_enrich_t_plot, known_go_function_enrich_d_plot, go_function2p_known, go_function2d_known, name2go_function, N_t_known, N_d)
		write_go_rank(query_type, known_go_process_outfile, known_go_process_enrich_t_plot, known_go_process_enrich_d_plot, go_process2p_known, go_process2d_known, name2go_process, N_t_known, N_d)
		write_go_rank(query_type, known_go_component_outfile, known_go_component_enrich_t_plot, known_go_component_enrich_d_plot, go_component2p_known, go_component2d_known, name2go_component, N_t_known, N_d)

		# merge        
		if d2t_predict != {}:
			in_file = merge_paths_outfile
			write_dt(pred_dt_outfile, d2t_predict,ind2drug,ind2pro)
			write_tar_rank(merge_tars_outfile, merge_tar_enrich_plot, t2d_merge, ind2pro, pro2path_merge, ind2path, M_d, N_d)
			write_path_rank(query_type, merge_paths_outfile, merge_path_enrich_t_plot, merge_path_enrich_d_plot, path2pro_merge, path2drug_merge, ind2path, N_t_merge, N_d)
			write_go_rank(query_type, merge_go_function_outfile, merge_go_function_enrich_t_plot, merge_go_function_enrich_d_plot, go_function2p_merge, go_function2d_merge, name2go_function, N_t_merge, N_d)
			write_go_rank(query_type, merge_go_process_outfile, merge_go_process_enrich_t_plot, merge_go_process_enrich_d_plot, go_process2p_merge, go_process2d_merge, name2go_process, N_t_merge, N_d)
			write_go_rank(query_type, merge_go_component_outfile, merge_go_component_enrich_t_plot, merge_go_component_enrich_d_plot, go_component2p_merge, go_component2d_merge, name2go_component, N_t_merge, N_d)

		# write all gene-pathway link
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