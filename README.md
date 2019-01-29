# Quantitative systems pharmacological analysis of drugs of abuse reveals the pleiotropy of their targets and the effector role of mTORC1 
Fen Pei†, Hongchun Li†, Bing Liu* and Ivet Bahar*
Department of Computational and Systems Biology, School of Medicine, University of Pittsburgh, PA, 15213, USA
* Correspondence: 
Bing Liu liubing@pitt.edu 
Ivet Bahar bahar@pitt.edu
† These authors made equal contributions

The data of drug-target interactions from DrugBank, chemical-protein interactions from STITCH and pathways from KEGG were pre-processed and stored on our server for future online server development. 'get_DrugBank_targets_and_pathways.py' and 'get_STITCH_targets_pathways.py' connect the database on our server, the password is not open to public at present. The codes posted here are for analysis methods demonstration purpose.

## Data processing
Before using the analysis code, readers can get access to the original data from the following links:

**[Download the latest DrugBank data](https://www.drugbank.ca/releases/latest#external-links)**  

Download drug information:     	

curl -Lfv -o Drugs_all.zip -u user:password https://www.drugbank.ca/releases/5-1-1/downloads/all-drug-links  
curl -Lfv -o Drugs_approved.csv -u user:password https://www.drugbank.ca/releases/5-1-1/downloads/approved-drug-links  

Download drug-target information:  
curl -Lfv -o dt_known_all.zip -u user:password https://www.drugbank.ca/releases/5-1-1/downloads/target-all-uniprot-links  
curl -Lfv -o dt_known_approved.zip -u user:password https://www.drugbank.ca/releases/5-1-1/downloads/target-approved-uniprot-links  

Download target identifiers:  
curl -Lfv -o pros_all.zip -u user:password https://www.drugbank.ca/releases/5-1-1/downloads/target-all-polypeptide-ids  
curl -Lfv -o pros_approved.zip -u user:password https://www.drugbank.ca/releases/5-1-1/downloads/target-approved-polypeptide-ids  

**[Download the latest STITCH data](http://stitch.embl.de/cgi/download.pl?UserId=PDvH3yDVUJM3&sessionId=FyYYPkIIkZKS)**  

choose organism: Homo sapiens  

Download [chemical-protein links](http://stitch.embl.de/download/protein_chemical.links.detailed.v5.0/9606.protein_chemical.links.detailed.v5.0.tsv.gz) (with detailed subscores)
	
**[Download Uniprot protein ID mapping](ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/)**  
Select 'HUMAN_9606_idmapping_selected.tab'  
[ReadMe file](ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/README) explains the detailed columns.  

**Download pathway links:**  
'path2gene.txt' # downloaded from http://rest.kegg.jp/list/pathway/hsa  
'path2name.txt' # downloaded from http://rest.kegg.jp/link/hsa/pathway  
'kegg_class.txt' # get this information from https://www.kegg.jp/kegg/pathway.html (the original version was downloadable)  

## Predicting drug-target interactions using probabilistic matrix factorization (PMF) method  

Click [here](http://balestra1.csb.pitt.edu/static/balestraweb.zip) for the code of applying PMF on DrugBank (small dataset)   
Click [here](http://bickson.blogspot.com/2012/12/collaborative-filtering-with-graphchi.html) for the method of applying PMF on STITCH (large scale computation)   

## Calculate 2D similarities for a input lis of drugs

2D_similarity.py  

## Get targets and pathways for a input list of drugs  

get_DrugBank_targets_and_pathways.py    
get_STITCH_targets_pathways.py    

## Perform pathway enrichment analysis given a input list of drug-target interaction pairs

pathway_enrichment.py (the enrichment analysis is based on the number of targets, same as input a list of targets)  
adjust_pval.py (FDR correction for pathway enrichment analysis)  


