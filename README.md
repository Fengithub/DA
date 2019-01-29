# Quantitative systems pharmacological analysis of drugs of abuse reveals the pleiotropy of their targets and the effector role of mTORC1 
Fen Pei†, Hongchun Li†, Bing Liu* and Ivet Bahar*
Department of Computational and Systems Biology, School of Medicine, University of Pittsburgh, PA, 15213, USA
* Correspondence: 
Bing Liu liubing@pitt.edu 
Ivet Bahar bahar@pitt.edu
† These authors made equal contributions

The data of drug-target interactions from DrugBank, chemical-protein interactions from STITCH and pathways from KEGG were pre-processed and stored on our server for future online server development. 'get_DrugBank_targets_and_pathways.py' and 'get_STITCH_targets_pathways.py' here are for analysis methods demonstration purpose, readers can get access to the original data from the following links:

# The latest DrugBank data is publically accessable: https://www.drugbank.ca/releases/latest#external-links.
Download drug information:
Command:	
curl -Lfv -o Drugs_all.zip -u user:password https://www.drugbank.ca/releases/5-1-1/downloads/all-drug-links
curl -Lfv -o Drugs_approved.csv -u user:password https://www.drugbank.ca/releases/5-1-1/downloads/approved-drug-links

Download drug-target information:
curl -Lfv -o dt_known_all.zip -u user:password https://www.drugbank.ca/releases/5-1-1/downloads/target-all-uniprot-links
curl -Lfv -o dt_known_approved.zip -u user:password https://www.drugbank.ca/releases/5-1-1/downloads/target-approved-uniprot-links

Download target identifiers:
curl -Lfv -o pros_all.zip -u user:password https://www.drugbank.ca/releases/5-1-1/downloads/target-all-polypeptide-ids
curl -Lfv -o pros_approved.zip -u user:password https://www.drugbank.ca/releases/5-1-1/downloads/target-approved-polypeptide-ids

# The latest STITCH data is publically accessable: http://stitch.embl.de/cgi/download.pl?UserId=PDvH3yDVUJM3&sessionId=FyYYPkIIkZKS

choose organism: Homo sapiens

Download chemical-protein links (with detailed subscores): http://stitch.embl.de/download/protein_chemical.links.detailed.v5.0/9606.protein_chemical.links.detailed.v5.0.tsv.gz
	
# Download uniprot protein ID mapping from: ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/
Select 'HUMAN_9606_idmapping_selected.tab'
ReadMe file explains the detailed columns:ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/README

# Download pathway links:
'path2gene.txt' # downloaded from http://rest.kegg.jp/list/pathway/hsa
'path2name.txt' # downloaded from http://rest.kegg.jp/link/hsa/pathway
'kegg_class.txt' # get this information from https://www.kegg.jp/kegg/pathway.html (the original version was downloadable)

# The code for applying PMF on DrugBank (small dataset) is accessible: http://balestra1.csb.pitt.edu/static/balestraweb.zip
# The method for applying PMF on STITCH (large scale computation): http://bickson.blogspot.com/2012/12/collaborative-filtering-with-graphchi.html





