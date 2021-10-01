import os
import pandas as pd
import numpy as np
from lxml import etree as etree_lxml
from biothings import config
from biothings.utils.dataload import dict_convert, dict_sweep
logging = config.logger

process_key = lambda k: k.replace(" ","_").lower()

#def setup_release(self):
    #release="2021-09-20"
    #return release

# Association Parser    
def load_hmdb_data(data_folder):
        
    # --- Set input XML file path ---
    protein_xml = os.path.join(data_folder, "hmdb_proteins.xml")
    meta_xml = os.path.join(data_folder, "hmdb_metabolites.xml")

    # --- Set Data ---
    def load_xml(opts):
        """Return the sample XML file as a string."""
        with open(meta_xml, opts) as xml:
            return xml.read()


    xml_as_bytes = load_xml('rb')
    tree = etree_lxml.fromstring(xml_as_bytes)

    # --- Upload XML Data (workaround for pandas >1.3, need to upgrade Biothings) --- 
    xml_data = open(protein_xml, 'r', encoding='UTF-8').read()  # Read file
    root = ET.XML(xml_data)  # Parse XML


    metabolite = tree.findall('{http://www.hmdb.ca}metabolite', {})

    metabolite_dict={}

    for meta in metabolite:
        #print(meta.tag)
        accession=meta.find('{http://www.hmdb.ca}accession')
        kegg=meta.find('{http://www.hmdb.ca}kegg_id')
        chemspider=meta.find('{http://www.hmdb.ca}chemspider_id')
        chebi=meta.find('{http://www.hmdb.ca}chebi_id')
        pubchem=meta.find('{http://www.hmdb.ca}pubchem_compound_id')

        metabolite_dict.setdefault(accession.text, {
                                                    "kegg_id":kegg.text,
                                                    "chemspider_id": chemspider.text,
                                                    "chebi_id": chebi.text,
                                                    "pubchem_compound_id": pubchem.text

        })

    #print(p)

    data_list=[] # final data holder 

    # iterate over the parsed tree
    for t in root[0].iter():
        metabolites = [] # holder for metabolite_associations.metabolite, the associations w/o references
        
        if t.tag == "{http://www.hmdb.ca}protein": # get the protein nodes only
            
            # we need the first accession number, this is main protein _id in our doc
            _id = t.find("{http://www.hmdb.ca}accession")
            _id = _id.text   
            protein_type=t.find("{http://www.hmdb.ca}protein_type")
            ct=1 # setup counter for the associations
            
            
            # ---------- Metabolite associations with references ------------  

            for m in t.findall("{http://www.hmdb.ca}metabolite_references"):
                
                id_list=[] # setup association id list
                
                for ref in m:
                    data={}
                    _id2=_id+"_%s"%ct
                    id_list.append(_id2)
                    ct+=1
                    data["_id"]=_id2 
                    data["pmid"]= None
                    data["subject"]={"protein_type": protein_type.text}
                    data["object"]={}
        
                    for met_ref in ref.findall("{http://www.hmdb.ca}reference"):
                        for refs in met_ref:
                            if "pubmed_id" in refs.tag:
                                data['pmid']=refs.text

                    for met in ref.findall("{http://www.hmdb.ca}metabolite"):
                        for info in met:
                            tag=info.tag.split("}")[1]
                            text=info.text
                            
                            data["object"][tag]=text

                            # get the extra IDs from the metabolite xml
                            #{'kegg_id': 'C01092', 'chemspider_id': '4578', 'chebi_id': '127029', 'pubchem_compound_id': '4740'}

                            if "accession" in tag:
                                #print(metabolite_dict[text])
                                data["object"]["kegg_id"]=metabolite_dict[text]["kegg_id"]
                                data["object"]["chemspider_id"]=metabolite_dict[text]["chemspider_id"]
                                data["object"]["chebi_id"]=metabolite_dict[text]["chebi_id"]
                                data["object"]["pubchem_compound_id"]=metabolite_dict[text]["pubchem_compound_id"]
                            
                    # setup subject info       
                    #uniprot_id, uniprot_name, genbank_protein_id, hgnc_id, genbank_gene_id, and gene_name.        
                    uniprot_id = t.find("{http://www.hmdb.ca}uniprot_id")
                    uniprot_id = uniprot_id.text
                    data["subject"]["uniprot_id"]=uniprot_id
                    
                    uniprot_name= t.find("{http://www.hmdb.ca}uniprot_name")
                    uniprot_name = uniprot_name.text
                    data["subject"]["uniprot_name"]=uniprot_name
                    
                    genbank_protein_id= t.find("{http://www.hmdb.ca}genbank_protein_id")
                    data["subject"]["genbank_protein_id"]=genbank_protein_id.text
                    
                    hgnc_id= t.find("{http://www.hmdb.ca}hgnc_id")
                    data["subject"]["hgnc_id"]=hgnc_id.text
                    
                    genbank_gene_id=t.find("{http://www.hmdb.ca}genbank_gene_id")
                    data["subject"]["genbank_gene_id"]=genbank_gene_id.text
                    
                    gene_name = t.find("{http://www.hmdb.ca}gene_name")
                    data["subject"]["gene_name"]=gene_name.text
            

                    data_list.append(data)

            # ---------- Metabolite associations without references ------------     
    
            for met_assc in t.findall("{http://www.hmdb.ca}metabolite_associations"):
                for met_assc_ in met_assc.findall("{http://www.hmdb.ca}metabolite"):
                    for met_assc_id in met_assc_.findall("{http://www.hmdb.ca}accession"):
                            
                            # if the metabolite association was already present above (in metabolite_refereces)
                            # we want to pass adding id to dict to avoid making a duplicate document 
                            pass_assc=False # set bool check for duplicates

                            # Check for duplicate in list, set pass_assc bool to True 
                            for elem in data_list:
                                if met_assc_id.text == elem['object']['accession']:
                                    #print('yes ', accession)
                                    pass_assc = True
                                    
                            # if bool is True pass making duplicate doc       
                            if pass_assc==True: 
                                pass
                            else:
                                print(pass_assc)
                                data={"_id": _id+"_%s"%ct, 'pmid': None, 'subject':{}, 'object':{'accession': met_assc_id.text} }
                                ct+=1
                                data["object"]["kegg_id"]=metabolite_dict[met_assc_id.text]["kegg_id"]
                                data["object"]["chemspider_id"]=metabolite_dict[met_assc_id.text]["chemspider_id"]
                                data["object"]["chebi_id"]=metabolite_dict[met_assc_id.text]["chebi_id"]
                                data["object"]["pubchem_compound_id"]=metabolite_dict[met_assc_id.text]["pubchem_compound_id"]

                                for met_assc_name in met_assc_.findall("{http://www.hmdb.ca}name"):
                                    data["object"]['name'] = met_assc_name.text

                                

                                #setup subject info       
                                #uniprot_id, uniprot_name, genbank_protein_id, hgnc_id, genbank_gene_id, and gene_name.        
                                uniprot_id = t.find("{http://www.hmdb.ca}uniprot_id")
                                uniprot_id = uniprot_id.text
                                data["subject"]["uniprot_id"]=uniprot_id
                                
                                uniprot_name= t.find("{http://www.hmdb.ca}uniprot_name")
                                uniprot_name = uniprot_name.text
                                data["subject"]["uniprot_name"]=uniprot_name
                                
                                genbank_protein_id= t.find("{http://www.hmdb.ca}genbank_protein_id")
                                data["subject"]["genbank_protein_id"]=genbank_protein_id.text
                                
                                hgnc_id= t.find("{http://www.hmdb.ca}hgnc_id")
                                data["subject"]["hgnc_id"]=hgnc_id.text
                                
                                genbank_gene_id=t.find("{http://www.hmdb.ca}genbank_gene_id")
                                data["subject"]["genbank_gene_id"]=genbank_gene_id.text
                                
                                gene_name = t.find("{http://www.hmdb.ca}gene_name")
                                data["subject"]["gene_name"]=gene_name.text
                        
                                data_list.append(data)

        


    for doc_ in data_list:
        yield doc_

    """
    # read in the input file to a dataframe
    #df=pd.read_xml(protein_file)
      
    # --- Clean Data --- 
    # get columns with all null values
    null_cols= list(df.loc[:, df.isna().all()].columns.values)
    data=df.drop(null_cols,axis=1) # remove the null cols
    





    # --- Push Data ---
    # loop through our data row and simply pass the row 
    for index, row in data[:10].iterrows():
        #print("\n ROW %s \n %s"%(index,row))  # check output
        yield row"""
        
    
