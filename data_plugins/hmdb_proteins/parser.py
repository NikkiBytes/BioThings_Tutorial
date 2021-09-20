import pandas as pd
import numpy as np
from biothings.utils.dataload import dict_convert, dict_sweep
from biothings import config
logging = config.logger

process_key = lambda k: k.replace(" ","_").lower()

#def setup_release(self):
    #release="2021-09-20"
    #return release

# Association Parser 
def load_hmdb_data():
    
    # --- Setup Data ---
    
    # set XML file path
    protein_file = os.path.join("c:\\Users\\19802\\Documents\\dev\scripps\\BioThings_SuLab\\data\\hmdb_proteins.xml")

    # we're reading in an XML
    df=pd.read_xml(protein_zipfolder)  
    results = {}
    # loop through our data rows and simply pass the row 
    for index, row in data[:4].iterrows():
        #print("\n ROW %s \n %s"%(index,row))  # check output
        yield row
        
        
        