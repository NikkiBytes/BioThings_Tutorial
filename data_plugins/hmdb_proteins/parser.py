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
def load_hmdb_data(data_folder):
    
    # --- Set Data ---
    
    # set input XML file path
    protein_file = os.path.join(data_folder, "hmdb_proteins.xml")

    # read in the input file to a dataframe
    df=pd.read_xml(protein_file)
      
    # --- Clean Data --- 
    # get columns with all null values
    null_cols= list(df.loc[:, df.isna().all()].columns.values)
    data=df.drop(null_cols,axis=1) # remove the null cols
    
    # --- Push Data ---
    # loop through our data row and simply pass the row 
    for index, row in data[:10].iterrows():
        #print("\n ROW %s \n %s"%(index,row))  # check output
        yield row
        
        
        