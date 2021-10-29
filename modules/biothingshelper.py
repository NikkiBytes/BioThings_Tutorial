"""
# Nichollette Acosta
# Package for the reusable methods used when developing for the
# Su and Wu lab. 
"""

from zipfile import ZipFile

class BioThingsHelper:
    
    def get_gene(gene_id):
        gene=gene_client.getgene(gene_id, fields='symbol,name')
        return gene;

    def zipload():
        print("hey")
        # open file
        #with ZipFile(input_folder, 'r') as zip:
         #   zip.printdir()

          #  print("[INFO] Extacting files now....")
           # zip.extractall()
            #print("[INFO] complete.")
