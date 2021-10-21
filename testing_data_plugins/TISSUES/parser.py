import csv, os


def load_tm_hu_data(file_path):
    """
    Load the textmining human data
    """
    with open(file_path) as file:
        json_docs=[]
        fieldnames = ['ensembl', 'symbol', 'doid',
                        'name', 'zscore', 'confidence', 'url']
        reader=csv.DictReader(file, fieldnames, delimiter='\t')
        for row in reader:
                row.pop('url')
                row['symbol'] = row['symbol']#fetch_symbol(row['symbol'])
                if row['symbol'] == None:
                    continue
                # convert string to float
                row['zscore'] = float(row['zscore'])
                row['confidence'] = float(row['confidence'])
                row['category'] = 'textmining'
                json_docs.append(dict(row))
    return json_docs


def load_data(data_folder):
    """
    main load function
    """
    tm_path=os.path.join(data_folder, "human_tissue_textmining_full.tsv")
    json_docs= load_tm_hu_data(tm_path)
    print("Document list length: %s \n%s"%(len(json_docs), json_docs[-1]))

load_data("/Users/nacosta/Documents")
print("\n[INFO] TISSUES data parser complete.\n")