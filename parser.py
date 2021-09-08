
import os, pandas, csv, re
import numpy as np
from biothings.utils.dataload import dict_convert, dict_sweep

from biothings import config
logging = config.logger

process_key = lambda k: k.replace(" ","_").lower()

def load_orthology(data, data_list):

    gene1_ids=data['Gene1ID']
    gene2_ids=data["Gene2ID"]

    assert len(gene1_ids) == len(gene2_ids)

    results = {}

    for rec in data_list:
        _id = rec['Gene1ID']
        rec = dict_convert(rec,keyfn=process_key)
        # remove NaN values, not indexable
        rec = dict_sweep(rec,vals=[np.nan])
        results.setdefault(_id,[]).append(rec)
        #print(rec)
    for _id,docs in results.items():
        doc = {"_id": _id, "orthology_data" : docs}
        yield doc
