import os, pandas, csv, re
import math

from biothings import config
from biothings.utils.dataload import dict_convert
logging = config.logger


def load_orthology(data_folder):
    # the path points to the target data
    infile = os.path.join(data_folder, "ORTHOLOGY-ALLIANCE_COMBINED.tsv")
    assert os.path.exists(infile)


