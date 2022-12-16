# -*- coding: utf-8 -*-
"""
GRDLOADZONE - Chargement depuis Copernic.

    python grdloadzone "data/Barro-Alto"

@author: Jérôme Lacaille
"""

__date__ = "2021-03-15"
__version__ = '0.1'


import os
import sys
import json
#from zipfile import ZipFile
from xml.etree import ElementTree 

#import numpy as np
#import pandas as pd

from sentinelsat import SentinelAPI, SentinelAPIError, InvalidChecksumError

# Description du travail.
assert len(sys.argv)==2, "geoloadzone <target-path>"

target = sys.argv[1]
print("Target directory:", target)
geofile = os.path.join(target,'geozone.json')
with open(geofile,'r') as fp:
    dx = json.load(fp)
rawdir = dx['GRD']
print("Downloading zone",dx['name'])
print("to directory:",rawdir)

bbox = dx['bbox']
query = '(footprint:"Contains(' +\
    'POLYGON(({1} {0}, {3} {0}, {3} {2}, {1} {2}, {1} {0}))'.format(*bbox) +\
    ')") AND  ((platformname:Sentinel-1 AND producttype:GRD))'

# Ouverture de l'API
api = SentinelAPI(None,None)

# query = '(footprint:"Contains(POLYGON((-49.06 -15.12,-48.92 -15.12,-48.92 -15.02,-49.06 -15.02,-49.06 -15.12)))" ) AND ( beginPosition:[2019-01-01T00:00:00.000Z TO 2021-02-21T23:59:59.999Z] AND endPosition:[2019-01-01T00:00:00.000Z TO 2021-02-21T23:59:59.999Z] ) AND ( (platformname:Sentinel-1 AND producttype:GRD))'

products = api.query(raw=query)
print("On va charger {} produits.".format(len(products)))

# convert to Pandas DataFrame
products_df = api.to_dataframe(products)

# sort and limit to first 5 sorted products
products_df_sorted = products_df.sort_values('ingestiondate', ascending=False)

errorcount = 0
for i in range(len(products_df_sorted)):
    product = products_df_sorted.iloc[i]
    filename = product.title + '.zip'
    pathname = os.path.join(rawdir,filename)
    if os.path.exists(pathname):
        print("{0:03d}: {1} --> already".format(i,product.beginposition.date()))
    else:
        print("{0:03d}: {1} ... loading".format(i,product.beginposition.date()))
        try:
            api.download(product.uuid,rawdir)
        except InvalidChecksumError:
            print("Invalid Checksum.")
        except SentinelAPIError as apierr:
            print("Sentinel error:",apierr)
            break
        except:
            print("Connection problem #{}".format(errorcount))
            errorcount += 1
            if errorcount == 3:
                break
            # On re-tente une connexion
            api = SentinelAPI(None,None)
            i -= 1
            
