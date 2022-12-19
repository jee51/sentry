# -*- coding: utf-8 -*-
"""
GRDLOADZONE - Chargement depuis Copernic.

        python grdloadzone.py data
    
    Les login et password sont dans ~/.netrc :

        machine apihub.copernicus.eu
        login *********
        password *********

    Par contre, de temps en temps cela ne marche ps et il faut alors remplir directement les informations dans l'appel de SentinelAPI().
    Il faut que les sous-répertoires <target> soient dans <data> 
    et qu'ils possèdent leur fichier geozone.json.

    Ces sous-répertoires peuvent être créés par dataset.buildzone
    depuis un fichier CSV.

    Pour éviter que le MacBook ne se mette en veille pendant
    les téléchargements, il suffit de lancer la commande suivante
    dans un terminal :
    
        caffeinate
    
    Un ^C suffit à arrêter ce mode.

@author: Jérôme Lacaille
"""

__date__ = "2021-04-30"
__version__ = '1.1'


import os
import netrc
import sys
import time
import json
from zipfile import ZipFile, BadZipFile
from xml.etree import ElementTree 

from esg_deforestation_radar.dataset import getzone, rviquery

from sentinelsat import SentinelAPI, SentinelAPIError, InvalidChecksumError

# Description du travail.
assert len(sys.argv)==2, "python grdloadzone.py data"

datadir = sys.argv[1]


# =========================================================================
def startSentinelAPI():
    """ Ouverture de l'API Sentinel HUB.

        Je ne sais pas pourquoi de temps en temps cela ne marche pas avec '.netrc'. Dans ce cas il faut recherger les login et mot de passe en direct.
    """
    #api = SentinelAPI(None,None)
    net = netrc.netrc()
    auth = net.authenticators('apihub.copernicus.eu')
    user = auth[0]
    password = auth[2]
    api = SentinelAPI(user,password,'https://apihub.copernicus.eu/apihub')

    return api

# =========================================================================
# Chargement des données.
# =========================================================================
def grdloadtarget(target):
    """ Chargement des données d'une zone. """

    print("Target directory:", target)
    zone = getzone(target)
    rawdir = os.path.join(zone['GRD'],'ZIP')
    print("Downloading zone",zone['name'])
    print("to directory:",rawdir)

    # Récupération de la requête pour la zone choisie.
    # Un exemple de requête.
    # query = '(footprint:"Contains(POLYGON((-49.06 -15.12,-48.92 -15.12,-48.92 -15.02,-49.06 -15.02,-49.06 -15.12)))" ) AND ( beginPosition:[2019-01-01T00:00:00.000Z TO 2021-02-21T23:59:59.999Z] AND endPosition:[2019-01-01T00:00:00.000Z TO 2021-02-21T23:59:59.999Z] ) AND ( (platformname:Sentinel-1 AND producttype:GRD))'
    query = rviquery(zone)

    # Ouverture de l'API
    api = startSentinelAPI()

    # Récupération de tous ls produits disponibles.
    products = api.query(raw=query)
    print("On va charger {} produits.".format(len(products)))

    # Convert to Pandas DataFrame
    products_df = api.to_dataframe(products)

    # Tris de produits du plus récent au plus ancien.
    products_df_sorted = products_df.sort_values('ingestiondate', ascending=False)

    errorcount = 0
    for i in range(len(products_df_sorted)):
        product = products_df_sorted.iloc[i]
        zipfilename = product.title + '.zip'
        pathname = os.path.join(rawdir,zipfilename)
        if os.path.exists(pathname):
            print("{0:03d}: {1} --> already in {2}".format(i,product.beginposition.date(),zipfilename))
        else:
            # Création des sous-répertoires.
            if not os.path.exists(zone['GRD']):
                os.mkdir(zone['GRD'])
            if not os.path.exists(rawdir):
                os.mkdir(rawdir)

            print("{0:03d}: {1} ... loading {2}".format(i,product.beginposition.date(),zipfilename))
            zf = None
            try:
                api.download(product.uuid,rawdir)
                # On teste si le fichier est bon.
                zf = ZipFile(pathname,'r')
                zf.close()
                errorcount = 0
            except InvalidChecksumError:
                print("Invalid Checksum.")
            except SentinelAPIError as apierr:
                print("Sentinel error:",apierr)
                break
            except BadZipFile:
                print('Bad archve.')
                baddir = os.path.join(rawdir,'bad')
                if not os.path.exists(baddir):
                    os.mkdir(baddir)
                os.replace(pathname,os.path.join(baddir,zipfilename))
            except KeyboardInterrupt:
                sys.exit(0)
            except:
                print("Connection problem #{}".format(errorcount))
                errorcount += 1
                if errorcount == 2:
                    break # On sort de la fonction.
                time.sleep(30)
                # On re-tente une connexion
                api = startSentinelAPI()
                i -= 1

# #########################################################################
# Exécution du code.
# #########################################################################

listdir = [os.path.join(datadir,d) for d in os.listdir(datadir) 
           if os.path.exists(os.path.join(datadir,d,'geozone.json'))]

# Test.
#grdloadtarget(listdir[0])

count = 3
while count>0:
    count -= 1
    for target in listdir:
        try:
            grdloadtarget(target)
        except KeyboardInterrupt:
            sys.exit(0)
        except:
            print('New tentative')
