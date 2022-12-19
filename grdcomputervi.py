# -*- coding: utf-8 -*-
"""
GRDCOMPUTERVI - Chargement des imagettes et calcul du RVI.
        python grdcomputervi.py data
    
    Pour éviter que le MacBook ne se mette en veille pendant
    les téléchargements, il suffit de lancer la commande suivante
    dans un terminal :
    
        caffeinate
    
    Un ^C suffit à arrêter ce mode.

    Note: revoir tdqm pour que cela marche hors notebook.
    
@author: Jérôme Lacaille
"""

__date__ = "2021-04-17"
__version__ = '0.1'

import os
import sys

from tqdm.notebook import tqdm

import esg_deforestation_radar as sentinel

# Description du travail.
assert len(sys.argv)==2, "geocomputervi.py <data-path>"

datadir = sys.argv[1]
listdir = [os.path.join(datadir,d) for d in os.listdir(datadir) 
           if os.path.exists(os.path.join(datadir,d,'geozone.json'))]

for i,target in tqdm(enumerate(listdir), desc="Zones", total=len(listdir)):
    zone = sentinel.getzone(target)
    sentinel.rvicompute(zone, renew=True)
    sentinel.dataset.rvicorrect(zone)
