# -*- coding: utf-8 -*-
"""
DATASET - Récupération de l'indicateur RVI depuis S1/SAR.

@author: Jérôme Lacaille
"""

__date__ = "2021-03-19"
__version__ = '0.1'

import os
import json
import numpy as np
import pandas as pd
from tqdm.notebook import tqdm
from .grdzone import GrdZone
from .ocli import get_ncfiles, get_ncmean

ZONEFILE = 'geozone.json'
RVIFILE = 'RVItable.csv'


# Récupération du nom du fichier CSV.
def getzone(tgtdir):
    """ Récupération d'une zone d'intérêt à partir du répertoire. """

    geofile = os.path.join(tgtdir,ZONEFILE)
    with open(geofile,'r') as fp:
        zone = json.load(fp)
    zone['target'] = tgtdir
    
    return zone

# -------------------------------------------------------------------------
def rvicompute(zone,renew=False):
    """ Recréation de la table des RVI. """
    
    # Récupération du nom du fichier CSV.
    
    tgtdir = zone['target']
    rawdir = zone['GRD']
    
    csvfile = os.path.join(tgtdir,RVIFILE)
    if os.path.exists(csvfile) and renew:
        # Destruction du fichier.
        os.remove(csvfile)
    
    if not os.path.exists(csvfile):
        with open(csvfile,'w') as fp:
            fp.write('DATE,RVI,RMSE,INCIDENCE,ELEVATION,' +
                     'DIRECTION,DLAT,DLON,ZIPFILE\n')
    
    # Récupération des éléments existants.
    df = pd.read_csv(csvfile).set_index('DATE')
    
    zipfilelist = [zipfile for zipfile in os.listdir(rawdir)
                   if zipfile.startswith('S1') 
                   and zipfile.endswith('.zip')
                   and not zipfile in list(df['ZIPFILE'])]
    
    for zipfile in tqdm(zipfilelist):
        geo = GrdZone(zone, os.path.join(rawdir,zipfile))
        rvi,std = geo.meanrvi()
        inc = geo.mask['incidence']
        ele = geo.mask['elevation']
        dv = geo.mask['dvec']
        if geo.data[0]['direction'] == "Descending":
            dr = -1
        else:
            dr = 1
        with open(csvfile,'a') as fp:
            fp.write("{},{},{},{},{},{},{},{},{}\n".format(geo.date,rvi,std,
                                                     inc,ele,dr,
                                                     dv[0],dv[1],
                                                     zipfile))
            

# -------------------------------------------------------------------------
def rvi(zone):
    """ Récupere la table de tous les RVI. """
        
    tgtdir = zone['target']
    csvfile = os.path.join(tgtdir,RVIFILE)
    assert os.path.exists(csvfile), "Need to create with 'rvicompute()' before."

    df = pd.read_csv(csvfile).set_index('DATE')
    df.index = pd.to_datetime(df.index)
    df.sort_index(inplace=True)
    return df


# -------------------------------------------------------------------------
def glscompute(zone,indice):
    """ Création du tableau de résultats
        On parcours chaque fichier et on ajoute les résultats calculés.
    """
    
    # Récupération de la liste des données.
    df = get_ncfiles(zone,indice)
    
    FF = []
    EE = []
    PP = []
    F = None
    E = None
    filelist = list(df['FILE'])
    for i,fname in tqdm(enumerate(filelist)):
        Fi,Ei = get_ncmean(zone,fname,indice)
        
        if E is None:
            E = Ei
            F = Fi
        
        F[~np.isnan(Fi)] = Fi[~np.isnan(Fi)]
        E[~np.isnan(Ei)] = Ei[~np.isnan(Ei)]
        
        m0 = np.nanmin(F-E)
        m1 = np.nanmedian(F)
        m2 = np.nanmax(F+E)
    
        n = np.prod(F.shape)
        n0 = np.sum(np.isnan(F))
        f = (np.nansum(F)+n0*m1)/n
        e = (np.nansum(E)+n0*(m2-m0)/2)/n
        n1 = np.sum(np.isnan(Fi))
        p = (n-n1)/n
        
        FF.append(f)
        EE.append(e)
        PP.append(p)
        
    df[indice] = FF
    df['RMSE'] = EE
    df['PR'] = PP
    
    df.to_csv(os.path.join(zone['target'],indice+"table.csv"))

# -------------------------------------------------------------------------
def gls(zone,indice):
    """ Récupération de l'indice si la table existe. """
    
    tgtdir = zone['target']
    csvfile = os.path.join(tgtdir,indice+"table.csv")
    assert os.path.exists(csvfile), "Need to create with 'glscompute()' before."

    df = pd.read_csv(csvfile).set_index('DATE')
    df.index = pd.to_datetime(df.index)
    df.sort_index(inplace=True)
    return df