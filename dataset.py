# -*- coding: utf-8 -*-
"""
DATASET - Utilisation des données calculées et outils de gestion.

Ce module gère les calculs et affichages des indicateurs.

Le nom du fichier 'geozone.json' est fixé par la variable constante ZONEFILE.

*Versions*

1.2) Rajout des autres sous-répertoires indices dans zone.

@author: Jérôme Lacaille
"""

__date__ = "2021-04-30"
__version__ = '1.2'

import os
import json
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from scipy.stats import kendalltau, pearsonr
from sklearn.cluster import KMeans
import statsmodels.formula.api as smf
from tqdm.notebook import tqdm

import matplotlib.pyplot as plt
import matplotlib.dates as mdt

from .grdzone import GrdZone, RVIFILE # pylint: disable=relative-beyond-top-level
from .ocli import Ocli, glslast # pylint: disable=relative-beyond-top-level

ZONEFILE = 'geozone.json'

# =========================================================================
# RGestion des zones.
# =========================================================================
def getzone(tgtdir):
    """ Récupération d'une zone d'intérêt à partir du répertoire. 

        * Le champ 'target' est rajouté au dictionnaire à partir de
          l'adresse fournie.
        * Le nom de la cible est rajoutée au chemin de stockage des
          données brutes (ZIP): GRD/<target>
    """

    if tgtdir.endswith('/'):
        tgtdir = tgtdir[:-1]
    geofile = os.path.join(tgtdir,ZONEFILE)
    with open(geofile,'r') as fp:
        zone = json.load(fp)
    zone['target'] = tgtdir
    
    # On rajoute le sous répertoire adéquat pour le stockage des GRD.
    zonename = os.path.basename(tgtdir)
    grd = zone['GRD']
    if grd.endswith('GRD') or grd.endswith('GRD/'):
        zone['GRD'] = os.path.join(grd,zonename)
    
    base = os.path.dirname(grd)
    for other in os.listdir(base):
        if other not in zone:
            zone[other] = os.path.join(base,other)

    return zone

# -------------------------------------------------------------------------
def buildzones(csvfile):
    """ Construit automatiquement un sous-répertoire de geozones
        à partir d'un fichier CSV.

        Ce code auxiliaire permet d'extraire des zones d'un fichier CSV
        disposant des champs :
        
            name, lat_min, lon_min, lat_max, lon_max
    """

    df = pd.read_csv(csvfile, sep=";", decimal='.')
    dirname = os.path.dirname(csvfile)
    with open(os.path.join(dirname,'_geozone.json'),'r') as fp:
        zone0 = json.load(fp)

    for i in range(len(df)):
        r = df.iloc[i]
        name = r['name']
        name = name.replace(' ','-')
        newdir = os.path.join(dirname,name)
        if not os.path.exists(newdir):
            # On construit un nouveau répertoire d'accueil.
            bbox = [r.lat_min, r.lon_min, r.lat_max, r.lon_max]
            os.mkdir(newdir)
            zone = zone0
            zone['bbox'] = bbox
            zone['name'] = name
            with open(os.path.join(newdir,'geozone.json'),'w') as fp:
                json.dump(zone,fp)
        

# =========================================================================
# Gestion des acquisitions Sentinel-1 / GRD
# =========================================================================
def rviquery(zone):
    """ Renvoie la requête à appliquer à  Sentinel Hub. 
        Convertie la bbox en une requête "Contains(POLYGON(...))"
    """
    
    bbox = zone['bbox']
    query = '(footprint:"Contains(' +\
        'POLYGON(({1} {0}, {3} {0}, {3} {2}, {1} {2}, {1} {0}))'.format(*bbox) +\
        ')") AND  ((platformname:Sentinel-1 AND producttype:GRD))'

    return query

# -------------------------------------------------------------------------
def rvi(zone):
    """ Récupere la table de tous les RVI. 

        Il s'agit des calculs stockés dans la table RVIFILE par la fonction
        grdzone.rvicompute().
    """
        
    tgtdir = zone['target']
    csvfile = os.path.join(tgtdir,RVIFILE)
    assert os.path.exists(csvfile), "Need to create with 'rvicompute()' before."

    df = pd.read_csv(csvfile).set_index('DATE')
    df.index = pd.to_datetime(df.index)
    df.sort_index(inplace=True)
    return df

# -------------------------------------------------------------------------
def rviplot(zone):
    """ Affichage de la courbe RVI. """
    
    df = rvi(zone)
    df['RVI'].plot(style='-*',legend=True)
    plt.title(zone['name']+" [RVI]")

# -------------------------------------------------------------------------
def rvicorrect(zone,plot=False,renew=False):
    """ Calcul des RVI corrigés et ajout à la table. 
    
        La table RVIFILE contient le calcul initial RVI et des données de
        localisation de la zone étudiée par rapport à la fenêtre GRD acquise.
        Notamment

        * DX donne la hauteur du centre de la zone.
        * DY donne son écart par rapport à la trajectoire du satellite.

        La correction présuppose que la zone est suffosemment petite pour que 
        l'on puisse se contenter d'un seul calcul moyen et pas une modification
        par pixel.

        La correction se passe en six temps :
        
        1. On calcule par un partitionnement simple (Kmean) le nombre de types
           de prises de vue différentes.
        2. On reconstruit par interpolation les points manquants pour chaque type
           de prise de vue. Ce sont des sortes d'enveloppes, en particulier s'il
           n'y a que deux prises.
        3. On construit une prise de vue virtuelle par une moyenne pondérée des
           enveloppes. Cette vue est un RVIM (moyen) utile juste pour les calculs.
        4. On construit l'écart entre le RVI et le RVIM, noté DRVI
        5. On estime une valeur de DRVI en cherchant des coefficients 
           d'interpolation linéaire (avec interaction d'ordre 1) depuis DX et DY.
           Puis on ajoute cet écart à RVI pour obtenir le RVIC (corrigé).
        6. Finalement on lisse RVIC par un filtre exponentiel de valeur 1/4.

        Le RVIC est rajouté à la table RVITABLE. Si cette valeur existait déjà
        et que l'on avait pas demandé de `plot` ou `renew`, elle n'est pas
        recalculée.

        Paramètres :
        
        * `renew`: permet de relancer le calcul de RVIC.
        * `plot`:  affiche les courbes enveloppes et le RVIC.
    """
    
    df = rvi(zone)
    if 'RVIC' in df.columns and not renew\
            and not plot and not df['RVIC'].isna().any():
        return df
    
    R = df.copy()
    
    # Calcul du nombre de clusters.
    DXY = R[['DX','DY']]
    N = 6
    X = np.arange(2,N)
    S = np.ones(N-2)
    for K in range(2,N):
        km = KMeans(K).fit(DXY)
        P = [np.sum(km.labels_==k)/len(DXY) for k in range(K)]
        L = -np.sum(P*np.log(P))
        S[K-2] = L/np.log(K)
    K = X[np.argmax(S)]

    # Calcul des enveloppes.
    km = KMeans(K).fit(DXY)
    P = [np.sum(km.labels_==k)/len(DXY) for k in range(K)]
    LR = []
    if plot:
        R['RVI'].plot(style='-*',legend=True)
    for k in range(K):
        rvi0 = R['RVI'].copy()
        rvi0.loc[km.labels_ != k] = np.nan
        rvi0.interpolate(method='index',inplace=True)
        if plot:
            rvi0.plot(style='k:')
        LR.append(rvi0)

    # Calcul du RVI moyen (qui sera oublié ensuite).
    R['DL'] = km.labels_
    R['RVIM'] = np.nanmean(LR,axis=0) # np.dot(P,LR) # 

    # Estimation de l'eccart.
    R['DRVI'] = R['RVIM']-R['RVI']
    
    # Calcul des coefficients de correction.
    res = smf.ols('DRVI ~ DX*DY', data=R).fit()
    R['RVIC'] = R['RVI']+res.predict(R)

    # Lissage du RVIC.
    df['RVIC'] = R['RVIC'].ewm(alpha=0.25).mean()
    if plot:
        df['RVIC'].plot(style='r-',legend=True)
        plt.title(zone['name']+" [RVIC]")

    # Sauvegarde du RVIC (corrigé)
    tgtdir = zone['target']
    csvfile = os.path.join(tgtdir,RVIFILE)
    df.to_csv(csvfile)
    
    return df

# =========================================================================
# Gestion des acquisitions Sentinel-3 / OCLI / GLS (FCOVER, LAI, NDVI...)
# =========================================================================
def gls(zone,ind):
    """ Récupération de l'indice si la table existe. 
    
        Les données ont été stockées dans la table <ind>table.csv par la fonction
        ocli.glscompute().
    """
    
    tgtdir = zone['target']
    csvfile = os.path.join(tgtdir,ind+"table.csv")
    assert os.path.exists(csvfile), "Need to create with 'glscompute()' before."

    df = pd.read_csv(csvfile).set_index('DATE')
    df.index = pd.to_datetime(df.index)
    df.sort_index(inplace=True)
    return df
     
# -------------------------------------------------------------------------
def glsplot(zone,ind,reg=False):
    """ Affichage des moyennes de l'indicateur avec 
        
        * le tube d'erreur en vert. 
        * un tube rouge un peu plus gros qui est augmenté par le taux de 
          valeus manquantes.
    """
    
    df = gls(zone,ind)
    x = df.index
    y = df[ind]
    e = df['RMSE']
    p = df['PR']
    sigma = np.sqrt(p*e*e + (1-p)*0.25)
    y.plot(style='*-.')
    plt.fill_between(x,y-e,y+e,color='g',alpha=0.2)
    plt.fill_between(x,y-sigma,y+sigma,color='g',alpha=0.2)
    if reg:
        d = mdt.date2num(x)
        m, b = np.polyfit(d, y.values, 1)
        plt.plot(x,m*d+b,':r')
        plt.annotate("{:.2f} %/an".format(m*100*365),
                     (np.mean(d),m*np.mean(d)+b),color='r')
    plt.title(zone['name']+" ["+ind+"]")
    
# =========================================================================
# Comparaison RVI/GLS
# =========================================================================
def rvicompare(zone,ind,plotind=False,plotrvi=False,
                        minrvi=None,maxrvi=None,
                        tube=True):
    """ Compare le RVI à un indicateur GLS.
        
        Il s'agit d'une fonction d'affichage uniquement. 
        Mais elle renvoie aussi un dictionnaire comprenant des calculs de 
        corrélation.

        Paramètres :

        * `plotind`: active l'affichage de l'indicateur GLS.
        * `plotrvi`: active l'affichage de RVIC.
        * `minrvi`et `maxrvi`: activent une normalisation par rapport à
           des valeurs min et max précalculées, par exemple sur un grand
           nombre de zones.
        * `tube`: active l'affichage des tubes autour de l'indicateur GLS.

        Quand RVI et GLS sont activés (`plotrvi` et `plotind`), le RVIC est 
        projeté sur l'indicateur GLS par une régression pondérée par l'erreur 
        en chaque point de l'indicateur. Il est alors possible de superposer
        les deux courbes qui ont la même échelle.

        Retourne un dictionnaire contenant les champs suivants :

        * NAME    : le nom de la zone.
        * PR      : le taux de points du RVIC qui tombent dans le tube rouge à
                    ±1 sigma.
        * R2      : le coefficient de détermination de la régression pondérée.
        * PEASON  : la corrélation de Pearson.
        * KENDALL : le Tau de Kendall (corrélation des rangs).
    """
    R = rvi(zone)
    F = gls(zone,ind)

    # Interpolations
    x = R.index
    Fc = interp1d(pd.to_numeric(F.index),F[ind].values,fill_value="extrapolate")
    fy = Fc(pd.to_numeric(x))
    R[ind] = fy
    Fe = interp1d(pd.to_numeric(F.index),F['RMSE'].values,fill_value="extrapolate")
    fe = Fe(pd.to_numeric(x))
    Fp = interp1d(pd.to_numeric(F.index),F['PR'].values,fill_value="extrapolate")
    fp = Fp(pd.to_numeric(x))
    
    sigma = np.sqrt(fp*fe*fe + (1.0-fp)*0.25)
    res = smf.wls('FCOVER ~ RVIC', weights=1.0/sigma, data=R).fit()
    y = res.fittedvalues
    out = (y >= fy-sigma) & (y <= fy+sigma)
    p = out.sum()/len(out)

    # Affichage
    if plotrvi:
        if plotind:
            y.plot(style='-*')
        else:
            if (minrvi is not None) and (maxrvi is not None):
                y = (R['RVIC']-minrvi)/(maxrvi-minrvi)
                plt.ylim((0.0,1.0))
            y.plot(style='-*')
            plt.title(zone['name']+" [RVI]")
    if plotind:
        plt.plot(x,fy,'-',color='darkgreen')
        if plotrvi and tube:
            plt.fill_between(x,fy-fe,fy+fe,color='g',alpha=0.2)
            plt.fill_between(x,fy-sigma,fy+sigma,color='r',alpha=0.2)
            plt.xlabel("Precision {:.1f}%".format(100.0*p))
            plt.legend(['RVIC', ind])
            plt.ylim((0.0,1.0))
        else:
            if (minrvi is not None) and (maxrvi is not None):
                plt.ylim((0.0,1.1))
            plt.title(zone['name']+ " ["+ind+"]")
    plt.xticks(rotation=30)
        
    return dict(NAME=zone['name'],
                PR=p,
                R2=res.rsquared,
                PEARSON=pearsonr(fy,y)[0],
                KENDALL=kendalltau(fy,y)[0])

# -------------------------------------------------------------------------
def zoneplot(zone,ind):
    """ Affichage des graphes d'une zone. 
    
        Cet affichage combine en un multiplot des graphiques précédents.
        Il est pratique utilisé avec une liste de zones interactive.
    """

    #plt.subplots_adjust(hspace=0.7,wspace=0)

    # Affichage de l'indice.
    plt.subplot(311)
    glsplot(zone,ind)
    plt.xlabel(None)
    
    # Affichage des images.
    plt.subplot(323)
    glslast(zone,ind) # Optique
    plt.title(None)
    plt.ylabel(ind)

    plt.subplot(325)
    geo = GrdZone(zone)
    geo[-1].plot() # Radar
    plt.title(None)
    plt.ylabel('RVI')

    # Affichage de la comparaison.
    plt.subplot(324)
    rvicompare(zone,'FCOVER',plotind=True, plotrvi=False)
    plt.title(None)
    plt.xlabel(None)
    plt.xticks([])

    # Affichage de l'indice RVI.
    plt.subplot(326)
    rvicompare(zone,'FCOVER',plotrvi=True)
    plt.title(None)
    plt.xlabel(None)
