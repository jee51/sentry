# -*- coding: utf-8 -*-
"""
GEOFCOVER - Outil d'extraction des données de S3/OCLI.

@author: Jérôme Lacaille
"""

__date__ = "2021-03-19"
__version__ = '0.1'

import os
import re
from netCDF4 import Dataset, num2date
import numpy as np
import pandas as pd
from tqdm.notebook import tqdm


# -------------------------------------------------------------------------
def get_ncfiles(zone,indice):
    """ Récupération de la liste des fichiers d'intérêt
        La fonction `get_ncfiles` récupère la liste des fichiers à traiter.
        Cette fonction ne conserve que les versions les plus récentes pour 
        une date donnée.
    """
    dirname = zone[indice]
    RT = []
    D = []
    N = []
    for froot,fdirs,fnames in os.walk(dirname):
        if len(fnames)>0:
            for fname in fnames:
                if fname.endswith('.nc'):
                    m = re.search('(RT[0-9])_([0-9]+)_',fname)
                    if m is None:
                        rt = 9
                        m = re.search('_([0-9]+)_',fname)
                        d = pd.to_datetime(m.group(1),format='%Y%m%d%H%M')
                    else:
                        rt = int(m.group(1)[2])
                        d = pd.to_datetime(m.group(2),format='%Y%m%d%H%M')
                    if d in D:
                        i = D.index(d)
                        if rt > RT[i]:
                            D[i] = d
                            RT[i] = rt
                            N[i] = os.path.join(froot,fname)
                    else:
                        N.append(os.path.join(froot,fname))
                        D.append(d)
                        RT.append(rt)
    df = pd.DataFrame({'FILE' : N, 'REVISION' : RT},index=pd.to_datetime(D))
    df.index.name = 'DATE'
    df.sort_index(inplace=True)
    return df

# -------------------------------------------------------------------------
def get_ncmean(zone,fname,colname):
    """ Lecture du contenu.
    Le contenu peut changer de nom...
    Les données Sentinel ont toujours une latitude et une longiture. 
    L'erreur est  RMSE mais la valeur étudiée est passée en paramètre.

    L'erreur RMSE est calculée sur chaque parcelle de 333m x 333m, en
    comptabilisant des observations indépendantes, éventuellement réajustées 
    à l'aide des passages précédents. Mais cela ne tiens pas compte des 
    cellules qui ne sont pas mesurées pour des raisons techniques ou
    météorologiques. J'estime dans ce cas l'indice à 1/2 et l'erreur de ±1/2.

    Pour la date, je prends celle du fichier car je ne comprends pas bien
    l'algorithme de calcul de la date des mesures stockée dans le fichier.
    """
    lat0,lon0,lat1,lon1 = zone['bbox']
    
    nc = Dataset(fname,'r')
    
    # Récupération de la longitude.
    lon = nc.variables['lon']
    x0 = np.min(np.where(lon[:]>=lon0))
    x1 = np.min(np.where(lon[:]>=lon1))
    
    # Récupération de la latitude.
    lat = nc.variables['lat']
    y1 = np.min(np.argwhere(lat[:]<=lat0))
    y0 = np.min(np.argwhere(lat[:]<=lat1))
    
    # Récupération de la date des données
    #t = pd.to_datetime(nc.variables['time'][:].data[0],unit='D')
    
    # Récupération des données.
    fc = nc.variables[colname]
    ec = None
    if 'time' in nc.variables:
        # Récupération des révisions de Sentinel 3
        fc = nc.variables[colname][0,y0:y1,x0:x1]
        if 'RMSE' in nc.variables:
            ec = nc.variables['RMSE'][0,y0:y1,x0:x1]
    else:
        # Récupération des collectes de PROBA-V
        fc = nc.variables[colname][y0:y1,x0:x1]
        if 'RMSE' in nc.variables:
            ec = nc.variables['RMSE'][y0:y1,x0:x1]
    F = fc.data
    F[fc.mask] = np.nan
    if ec is not None:
        E = ec.data
        E[ec.mask] = np.nan
    else:
        E = fc.data*0.0
        E[fc.mask] = np.nan
    
    nc.close()
    return F,E

# -------------------------------------------------------------------------
def plot_ncmean(df,colname,title=None,scale='Proportion',reg=False):
    if title is None:
        title = colname
        scale = 'Proportien de' + colname
    fig,ax = plt.subplots()
    x = df.index
    y = df[colname]
    e = df['RMSE']
    df[colname].plot(style='*-.')
    ax.fill_between(x,y-e,y+e,color='g',alpha=0.2)
    if reg:
        d = mdt.date2num(x)
        m, b = np.polyfit(d, y.values, 1)
        plt.plot(x,m*d+b,':r')
        plt.annotate("{:.2f} %/an".format(m*100*365),
                     (np.mean(d),m*np.mean(d)+b),color='r')
    plt.ylabel(scale,fontsize=12)
    plt.title(title,fontsize=15)
    
