# -*- coding: utf-8 -*-
"""
OCLI - Outil d'extraction des données de Sentinel-3/OCLI.

Ce fichier est composé de trois groupes de codes.

* Le premier groupe comporte deux fonctions pour l'accès au fichier 
  et l'extraction de la zone d'intérêt de l'ensemble du globe.
* Le second groupe de codes définit l'objet Ocli qui permet de 
  manipuler simplement des zone géographiques téléchargées à 
  différentes dates.
* Finalement le troisième groupe correspond à des codes de calcul de
  l'indicateur téléchargé. On produit une moyenne, un écart-type et une
  proportion de valeurs manquantes.

@author: Jérôme Lacaille
"""

__date__ = "2021-04-23"
__version__ = '1.1'

import os
import re
import pickle
from netCDF4 import Dataset, num2date # pylint: disable=no-name-in-module
import numpy as np
import pandas as pd
from tqdm.notebook import tqdm

import matplotlib.pyplot as plt
from tqdm.std import TqdmDeprecationWarning

# Elargissement de la zone de température.
TDEG = 5.0 # Nombre de degrés pris autour du centre de la zone.

# -------------------------------------------------------------------------
def get_ncfiles(rawdir,ind):
    """ Récupération de la liste des fichiers d'intérêt
        La fonction `get_ncfiles` récupère la liste des fichiers à traiter.
        Cette fonction ne conserve que les versions les plus récentes pour 
        une date donnée.

        Chaque téléchargement possède un numéro de révision allant de 0 à 6.
        Les fichiers les plus anciens téléchargés de PROBA-V sont associé au 
        numéro de révision 9.

        Les anciennes révisions ne sont pas supprimées, on pourrait le faire
        pour gagner de la place, mais cela permet éventuellement d'étudier la 
        méthode d'amélioration des acquisitions.

        La date associée à l'acquisition est issue du nom du fichier.
    """

    RT = []
    D = []
    N = []
    for froot,fdirs,fnames in os.walk(rawdir): # pylint: disable=unused-variable
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
    
        L'indicateur recherché se trouve stocké dans une colonne 'colname'
        qui correspond à 'ind' dans la plupart des autres fonctions.
        La distinction est faite ici pour rappeler que le codage peut varier
        d'un indicateur à l'autre. Pour l'isntant ce code a été testé sur les
        indicateurs FCOVER, LAI et NDVI.

        Ce dernier (NDVI), est stocké de manière différente puisque l'erreur
        RMSE n'est pas enregistrée.
        Les données Sentinel ont toujours une latitude et une longiture. 
        L'erreur est RMSE (quand elle est présente).

        L'erreur RMSE est calculée sur chaque parcelle de 333m x 333m, en
        comptabilisant des observations indépendantes, éventuellement réajustées 
        à l'aide des passages précédents. Mais cela ne tiens pas compte des 
        cellules qui ne sont pas mesurées pour des raisons techniques ou
        météorologiques.

        :todo:  Il est possible que des informations horaires soeint ajoutées 
                pour certains indices à l'aide de la dimension `time`. On prendra
                par exemple une moyenne sur l'ensemble des dimensions plutôt que
                se contenter de la première.
    """
    lat0,lon0,lat1,lon1 = zone['bbox']
    
    # Dans le cas de la température on élargie la zone
    if colname == 'TCI':
        lat = (lat0 + lat1)/2.0
        lon = (lon0 + lon1)/2.0
        lat0 = lat - TDEG
        lat1 = lat + TDEG
        lon0 = lon - TDEG
        lon1 = lon + TDEG
        
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
    #fc = nc.variables[colname]
    ec = None
    if 'time' in nc.variables:
        # Récupération des révisions de Sentinel 3
        # :todo: revoir la dimension temporelle.
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


# =========================================================================
# L'objet Ocli permet de gérer les acquisitions de GLS.
# =========================================================================
class Ocli:
    """ L'objet décrivant la zone.
        
        Cet objet peut être itéré automatiquement sur la liste
        des données disponibles (téléchargées).
    """
    
    def __init__(self, zone, ind='FCOVER'):
        """ Initialisation 
        
            L'initialisation a besoin de la zone géographique 
            et du nom de l'indicateur calculé.
        """
        
        self.zone = zone
        self.ind = ind
        self.name = zone['name']
        
        rawdir = zone[ind]
        self.datalist = get_ncfiles(rawdir,ind)
        
        self.F = None
        self.E = None
        self.rev = None
        self.date = None
        self.pos = None
        
        self.datemin = self.datalist.index[0].date()
        self.datemax = self.datalist.index[-1].date()
        
    def __len__(self):
        """ Le nombre de données. """
        return len(self.datalist)
    
    def __getitem__(self, i):
        """ Récupère le  i-ième élément. """
        if i<0:
            i = len(self)+i
        return self.load(i)
        
    def __iter__(self):
        return self[0]
    
    def __next__(self):
        pos = self.pos+1
        if pos==len(self.datalist):
            raise StopIteration
        else:
            return self[pos]
    
    # ---------------------------------------------------------------------   
    def load(self,d):
        """ Chargement d'un fichier GRD. 
        
            Il est possible de retrouver une image depuis son numéro dans
            la liste des acquisitions téléchargées, ou en précisant une date,
            dans ce cas c'est la dernière image acquise avant cette date qui
            est chargée.
        """
        
        df = self.datalist
        if isinstance(d,int):
            # A partir du numéro dans la liste triée par date.
            assert d>=0 and d<len(df),\
                "Integer value must be in range 0 to {}."\
                .format(len(df)) 
            pos = d
        else:
            # A partir de la date.
            d = pd.to_datetime(d)
            assert d >= df.index[0], "Date is too early, min is {}."\
                .format(df.index[0])
            ind = df.index<=d
            pos = sum(ind)-1
        
        if self.pos == pos:
            return self
        
        self.pos = pos
        ncfile = df.iloc[pos].FILE
        self.rev = df.iloc[pos].REVISION
        self.data = get_ncmean(self.zone, ncfile, self.ind)
        self.date = df.index[pos].date()
        
        return self
    
    # -------------------------------------------------------------------------
    def mean(self):
        """ Récupération des stats moyennes. 
        
            Les retours sont :
                'f' : l'indicateur moyen.
                'e' : l'écart-type.
                'p' : la proportion de données non manquantes.

            L'erreur est corrigée en utilisant l'eccart maximal sur les données
            connues quand on a des valeurs manquantes.
            On transmet aussi un indicateur correspondant au taux de valeurs
            non-manquantes.
        """
        
        F,E = self.data
                
        m0 = np.nanmin(F-E)
        m1 = np.nanmedian(F)
        m2 = np.nanmax(F+E)
    
        n = np.prod(F.shape)
        n0 = np.sum(np.isnan(F))
        f = (np.nansum(F)+n0*m1)/n
        e = (np.nansum(E)+n0*(m2-m0)/2)/n
        n1 = np.sum(np.isnan(F))
        p = (n-n1)/n

        return f, e, p
        
    # -------------------------------------------------------------------------
    def plot(self):
        """ Affichage de l'image chargée. """

        assert self.pos is not None, "Nothing to plot."
        
        f = plt.imshow(self.data[0])
        ax = f.axes
        ax.set_title(self.zone['name'] + " [" + self.ind + "] " +
                     str(self.date) + " (rev. {})".format(self.rev),
                     fontsize=12)

# #########################################################################
# Fonction de création de la liste des indicateurs.
# #########################################################################
def glscompute(zone,ind):
    """ Création du tableau de résultats
        On parcours chaque fichier et on ajoute les résultats calculés.
        
        Le résultat est stocké dans une table CSV qui de nom 
            <indice>table.csv.

        Le même calcul sur l'erreur est effectué que pour la méthode mean()
        mais on améliore encore l'image en exploitant les données précédentes.
        Cela peut porter à confusion car on lisse le signal, et on réduit
        l'erreur, mais on garde quand même l'indicateur de précision qui renvoie
        la proportion de points visibles sur l'image.
    """
            
    # Récupération de la liste des données.
    ocli = Ocli(zone,ind)
    df = ocli.datalist.copy()
    
    FF = []
    EE = []
    PP = []
    F = None
    E = None
    name = os.path.basename(zone['target'])
    for i in tqdm(range(len(ocli)),desc=name):
        try:
            Fi,Ei = ocli[i].data
        except:
            # Normalement cette erreur arrive en fin de liste.
            Ei = E
            Fi = F
        
        if E is None:
            E = Ei
            F = Fi
        
        # On garde les ancien points visibles (attention !)
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
        
    df[ind] = FF
    df['RMSE'] = EE
    df['PR'] = PP
    
    # On sauvegarde la table.
    df.to_csv(os.path.join(zone['target'],ind+"table.csv"))
    
    # On sauvegarde la dernière vue.
    imgfile = os.path.join(zone['target'],ind+"last.pkl")
    with open(imgfile,'wb') as fp:
        pickle.dump(dict(FCOVER=F,RMSE=E,date=ocli.date),fp)

# -------------------------------------------------------------------------
def glslast(zone,ind):
    """ Affiche la dernière image lissée. 

        Comme on a pris soin de sauvegarder la dernière image lissée
        il sera possible de l'afficher.

        En fait cette sauvegarde devrait permettre de reprendre le
        lissage quand on rajoute des données, mais cela ne marche pas vraiment
        à cause du système de révision. Le mieux serait de garder aussi la
        dernière révision 6.
    """
    
    imgfile = os.path.join(zone['target'],ind+"last.pkl")
    assert os.path.exists(imgfile), "Need to compute before."
    
    with open(imgfile,'rb') as fp:
        d = pickle.load(fp)
        
    plt.imshow(d[ind])
    plt.title(zone['name'] + " [" + ind + "] last:" +
                 str(d['date']))
