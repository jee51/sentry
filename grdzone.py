# -*- coding: utf-8 -*-
"""
GRDZONE - Un objet contenant un réctangle géographique dans une zones.

Le rectangle est issu d'une image de balayage rectangulaire définie 
par une acquisition GRD SAFE. Le rectangle géographique est construit
à partir de bornes en latitude et longitude.

Comme l'orbite du satellite est en biais par rapport au méridiens, la
zone recherchée n'est pas un rectangle parallèle au bords de l'image
acquise. On va donc créer un amsque booléen pour préciser les pixels
appartenant à notre zone.

La résolution de l'image radar étant importante, on utilise un
échantillonnage défini par un facteur de réduction 'STEP'. Finalement
on ne conservera qu'une sous-image contenant notre zone géographique
à STEP pixels près.

Ce fichier est décomposé en trois parties :

* D'abord des fonctions de manipulation des métadonnées XML et des
  images TIFF.
* Puis un objet GRdZone qui se charge de gérer les acquisitions et 
  propose un itérateur sur les images successives.
* Finalement des codes de calcul et de sauvegarde de la série temporelle
  RVI.

Les indicateurs RVI finalement calculés sont stockés dans une table CSV
dont le nom est donné par la constante RVIFILE.

@author: Jérôme Lacaille
"""

__date__ = "2021-04-23"
__version__ = '1.1'

import os
import re
import shutil
import json
import pickle
from zipfile import ZipFile
from xml.etree import ElementTree 

import datetime as dt
import numpy as np
import pandas as pd
from scipy.interpolate import griddata
import statsmodels.formula.api as smf
from tqdm.notebook import tqdm

from osgeo import gdal

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

RVIFILE = 'RVItable.csv'
KEEPSAFE = False

# #########################################################################
# Des fonctions pratiques.
# #########################################################################
def get_grdfiles(rawdir):
        """ Renvoie une liste de données disponibles 
            sous la forme d'une table.

            On scanne toutes les acquisitions effectuées qui se trouvent
            dans des archives ZIP du répertoire 'rawdir'.

            La fonction renvoie une table des noms d'archives indexée
            par la date d'acquisition.

            TODO: J'ai considéré qu'une seule image ne pouvait être acquise
                  par jour et je ne stocke que la date sans l'heure, mais 
                  avec 2 satellites, on peut avoir deux images, du coup il
                  faudra corriger cette bêtise.
        """
        
        N = []
        D = []
        zipdir = os.path.join(rawdir,'ZIP')
        zipfilelist = [zipfile for zipfile in os.listdir(zipdir)
                       if zipfile.startswith('S1') 
                       and zipfile.endswith('.zip')]
        for zipfile in zipfilelist:
            m = re.search('_([0-9]+T[0-9]+)_',zipfile)
            d = pd.to_datetime(m.group(1))
            D.append(d.date())
            N.append(os.path.join(zipdir,zipfile))
            
        df = pd.DataFrame({'FILE': N},index=pd.to_datetime(D))
        df.index.name = 'DATE'
        df.sort_index(inplace=True)
        return df

# -------------------------------------------------------------------------
def tif_to_xmlname(tifname):
    """ Récupère un XML associé au TIF. """

    root,name = os.path.split(tifname)
    root,_local = os.path.split(root)
    radical,_ext = os.path.splitext(name)
    xmlname = os.path.join(root,'annotation',radical+'.xml')
    return xmlname

# ---------------------------------------------------------------------
def xmlinfos(xmlfilename):
    """ Lecture des informations d'une image. 

        C'est dans cette fonction qu'il est possible d'ajouter des 
        éléments supplémentaires décrivant l'image à partir de ses
        métadonnées.
    """
    
    root = ElementTree.parse(xmlfilename)
    
    # La polarisation.
    for e in root.iter('polarisation'):
        pl = e.text
        break
    # La date d'acquisition.
    for e in root.iter('startTime'):
        date = dt.datetime.fromisoformat(e.text)
        break
    # La direction.
    for e in root.iter('pass'):
        dr = e.text
        break
        
    # On récupère l'imageInformation, n'apparait qu'une fois.
    # On s'ensert ensuite comme racine, c'est plus rapide.
    for imageinfo in root.iter('imageInformation'):
        break
    # Lignes
    for e in imageinfo.iter('numberOfLines'):
        h = int(e.text)
        break
    # Pixels
    for e in imageinfo.iter('numberOfSamples'):
        w = int(e.text)
        break
        
    M = [(int(e.find('line').text),int(e.find('pixel').text),
          float(e.find('latitude').text), float(e.find('longitude').text),
          float(e.find('incidenceAngle').text),
          float(e.find('elevationAngle').text) 
         )
         for e in root.iter('geolocationGridPoint')]
    code = pd.DataFrame(M,columns=['line','pixel','lat','lon',
                                   'incidence','elevation'])
    desc = dict(date=date, polarisation=pl, direction=dr, 
                heigth=h, width=w)
    return desc, code
    
# ---------------------------------------------------------------------
def llmask(zone,desc,code):
    """ Création du masque des pixels de la zone sur la sous-image. """
    
    # Calcul de la bbox minimale.
    heigth = desc['heigth']
    width = desc['width']

    # On recherche sur un box plus fruste.
    STEP = zone['step']
    grid_x, grid_y = np.mgrid[0:heigth-1:STEP, 0:width-1:STEP]
    points = code[['line','pixel']].values
    latitude = code['lat'].values
    Lat = griddata(points, latitude, (grid_x, grid_y), method='linear')
    longitude = code['lon'].values
    Lon = griddata(points, longitude, (grid_x, grid_y), method='linear')

    # Le masque
    lat0,lon0,lat1,lon1 = zone['bbox']
    Mask = (Lat>=lat0) & (Lat<=lat1) & (Lon>=lon0) & (Lon<=lon1)

    dlat = (lat1+lat0)/2 - (latitude.max()+latitude.min())/2
    dlon = (lon1+lon0)/2 - (longitude.max()+longitude.min())/2
    
    ux = np.argwhere(Mask.any(axis=0))
    minx = int(max(STEP*(ux.min()-STEP),0))
    maxx = int(min(STEP*(ux.max()+STEP),width))

    uy = np.argwhere(Mask.any(axis=1))
    miny = int(max(STEP*(uy.min()-STEP),0))
    maxy = int(min(STEP*(uy.max()+STEP),heigth))

    # Le masque sur l'extraction
    grid_x, grid_y = np.mgrid[miny:maxy, minx:maxx]
    lat = griddata(points, latitude, (grid_x, grid_y), method='linear')
    lon = griddata(points, longitude, (grid_x, grid_y), method='linear')
    belong = (lat>=lat0) & (lat<=lat1) & (lon>=lon0) & (lon<=lon1)
    
    # Calcul de l'angle moyen d'incidence.
    incidence = code['incidence'].values
    incgrid = griddata(points, incidence, (grid_x, grid_y), method='linear')
    inc = np.nanmean(incgrid[belong])
    
    # Calcul de l'angle moyen d'élévation.
    elevation = code['elevation'].values
    elegrid = griddata(points, elevation, (grid_x, grid_y), method='linear')
    ele = np.nanmean(elegrid[belong])
    
    dx = (maxx+minx)/2 - width/2
    dy = (maxy+miny)/2 - heigth/2
    
    mask = dict(bbox=(minx, miny, maxx, maxy),
                incidence=inc, elevation=ele,
                dvec=(dlat,dlon), dxy=(dx,dy),
                lat=lat, lon=lon, belong=belong)

    return mask
    
# ---------------------------------------------------------------------
def extractdata(zone, zipfile):
    """ Extraction des données du fichier ZIP. """

    # Sauvegarde temporaire.
    tgtdir = os.path.join(zone['target'],'SAFE')
    
    # Création du sous-répertoire SAFE.
    if not os.path.exists(tgtdir):
        os.mkdir(tgtdir)

    try:
        zf = ZipFile(zipfile,'r')
    except:
        print(zipfile)
        raise
    
    basename = os.path.basename(zipfile)[:-4] # On enlève '.zip'
    pathdir = os.path.join(tgtdir,basename+'.SAFE')
    if not os.path.exists(pathdir):
        os.mkdir(pathdir)

    tifs = [name for name in zf.namelist() if name.endswith('.tiff')]

    data = []
    mask = None
    for tiffile in tifs:
        # Annotations.
        xmlname = tif_to_xmlname(tiffile)
        xmlfilename = os.path.join(tgtdir,xmlname)
        if not os.path.exists(xmlfilename):
            xmlfilename = zf.extract(xmlname,tgtdir)
        desc,code = xmlinfos(xmlfilename)

        # Suppression du fichier XML.
        if not KEEPSAFE:
            os.remove(xmlfilename)
        
        # Analyse des coordonnées.
        if mask is None:
            # Les deux polarisations ont le même masque.
            mask = llmask(zone, desc, code) 
            
        # Gestion des données numériques
        tiffilename = os.path.join(tgtdir,tiffile)
        if not os.path.exists(tiffilename):
            tiffilename = zf.extract(tiffile,tgtdir)

        # Taille de l'image complète.
        width = desc['width']
        heigth = desc['heigth']
        
        # Récupération de la bande.
        ds = gdal.Open(tiffilename,gdal.GA_ReadOnly)
        assert ds.RasterCount == 1, "Une seule bande sur GRD"
        assert ds.RasterXSize == width,\
            "Largeur d'image lue {} vs. {} attendue.".\
            format(ds.RasterXSize,width)
        assert ds.RasterYSize == heigth,\
            "Hauteur d'image lue {} vs. {} attendue.".\
            format(ds.RasterYSize,heigth)
        band = ds.GetRasterBand(1)
        
        # Position de la zone.
        minx,miny,maxx,maxy = mask['bbox'] 
        img = band.ReadAsArray(minx,miny,maxx-minx,maxy-miny)
        desc['image'] = img
        
        
        # Sauvegarde du masque
        data.append(desc)

        # Supression du fichier TIFF.
        if not KEEPSAFE:
            os.remove(tiffilename)
            
    zf.close()

    # Nettoyage intermédiaire.
    if not KEEPSAFE:
        try:
            shutil.rmtree(pathdir)
        except OSError as e:
            print("Cannot clear SAFE directory {}:{}".format(tgtdir,e.strerror))

    return code,mask,data
    
# ---------------------------------------------------------------------
def imgnorm(img, p):
    """ Normalise une image entre deux quantiles. """
    
    if p>0:
        q0,q1 = np.quantile(img,[p,1-p])
    else:
        q0,q1 = img.min(), img.max()
    delta = float(q1-q0)
    fimg = (1.0*img-q0)/delta
    fimg[fimg>1.0] = 1.0
    fimg[fimg<0.0] = 0
    return fimg

def rviband(data):
    """ Convertion de vv et vh en rvi. 
        Sur le S1B les bandes sont hh et hv.
    """
    
    bands = [desc['polarisation'] for desc in data]
    if 'VV' in bands:
        ixy = bands.index('VH')
        ixx = bands.index('VV')
    else:
        ixy = bands.index('HV')
        ixx = bands.index('HH')
 
    xy = data[ixy]['image']
    xx = data[ixx]['image']
    
    # Passer en float et éviter les erreurs d'arrondis.
    xx = xx/10000.0
    xy = xy/10000.0

    # Passage en linéaire.
    xy = 10**(-xy)
    xx = 10**(-xx)
    
    rvi = 4.0*xy/(xx+xy)
 
    # Correction vue dans les codes en lignes.
    #dop = xx/(xx+xy)
    #rvi = np.sqrt(dop)*rvi

    return rvi
    
# =========================================================================
# L'objet GrdZone permet de gérer les acquisitions GRD de Sentinel-1.
# =========================================================================
class GrdZone:
    """ L'objet décrivant la zone."""
    
    def __init__(self, zone):
        """ Initialisation 
        """
        
        self.zone = zone
        self.name = zone['name']
        
        rawdir = zone['GRD']
        self.datalist = get_grdfiles(rawdir)
        
        self.code = None
        self.mask = None
        self.data = None
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
        """ Chargement d'un fichier GRD. """
        
        df = self.datalist
        if isinstance(d,int):
            assert d>=0 and d<len(df),\
                "Integer value must be in range 0 to {}."\
                .format(len(df)) 
            pos = d
        else:
            d = pd.to_datetime(d)
            assert d >= df.index[0], "Date is too early, min is {}."\
                .format(df.index[0])
            ind = df.index<=d
            pos = sum(ind)-1
        
        if self.pos == pos:
            return self
        
        if not self.restore(pos):
            zipfile = df.iloc[pos].FILE
            code, mask, data = extractdata(self.zone, zipfile)
            self.pos = pos
            self.code = code
            self.mask = mask
            self.data = data
            self.date = data[0]['date'].date()
            self.dump()

        return self
        
    # ---------------------------------------------------------------------   
    def dump(self):
        """ Sauvegarde d'une image geo. """

        tgtdir = self.zone['target']
        savedir = os.path.join(tgtdir,'GEO')
        if not os.path.exists(savedir):
            os.mkdir(savedir)
        filename = os.path.join(savedir,"GEO"+str(self.date)+".pkl")

        with open(filename,'wb') as fp:
            pickle.dump(dict(pos=self.pos,
                             date=self.date,
                             code=self.code,
                             mask=self.mask,
                             data=self.data),fp)
    
    # ---------------------------------------------------------------------
    def restore(self,pos):
        """ Récupération d'un contenu. """
        
        date = self.datalist.index[pos].date()
        tgtdir = self.zone['target']
        savedir = os.path.join(tgtdir,'GEO')
        filename = os.path.join(savedir,"GEO"+str(date)+".pkl")
        ret = False
        if os.path.exists(filename):
            with open(filename,'rb') as fp:
                g = pickle.load(fp)
                self.data = g['data']
                self.mask = g['mask']
                self.code = g['code']
                self.date = g['date']
                self.pos = pos
                ret = True
        
        return ret

    # ---------------------------------------------------------------------   
    def band(self,band='RVI'):
        """ Récupération d'une matrice de bande VV ou VH. 
            On peut aussi passer RVI pour récupérer l'indicateur de végétation.
            
            :param band: VH ou VV ou RVI.
            :return img: une image rectangulaire encadrant la zone.
        """
        
        assert self.pos is not None, "Need to load before."
        
        if band == 'RVI':
            img = rviband(self.data)
        else:
            bands = [desc['polarisation'] for desc in self.data]
            assert band in bands, "Band {} does not exist.".format(band)
            
            i = bands.index(band)
            img = self.data[i]['image']

        return img
        
    # ---------------------------------------------------------------------
    def bandlist(self):
        """ Renvoie la liste des bandes disponibles. """
        return [desc['polarisation'] for desc in self.data] + ['RVI']

    # ---------------------------------------------------------------------   
    def plot(self,band='RVI', p=0.01):
        """ Affichage de l'image. 
        
            :param band: VH ou VV, HH ou HV ou RVI.
        """
        
        assert self.pos is not None, "Nothing to plot."
            
        img = self.band(band)
        
        fimg = imgnorm(img, p)
        #fig, ax = plt.subplots()
        f = plt.imshow(fimg[:,::-1])
        ax = f.axes
        
        belong = self.mask['belong']
        ax.imshow(belong[:,::-1], alpha=0.2)
        ax.set_title(self.name + " [" + band + "] " + str(self.date),
                     fontsize=12)
        

    # ---------------------------------------------------------------------   
    def meanrvi(self):
        """ Calcule le RVI moyen et l'écart-type """
        
        assert self.pos is not None, "Nothing to compute."
        
        img = rviband(self.data)
        belong = self.mask['belong']
        values = img[belong]
        rvi = values.mean()
        std = values.std()
        
        return rvi,std
    
    # ---------------------------------------------------------------------   
    def plotmask(self,add=None,alpha=0.1):
        """ Affichage du masque. """
        
        # Calcul de la bbox minimale.
        heigth = self.data[0]['heigth']
        width = self.data[0]['width']

        # On recherche sur un box plus fruste.
        STEP = self.zone['step']
        grid_x, grid_y = np.mgrid[0:heigth-1:STEP, 0:width-1:STEP]
        points = self.code[['line','pixel']].values
        latitude = self.code['lat'].values
        Lat = griddata(points, latitude, (grid_x, grid_y), method='linear')
        longitude = self.code['lon'].values
        Lon = griddata(points, longitude, (grid_x, grid_y), method='linear')

        # Le masque
        lat0,lon0,lat1,lon1 = self.zone['bbox']
        Mask = (Lat>=lat0) & (Lat<=lat1) & (Lon>=lon0) & (Lon<=lon1)

        ux = np.argwhere(Mask.any(axis=0))
        minx = int(max(STEP*(ux.min()-STEP),0))
        maxx = int(min(STEP*(ux.max()+STEP),width))

        uy = np.argwhere(Mask.any(axis=1))
        miny = int(max(STEP*(uy.min()-STEP),0))
        maxy = int(min(STEP*(uy.max()+STEP),heigth))

        wx = maxx-minx
        wy = maxy-miny
        
        _fig, ax = plt.subplots()
        cb = None
        if add is 'Lon':
            cb = ax.imshow(Lon[:,::-1])
            ax.imshow(Mask[:,::-1],alpha=alpha)
        elif add is 'Lat':
            cb = ax.imshow(Lat[:,::-1])
            ax.imshow(Mask[:,::-1],alpha=alpha)
        else:
            add = 'Mask'
            ax.imshow(Mask[:,::-1])
        rect = Rectangle(((width-maxx)/STEP,miny/STEP),wx/STEP, wy/STEP,
                         fill=False,edgecolor='r')
        ax.add_patch(rect)
        if cb is not None:
            plt.colorbar(cb)
        ax.set_xlabel("x {} (m)".format(STEP))
        ax.set_ylabel("x {} (m)".format(STEP))
        ax.set_title(self.name + " [" + add + "] " + str(self.date),
                     fontsize=12)


# #########################################################################
# Des fonctions de gestion des listes d'images radar.
# #########################################################################
def rvicompute(zone,renew=False):
    """ Recréation de la table des RVI. 
    
        Le résultat est stocké dans le fichier RVItable.csv dont le nom
        est donné par la constante RVIFILE.

        J'ai choisi d'utiliser un format CSV (;.) plutôt que le format
        standard car du coup il est automatiquement converti quand on 
        l'ouvre sous EXCEL.
    """
    
    # Récupération du nom du fichier CSV.
    
    tgtdir = zone['target']
    
    csvfile = os.path.join(tgtdir,RVIFILE)
    if os.path.exists(csvfile) and renew:
        # Destruction du fichier.
        os.remove(csvfile)
    
    if not os.path.exists(csvfile):
        with open(csvfile,'w') as fp:
            fp.write('DATE,RVI,RMSE,INCIDENCE,ELEVATION,' +
                     'DIRECTION,DLAT,DLON,DX,DY,ZIPFILE\n')
    
    # Récupération des éléments existants.
    df = pd.read_csv(csvfile).set_index('DATE')
    
    geo = GrdZone(zone)
    name = os.path.basename(tgtdir)
    for i in tqdm(range(len(geo)), desc=name):
        filename = geo.datalist.iloc[i].FILE
        zipfile = os.path.basename(filename)
        if zipfile not in list(df['ZIPFILE']):
            geo.load(i)
            rvi,std = geo.meanrvi()
            inc = geo.mask['incidence']
            ele = geo.mask['elevation']
            dv = geo.mask['dvec']
            dxy = geo.mask['dxy']
            if geo.data[0]['direction'] == "Descending":
                dr = -1
            else:
                dr = 1
            with open(csvfile,'a') as fp:
                fp.write("{},{},{},{},{},{},{},{},{},{},{}\n".format(geo.date,
                                            rvi,std,inc,ele,dr,
                                            dv[0],dv[1],dxy[0],dxy[1],
                                            zipfile))
