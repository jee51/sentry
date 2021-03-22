# -*- coding: utf-8 -*-
"""
GEOZONE - Un objet contenant un réctangle géographique dans une zones.

Le rectangle est issu d'une image de balayage rectangulaire définie 
par une acquisition GRD SAFE. Le rectangle géographique est construit
à partir de bornes en latitude et longitude.

@author: Jérôme Lacaille
"""

__date__ = "2021-03-14"
__version__ = '0.1'

import os
import json
from zipfile import ZipFile
from xml.etree import ElementTree 

import datetime as dt
import numpy as np
import pandas as pd
from scipy.interpolate import griddata

from osgeo import gdal

import matplotlib.pyplot as plt
import matplotlib.dates as mdt
from matplotlib.patches import Rectangle


###########################################################################
#%% Des fonctions pratiques.

# -------------------------------------------------------------------------
def tif_to_xmlname(tifname):
    """ Récupère un XML associé au TIF """

    root,name = os.path.split(tifname)
    root,local = os.path.split(root)
    radical,ext = os.path.splitext(name)
    xmlname = os.path.join(root,'annotation',radical+'.xml')
    return xmlname

# ---------------------------------------------------------------------
def xmlinfos(xmlfilename):
    """ Lecture des informations d'une image. """
    
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
    """ Revoie les matrices de latitude, longitude et le masque. """
    
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
    inc = np.nanmean(incgrid*belong)
    
    # Calcul de l'angle moyen d'élévation.
    elevation = code['elevation'].values
    elegrid = griddata(points, elevation, (grid_x, grid_y), method='linear')
    ele = np.nanmean(elegrid*belong)
    
    mask = dict(bbox=(minx, miny, maxx, maxy),
                incidence=inc, elevation=ele,
                dvec=(dlat,dlon),
                lat=lat, lon=lon, belong=belong)

    return mask
    
# ---------------------------------------------------------------------
def extractdata(zone, zipfile):
    """ Extraction des données du fichier ZIP. """

    # Sauvegarde temporaire.
    tgtdir = os.path.join(zone['target'],'GRD')
    
    try:
        zf = ZipFile(zipfile,'r')
    except:
        print(zipfile)
        raise
    
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
            
    zf.close()

    return mask,data
    
    
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

def rviband(vh,vv):
    """ Convertion de vv et vh en rvi. """
      
    vh = 10**(-vh/10000)
    vv = 10**(-vv/10000)
    #rvi = 4.0*vh/(vv+vh)
    
    #dop = vv/(vv+vh)
    rvi = 4.0*vh/(vv+vh)
    
    return rvi
    
###########################################################################
#%% Un affichage interactif permettant de visualiser le contenu de la liste.
class GrdZone:
    """ L'objet décrivant la zone."""
    
    def __init__(self, zone, zipfile):
        """ Initialisation 
        
            :param bbox: (lat0, lon0, lat1, lon1)
            :param filename: zip filename ou SAFE dirname
            :param tgtrep: répertoire de copie des données décompressées
            
            Par défaut s'il n'est pas précisé, le répertoire de sortie 
            est le même que le répertoire des données comprimées.
        """
        
        self.zipfile = zipfile # Fichier ZIP.
        self.zone = zone
        self.name = zone['name']
        
        mask, data = extractdata(zone, zipfile)
        self.mask = mask
        self.data = data
        self.date = data[0]['date'].date()
        
        
    def getband(self,band):
        """ récupération d'une matrice de bande VV ou VH. """
        
        bands = [desc['polarisation'] for desc in self.data]
        ivh = bands.index('VH')
        ivv = bands.index('VV')
        vh = self.data[ivh]['image']
        vv = self.data[ivv]['image']
        if band == 'VH':
            img = vh
        elif band == 'VV':
            img = vv
        elif band == 'RVI':
            img = rviband(vh,vv)
            
        return img
        
    def plot(self,band='RVI', p=0.01):
        """ Affichage de l'image. """
         
        img = self.getband('RVI')
        
        fimg = imgnorm(img, p)
        #fig, ax = plt.subplots()
        f = plt.imshow(fimg[:,::-1])
        ax = f.axes
        
        belong = self.mask['belong']
        ax.imshow(belong[:,::-1], alpha=0.2)
        ax.set_title(band,fontsize=12)
        
    def meanrvi(self):
        """ Calcule le RVI moyen. """
        
        bands = [desc['polarisation'] for desc in self.data]
        ivh = bands.index('VH')
        ivv = bands.index('VV')
        vh = self.data[ivh]['image']
        vv = self.data[ivv]['image']
        
        # Passage db->linéaire
        
        
        img = rviband(vh,vv)
        belong = self.mask['belong']
        values = img[belong]
        rvi = values.mean()
        std = values.std()
        
        return rvi,std