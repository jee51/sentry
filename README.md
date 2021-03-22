# sentry
Des codes de traitement de données satellites Sentinel.

Ces codes permettent de récupérer les données GLS (Global Land Service) du capteur optique OCLI de Sentinel-3 et les données radar GRD de Sentinel-1.
Le code permet de calculer le Radar végétation Index (RVI).


Pour l’instant ce sont plutôt des brouillons, mais ils peuvent servir pour accélérer le travail de récupération et l’automatisation des traitements.

Pour analyser une zone donnée, il faut créer un fichier geozone.json qui contient des références géographiques et les répertoires de stockage et de travail.

L’utilisation d’un répertoire de travail différent du répertoire de stockage permet d’utiliser un disque plus rapide pour les données temporaires.

*	Le nom va de soi. Il sert uniquement aux affichages.
*	La « bbox » correspond aux bornes en latitudes et longitude : (lat0, lon0, lat1, lon1), le point 0 étant le coin supérieur gauche, de latitude et longitude minimales, et le point 1, le coins inférieur droit.
*	« step » est un facteur de décimation de la fenêtre pour le calcul de l’imagette, le rectangle rouge de l’image (Figure 11). On déborde forcément de la zone d’intérêt, mais d’au plus « step » pixels lors du chargement en mémoire.
*	« GRD » est le répertoire dans lequel se trouvent, ou doivent se trouver les archives .zip GRD de Sentinel-1.
*	« FCOVER » donne l’adresse du répertoire dans lequel se situent les fichiers netCDF téléchargés depuis GLS pour l’indicateur FCOVER. Même chose pour les indicateurs « LAI » ou « NDVI » s’ils sont téléchargés.

Pour récupérer des données GRD du serveur « SentielHub » il faut juste exécuter le code

    python grdloadzone.py <adresse du répertoire contenant le fichier json>

Par exemple

    python grdloadzone.py data/Barro-Alto
 
Le chargement s’effectuera automatiquement.

Les notebooks explicatifs sont dans le dossier « notebooks ».
Le fichier « TraceRVI.ipynb » présente des calculs ayant abouti aux graphes présentés dans ce document. Des données CSV sont disponibles dans le répertoire « data/Barro-Alto », ce qui évite d’avoir à relancer les traitements depuis les données brutes. Si vous avez téléchargé des données personnelles il faut juste supprimer le commentaire des lignes « compute… »
