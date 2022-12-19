# ESG_DEFORESTATION_RADAR

Des codes de traitement de données satellites Sentinel.

Ces codes permettent de récupérer les données GLS (Global Land Service) du capteur optique OCLI de Sentinel-3 et les données radar GRD de Sentinel-1.
Le code permet de calculer le Radar végétation Index (RVI).


Pour l’instant ce sont plutôt des brouillons, mais ils peuvent servir pour accélérer le travail de récupération et l’automatisation des traitements.

## Travail à effectuer (issues)

- [x] Lisser l'indicateur RVI.
- [x] Construire une validation dépendant de la qualité des extractions.
- [x] Récupérer des données de température.
- [x] Exploiter les saisonnalités.
- [ ] Utiliser geojson.io pour récupérer les coordonnées géographiques d'une zone.
- [ ] Lisser progressivement les images GRD plutôt que de lisser l'indicateur.
- [ ] Extraire des images d'une zone au lieu d'un rectangle et d'un masque.
- [ ] Cumuler les extraction de zones pour pouvoir créer une image de zone incomplète à partir de plusieurs images combinées.

## Utilisation
Pour analyser une zone donnée, il faut créer un fichier geozone.json qui contient des références géographiques et les répertoires de stockage et de travail.

L’utilisation d’un répertoire de travail différent du répertoire de stockage permet d’utiliser un disque plus rapide pour les données temporaires.

*	Le nom va de soi. Il sert uniquement aux affichages.
*	La « bbox » correspond aux bornes en latitudes et longitude : (lat0, lon0, lat1, lon1), le point 0 étant le coin supérieur gauche, de latitude et longitude minimales, et le point 1, le coins inférieur droit.
*	« step » est un facteur de décimation de la fenêtre pour le calcul de l’imagette, le rectangle rouge de l’image (Figure 11). On déborde forcément de la zone d’intérêt, mais d’au plus « step » pixels lors du chargement en mémoire.
*	« GRD » est le répertoire dans lequel se trouvent, ou doivent se trouver les archives .zip GRD de Sentinel-1.

Les autres répertoires contiennet des données de Copernicus Global Land Services (GLS). On peut entrer les répertoire dans le jichier JSON, mais sinon ils seront retrouvés par défaut dans le même dossier que celui où se trouve GRD.

*	« FCOVER » donne l’adresse du répertoire dans lequel se situent les fichiers netCDF téléchargés depuis GLS pour l’indicateur FCOVER. Même chose pour les indicateurs « LAI » ou « NDVI » s’ils sont téléchargés.
*   « TCI » est l'indicateur de température (Thermal Condition Index). Pour l'instant ces données sont d'asez mauvaise qualité. Il va falloir truver d'autres sources.

### Récupération des données RADAR de Sentinel-1
Pour récupérer des données GRD du serveur « SentielHub » il faut juste exécuter le code

    python grdloadzone.py <adresse du dossier contenant les zones>

Chaque zone est un sous-dossier qui contien son sous-répertoir JSON.

Par exemple

    python grdloadzone.py data
 
Le chargement s’effectuera automatiquement zone après zone. Il ne s'arrêtera pas. Il faudra un ^C.
Les données sont recherchées des plus récentes aux plus anciennes. On retourne sur les données anciennes car on espère progressivement les faires passer de zone froide en zone chaude sur le serveur.
L'accès à Sentinel Hub nécessite un login et un mot de passe qui peuvent être stockés dans le fichier ~/.netrc :

    machine apihub.copernicus.eu login ********* password **********


## Des exemples
Les notebooks explicatifs sont dans le dossier « notebooks ».

* Chaque fichier à son notebook associé : par exemple dataset_doc.ipynb donne des explications du fonctionnement du code du module dataset.

Le fichier « TraceRVI.ipynb » présente des calculs ayant abouti aux graphes présentés dans le document #7. Des données CSV sont disponibles dans le répertoire « data/target » pour chaque zone, ce qui évite d’avoir à relancer les traitements depuis les données brutes.

Le fichier « TrendAnalysis.ipynb » contient des exemples de calculs permettant de prendre en compte les variations saisonnières pour interpréter l'évolution de la végétation. On rajoute aussi une information de température pour ne garder en résidus que ce qui ne dépend plus ni de la saison, ni de la température.

Pour l'instant les calculs sont des stats linéaires élémentaires permettant d'obtenir un résultat rapide et ouvant certaines perspectives.

La saisonnalité est étudiée sur les données optiques GLS et pourra être prise en compte ensuite avec les donnés radar.