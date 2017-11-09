#!/usr/bin/python
# Author : Alexandra Benamar, Robin Droit

import os                      # gestion de fichiers et de dossiers
import sys                     # gestion des erreurs et des arguments
import string
from math import *
import numpy as np
import matplotlib.pyplot as plt
from pylab import *
from fonctions import *

try:
    # on veut que le fichier .pdb soit en argument juste apres -pdb
    pdb_file = sys.argv[1]
    # affichage
    print "Fichier pdb a traiter : ", pdb_file

except:
    usage()
    print "ERREUR: veuillez entrer le nom du fichier .pdb"
    sys.exit()


dico = parserPDB(pdb_file)
centerMassOfConf(dico)
centerMassOfRes(dico)
distance(dico)
rayon_giration(dico)
rmsd_local(dico)
dist_moy = distance_moyenne(dico)
rmsd_moy = rmsd_moyen(dico)
write_PDB(dico)
write_Gira(dico)
write_DMoy(dist_moy)
write_rmsd_moy(rmsd_moy)


dico2 = parserPDB_CA(pdb_file)
rmsd_global(dico2)
rmsd_local(dico2)
write_PDB_CA(dico2)
write_RMSD(dico2)

#------------------------------------------------------------------
# Graphiques

# LISTE : RMSD DE CHAQUE CONFORMATION
list_RMSD = []
for model in dico2.keys():
    list_RMSD.append(dico2[model]["rmsd"])


# LISTE : RAYON DE GIRATION DE CHAQUE CONFORMATION
list_Gira = []
for model in dico.keys():
    list_Gira.append(dico[model]["giration"])


# GRAPHIQUES DU RMSD GLOBAL ET DU RAYON GIRATOIRE DE CHAQUE CONFORMATION
fig, ax = plt.subplots(2, sharex=True)
fig.subplots_adjust(hspace=0.3)
ax[0].plot(range(1,len(list_RMSD)),list_RMSD[1:])
ax[0].set_ylabel("RMSD (en Angstrom)")
ax[1].plot(range(1,len(list_RMSD)),list_Gira[1:],range(1,len(list_RMSD)),[list_Gira[0] for x in range(1,len(list_RMSD))],'r-')
ax[1].set_ylabel("Rayon giratoire (en Angstrom)"); ax[1].set_xlabel("Modeles")
fig.suptitle('Courbe RMSD', fontsize=12)
fig.text(.5,.5,'Courbe rayon giratoire',fontsize=12,ha='center')
#plt.show()

# LISTE : RMSD MOYENS POUR CHAQUE AA
rmsd_list = []
reslist = rmsd_moy["reslist"]
for res in reslist :
    rmsd_list.append(rmsd_moy[res])

# LISTE : DISTANCES MOYENNES AU CENTRE DE MASSE
moy_list = []
reslist = dist_moy["reslist"]
for res in reslist :
    moy_list.append(dist_moy[res])


# GRAPHE DONNANT LES RMSD MOYENS POUR CHAQUE ACIDE AMINE
# ET LES DISTANCES MOYENNES AU CENTRE DE MASSE
fig, ax = plt.subplots(2, sharex=True)
fig.subplots_adjust(hspace=0.3)
ax[0].plot(range(0,len(rmsd_list)),rmsd_list)
ax[0].set_ylabel("RMSD (en Angstrom)")
ax[1].plot(range(0,len(rmsd_list)),moy_list,'r-')
ax[1].set_ylabel("Distance (en Angstrom)"); ax[1].set_xlabel("AA")
fig.suptitle('RMSD moyen', fontsize=12)
fig.text(.5,.5,'Enfouissement',fontsize=12,ha='center')
#plt.show()

#recuperation des donnes de rmsd des residus 39 et 76

def recup_39(dPDB):
    vect39 = []
    for i in sorted(dPDB.keys()):

        for y in dPDB[i]["listChains"]:
            if y == "39 ":
                vect39.append(dPDB[i][y]["rmsd"])

    return vect39

def recup_76(dPDB):

    vect76 = []
    for i in sorted(dPDB.keys()):

        for y in dPDB[i]["listChains"]:
            if y == "76 ":
                vect76.append(dPDB[i][y]["rmsd"])

    return vect76

#------------------------------------------------------------------
# Ecriture des resultats dans des fichiers de sortie

def write_PDB(dPDB, filout="output_pdb.pdb"):
 fig.suptitle('RMSD moyen', fontsize=12)
 fig.text(.5,.5,'Enfouissement',fontsize=12,ha='center')
 #plt.show()

#GRAPHE DONNANT LES RMSD DES RESIDUS 39 ET 76
vect39 = recup_39(dico)
vect76 = recup_76(dico)

fig, ax = plt.subplots(2, sharex=True)
fig.subplots_adjust(hspace=0.3)
ax[0].plot(range(0,len(vect39)),vect39)
ax[0].set_ylabel("RMSD 39 (en Angstrom)")
ax[1].plot(range(0,len(vect76)),vect76)
ax[1].set_ylabel("RMSD 76 (en Angstrom)")
fig.suptitle('Variation du RMSD au cours du temps - 39 ASP', fontsize=12)
fig.text(.5,.5,'Variation du RMSD au cours du temps - 76 GLU',fontsize=12,ha='center')
plt.show()
