#!/usr/bin/python
# Author : Alexandra Benamar, Robin Droit

import os                      # gestion de fichiers et de dossiers
import sys                     # gestion des erreurs et des arguments
import string                  # gestion des chaines de caractere
from math import *             # necessaire a l'utilisation de sqrt
import numpy as np

#------------------------------------------------------------------
# Parseur general

def parserPDB(infile):
    """ Cette fonction a pour but de parser un fichier de type pdb afin
    d'en recuperer les informations sur les differents atomes qui composent
    les acides amines.

    Input : fichier pdb a parser.
    Output : dictionnaire contenant les informations sur le fichier.

    """

    # Test d'ouverture du fichier
    try:
		f = open(infile,'r')
    except:
		print "Erreur, le fichier n'a pas pu s'ouvrir"
		sys.exit(1)

    lines = f.readlines()
    number_of_lines=len(lines)


    # Test fichier vide
    if number_of_lines==0:
		print "Erreur, le fichier est vide"
		sys.exit(1)

    # Initialisation d'un dictionnaire
    dicPDB = {}
    # Initialisation d'un flag pour remplir la liste des atomes
    Flag = False

    for line in lines:

		if line[0:5] == "MODEL":
			modelnumber = int(string.strip(line[10:14]))
			dicPDB[modelnumber] = {}
			modelnumber = int(modelnumber)
			dicPDB[modelnumber]["listChains"] = []  # numero du residu
			dicPDB[modelnumber]["listRes"] = []     # noms des residus
			# le dictionnaire a la cle "listChains" qui prend une liste

												   # Pour toutes les lignes qui commencent par ATOM (celles qui ont des atomes)
		elif line[0:4] == "ATOM":
			chain = line[24:27]
												   # on ne selectionne que les lignes qui contiennent des ATOM
			if chain not in dicPDB[modelnumber]["listChains"]:
				dicPDB[modelnumber]["listChains"].append(chain)
				dicPDB[modelnumber]["listRes"].append(string.strip(line[17:20]))
				dicPDB[modelnumber][chain] = {}
                                                    # pour la cle number ayant pour cle "resname"
			if dicPDB[modelnumber][chain].has_key("resname") == False:
				dicPDB[modelnumber][chain]["resname"] = string.strip(line[17:20])
				dicPDB[modelnumber][chain]["atomlist"] = []  # a pour cle atomlist et prend une liste

			atomtype = string.strip(line[12:16])
			#print atomtype

			dicPDB[modelnumber][chain]["atomlist"].append(atomtype) # ajout de l'atome a la liste

			dicPDB[modelnumber][chain][atomtype] = {}    # cree un dictionnaire dans dicPBD[chain][number]

			dicPDB[modelnumber][chain][atomtype]["x"] = float(line[30:38])
			dicPDB[modelnumber][chain][atomtype]["y"] = float(line[38:46])
			dicPDB[modelnumber][chain][atomtype]["z"] = float(line[46:54])
			dicPDB[modelnumber][chain][atomtype]["id"] = line[7:11].strip()
			dicPDB[modelnumber][chain][atomtype]["lyst"] = [float(line[30:38]),float(line[38:46]),float(line[46:54]),line[6:11].strip()]



	# Test presence d'ATOM
    if dicPDB==0:
		print "Le fichier ne contient pas d'ATOM"
		sys.exit(1)

    # Fermeture du fichier
    f.close()
    return dicPDB


#------------------------------------------------------------------
# Parseur avec carbones alpha uniquement

def parserPDB_CA(infile):
    """ Cette fonction a pour but de parser un fichier de type pdb afin
    d'en recuperer les informations sur les carbones alpha.

    Input : fichier pdb a parser.
    Output : dictionnaire contenant les informations sur les carbones alpha du fichier.

    """

    # Test d'ouverture du fichier
    try:
		f = open(infile,'r')
    except:
		print "Erreur, le fichier n'a pas pu s'ouvrir"
		sys.exit(1)

    lines  = f.readlines()
    number_of_lines=len(lines)


    # Test fichier vide
    if number_of_lines==0:
		print "Erreur, le fichier est vide"
		sys.exit(1)

    # Initialisation d'un dictionnaire
    dicPDB = {}
    # Initialisation d'un flag pour remplir la liste des atomes
    Flag = False

    for line in lines:

		if line[0:5] == "MODEL":
			modelnumber = int(string.strip(line[10:14]))
			dicPDB[modelnumber] = {}
			modelnumber = int(modelnumber)
			dicPDB[modelnumber]["listChains"] = []  # numero du residu
			dicPDB[modelnumber]["listRes"] = []     # noms des residus
			# le dictionnaire a la cle "listChains" qui prend une liste

												   # Pour toutes les lignes qui commencent par ATOM (celles qui ont des atomes)
		elif line[0:4] == "ATOM" and line[13:15] == "CA": # pour tous les atomes etant des CA
			chain = line[24:27]

												   # on ne selectionne que les lignes qui contiennent des ATOM
			if chain not in dicPDB[modelnumber]["listChains"]:
				dicPDB[modelnumber]["listChains"].append(chain)
				dicPDB[modelnumber]["listRes"].append(string.strip(line[17:20]))
				dicPDB[modelnumber][chain] = {}
                                                    # pour la cle number ayant pour cle "resname"
			if dicPDB[modelnumber][chain].has_key("resname") == False:
				dicPDB[modelnumber][chain]["resname"] = string.strip(line[17:20])

				dicPDB[modelnumber][chain]["atomlist"] = []  # a pour cle atomlist et prend une liste

			atomtype = string.strip(line[13:16])

			dicPDB[modelnumber][chain]["atomlist"].append(atomtype) # ajout de l'atome a la liste
			dicPDB[modelnumber][chain][atomtype] = {}    # cree un dictionnaire dans dicPBD[chain][number]

			dicPDB[modelnumber][chain][atomtype]["x"] = float(line[30:38])
			dicPDB[modelnumber][chain][atomtype]["y"] = float(line[38:46])
			dicPDB[modelnumber][chain][atomtype]["z"] = float(line[46:54])
			dicPDB[modelnumber][chain][atomtype]["id"] = line[6:11].strip()
			dicPDB[modelnumber][chain][atomtype]["lyst"] = [float(line[30:38]),float(line[38:46]),float(line[46:54]),line[6:11].strip()]

	# Test presence d'ATOM
    if dicPDB==0:
		print "Le fichier ne contient pas d'ATOM"
		sys.exit(1)

    # Fermeture du fichier
    f.close()
    return dicPDB


#------------------------------------------------------------------
# RMSD global et RMSD local

def rmsd_global(dCA):
    """ Calcul du RMSD global de chaque modele
        Input : dictionnaire des carbones alpha
        Output : RMSD entre le modele pris en compte et le modele de reference
    """

    for model in dCA.keys():      # numeros de conformation
        n=0
        somme=0
        chainlist = dCA[model]["listChains"]
        for chain in chainlist:    # numeros des residus
            atomlist = dCA[model][chain]["atomlist"]
            for atom in atomlist: # numeros des atomes
                distance=(dCA[model][chain][atom]["x"]-dCA[0][chain][atom]["x"])**2 + (dCA[model][chain][atom]["y"]-dCA[0][chain][atom]["y"])**2 + (dCA[model][chain][atom]["z"]-dCA[0][chain][atom]["z"])**2
                somme += distance
                n+=1
        rmsd = sqrt(somme/n)
        dCA[model]["rmsd"] = rmsd  #ajout d'une liste des rmsd de chaque conformation


def rmsd_local(dPDB): #on donne soit un dico specifique soit le dico general et les numeros des modeles et on renvoi une liste de tous les rmsd de tous les residus
    """
        Permet de calculer l'ensemble des distances entre les differents carbones alpha
        entre deux modeles.
        Input : dictionnaire (specifique ou general).
        Output : liste du RMSD de chaque residu.
    """
    for model in dPDB.keys():      # numeros de conformation
        chainlist = dPDB[model]["listChains"]
        for chain in chainlist:    # numeros des residus
            n = 0
            somme = 0
            atomlist = dPDB[model][chain]["atomlist"]
            for atom in atomlist:  # numeros des atomes
                distance=(dPDB[model][chain][atom]["x"]-dPDB[0][chain][atom]["x"])**2 + (dPDB[model][chain][atom]["y"]-dPDB[0][chain][atom]["y"])**2 + (dPDB[model][chain][atom]["z"]-dPDB[0][chain][atom]["z"])**2
                somme += distance
                n += 1
            dPDB[model][chain]["rmsd"] = sqrt(somme/n)  # liste de rmsd pour chaque residu

#------------------------------------------------------------------
# Rayon de giration

def rayon_giration(dPDB) :
    """ Calcul du rayon de giration de chaque conformation du dictionnaire.
        Input : dictionnaire global contenant les coordonnees des centres
        de masses des conformations de chaque atome.
        Output : ajout dans le dictionnaire d'une liste de rayon de giration
        pour chaque conformation.
    """
    for model in range(0,len(dPDB)):    # numeros de conformation
        chainlist = dPDB[model]["listChains"]
        raymax=0
        for chain in chainlist:         # numeros des residus
            atomlist = dPDB[model][chain]["atomlist"]
            for atom in atomlist:       # numeros des atomes
                dist = math.sqrt((dPDB[model][chain][atom]["x"]-dPDB[model]["XCM"])**2
                +(dPDB[model][chain][atom]["y"]-dPDB[model]["YCM"])**2
                +(dPDB[model][chain][atom]["z"]-dPDB[model]["ZCM"])**2)
                if (dist>=raymax):
                    raymax = dist
        dPDB[model]["giration"] = raymax

#------------------------------------------------------------------
# Distance entre le centre de masse d'un residu et la conformation
def distance(dPDB) :
    """
        Calcul de la distance entre le centre de masse
        d'un residu et la conformation.
        Input : dictionnaire contenant les informations de centre de masse.
    """
    for model in range(0,len(dPDB)):        # numeros de conformation
        chainlist = dPDB[model]["listChains"]
        for chain in chainlist:             # numeros des residus
            distance = sqrt((dPDB[model][chain]["XCM"]-dPDB[model]["XCM"])**2+(dPDB[model][chain]["YCM"]-dPDB[model]["YCM"])**2+(dPDB[model][chain]["ZCM"]-dPDB[model]["ZCM"])**2)
            dPDB[model][chain]["dist_CM"] = distance

#------------------------------------------------------------------
# Centre de masse

def centerOfMassConf(dPDB):
    """ Calcul du centre de masse de chaque conformation.
        Tient compte de tous les atomes de chaque residus.

        Input : dictionnaire provenant du parseur global.
        Output : dictionnaire contenant le centre de masse
        de chaque conformation.
    """
    for model in range(0,len(dPDB)):
        x=y=z=cpt=0.0
        chainlist = dPDB[model]["listChains"]
        for chain in chainlist :
            atomlt = dPDB[model][chain]["atomlist"]
            for atom in atomlt:
                x += dPDB[model][chain][atom]['x']
                y += dPDB[model][chain][atom]['y']
                z += dPDB[model][chain][atom]['z']
                cpt += 1
        Xcm = float(x)/cpt
        Ycm = float(y)/cpt
        Zcm = float(z)/cpt
        dPDB[model]["XCM"] = Xcm
        dPDB[model]["YCM"] = Ycm
        dPDB[model]["ZCM"] = Zcm

def centerOfMassRes(dPDB):
    """ Calcul du centre de masse d'un residu en tenant compte
        de tous ses atomes.

        Input : dictionnaire global
        Output : dictionnaire contenant le centre de masse
        de chaque residus.
    """

    for model in range(0,len(dPDB)):
        chainlist = dPDB[model]["listChains"]
        for chain in chainlist :
            cpt=x=y=z=0.0
            atomlt = dPDB[model][chain]["atomlist"]
            for atom in atomlt:
                x += dPDB[model][chain][atom]['x']
                y += dPDB[model][chain][atom]['y']
                z += dPDB[model][chain][atom]['z']
                cpt += 1
            Xcm = float(x)/cpt
            Ycm = float(y)/cpt
            Zcm = float(z)/cpt
            dPDB[model][chain]["XCM"] = Xcm
            dPDB[model][chain]["YCM"] = Ycm
            dPDB[model][chain]["ZCM"] = Zcm

#------------------------------------------------------------------
# Calcul des distances moyennes

def mean_dist(dPDB):
    """ Calcul de la distance moyenne au cours du temps entre le
        centre de masse d'un residu et le centre de masse de la proteine.
        Input : dictionnaire global contenant les distances entre le centre
        de masse des residus et la conformation.
        Output : ajout dans le dictionnaire des distances moyennes.
    """
    mean_dist = {}
    mean_dist["reslist"] = []
    n=0

    for model in range(0,len(dPDB)):
        chainlist = dPDB[model]["listChains"]
        for chain in chainlist:
            if not chain in mean_dist["reslist"]:
                mean_dist["reslist"].append(chain)
                mean_dist[chain]=dPDB[model][chain]["dist_CM"]
            else:
                mean_dist[chain]+=dPDB[model][chain]["dist_CM"]
        n+=1
    reslist=mean_dist["reslist"]
    for res in reslist:
        mean_dist[res]=mean_dist[res]/n

    return mean_dist


def rmsd_moyen(dPDB):
    """
        Calcul du RMSD local moyen au cours du temps pour chaque acide amine.
        Input : dictionnaire global contenant les RMSD locaux pour chaque residu.
        Output : ajout dans le dictionnaire des RMSD moyens.
    """
    dRMSDMean={}
    dRMSDMean["reslist"]=[]
    n=0
    for model in range(0,len(dPDB)):
        chainlist = dPDB[model]["listChains"]
        for chain in chainlist:
            if not chain in dRMSDMean["reslist"]:
                dRMSDMean["reslist"].append(chain)
                dRMSDMean[chain]=dPDB[model][chain]["rmsd"]
            else:
                dRMSDMean[chain]+=dPDB[model][chain]["rmsd"]
        n+=1
    reslist=dRMSDMean["reslist"]
    for res in reslist:
        dRMSDMean[res]=dRMSDMean[res]/n

    return dRMSDMean

#------------------------------------------------------------------
# Ecriture des resultats dans des fichiers de sortie

def write_PDB(dPDB, filout="output_pdb.pdb"):
    """
        Ecriture des informations contenues dans le dictionnaire
        global dans un fichier de sortie au format pdb.
    """
    fout = open(filout, "w")

    for model in range(0,len(dPDB)):
        fout.write("CEMA  %4s                          %8.3f%8.3f%8.3f\n"
        %(model, dPDB[model]["XCM"],dPDB[model]["YCM"],dPDB[model]["ZCM"]))
        fout.write("RAGI  %4s                                                   %8.3f\n"
        %(model, dPDB[model]["giration"]))
        for chain in dPDB[model]["listChains"]:
            fout.write("CMRE   %4s                  %3s     %8.3f%8.3f%8.3f\n"
            %(model, chain, dPDB[model][chain]["XCM"], dPDB[model][chain]["YCM"], dPDB[model][chain]["ZCM"]))
            fout.write("D_CM   %4s                  %3s                     %8.3f\n"
            %(model, chain, dPDB[model][chain]["dist_CM"]))
            fout.write("RMSD   %4s                  %3s                     %8.3f\n"
            %(model, chain, dPDB[model][chain]["rmsd"]))
            for atom in dPDB[model][chain]["atomlist"]:
                fout.write("ATOM  %4s       %4s      %4s     %8.3f%8.3f%8.3f\n"
                %(model, chain, atom, dPDB[model][chain][atom]["x"], dPDB[model][chain][atom]["y"], dPDB[model][chain][atom]["z"]))
    fout.close()

def write_PDB_CA(dPDB, filout="output2_pdb.pdb"):
    """
        Ecriture des informations du dictionnaire de carbones alpha
        dans un fichier de sortie au format pdb.
    """
    fout = open(filout, "w")

    for model in range(0,len(dPDB)):
        for chain in dPDB[model]["listChains"]:
            fout.write("RMSD        %8.3f\n"%(dPDB[model][chain]["rmsd"]))
            for atom in dPDB[model][chain]["atomlist"]:
                fout.write("ATOM   %4s      %4s      %8.3f%8.3f%8.3f\n"%(model, chain, dPDB[model][chain][atom]["x"], dPDB[model][chain][atom]["y"], dPDB[model][chain][atom]["z"]))

    fout.close()


def write_Gira(dPDB, fileout = "giration.pdb"):
    """ Ecriture des rayons de giration dans un fichier de sortie
        au format pdb.
    """

    fout = open(fileout, "w")
    for model in range(0,len(dPDB)):
        fout.write("CONFO   %4s   RAYON_GIRATION   %8.3f\n"
        %(model, dPDB[model]["giration"]))
    fout.close()


def write_RMSD(dPDB, fileout = "rmsd.pdb"):
    """
        Ecriture des RMSD globaux dans un fichier de sortie au format pdb
    """

    fout = open(fileout, "w")
    for model in range(0,len(dPDB)):
        fout.write("CONFO   %4s   RMSD   %8.3f\n"
        %(model, dPDB[model]["rmsd"]))
    fout.close()


def write_DMoy(dist_moy, fileout="dist_moy.pdb"):
    """
        Ecriture des distances moyennes de chaque residu
        dans un fichier de sortie au format pdb.
    """

    fout=open(fileout,"w")
    reslist=dist_moy["reslist"]
    for res in reslist:
        fout.write("RES   %4s   DIST_MOY   %8.3f\n"
        %(res, dist_moy[res]))
    fout.close()


def write_rmsd_moy(dRMSD, fileout="rmsd_moy.pdb"):
    """ Ecriture des RMSD moyens de chaque residu
        dans un fichier de sortie au format pdb.
    """

    fout=open(fileout,"w")
    reslist=dRMSD["reslist"]
    for res in reslist:
        fout.write("RES   %4s   RMSD_MOY   %8.3f\n"
        %(res, dRMSD[res]))
    fout.close()


#------------------------------------------------------------------
# Informations sur le projet

def usage() :
    """
        Affichage des informations concernant le projet et son utilisation
    """
    print("\n\nProjet d'Alexandra BENAMAR et Robin DROIT\nCe programme permet de parser un fichier .pdb\net calcule les distances globales et locales de la proteine.\nUsage : ./projet_barstar.py -pdb <nom du fichier PDB>\n")

#------------------------------------------------------------------
