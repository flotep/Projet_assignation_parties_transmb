import os
import sys
import re as reg
import pandas as pd
import numpy as np
import biotite.structure as struc
import biotite.structure.io as strucio
from Bio import PDB
from Bio.PDB.DSSP import *
from Bio.PDB import *
from sympy import *

script_dir = os.path.dirname( __file__ )
mymodule_dir = os.path.join( script_dir, '..', 'src' )
sys.path.append( mymodule_dir )

from membrane_function import *
from math_function import *
from random_geometry_points.plane import Plane


p = PDBParser()
structure = p.get_structure("2n90", "2n90.pdb")
model_protein = structure[0]
dssp = DSSP(model_protein, "2n90.pdb", dssp='mkdssp')


dssp_tuple=dssp_dict_from_pdb_file("2n90.pdb") 
dssp_dict = dssp_tuple[1]

id_tuple = [x[1] for x in dssp_dict] 
id=[y[1] for y in id_tuple]

df_DSSP = pd.DataFrame(dssp) 
column=[0,1,3]
df_DSSP = df_DSSP.iloc[:,column] 
df_DSSP.columns=['#','amino_acid','Relative ASA']
df_DSSP['#']=id  


threshold = 0.20 #Seuil minimum de la relative ASA pour considérer qu'un résidu est accessible au solvant.
residue_hydrophobic=['F','G','I','L','M','V','W','A','C','S','H'] #Classification binaire des résidus selon Postic et al., 2015
residue_hydrophilic=['D','E','K','N','P','Q','R']


#Conserve seulement les CA si leur résidu est accessible au solvant et classe les hydrophobes et les hydrophiles dans deux listes séparées
accessible_ca = df_DSSP.loc[(df_DSSP['Relative ASA']>threshold)] 
accessible_ca_hydrophobic = accessible_ca.loc[accessible_ca['amino_acid'].str.contains('|'.join(residue_hydrophobic))]
accessible_ca_hydrophilic = accessible_ca.loc[accessible_ca['amino_acid'].str.contains('|'.join(residue_hydrophilic))] 

#Récupère les résidus IDs des CA accessibles au solvant pour les résidus hydrophobes et hydrophiles.
accessible_residue_hydrophobic= accessible_ca_hydrophobic['#'].tolist() 
accessible_residue_hydrophilic= accessible_ca_hydrophilic['#'].tolist() 



#Lecture du fichier PDB en entrée pour récupérer les coordonnées de tous les CA accessibles au solvant.
#Ces coordonnées ne sont pas accessible depuis le fichier DSSP ou le dictionnaire DSSP.
io = PDB.PDBIO
struct = p.get_structure("2n90", "2n90.pdb") 
model = structure[0] #Ne parcours que le modèle de la protéine souhaité, par défaut le 1er modèle.

#Liste les coordonnées x,y,z des CA des résidus hydrophobes et hydrophiles.
ca_hydrophobic_coord_list = [] 
ca_hydrophilic_coord_list = []

for chain in model:
    for residue in chain: 
        res_id=residue.get_full_id()[3][1] #On ne récupère que le ID du résidu parmi tous les IDs obtenus avec le get_full_id
        if res_id in accessible_residue_hydrophobic :
            ca=residue["CA"]
            ca_hydrophobic_coord=ca.get_coord() #Récupère uniquement les coordonnées des CA.
            ca_hydrophobic_coord_list.append(ca_hydrophobic_coord) 
            
        elif res_id in accessible_residue_hydrophilic : 
            ca=residue["CA"]
            ca_hydrophilic_coord=ca.get_coord() 
            ca_hydrophilic_coord_list.append(ca_hydrophilic_coord) 

#Calcul le centre de masse de la protéine avec les coordonnées de tous les CA des résidus accessibles au solvant.
ca_coord_list = np.array(ca_hydrophilic_coord_list + ca_hydrophobic_coord_list)
com=center_of_mass(ca_coord_list)
   

#Transforme les coordonnées de tous les CA par rapport aux coordonnées du centre de masse (qui représentera le point (0,0,0))
transformed_coord_hydrophilic=transform_coordinates(ca_hydrophilic_coord_list,com)
transformed_coord_hydrophobic=transform_coordinates(ca_hydrophobic_coord_list,com)

#Détermine les vecteurs (a,b,c) à partir de l'échantillonnage de points sur sphère.
vecteurs=sample_spherical(30,3)
#Création d'un plan orthogonal au vecteur passant par le point (0,0,0).
orthogonal_planes=get_orthogonal_planes(vecteurs,x,y,z,d)


#Création des listes qui vont stocker la position et le score maximum obtenus dans chaque direction.
optimal_position = []
optimal_score_direction =[]

#Boucle "for" pour trouver le meilleur score obtenu dans chaque direction, et à quelle position il a été obtenu.
for i in orthogonal_planes : 
    list_score=[]
    position=0 #The initial position is always 0

    #Donne les équations des deux plans de la membrane à la position 0.
    membrane1, membrane2 = get_membrane_planes(i) 


    #Calcul le score obtenu avec la membrane à la position 0.
    score, nb_hydrophobic_atom_inside=objective_function(membrane1, membrane2, transformed_coord_hydrophilic, transformed_coord_hydrophobic)
    list_score.append(score)
    #print (atoms_inside)

    #Déplace la membrane dans la direction du vecteur par tranche de 2 A tant que le nombre de CA hydrophobes entre les deux plans > 0.
    while nb_hydrophobic_atom_inside > 0:
        position += 2

        #Donne les équations des plans de la membrane à la position donnée.
        membrane1, membrane2 = get_membrane_planes(i,pos=position)

        #Calcul le score obtenu avec cette membrane à la position donnée.
        score, nb_hydrophobic_atom_inside=objective_function(membrane1, membrane2, transformed_coord_hydrophilic, transformed_coord_hydrophobic)
        list_score.append(score)

    #Trouve le score maximum obtenu pour une membrane dans cette direction, et à quel index de la liste il se trouve.
    max_score_direction=max(list_score)
    idx_max_score = list_score.index(max_score_direction) 

    
    #Récupère le score maximal et l'index associé pour chaque direction.
    optimal_position.append(idx_max_score)
    optimal_score_direction.append(max_score_direction)

#print (optimal_position)
 
#Liste les équations des plans optimales après épaississement dans chaque direction, et leur score associé.
optimal_membrane1_list=[]
optimal_membrane2_list=[]
optimal_score=[]

#Boucle "for" pour trouver le meilleur score obtenu après épaississement à la position optimale de chaque direction.
for plan,idx_position,score in zip(orthogonal_planes,optimal_position, optimal_score_direction)  : #for each membrane with the highest number of atoms, we increase their thickness 1 A by 1 A as long as the score increase
    add_to_pos=2 #correspond to the size of the pas during the sliding of the membrane
    thickness = 14
    max_score=score 

    list_mb1=[] 
    list_mb2=[]
    list_score=[]
    

    thickness +=1 #Incrémentation de l'épaisseur de la membrane de 1 A par 1 A.
    position = 0 + add_to_pos*idx_position #Retrouve la position à laquelle la meilleur membrane a été trouvée dans cette direction.
    

    #Calcul les nouvelles équations et le score des plans de la meilleur position après augmentation de son épaisseur.
    membrane1, membrane2 = get_membrane_planes(plan, thickness= thickness,pos= position)
    score,nb_hydrophobic_atom_inside=objective_function(membrane1, membrane2, transformed_coord_hydrophilic, transformed_coord_hydrophobic)

    list_mb1.append(membrane1)
    list_mb2.append(membrane2)
    list_score.append(score)


    #Si le score de la membrane après épaississement est supérieur au score précédent, continue d'épaissir la membrane tant que c'est le cas.
    if score>=max_score: 
        while score >= max_score : 
            max_score=score
            thickness += 1
            membrane1, membrane2 = get_membrane_planes(plan, thickness=thickness,pos= position)
            new_score,nb_hydrophobic_atom_inside=objective_function(membrane1, membrane2, transformed_coord_hydrophilic, transformed_coord_hydrophobic)
        
            list_mb1.append(membrane1)
            list_mb2.append(membrane2)
            list_score.append(new_score)

            if new_score >= score : 
                max_score = score
                score = new_score
            
            #Si le nouveau score est inférieur au précédent, récupère les équations et le score de la membrane précédente. 
            else : 
                membrane1, membrane2 = get_membrane_planes(plan, thickness=thickness,pos= position) #if the score doesnt increase for the new membrane, get the equation of the previous membrane + its score 
                optimal_membrane1_list.append(list_mb1[-2])
                optimal_membrane2_list.append(list_mb2[-2])
                optimal_score.append(list_score[-2])
                break

    #Si le premier épaissisement n'a pas amélioré le score, garde les équations et le score obtenu sans épaississement.
    else : 
        membrane1, membrane2 = get_membrane_planes(plan,pos= position) #keep the original thickness if the incremation doesn't increase the score from the beginning
        score, nb_hydrophobic_atom_inside=objective_function(membrane1, membrane2, transformed_coord_hydrophilic, transformed_coord_hydrophobic)
        optimal_membrane1_list.append(membrane1)
        optimal_membrane2_list.append(membrane2)
        optimal_score.append(score)
  

#Trouve le score maximal et l'index associé obtenus après épaississement de toutes les membranes à la meilleure position de chaque direction. 
max_score_membrane=max(optimal_score) 
idx_score_optimal = optimal_score.index(max_score_membrane) 


#Retourne les équations des plans de la membrane optimale pour la protéine. 
optimal_membrane1=optimal_membrane1_list[idx_score_optimal] 
optimal_membrane2=optimal_membrane2_list[idx_score_optimal]
print (optimal_membrane1, optimal_membrane2)


#Récupération de chaque argument au sein des équations.
#Doit utiliser findall car les équations sont composées de symbole sympy, l'ordre des arguments retournés par .arg[] change entre les équations.
str_mb1_equation=str(optimal_membrane1)
str_mb2_equation=str(optimal_membrane2)
str_mb1_parameters= reg.findall(r"(?:\+\s*|\-\s*)?[-+]?\d*\.*\d+", str_mb1_equation)
str_mb2_parameters= reg.findall(r"(?:\+\s*|\-\s*)?[-+]?\d*\.*\d+", str_mb2_equation)
mb1_parameters=[float(str(i).replace(' ','')) for i in str_mb1_parameters]
mb2_parameters=[float(str(i).replace(' ','')) for i in str_mb2_parameters]


a,b,c=mb1_parameters[0],mb1_parameters[1],mb1_parameters[2]
d1=mb1_parameters[3]
d2=mb2_parameters[3]

normal_vec=(a,b,c)

#Génération de points (x,y,z) random situés dans le plan de la 1ere couche de la membrane.
r=d1/c
refpoint=(0.0,0.0,r)
plan1=Plane(normal_vec,d1,refpoint,30)
random_plane_points = plan1.create_random_points(1000) #https://pypi.org/project/random-geometry-points/

#Génération de points random (x,y,z) situés dans le plan de la 2ere couche de la membrane.
r2=d2/c
refpoint2=(0.0,0.0,r2)
plan2=Plane(normal_vec,d2,refpoint2,30)
random_plane_points_2 = plan2.create_random_points(1000) 

#Lecture du fichier PDB en entrée puis transformation des coordonnées de chaque atome par rapport au centre de masse.
parser = PDB.PDBParser()
io = PDB.PDBIO()
struct_transform = parser.get_structure("2n90", "2n90.pdb")
model_transform =struct_transform[0]
rotation_matrix = PDB.rotmat(PDB.Vector([0, 0, 0]), PDB.Vector([0, 0, 0])) #Matrice neutre, pas de rotation des coordonnées souhaitée.
transform=(-com[0],-com[1], -com[2])  

for chain in model_transform:
    for residue in chain:
        for atom in residue:
            atom.transform(rotation_matrix, transform)

#Retourne un fichier PDB avec les coordonnées de tous les atomes transformées.
io.set_structure(model_transform)
io.save("2n90_transform.pdb")


#Lis le fichier contenant les nouvelles coordonnées avec le module Biotite.
atom_array = strucio.load_structure("2n90_transform.pdb")

#Pour chaque point généré aléatoirement sur la 1ere couche de la membrane , créé un hétéroatome possédant ses coordonnées.
for i in random_plane_points : 
    atom = struc.Atom(
        coord = [i[0],i[1],i[2]],
        chain_id = "N",
        res_id = atom_array.res_id[-1] + 1, #Le résidu id de l'atome est égale au dernier résidu ID du fichier +1
        res_name = "MB1", #Chaque atome créée pour la 1ère couche possède le nom de résidu MB1.
        hetero = True,
        atom_name = "CA",
        element = "C"
    )
    atom_array += struc.array([atom]) #Ajout de l'hétéroatome à la suite des atomes présent dans le fichier PDB.

#Pour chaque point généré aléatoirement sur la 2eme couche de la membrane, créé un hétéroatome possédant ses coordonnées.
for j in random_plane_points_2 : 
    atom = struc.Atom(
        coord = [j[0],j[1],j[2]],
        chain_id = "N",
        res_id = atom_array.res_id[-1] + 1,
        res_name = "MB2", #Chaque atome créée pour la 1ère couche possède le nom de résidu MB2.
        hetero = True,
        atom_name = "CA",
        element = "C"
    )
    atom_array += struc.array([atom]) #Ajout de l'hétéroatome à la suite des atomes présent dans le fichier PDB.

#Retourne un fichier PDB édité avec les hétéroatomes représentant chaque couche de la membrane ajoutés à la structure de la protéine.
strucio.save_structure("2n90_edited.pdb", atom_array)


