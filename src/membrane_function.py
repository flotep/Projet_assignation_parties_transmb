#!/usr/bin/env python

"""
Fonctions qui permettent de définir l'équation des plans des membranes et de calculer leur score.
La formule mathématique utilisée pour le scoring est issue de l'algorithme ANVIL.
"""
__author__ = "Florian TEP"
__email__ = "tepflorian@gmail.com"

import math 
from sympy import *

x = Symbol("x") # Définit un symbole sympy pour chaque paramètre inconnu de l'équation d'un plan.
y = Symbol("y")
z = Symbol("z")
d = Symbol("d")

def get_orthogonal_planes(vector_list,x,y,z,d): 
    """
    Permet de calculer un plan orthogonal à un vecteur passant par le point (0,0,0).
    Parametres 
    ----------
    vector_list : liste
        Une liste de vecteur sous la forme [(a1,b1,c1),(a2,b2,c2),...]

    x, y, z, d : symboles

    Return 
    ------
    liste
        Une liste avec un plan orthogonal calculé pour chaque vecteur 
        sous la forme [((a1*x)+(b1*y)+(c1*z)+d), (a2*x)+(b2*y)+(c2*z)+d]
    """
    all_orthogonal_planes=[]
    for i in vector_list:
        a=i[0]
        b=i[1]
        c=i[2]
        orthogonal_plane=((a*x) + (b*y) + (c*z) + d) 
        all_orthogonal_planes.append(orthogonal_plane)
    return all_orthogonal_planes

def get_membrane_planes(plane, thickness =14, pos=0): # get the 2 initial equation plane of the membrane for each vector, can choose the thickness
    """
    Permet d'obtenir les équations des deux plans d'une membrane à partir d'un plan orthogonal, 
    de l'épaisseur de la membrane voulue, et d'une position précise par rapport au vecteur (0 étant la position d'origine).
    Parametres 
    ----------
    plane : equation sympy
        L'équation d'un plan.

    thickness : int
        L'épaisseur de la membrane souhaitée.
    
    pos : int
        La position voulue de la membrane par rapport au vecteur normal.

    Return 
    ------
    mb1, mb2 : equation sympy
        Les équations des plans des deux couches de la membrane avec une épaisseur et une position donnée.
    """
    mb1= plane.subs(d,(thickness/2)+pos)
    mb2= plane.subs(d,(-(thickness/2))+pos)
    return mb1, mb2 

def objective_function (mb1, mb2, hydrophilic_atom_coordinates, hydrophobic_atom_coordinates): #get the number of atom inside each membrane
    """
    Calcul le score d'une membrane selon la formule mathématique contenue dans la fonction.
    Parametres 
    ----------
    mb1, mb2 : equation sympy
        Les équations des deux plans de la membrane

    hydrophilic_atom_coordinates, hydrophobic_atom_coordinates : listes
        Listes des coordonnées des atomes hydrophiles et hydrophobes accessibles au solvant.

    Return 
    ------
    score : int
        Le score de la membrane étudié.
    mi : int 
        Le nombre d'atomes issus de résidus hydrophobes situés dans la membrane.
    """
    mi=0 #atomes issus de résidus hydrophobes à l'intérieur de la membrane 
    me=0 #atomes issus de résidus hydrophobes à l'extérieur de la membrane 
    si=0 #atomes issus de résidus hydrophiles à l'intérieur de la membrane 
    se=0 #atomes issus de résidus hydrophiles à l'extérieur de la membrane 

    for i in hydrophobic_atom_coordinates:
        result_upper = mb1.subs([(x, i[0]),(y, i[1]), (z,i[2])]) #change les valeurs de x, y et z de l'équation par les valeurs des coordonées x,y,z  de chaque atome (coord transformée par rapport au centre de masse)
        result_under = mb2.subs([(x, i[0]),(y, i[1]), (z,i[2])]) 
        if (result_upper<0 and result_under>0) or (result_upper>0 and result_under<0): #see if the atom is in between the 2 planes
            mi += 1
        else : 
            me +=1

    for j in hydrophilic_atom_coordinates:
        result_upper = mb1.subs([(x, j[0]),(y, j[1]), (z,j[2])]) #change les valeurs de x, y et z de l'équation par les valeurs des coordonées x,y,z  de chaque atome (coord transformée par rapport au centre de masse)
        result_under = mb2.subs([(x, j[0]),(y, j[1]), (z,j[2])])
        if (result_upper<0 and result_under>0) or (result_upper>0 and result_under<0): #see if the atom is in between the 2 planes
            si += 1
        else :
            se +=1
        
    try :
        score = (mi * se - si * me) / (math.sqrt((mi + si) * (mi + me) * (si + se) * (se + me)))

    except ZeroDivisionError:
        score = 0

    return score, mi 
