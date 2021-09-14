#!/usr/bin/env python

"""
Fonctions qui permettent de réaliser des opérations mathématiques tels que 
le calcul du centre de masse, la transformation de coordonées par rapport à un centre,
ainsi que la génération aléatoire de points sur la surface d'une sphère.
"""
__author__ = "Florian TEP"
__email__ = "tepflorian@gmail.com"


import numpy as np

def center_of_mass(list_coordinates): #calculate the center of mass of the CA
    x,y,z=0,0,0
    num_atoms=list_coordinates.shape[0]
    for i in list_coordinates :
        x += i[0]
        y += i[1]
        z += i[2]
    center= [x / num_atoms, y / num_atoms, z / num_atoms]
    return center


def transform_coordinates(atom_coordinates,center_coordinates):
    for atom in atom_coordinates :
        atom[0]=atom[0]- center_coordinates[0]
        atom[1]=atom[1]- center_coordinates[1]
        atom[2]=atom[2]- center_coordinates[2]
    return atom_coordinates


def sample_spherical(npoints,ndim=3):
    """
    Génère un vecteur constitué d'échantillons indépendants provenant de trois distributions normales standards.
    """ 
    vecteurs = np.random.randn(ndim, npoints)
    vecteurs /= np.linalg.norm(vecteurs, axis=0)
    vecteurs = vecteurs.transpose()  
    return vecteurs



