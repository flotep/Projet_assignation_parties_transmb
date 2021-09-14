# Assignation et détection des parties transmembranaires d'une protéine.

## projet_transmb.py

Détermine les zones transmembranaires d'une protéine à partir d'un fichier PDB.

- Usage:
    python projet_transmb.py infile outfile --model --vector_nb

- Example:
    python projet_transmb.py /home/ftep/Documents/Projet_court/data/2n90.pdb 
    /home/ftep/Documents/Projet_court/results --model=1 --vector_nb=30

Ce programme retourne : 

-Les deux équations des plans de la membrane optimale pour la protéine.

-Un fichier PDB avec les coordonnées transformées par rapport au centre de masse du fichier PDB en entrée ('...'_transform.pdb).

-Un fichier PDB édité avec des hétéroatomes ajoutés à la structure de la protéine pour représenter la membrane ('...'_edited.pdb). 

En ouvrant ce fichier PDB édité depuis PyMOL, les couches de la membrane peuvent être aisément mises en évidence en entrant les commandes suivantes:

```

select mb1, resn mb1

select mb2, resn mb2

```


Il est alors possible de changer la couleur des hétéroatomes séléctionnés pour mieux visualiser les deux plans.

Cet outil est inspiré de la méthode proposée par l'algorithme TMDET. 


## math_function.py
    
Fonctions qui permettent de réaliser des opérations mathématiques tels que 
le calcul du centre de masse, la transformation de coordonées par rapport à un centre,
ainsi que la génération aléatoire de points sur la surface d'une sphère.
    
def center_of_mass(list_coordinates): 
  
    Calcul le centre de masse en fonction des coordonnées des atomes contenus dans une liste.
    
    Parametres 
    ----------
    list_coordinates : liste
        Liste contenant les coordonées des atomes.
    
    Return
    ------
    center : liste
        Liste contenant les coordonées x,y et z du centre de masse calculé.
 
 
 def transform_coordinates(atom_coordinates,center_coordinates):
 
    Transforme des coordonnées atomiques par rapport au centre des coordonnées indiqués.
    
    Parametres 
    ----------
    atom_coordinates : list
        Liste contenant les coordonées des atomes [(x1,y1,z1),(x2,y2,z2)].
    center_coordinates : list
        Liste contenant les coordonnées du centre (x,y,z).
    
    Return
    ------
    atom_coordinates : list
        Listes contenant les coordonées transformées des atomes.
    

  def sample_spherical(npoints,ndim=3):
   
    Génère des vecteurs constitués d'échantillons indépendants provenant de trois distributions normales standards.
    
    Parametres 
    ----------
    npoints : int
        Nombre de points qui vont être générés aléatoirement
         sur la surface de la sphère pour créer les vecteurs.

    ndim : int
        nombres de distributions normales standards (dimensions)

    Return 
    ------
    vecteurs : liste
        Liste de vecteurs issus de l'échantillonnage de points sur la sphère.
  
  
  
  
  
## membrane_function.py

Fonctions qui permettent de définir l'équation des plans des membranes et de calculer leur score.
La formule mathématique utilisée pour le scoring est issue de l'algorithme ANVIL.

  
def get_orthogonal_planes(vector_list,x,y,z,d): 
   
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
   
  
def get_membrane_planes(plane, thickness =14, pos=0): 
   
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
  
    
  def objective_function (mb1, mb2, hydrophilic_atom_coordinates, hydrophobic_atom_coordinates): 
  
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
  
 

