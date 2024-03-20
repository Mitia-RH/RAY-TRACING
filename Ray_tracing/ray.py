from math import dist, sqrt
from csv import DictReader
from sys import argv, exit
import corrige_ray


NOM = "RANAIVO-HARISOA"
PRENOM = "Mitia"
EMAIL = "mitia.ranaivo-harisoa@u-psud.fr"

# Opérations sur les vecteurs
# Un vecteur est simplement représenté par un couple (x, y, z)

def add(v1, v2):
  """Addition de deux vecteurs"""
  return (v1[0]+v2[0], v1[1]+v2[1], v1[2]+v2[2])

def sub(v1, v2):
  """Différence entre deux vecteurs"""
  return (v1[0]-v2[0], v1[1]-v2[1], v1[2]-v2[2])

def mul(k, v):
  """Multiplication d'un vecteurs par une constante"""
  return (k*v[0], k*v[1], k*v[2])

def dot(v1, v2):
  """Produit scalaire de deux vecteurs"""
  return float((v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]))

def norm(v):
  """Norme d'un vecteur """
  return float(sqrt(v[0]**2 + v[1]**2 + v[2]**2))

def normalize(v):
  #note: nmz = normalize
  normv = norm(v)
  try:
    nmz = mul(1/(normv), v)
    return nmz
  except:
    raise ZeroDivisionError("La norme du vecteur est 0, impossible de le normaliser")

def test_add():
  """Fonction de test pour add"""
  assert(add((0,2.2,4), (1,2,3)) == (1,4.2,7))
  assert(add((0,0,0), (0,0,0)) == (0,0,0))
  assert(add((-2,1,-1), (1,2,3)) == (-1,3,2))


test_add()

def test_sub():
  """Fonction de test pour sub"""
  assert(sub((0,2.2,4), (1,2,3)) == (-1,2.2-2,1))
  assert(sub((0,0,0), (0,0,0)) == (0,0,0))
  assert(sub((-2,1,-1), (1,2,3)) == (-3,-1,-4))

test_sub()

def test_mul():
  """Fonction de test pour mul"""
  assert(mul(0, (2,4,6)) == (0,0,0))
  assert(mul(-1, (2,-4,6)) == (-2,4,-6))
  assert(mul(4, (2.1,4,6)) == (8.4,16,24))


test_mul()

def test_dot():
  """Fonction de test pour dot"""
  assert(dot((0,0,0), (0,0,0))== 0)
  assert(dot((1,2,4.5), (0,-2,3)) == 0-4+13.5)

test_dot()

def test_norm():
  """Fonction de test pour norm"""
  assert(norm((0,0,0)) == 0)
  assert(norm((0,-2,4.1)) == sqrt(0**2 + (-2)**2 + 4.1**2))
  

test_norm()

def test_normalize():
  """Fonction de test pour normalize"""
  assert(normalize((0,1,0)) == (0,1,0))
  assert(normalize((1,2,2)) == (((1/3) * 1), ((1/3) *2), ((1/3)*2)))

test_normalize()



### Opérations sur les images
### Une image est représentée par un triplet (i, w, h) où :
### - i est un bytearray de taille w * h * 3
### - w est la largeur de l'image (en pixels)
### - h est la hauteur de l'image (en pixels)

def init_image(w, h):
  """Initialise une image de w pixels de large et h pixel de haut"""
  return (bytearray(w*h*3), w, h)

def set_pixel(img, x, y, c):
  """Met le pixel au coordonnées (x, y) à la couleur c. C'est un est
  triplet (r, v, b) de valeurs. Les valeurs supérieures à 1 (resp. inférieures à
  0) sont mises à 1 (resp. 0).
  """
  for i in range(len(c)):
    px = c[i]
    if px < 0:
      img[0][(img[1]*(y)+x)*3+i] = 0
    elif px > 1:
      img[0][(img[1]*(y)+x)*3+i] = 255
    else:
      img[0][(img[1]*(y)+x)*3+i] = int(px * 255)

def save_image(chemin, img):
  """Écrit l'image img dans le fichier dont le chemin est donné. Si
  le fichier existe, il est supprimé. L'image est stockée au format PPM"""
  buff, w, h = img
  with open(chemin, "wb") as f:
    f.write(b'P6\n')
    f.write(str(w).encode())
    f.write(b' ')
    f.write(str(h).encode())
    f.write(b'\n255\n')
    f.write(buff)

def test_img():
  """Test des fonctions set_pixel et init_image"""
  save_image('black100.ppm', init_image(100,100))
  #affiche un écran noire de dimension 100x100
  #########################################################
  img = init_image(600,400)
  for j in range(200):
    for i in range(200):
      set_pixel(img, 200 + i, 100+j, (246, 203, 100))
  save_image('carre.ppm', img)
  #affiche un carré de couleur blanc sur fond noire et de dimension: 200 x 200 px
  #note: la taille de la fenetre est 600x400 px et le carré se situe a peu prêt au centre de la fenêtre


test_img()

### Fonctions de ray tracing

def pixel_to_point(w, h, xmin, xmax, ymin, ymax, px, py):
  """Convertit un pixel (px, py) en un point du plan."""
  return (((xmax-xmin)*px)/w + xmin, ((ymax-ymin)*py)/h +ymin)


def sphere_intersect (c, r, v, d):
  """Calcule l'intersection entre une sphere de centre c (vecteur) et de rayon r
     (flottant) et une droite passant par v (vecteur) et de direction d (vecteur)
  """
  # cste = constante
  # k2 est la solution  
  a = 1
  b = 2 * dot(d, sub(v, c))
  cste = norm(sub(v, c))*norm(sub(v, c)) - (r*r)
  delta = (b*b) - (4 * a * cste)
  if delta <= 0:
    return None
  else:
    k2 = (-1*b - sqrt(delta))/2 
    if k2 >= 0:
      return k2
    else:
      return None
 
INF = float('inf')

def nearset_intersection(objs, o, d):
  """Renvoie la sphère la plus proche qui intersecte la droite partant de o dans
  la direction d, ainsi que la distance d'intersection depuis o. S'il n'y a pas
  d'intersection, renvoie (None, INF)"""

  # distances est un tableau de tableau avec distances[i][0] l'indice de l'objet, et distances[i][1] la distance
  # m est la distance minimal

  distances = [(i, sphere_intersect(objs[i]['center'], objs[i]['radius'], o, d)) for i in range(len(objs))]
  m = distances[0][1]
  indexMin = distances[0][0]
  
  for i in range(len(distances)):
    if distances[i][1] == None:
      vide = True
    else:
      vide = False
      break

  if vide == False:
    for i in range(len(distances)):
      if m == None and distances[i][1] == None or m != None and distances[i][1] == None:
        pass
      elif (m == None and distances[i][1] != None) or m > distances[i][1]:
        m = distances[i][1]
        indexMin = i
    return (objs[indexMin], m)
  else:
    return (None, INF)

def compute_color (obj, v, n, l):
  """calcule la couleur du point v se trouvant à la surface de l'objet obj.
  n est le vecteur normal au point d'intersection et l le vecteur unitaire dans
  la direction de la source de lumière.
  """
  # note: j'ai du changé la formule du sujet car il y manquait une mise a la valeur absolue 
  # et un symbole "*" a été remplacée par un "+"
  a = obj["ambiant"]
  d = obj["diffuse"] 
  s = obj["specular"]
  alpha = obj["shininess"]

  v1 = add(a, mul(dot(l, n) ,d))
  v2 = mul(abs(dot(n, normalize(add(v, l)))**(alpha/4)), s)
  
  return add(v1, v2)
  
def trace(w, h, xmin, xmax, ymin, ymax, camera, light,objs):

  img = init_image(w, h)

  for py in range(h):
    for px in range(w):
      x, y = pixel_to_point(w, h, xmin, xmax, ymin, ymax, px, py)
      p = (x, y, 0)
      vp = sub(p, camera)
      d = normalize(vp)

      obj, dist = nearset_intersection(objs, camera, d)
      if obj is None:
        couleur = (0, 0, 0)
      else:
        x_point = add(camera, mul(dist, d))
        l = normalize(sub(light, x_point))
        obstacle, dist_obst = nearset_intersection(objs, x_point, l)
        if dist_obst < norm(sub(light, x_point)):
          couleur = (0, 0, 0)
        else:
          n = normalize (sub(x_point, obj["center"]))
          couleur = compute_color(obj, camera, n, l)

      set_pixel(img, px, h - py - 1, couleur)

  return img



### Lecture de description de scene
###

def read_vector(s):
  fields = s.split(",")
  if len(fields) != 3:
    raise ValueError("Erreur de chargement")
  return [ float (n) for n in fields ]

def load_scene (chemin):
  """Charge un fichier de description de scène. En cas d'erreur, la fonction
     lève une exception 'Exception("Erreur de chargement")'
  """
  try:
    with open (chemin, "r") as f:
      w = int (f.readline())
      h = int (f.readline())
      xmin = float(f.readline())
      xmax = float(f.readline())
      ymin = float(f.readline())
      ymax = float(f.readline())
      camera = read_vector(f.readline())
      light = read_vector(f.readline())
      objects = list(DictReader(f, delimiter=";"))
      for obj in objects:
        obj['center'] = read_vector(obj['center'])
        obj['radius'] = float(obj['radius'])
        obj['ambiant'] = read_vector(obj['ambiant'])
        obj['diffuse'] = read_vector(obj['diffuse'])
        obj['specular'] = read_vector(obj['specular'])
        obj['shininess'] = min(100, max(0, float(obj['shininess'])))
        obj['reflection'] = min(1, max(0, float(obj['reflection'])))
    return (w, h, xmin, xmax, ymin, ymax, camera, light, objects)
  except:
      raise ValueError("Erreur de chargement")

def usage():
  print(f"Usage: {argv[0]} <fichier.scene>")
  exit (1)

if __name__ == "__main__":
  if len (argv) != 2:
    usage()

  fichier = argv[1]
  if len(fichier) < 6 or fichier[-6:] != ".scene":
    usage()

  out = fichier[0:-6] + ".ppm"

  w, h, xmin, xmax, ymin, ymax, camera, lum, objs = load_scene(fichier)
  img = trace(w, h, xmin, xmax, ymin, ymax, camera, lum, objs)
  save_image(out, img)
