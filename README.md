# Projet CHP — Advection 2D par différences finies avec MPI

**Cours** :Calcul Haute Performance  
**École** : Enseirb-Matmeca, Master 2 – 2025/2026  


---

## Table des matières

1. [Présentation du problème](#1-présentation-du-problème)
2. [Schémas numériques implémentés](#2-schémas-numériques-implémentés)
3. [Structure du projet](#3-structure-du-projet)
4. [Compilation et exécution](#4-compilation-et-exécution)
5. [Le fichier de paramètres `para.dat`](#5-le-fichier-de-paramètres-paradat)
6. [Description détaillée de chaque fichier](#6-description-détaillée-de-chaque-fichier)
   - [init.h / init.cpp](#61-inith--initcpp)
   - [fonctions.h / fonctions.cpp](#62-fonctionsh--fonctionscpp)
   - [Prod.h / Prod.cpp](#63-prodh--prodcpp)
   - [BICgSTAB.h / BICgSTAB.cpp](#64-bicgstabh--bicgstabcpp)
   - [main.cpp](#65-maincpp)
7. [Classes et types C++ utilisés — guide pour débutants](#7-classes-et-types-c-utilisés--guide-pour-débutants)
8. [Visualisation — génération des GIFs animés](#8-visualisation--génération-des-gifs-animés)
9. [Cas tests](#9-cas-tests)

---

## 1. Présentation du problème

On résout l'**équation d'advection scalaire 2D** sur un domaine rectangulaire Ω = [xmin, xmax] × [ymin, ymax] avec des **conditions aux bords périodiques** :

```
∂u/∂t + v · ∇u = 0
u(xmin, y, t) = u(xmax, y, t)   ∀ y, t
u(x, ymin, t) = u(x, ymax, t)   ∀ x, t
```

Le champ de vitesse `v = (vx, vy)` transporte la quantité scalaire `u(x, y, t)`.  
La solution exacte est connue pour les cas à vitesse constante : la condition initiale se déplace sans se déformer, et réapparaît de l'autre côté du domaine grâce à la périodicité.

---

## 2. Schémas numériques implémentés

Le domaine spatial est discrétisé avec **Nx × Ny nœuds**.  
La numérotation est **en C (0-indexée)** : le nœud `(i, j)` a l'indice global `I = j * Nx + i`.

> ⚠️ Les bords gauche (`i=0`) et bas (`j=0`) **ne sont pas stockés** car les conditions périodiques imposent `u[i=0] = u[i=Nx]`. On stocke donc les indices `i ∈ [0, Nx-1]` et `j ∈ [0, Ny-1]`.

Les nombres de CFL sont définis par :
```
νx = vx · dt / dx      νy = vy · dt / dy
```

### Schémas en espace × temps

| Code | Espace | Temps | Stabilité |
|------|--------|-------|-----------|
| `space_scheme=2, time_scheme=1` | Upwind | Explicite (EU) | Stable si νx + νy ≤ 1 |
| `space_scheme=1, time_scheme=1` | Centré | Explicite (EC) | Inconditionnellement **instable** |
| `space_scheme=2, time_scheme=2` | Upwind | Implicite (IU) | Inconditionnellement stable |
| `space_scheme=1, time_scheme=2` | Centré | Implicite (IC) | Inconditionnellement stable |

### Schéma Upwind Explicite (EU) — recommandé

Pour chaque nœud `(i, j)` :

```
u^{n+1}_{i,j} = (1 - |νx| - |νy|) · u^n_{i,j}
              + |νx| · u^n_{i-1,j}   (si vx ≥ 0, amont gauche)
              + |νy| · u^n_{i,j-1}   (si vy ≥ 0, amont bas)
```

Condition de stabilité CFL : `νx + νy ≤ 1` ⟹ choisir `CFL ≤ 0.45` en pratique.

### Schéma Implicite — résolution par BiCGSTAB

Pour le schéma implicite, on résout à chaque pas de temps le système linéaire :

```
A · U^{n+1} = U^n
```

La matrice `A` est **creuse** (3 à 5 coefficients non nuls par ligne). Elle n'est **jamais stockée en entier** : on ne stocke que ses coefficients non nuls au format **CSR** (voir section Prod).

---

## 3. Structure du projet

```
CHP_2026/
│
├── Include/                   ← Fichiers d'en-tête (.h)
│   ├── init.h                 ← Structure Params + déclarations init
│   ├── fonctions.h            ← Déclarations vx, vy, u0, u_exact, ...
│   ├── Prod.h                 ← Déclarations CSR + algèbre linéaire
│   └── BICgSTAB.h             ← Déclaration du solveur BiCGSTAB
│
├── src/                       ← Fichiers source (.cpp)
│   ├── init.cpp               ← Lecture de para.dat, paramètres dérivés
│   ├── fonctions.cpp          ← Vitesses, conditions initiales, erreur L2
│   ├── Prod.cpp               ← Construction CSR, matvec, dot, axpy, ...
│   ├── BICgSTAB.cpp           ← Algorithme BiCGSTAB complet
│   └── main.cpp               ← Boucle en temps, écriture des sorties
│
├── Sortie/                    ← Fichiers de résultats (sol.XXXX.dat)
├── Gnuplot/                  ← Figures et GIFs animés
│
├── para.dat                   ← Fichier de paramètres (éditable)
├── Makefile                   ← Compilation automatisée
├── gnuplot.gp                 ← Script gnuplot pour visualisation
```

---

## 4. Compilation et exécution

### Prérequis

- Compilateur C++11 ou supérieur : `g++`
- `make`
- Python 3 + matplotlib (pour les GIFs)

### Compiler

```bash
make
```

Cela compile tous les `.cpp` et produit l'exécutable `./advection`.

### Lancer une simulation

```bash
./exec
```

Le programme lit automatiquement `para.dat` dans le répertoire courant, lance la simulation et écrit les résultats dans `Sortie/sol.XXXX.dat`.

### Autres commandes make

```bash
make clean       # supprime les fichiers objets (.o) et l'exécutable
make distclean   # idem + supprime les fichiers Sortie/*.dat
make help        # affiche l'aide
```

## 5. Le fichier de paramètres `para.dat`

Le fichier `para.dat` est lu par la fonction `Lecture()`.  
Il accepte **n'importe quelle mise en forme** : valeurs sur une ligne, une par ligne, avec des commentaires `#` n'importe où (en début de ligne, en fin de ligne, ou sur une ligne entière). L'ordre de lecture est fixe.

### Exemple — format commenté, une valeur par ligne

```
# cas (1=gaussienne, 2=cylindre, 3=rotation, 4=sinus manufacture)
1
# xmin   xmax
-1.0     1.0
# ymin   ymax
-0.5     0.5
# Tf (temps final)
2.0
# Nx     Ny
100      50
# CFL
0.45
# space_scheme (1=centre  2=upwind)
2
# time_scheme  (1=explicite  2=implicite)
1
```

### Exemple — format compact sur une ligne

```
# cas xmin xmax ymin ymax Tf Nx Ny CFL space time
1 -1.0 1.0 -0.5 0.5 2.0 100 50 0.45 2 1
```

### Exemple — valeur puis commentaire sur la même ligne

```
1       # cas
-1.0    # xmin
1.0     # xmax
-0.5    # ymin
0.5     # ymax
2.0     # Tf
100     # Nx
50      # Ny
0.45    # CFL
2       # space_scheme
1       # time_scheme
```

### Description de chaque paramètre

| Paramètre | Type | Description |
|-----------|------|-------------|
| `cas` | entier | Numéro du cas test (1 à 4) |
| `xmin`, `xmax` | réels | Bornes du domaine en x |
| `ymin`, `ymax` | réels | Bornes du domaine en y |
| `Tf` | réel | Temps final de la simulation |
| `Nx`, `Ny` | entiers | Nombre de nœuds dans chaque direction |
| `CFL` | réel | Nombre CFL cible (< 1 pour upwind explicite) |
| `space_scheme` | entier | Schéma en espace : 1 = centré, 2 = upwind |
| `time_scheme` | entier | Schéma en temps : 1 = explicite, 2 = implicite |

---

## 6. Description détaillée de chaque fichier

---

### 6.1 `init.h` / `init.cpp`

Ce module gère la **lecture du fichier de paramètres** et le calcul des quantités dérivées.

#### La structure `Params`

`Params` est une structure C++ qui regroupe tous les paramètres du problème dans un seul objet. On peut penser à une structure comme à une "boîte" qui contient plusieurs variables nommées.

```cpp
struct Params {
    int    cas;            // numéro du cas test
    double xmin, xmax;    // bornes du domaine en x
    double ymin, ymax;    // bornes du domaine en y
    double Tf;            // temps final
    int    Nx, Ny;        // nombre de nœuds
    double CFL;           // nombre CFL
    int    space_scheme;  // 1=centré, 2=upwind
    int    time_scheme;   // 1=explicite, 2=implicite

    // Calculés automatiquement après lecture :
    double dx, dy;        // pas d'espace
    double dt;            // pas de temps
    double vx_cst, vy_cst; // vitesses constantes (cas 1,2,4)
    int    N;             // nombre total d'inconnues = Nx * Ny
};
```

#### Fonctions

```cpp
void Lecture(const std::string& filename,
             int& cas, double& xmin, double& xmax,
             double& ymin, double& ymax, double& Tf,
             int& Nx, int& Ny, double& CFL,
             int& space_scheme, int& time_scheme);
```
Ouvre `filename`, ignore tous les commentaires (`#`) et lit les 11 paramètres dans l'ordre, qu'ils soient sur une ou plusieurs lignes. Les paramètres sont passés **par référence** (`&`) : cela signifie que la fonction modifie directement les variables de l'appelant, sans faire de copie. Lance une exception `std::runtime_error` si le fichier est introuvable ou si des valeurs manquent.

```cpp
void calculer_parametres_derives(Params& p);
```
À partir des paramètres lus, calcule `dx`, `dy`, `N`, `vx_cst`, `vy_cst`, et surtout `dt` via la formule CFL : `dt = CFL × min(dx, dy) / vmax`.

```cpp
void afficher_parametres(const Params& p);
```
Affiche un tableau récapitulatif de tous les paramètres dans le terminal.

---

### 6.2 `fonctions.h` / `fonctions.cpp`

Ce module contient les **données physiques** du problème : champ de vitesse, condition initiale, solution exacte.

#### Les 4 cas tests

| Cas | Description | v(x,y) | u₀(x,y) |
|-----|-------------|---------|---------|
| 1 | Translation d'une gaussienne | (1, 0) | exp(-(x²+y²)/0.0075) |
| 2 | Translation d'un cylindre | (0.5, 0.5) | 1 si ‖(x,y)‖ ≤ 0.4, sinon 0 |
| 3 | Rotation d'un cylindre | (-y, x) | 1 si ‖(x,y)‖ ≤ 0.4, sinon 0 |
| 4 | Sinus manufacturé (convergence) | (1, 1) | sin(πx)·sin(πy) |

Le cas 4 a une solution exacte fermée `u(x,y,t) = sin(π(x-t))·sin(π(y-t))` qui permet de mesurer l'ordre de convergence du schéma.

#### Fonctions

```cpp
double vx(double x, double y, double t, const Params& p);
double vy(double x, double y, double t, const Params& p);
```
Retournent les composantes de la vitesse au point `(x, y)` à l'instant `t`. Pour les cas 1, 2, 4 la vitesse est constante. Pour le cas 3 elle dépend de la position : `vx = -y`, `vy = x`.

```cpp
double u0(double x, double y, const Params& p);
```
Retourne la condition initiale `u(x, y, t=0)` selon le cas.

```cpp
double u_exact(double x, double y, double t, const Params& p);
```
Retourne la solution analytique exacte `u(x, y, t)`. Pour les cas 1 et 2 (transport à vitesse constante), on remonte au point d'origine en soustrayant `v·t` et on applique la périodicité. Pour le cas 3 on applique une rotation inverse.

```cpp
void initialiser_U(std::vector<double>& U, const Params& p);
```
Remplit le vecteur solution `U` avec les valeurs de `u0` sur toute la grille. Le nœud `(i, j)` est stocké à l'indice `I = j * Nx + i`.

```cpp
double erreur_L2(const std::vector<double>& U, double t, const Params& p);
```
Calcule l'erreur L2 discrète entre la solution numérique `U` et la solution exacte à l'instant `t` :  
`‖e‖_L₂ = sqrt(dx·dy · Σ (U_I - u_exact(x_i, y_j, t))²)`

---

### 6.3 `Prod.h` / `Prod.cpp`

C'est le module le plus important du projet. Il implémente le **stockage creux de la matrice** et toutes les **opérations d'algèbre linéaire**.

#### Pourquoi ne pas stocker la matrice entière ?

Pour un problème de taille N = Nx × Ny, la matrice A est de taille N×N. Si N = 10 000, cela représente 10 000² = 100 millions de valeurs double, soit **800 Mo de mémoire**, dont la grande majorité serait des zéros. C'est inacceptable.

La matrice est **creuse** : chaque ligne ne contient que 3 coefficients non nuls (upwind) ou 5 (centré). On ne stocke que ceux-là.

#### Le format CSR (Compressed Sparse Row)

Le format CSR représente une matrice creuse avec trois tableaux :

```
val     : toutes les valeurs non nulles, rangées ligne par ligne
col_ind : l'indice de colonne correspondant à chaque valeur dans val
row_ptr : row_ptr[i] = position dans val du premier élément de la ligne i
          row_ptr[N] = nombre total de non-zéros (nnz)
```

**Exemple** pour une matrice 4×4 avec des coefficients sur la diagonale et la sous-diagonale :

```
Matrice :          val     = [2, 1, 2, 1, 2, 1, 2]
2 0 0 0            col_ind = [0, 0, 1, 1, 2, 2, 3]
1 2 0 0    →       row_ptr = [0, 1, 3, 5, 7]
0 1 2 0
0 0 1 2
```

Pour lire la ligne `i`, on accède aux éléments de `val[row_ptr[i]]` jusqu'à `val[row_ptr[i+1] - 1]`.

**Avantage** : pour N = 10 000 avec upwind (3 nnz/ligne), on stocke `3 × 10 000 = 30 000` valeurs au lieu de 100 millions. Facteur **3 300** de réduction mémoire.

#### Fonctions

```cpp
void construire_matrice(double t, const Params& p,
                        std::vector<double>& val,
                        std::vector<int>&    col_ind,
                        std::vector<int>&    row_ptr);
```
Construit la matrice CSR en parcourant tous les nœuds `(i, j)`. Pour chaque nœud, calcule la vitesse locale, les CFL locaux, et ajoute les coefficients non nuls correspondant au schéma choisi (explicite ou implicite, upwind ou centré). La matrice est reconstruite à chaque pas de temps (nécessaire pour le cas 3 où la vitesse dépend de la position).

```cpp
void matvec(const std::vector<double>& val,
            const std::vector<int>&    col_ind,
            const std::vector<int>&    row_ptr,
            const std::vector<double>& x,
                  std::vector<double>& y);
```
Calcule le produit `y = A · x` en n'utilisant que les coefficients non nuls. Pour chaque ligne `i`, parcourt uniquement les `row_ptr[i+1] - row_ptr[i]` éléments non nuls : `y[i] = Σ val[k] * x[col_ind[k]]`.

```cpp
double dot(const std::vector<double>& u, const std::vector<double>& v);
```
Produit scalaire : retourne `Σ u[i] * v[i]`.

```cpp
double norme2(const std::vector<double>& u);
```
Norme Euclidienne : retourne `sqrt(Σ u[i]²)`. Utilise `dot(u, u)` en interne.

```cpp
void axpy(double alpha, const std::vector<double>& x, std::vector<double>& y);
```
Opération BLAS de base : `y = y + alpha * x`. Le nom vient de "Alpha X Plus Y".

```cpp
void axpby(double alpha, const std::vector<double>& x,
           double beta,  std::vector<double>& y);
```
Opération généralisée : `y = alpha * x + beta * y`.

---

### 6.4 `BICgSTAB.h` / `BICgSTAB.cpp`

Ce module implémente le solveur itératif **BiCGSTAB** (BiConjugate Gradient Stabilized).

#### Pourquoi BiCGSTAB ?

Pour le schéma implicite, on doit résoudre `A · U^{n+1} = U^n` à chaque pas de temps. La matrice A est creuse et non symétrique. BiCGSTAB est un algorithme itératif adapté à ce type de système : il construit progressivement une solution approchée en utilisant uniquement des produits matrice-vecteur, sans jamais factoriser la matrice.

#### Algorithme (référence : Wikipedia)

À partir d'un vecteur initial `x₀`, BiCGSTAB produit une suite de vecteurs qui convergent vers la solution exacte. À chaque itération :
1. On calcule le résidu `r = b - A·x`
2. On construit des directions de descente `p`
3. On effectue des pas en minimisant le résidu

L'algorithme s'arrête quand `‖r‖ / ‖b‖ < tol` (la tolérance, par défaut `1e-10`).

#### Fonction

```cpp
int bicgstab(const std::vector<double>& val,
             const std::vector<int>&    col_ind,
             const std::vector<int>&    row_ptr,
             const std::vector<double>& b,
                   std::vector<double>& x,
             double tol      = 1e-10,
             int    max_iter = 1000);
```

- **Entrée** : la matrice `A` au format CSR (`val`, `col_ind`, `row_ptr`), le membre droit `b`, et un vecteur `x` servant de point de départ (on passe `U^n` comme point de départ).
- **Sortie** : `x` est modifié en place pour contenir la solution `U^{n+1}`.
- **Retour** : nombre d'itérations effectuées.
- Les paramètres `tol` et `max_iter` ont des **valeurs par défaut** (`= 1e-10` et `= 1000`), ce qui signifie qu'on n'est pas obligé de les passer lors de l'appel.

---

### 6.5 `main.cpp`

Programme principal qui orchestre toute la simulation.

#### Fonctions

```cpp
void ecrire_solution(const std::vector<double>& U, int ite, const Params& p);
```
Écrit la solution à l'itération `ite` dans le fichier `Sortie/sol.XXXX.dat` avec trois colonnes : `x  y  U`. Une ligne vide sépare chaque bloc de `j` constant (format compatible gnuplot `pm3d`).

```cpp
void verifier_cfl(const Params& p);
```
Pour le schéma upwind explicite, vérifie que `νx + νy ≤ 1`. Si ce n'est pas le cas, affiche un avertissement. Pour le schéma centré explicite, affiche systématiquement une info (ce schéma est inconditionnellement instable).

#### Boucle en temps

```
Pour chaque pas de temps n = 0, ..., Nstep-1 :
  1. Construire la matrice CSR au temps t  →  construire_matrice(t, p, val, col_ind, row_ptr)
  2a. Si EXPLICITE : calculer U^{n+1} = A · U^n  →  matvec(...)
  2b. Si IMPLICITE : résoudre A · U^{n+1} = U^n  →  bicgstab(...)
  3. Avancer le temps : t += dt
  4. Échanger U et U1  →  std::swap(U, U1)
  5. Écrire la solution dans Sortie/sol.XXXX.dat
  6. Calculer et afficher l'erreur L2 (cas 1, 2, 4)
```

L'échange à l'étape 4 utilise `std::swap(U, U1)` qui échange les contenus des deux vecteurs en temps constant (O(1)) — sans copier les données.

---

## 7. Classes et types C++ utilisés:

Cette section explique chaque classe et mécanisme C++ utilisé dans le projet.

---

### `std::string`

**Fichier d'inclusion** : `#include <string>`

Une chaîne de caractères. Contrairement aux tableaux `char[]` du C, `std::string` gère automatiquement la mémoire et offre des opérations pratiques.

```cpp
std::string filename = "para.dat";
std::string message  = "Erreur : " + filename;  // concaténation avec +
filename.find('#');    // cherche le caractère '#', retourne sa position
filename.substr(0, 5); // sous-chaîne de la position 0, de longueur 5
filename.empty();      // true si la chaîne est vide
```

---

### `std::vector<T>`

**Fichier d'inclusion** : `#include <vector>`

Un tableau de taille **dynamique** (qui peut grandir ou rétrécir). C'est le conteneur le plus utilisé en C++ moderne. `T` est le type des éléments : `std::vector<double>` pour des réels, `std::vector<int>` pour des entiers.

```cpp
std::vector<double> U(1000);       // vecteur de 1000 doubles initialisés à 0
std::vector<double> v(50, 3.14);   // vecteur de 50 éléments valant 3.14
U[42] = 1.5;                       // accès par indice comme un tableau
U.size();                          // nombre d'éléments
U.resize(2000);                    // redimensionner
U.assign(500, 0.0);                // remettre à zéro avec 500 éléments
U.push_back(9.9);                  // ajouter un élément à la fin
U.clear();                         // vider le vecteur (size devient 0)
U.reserve(5000);                   // réserver de la mémoire sans changer size
U.begin(), U.end();                // itérateurs (pour les boucles, insert, etc.)
U.insert(U.end(), V.begin(), V.end()); // copier V à la fin de U
```

**Pourquoi utiliser `std::vector` plutôt que `double tab[1000]` ?**
- Taille connue uniquement à l'exécution (pas besoin de la connaître à la compilation)
- Gestion automatique de la mémoire (pas de `new`/`delete`)
- Copie simple avec `=` : `std::vector<double> copie = original;`
- Échange rapide avec `std::swap`

---

### `std::swap`

**Fichier d'inclusion** : `#include <vector>` (inclus avec vector)

Échange le contenu de deux variables **sans copier les données**.

```cpp
std::swap(U, U1);  // échange les contenus de U et U1 en O(1)
```

Pour deux `std::vector`, `swap` échange simplement les pointeurs internes — c'est instantané même pour des millions d'éléments. C'est pour ça qu'on l'utilise à chaque pas de temps au lieu de copier U dans U1.

---

### `std::ifstream`

**Fichier d'inclusion** : `#include <fstream>`

Flux de lecture depuis un fichier ("input file stream").

```cpp
std::ifstream fichier("para.dat");  // ouvre para.dat en lecture
if (!fichier.is_open())             // vérifie si l'ouverture a réussi
    throw ...;
std::string ligne;
std::getline(fichier, ligne);       // lit une ligne entière
fichier >> valeur;                  // lit un token (mot ou nombre)
// Le fichier se ferme automatiquement quand fichier sort de scope
```

---

### `std::ofstream`

**Fichier d'inclusion** : `#include <fstream>`

Flux d'écriture dans un fichier ("output file stream").

```cpp
std::ofstream f("Sortie/sol.0001.dat");  // crée ou écrase le fichier
if (!f.is_open()) throw ...;
f << std::scientific << std::setprecision(10);  // format scientifique
f << x << "  " << y << "  " << U[I] << "\n";   // écriture ligne par ligne
// Le fichier se ferme automatiquement
```

---

### `std::istringstream`

**Fichier d'inclusion** : `#include <sstream>`

Flux de lecture depuis une **chaîne de caractères** (au lieu d'un fichier). Très pratique pour parser du texte.

```cpp
std::string ligne = "3 -1.0 1.0 100";
std::istringstream iss(ligne);
int cas; double xmin, xmax; int Nx;
iss >> cas >> xmin >> xmax >> Nx;  // lit depuis la chaîne comme depuis cin
```

Dans `init.cpp`, on concatène toutes les lignes non-commentées en une seule chaîne, puis on la parse avec `istringstream` — ce qui permet de lire les paramètres quelle que soit leur disposition dans le fichier.

---

### `std::cout` et `std::cerr`

**Fichier d'inclusion** : `#include <iostream>`

Flux de sortie standard :
- `std::cout` : sortie normale (terminal)
- `std::cerr` : sortie d'erreur (terminal, non bufferisée, pour les avertissements)

```cpp
std::cout << "Résultat : " << valeur << "\n";
std::cerr << "[ERREUR] Fichier introuvable\n";
```

L'opérateur `<<` s'appelle l'opérateur d'insertion. On peut enchaîner plusieurs `<<` sur la même ligne.

**Modificateurs de format** (`#include <iomanip>`) :
```cpp
std::cout << std::scientific       // notation scientifique : 1.23e+04
          << std::setprecision(10) // 10 chiffres après la virgule
          << std::setw(8)          // largeur minimale de 8 caractères
          << valeur << "\n";
```

---

### `std::runtime_error`

**Fichier d'inclusion** : `#include <stdexcept>`

Une **exception** C++ qui représente une erreur survenue à l'exécution.

```cpp
throw std::runtime_error("Fichier 'para.dat' introuvable");
```

`throw` interrompt l'exécution de la fonction courante et remonte la pile d'appels jusqu'à un bloc `try/catch`. Si personne n'attrape l'exception, le programme s'arrête avec le message d'erreur.

---

### `try` / `catch`

Mécanisme de gestion d'erreurs en C++. Dans `main.cpp` :

```cpp
try {
    // Tout le code de la simulation
    Lecture("para.dat", ...);
    ...
} catch (const std::exception& e) {
    // Si une exception est lancée n'importe où dans le try :
    std::cerr << "\n[ERREUR] " << e.what() << "\n";
    return EXIT_FAILURE;
}
```

`e.what()` retourne le message d'erreur de l'exception.  
`EXIT_FAILURE` est une constante (valeur 1) qui indique au système d'exploitation que le programme a échoué.

---

### `const` et `&` dans les paramètres

Deux concepts clés en C++ pour passer des arguments à une fonction :

**Par référence (`&`)** : la fonction reçoit une "adresse" vers la variable originale et peut la modifier.
```cpp
void remplir(std::vector<double>& U);  // modifie U dans l'appelant
```

**Par référence constante (`const &`)** : la fonction accède à la variable mais ne peut pas la modifier. Cela évite une copie coûteuse.
```cpp
void lire(const std::vector<double>& U);  // lit U sans le modifier, sans copie
```

**Sans référence** : la fonction reçoit une copie. Toute modification est locale.
```cpp
void f(double x);  // x est une copie, l'original n'est pas modifié
```

En règle générale :
- Grands objets (`std::vector`, `std::string`, `Params`) → toujours passer par `const&` en lecture, par `&` en écriture
- Petits types (`int`, `double`) → passer par valeur (copie)

---

### Valeurs par défaut dans les fonctions

```cpp
int bicgstab(..., double tol = 1e-10, int max_iter = 1000);
```

Les paramètres avec `=` ont une valeur par défaut. On peut les omettre lors de l'appel :

```cpp
bicgstab(val, col_ind, row_ptr, b, x);          // utilise tol=1e-10, max_iter=1000
bicgstab(val, col_ind, row_ptr, b, x, 1e-8);    // tol=1e-8, max_iter=1000
bicgstab(val, col_ind, row_ptr, b, x, 1e-8, 500); // tol=1e-8, max_iter=500
```

---

### `std::to_string`

Convertit un entier ou un réel en `std::string`.

```cpp
int ite = 42;
std::string s = std::to_string(ite);  // s = "42"
// Utilisé dans main.cpp pour construire le nom de fichier sol.0042.dat
```

---

### `std::min` / `std::max`

**Fichier d'inclusion** : `#include <algorithm>` (ou souvent déjà inclus)

```cpp
double h = std::min(p.dx, p.dy);  // retourne le plus petit des deux
```

---

### `std::ceil`, `std::fabs`, `std::sqrt`, `std::floor`

**Fichier d'inclusion** : `#include <cmath>`

Fonctions mathématiques standard :
```cpp
std::ceil(2.3)   // → 3.0  (arrondi supérieur)
std::floor(2.7)  // → 2.0  (arrondi inférieur)
std::fabs(-3.5)  // → 3.5  (valeur absolue pour double)
std::sqrt(4.0)   // → 2.0  (racine carrée)
std::sin(M_PI)   // → 0.0  (sinus, M_PI ≈ 3.14159...)
std::cos(0.0)    // → 1.0  (cosinus)
std::exp(-1.0)   // → 0.368 (exponentielle)
```

---

## 8. Visualisation — génération des GIFs animés

Après une simulation, les fichiers `Sortie/sol.XXXX.dat` contiennent la solution à chaque instant.

### Avec gnuplot

```bash
gnuplot gnuplot.gp
```

### Format des fichiers de sortie

Chaque fichier `Sortie/sol.XXXX.dat` contient trois colonnes :

```
# Iteration 42
# x              y              U
-1.0000000000e+00  -5.0000000000e-01  0.0000000000e+00
-9.8000000000e-01  -5.0000000000e-01  3.4512e-05
...
                                       ← ligne vide entre blocs j
```

La ligne vide entre chaque bloc de `j` constant est nécessaire pour que gnuplot puisse dessiner une surface `pm3d`.

---

## 9. Cas tests

### Cas 1 — Translation d'une gaussienne

```
# cas
1
# xmin xmax
-1.0  1.0
# ymin ymax
-0.5  0.5
# Tf
2.0
# Nx Ny
100  50
# CFL
0.45
# space time
2  1
```

La bosse gaussienne se translate vers la droite à vitesse 1, disparaît à droite et réapparaît à gauche par périodicité.

### Cas 2 — Translation d'un cylindre

```
2  -1.0  1.0  -1.0  1.0  4.0  80  80  0.4  2  1
```

Le cylindre (bosse plate de hauteur 1) se déplace en diagonale.

### Cas 3 — Rotation d'un cylindre

```
# cas
3
# xmin xmax ymin ymax Tf
-1.0  1.0  -1.0  1.0  6.2832
# Nx Ny CFL space time
80   80  0.35  2   1
```

`Tf = 2π` correspond à un tour complet. On peut vérifier que le cylindre revient à sa position initiale (avec la diffusion numérique du schéma upwind).

### Cas 4 — Solution manufacturée (étude de convergence)

```
4  -1.0  1.0  -1.0  1.0  0.5  60  60  0.8  1  2
```

Solution exacte connue `u = sin(π(x-t))·sin(π(y-t))`. Permet de calculer l'erreur L2 et de vérifier l'ordre de convergence du schéma.

---

*Projet réalisé dans le cadre du cours de Calcul Haute Performance, Enseirb-Matmeca M2 2025-2026.*
