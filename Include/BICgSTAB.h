#ifndef BICGSTAB_H
#define BICGSTAB_H

// ============================================================
//  BICgSTAB.h  –  Solveur iteratif BiConjugate Gradient Stabilise
//  Projet CHP – Advection 2D par differences finies
//
//  Reference : https://en.wikipedia.org/wiki/Biconjugate_gradient_stabilized_method
//
//  Resout A * x = b via des produits matrice-vecteur creux.
//  La matrice est representee par ses trois vecteurs CSR.
// ----------------------------------------------------------
//  Parametres :
//    val, col_ind, row_ptr : matrice A au format CSR
//    b        : membre droit
//    x        : en entree point de depart, en sortie la solution
//    tol      : tolerance ||r||/||b||
//    max_iter : nombre maximal d'iterations
//
//  Retourne le nombre d'iterations effectuees.
// ============================================================

#include <vector>

int bicgstab(const std::vector<double>& val,
             const std::vector<int>&    col_ind,
             const std::vector<int>&    row_ptr,
             const std::vector<double>& b,
                   std::vector<double>& x,
             double tol      = 1e-10,
             int    max_iter = 1000);

#endif // BICGSTAB_H
