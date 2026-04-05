#ifndef PROD_H
#define PROD_H

// ============================================================
//  Prod.h  –  Stockage creux CSR et algebre lineaire
//  Projet CHP – Advection 2D par differences finies
//
//  La matrice A n'est pas stockee en plein.
//  On stocke uniquement les coefficients non nuls (val) et leurs
//  positions (col_ind, row_ptr) dans trois std::vector<double/int>.
//
//  Format CSR :
//    val     [k] : valeur du k-ieme element non nul
//    col_ind [k] : indice de colonne du k-ieme element non nul
//    row_ptr [i] : indice dans val du premier non-nul de la ligne i
//                  row_ptr[N] = nnz total
// ============================================================

#include <vector>
#include "init.h"

// ----------------------------------------------------------
//  Construction de la matrice creuse CSR
//
//  Remplit les trois vecteurs val, col_ind, row_ptr selon le
//  schema (espace + temps) et le temps t (vitesse locale).
//
//  Schema EXPLICITE : Un+1 = A * Un   (coefficients du membre droit)
//  Schema IMPLICITE : A * Un+1 = Un   (coefficients du membre gauche)
// ----------------------------------------------------------
void construire_matrice(double t,
                        const Params& p,
                        std::vector<double>& val,
                        std::vector<int>&    col_ind,
                        std::vector<int>&    row_ptr);

// ----------------------------------------------------------
//  Produit matrice-vecteur creux : y = A * x
//  Exploite uniquement les nnz coefficients non nuls.
// ----------------------------------------------------------
void matvec(const std::vector<double>& val,
            const std::vector<int>&    col_ind,
            const std::vector<int>&    row_ptr,
            const std::vector<double>& x,
                  std::vector<double>& y);

// ----------------------------------------------------------
//  Algebre lineaire sur std::vector<double>
// ----------------------------------------------------------

double dot   (const std::vector<double>& u,
              const std::vector<double>& v);

double norme2(const std::vector<double>& u);

void   axpy  (double alpha,
              const std::vector<double>& x,
                    std::vector<double>& y);

void   axpby (double alpha,
              const std::vector<double>& x,
              double beta,
                    std::vector<double>& y);

#endif // PROD_H
