// ============================================================
//  Prod.cpp  –  Construction CSR, produit A*x et algebre lineaire
//  Projet CHP – Advection 2D par differences finies
//
//  La matrice A n'est pas stockee en plein.
//  On stocke uniquement ses coefficients non nuls et leurs
//  positions via trois std::vector :
//    val     : valeurs non nulles
//    col_ind : indices de colonne correspondants
//    row_ptr : debut de chaque ligne dans val/col_ind
//
//  Nombre de non-zeros theorique (nnz) :
//    Upwind (EU/IU) : 3 coefficients par noeud => nnz = 3*N
//    Centre (EC/IC) : 5 coefficients par noeud => nnz = 5*N
//
//  Numerotation C (0-indexee) : I = j*Nx + i
//  Periodicite : voisin gauche de i=0 est i=Nx-1, etc.
// ============================================================

#include <cmath>
#include <vector>
#include "../Include/Prod.h"
#include "../Include/fonctions.h"
#include "../Include/init.h"

// ----------------------------------------------------------
//  Ajoute un coefficient non nul dans les listes temporaires
//  val_tmp et col_tmp pour la ligne courante.
// ----------------------------------------------------------
void ajouter_nonzero(std::vector<double>& val_tmp,
                     std::vector<int>&    col_tmp,
                     int col, double valeur)
{
    if (std::fabs(valeur) > 0.0) {
        val_tmp.push_back(valeur);
        col_tmp.push_back(col);
    }
}

// ----------------------------------------------------------
//  Construction de la matrice creuse CSR
// ----------------------------------------------------------
void construire_matrice(double t,
                        const Params& p,
                        std::vector<double>& val,
                        std::vector<int>&    col_ind,
                        std::vector<int>&    row_ptr)
{
    const int Nx = p.Nx;
    const int Ny = p.Ny;
    const int N  = p.N;

    val    .clear();
    col_ind.clear();
    row_ptr.clear();
    row_ptr.reserve(N + 1);
    row_ptr.push_back(0);

    // Listes temporaires pour la ligne courante
    std::vector<double> val_tmp;
    std::vector<int>    col_tmp;
    val_tmp.reserve(5);
    col_tmp.reserve(5);

    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {

            val_tmp.clear();
            col_tmp.clear();

            const double xi   = p.xmin + i * p.dx;
            const double yj   = p.ymin + j * p.dy;

            // CFL locaux signes
            const double nu_x = vx(xi, yj, t, p) * p.dt / p.dx;
            const double nu_y = vy(xi, yj, t, p) * p.dt / p.dy;

            // Indices globaux periodiques
            const int I  = j * Nx + i;
            const int Im = j * Nx + (i - 1 + Nx) % Nx;   // gauche
            const int Ip = j * Nx + (i + 1)       % Nx;   // droite
            const int Jm = ((j - 1 + Ny) % Ny) * Nx + i;  // bas
            const int Jp = ((j + 1)       % Ny) * Nx + i;  // haut

            if (p.time_scheme == 1) {
                // ------------------------------------------------
                //  EXPLICITE  :  Un+1 = A * Un
                // ------------------------------------------------
                if (p.space_scheme == 2) {
                    // Upwind explicite
                    // u^{n+1}_I = (1-|nu_x|-|nu_y|)*u^n_I
                    //           + |nu_x|*u^n_amont_x
                    //           + |nu_y|*u^n_amont_y
                    const double anx = std::fabs(nu_x);
                    const double any = std::fabs(nu_y);
                    ajouter_nonzero(val_tmp, col_tmp, I, 1.0 - anx - any);
                    ajouter_nonzero(val_tmp, col_tmp,
                                    (nu_x >= 0.0) ? Im : Ip, anx);
                    ajouter_nonzero(val_tmp, col_tmp,
                                    (nu_y >= 0.0) ? Jm : Jp, any);
                } else {
                    // Centre explicite
                    // u^{n+1}_I = u^n_I
                    //           - (nu_x/2)*(u^n_{i+1,j} - u^n_{i-1,j})
                    //           - (nu_y/2)*(u^n_{i,j+1} - u^n_{i,j-1})
                    ajouter_nonzero(val_tmp, col_tmp, Im,  nu_x / 2.0);
                    ajouter_nonzero(val_tmp, col_tmp, I,   1.0);
                    ajouter_nonzero(val_tmp, col_tmp, Ip, -nu_x / 2.0);
                    ajouter_nonzero(val_tmp, col_tmp, Jm,  nu_y / 2.0);
                    ajouter_nonzero(val_tmp, col_tmp, Jp, -nu_y / 2.0);
                }
            } else {
                // ------------------------------------------------
                //  IMPLICITE  :  A * Un+1 = Un
                // ------------------------------------------------
                if (p.space_scheme == 2) {
                    // Upwind implicite
                    // (1+|nu_x|+|nu_y|)*u^{n+1}_I
                    //   - |nu_x|*u^{n+1}_amont_x
                    //   - |nu_y|*u^{n+1}_amont_y  = u^n_I
                    const double anx = std::fabs(nu_x);
                    const double any = std::fabs(nu_y);
                    ajouter_nonzero(val_tmp, col_tmp,
                                    (nu_x >= 0.0) ? Im : Ip, -anx);
                    ajouter_nonzero(val_tmp, col_tmp,
                                    (nu_y >= 0.0) ? Jm : Jp, -any);
                    ajouter_nonzero(val_tmp, col_tmp, I, 1.0 + anx + any);
                } else {
                    // Centre implicite
                    // -(nu_x/2)*u^{n+1}_{i-1,j} + u^{n+1}_{i,j}
                    // +(nu_x/2)*u^{n+1}_{i+1,j}
                    // -(nu_y/2)*u^{n+1}_{i,j-1}
                    // +(nu_y/2)*u^{n+1}_{i,j+1}  = u^n_{i,j}
                    ajouter_nonzero(val_tmp, col_tmp, Im, -nu_x / 2.0);
                    ajouter_nonzero(val_tmp, col_tmp, Jm, -nu_y / 2.0);
                    ajouter_nonzero(val_tmp, col_tmp, I,   1.0);
                    ajouter_nonzero(val_tmp, col_tmp, Ip,  nu_x / 2.0);
                    ajouter_nonzero(val_tmp, col_tmp, Jp,  nu_y / 2.0);
                }
            }

            // Ajout des coefficients de la ligne I dans val et col_ind
            val    .insert(val    .end(), val_tmp.begin(), val_tmp.end());
            col_ind.insert(col_ind.end(), col_tmp.begin(), col_tmp.end());
            row_ptr.push_back(static_cast<int>(val.size()));
        }
    }
}

// ----------------------------------------------------------
//  Produit matrice-vecteur creux : y = A * x
//  On ne parcourt que les nnz coefficients non nuls.
// ----------------------------------------------------------
void matvec(const std::vector<double>& val,
            const std::vector<int>&    col_ind,
            const std::vector<int>&    row_ptr,
            const std::vector<double>& x,
                  std::vector<double>& y)
{
    const int N = static_cast<int>(row_ptr.size()) - 1;
    y.assign(N, 0.0);
    for (int i = 0; i < N; ++i) {
        double sum = 0.0;
        for (int k = row_ptr[i]; k < row_ptr[i + 1]; ++k)
            sum += val[k] * x[col_ind[k]];
        y[i] = sum;
    }
}

// ----------------------------------------------------------
//  Algebre lineaire
// ----------------------------------------------------------

double dot(const std::vector<double>& u, const std::vector<double>& v)
{
    double s = 0.0;
    const int n = static_cast<int>(u.size());
    for (int i = 0; i < n; ++i) s += u[i] * v[i];
    return s;
}

double norme2(const std::vector<double>& u)
{
    return std::sqrt(dot(u, u));
}

void axpy(double alpha,
          const std::vector<double>& x,
                std::vector<double>& y)
{
    const int n = static_cast<int>(x.size());
    for (int i = 0; i < n; ++i) y[i] += alpha * x[i];
}

void axpby(double alpha,
           const std::vector<double>& x,
           double beta,
                 std::vector<double>& y)
{
    const int n = static_cast<int>(x.size());
    for (int i = 0; i < n; ++i) y[i] = alpha * x[i] + beta * y[i];
}
