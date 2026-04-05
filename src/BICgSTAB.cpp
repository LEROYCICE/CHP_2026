// ============================================================
//  BICgSTAB.cpp  –  Solveur BiConjugate Gradient Stabilise
//  Projet CHP – Advection 2D par differences finies
//
//  Reference : https://en.wikipedia.org/wiki/Biconjugate_gradient_stabilized_method
//
//  Resout A * x = b via des produits matrice-vecteur creux.
//  La matrice est representee par ses trois vecteurs CSR
//  (val, col_ind, row_ptr), sans aucune structure intermediaire.
// ============================================================

#include <cmath>
#include <iostream>
#include <vector>
#include "../Include/BICgSTAB.h"
#include "../Include/Prod.h"

int bicgstab(const std::vector<double>& val,
             const std::vector<int>&    col_ind,
             const std::vector<int>&    row_ptr,
             const std::vector<double>& b,
                   std::vector<double>& x,
             double tol,
             int    max_iter)
{
    const int n = static_cast<int>(b.size());

    // Vecteurs de travail
    std::vector<double> r(n), r_hat(n),
                        p_vec(n, 0.0), v(n, 0.0),
                        s(n), t_vec(n), Ax(n);

    // r0 = b - A*x0
    matvec(val, col_ind, row_ptr, x, Ax);
    for (int i = 0; i < n; ++i)
        r[i] = b[i] - Ax[i];

    r_hat = r;   // choix standard : r_hat = r0

    double rho_old = 1.0;
    double alpha   = 1.0;
    double omega   = 1.0;

    double norm_b = norme2(b);
    if (norm_b < 1e-30) norm_b = 1.0;

    int iter = 0;

    while (iter < max_iter) {

        if (norme2(r) / norm_b < tol) break;

        const double rho_new = dot(r_hat, r);

        if (std::fabs(rho_new) < 1e-300) {
            std::cerr << "[BiCGSTAB] rho -> 0 a l'iteration " << iter << "\n";
            break;
        }

        const double beta = (rho_new / rho_old) * (alpha / omega);

        // p = r + beta*(p - omega*v)
        for (int i = 0; i < n; ++i)
            p_vec[i] = r[i] + beta * (p_vec[i] - omega * v[i]);

        // v = A * p
        matvec(val, col_ind, row_ptr, p_vec, v);

        const double r_hat_v = dot(r_hat, v);
        if (std::fabs(r_hat_v) < 1e-300) {
            std::cerr << "[BiCGSTAB] <r_hat,v> -> 0 a l'iteration " << iter << "\n";
            break;
        }

        alpha = rho_new / r_hat_v;

        // s = r - alpha*v
        for (int i = 0; i < n; ++i)
            s[i] = r[i] - alpha * v[i];

        if (norme2(s) / norm_b < tol) {
            axpy(alpha, p_vec, x);
            ++iter;
            break;
        }

        // t = A * s
        matvec(val, col_ind, row_ptr, s, t_vec);

        const double t_t = dot(t_vec, t_vec);
        if (std::fabs(t_t) < 1e-300) {
            std::cerr << "[BiCGSTAB] ||t||^2 -> 0 a l'iteration " << iter << "\n";
            break;
        }

        omega = dot(t_vec, s) / t_t;

        // x = x + alpha*p + omega*s
        for (int i = 0; i < n; ++i)
            x[i] += alpha * p_vec[i] + omega * s[i];

        // r = s - omega*t
        for (int i = 0; i < n; ++i)
            r[i] = s[i] - omega * t_vec[i];

        rho_old = rho_new;
        ++iter;

        if (std::fabs(omega) < 1e-300) {
            std::cerr << "[BiCGSTAB] omega -> 0 a l'iteration " << iter << "\n";
            break;
        }
    }

    if (iter == max_iter)
        std::cerr << "[BiCGSTAB] non convergence apres " << max_iter
                  << " iterations (||r||/||b|| = "
                  << norme2(r) / norm_b << ")\n";

    return iter;
}
