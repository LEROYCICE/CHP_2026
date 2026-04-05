// ============================================================
//  main.cpp  –  Programme principal
//  Projet CHP – Advection 2D par differences finies (sequentiel)
//
//  Sorties : Sortie/sol.NNNN.dat  (colonnes : x  y  U)
// ============================================================

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "../Include/init.h"
#include "../Include/fonctions.h"
#include "../Include/Prod.h"
#include "../Include/BICgSTAB.h"

// Frequence d'ecriture : 1 = a chaque pas de temps
const int WRITE_FREQ = 1;

// ----------------------------------------------------------
//  Ecriture de la solution dans Sortie/sol.NNNN.dat
//  Format : x  y  U  (trois colonnes)
//  Ligne vide entre chaque bloc j -> compatible gnuplot pm3d
// ----------------------------------------------------------
void ecrire_solution(const std::vector<double>& U, int ite, const Params& p)
{
    // Construction du nom avec zero-padding sur 4 chiffres
    const std::string num  = std::to_string(ite);
    const std::string nom  = "Sortie/sol."
                           + std::string(4 - num.size(), '0')
                           + num + ".dat";

    std::ofstream f(nom);
    if (!f.is_open())
        throw std::runtime_error("Impossible d'ouvrir " + nom);

    f << std::scientific << std::setprecision(10);
    f << "# Iteration " << ite << "\n# x              y              U\n";

    for (int j = 0; j < p.Ny; ++j) {
        const double yj = p.ymin + j * p.dy;
        for (int i = 0; i < p.Nx; ++i) {
            const double xi = p.xmin + i * p.dx;
            f << xi << "  " << yj << "  " << U[j * p.Nx + i] << "\n";
        }
        f << "\n";   // separateur de bloc pour gnuplot pm3d
    }
}

// ----------------------------------------------------------
//  Verification de stabilite CFL (schemas explicites)
// ----------------------------------------------------------
void verifier_cfl(const Params& p)
{
    if (p.time_scheme != 1) return;

    const double nu_x = std::fabs(p.vx_cst) * p.dt / p.dx;
    const double nu_y = std::fabs(p.vy_cst) * p.dt / p.dy;

    if (p.space_scheme == 2 && nu_x + nu_y > 1.0 + 1e-10)
        std::cerr << "\n[main] ATTENTION : schema EU instable !"
                  << "  nu_x + nu_y = " << nu_x + nu_y << " > 1\n\n";

    if (p.space_scheme == 1)
        std::cerr << "\n[main] INFO : schema Centre Explicite "
                  << "inconditionnellement instable.\n\n";
}

// ----------------------------------------------------------
//  main
// ----------------------------------------------------------
int main()
{
    try {
        // ---- Lecture de para.dat via Lecture() ----
        Params p;
        Lecture("para.dat",
                p.cas,
                p.xmin, p.xmax,
                p.ymin, p.ymax,
                p.Tf,
                p.Nx, p.Ny,
                p.CFL,
                p.space_scheme,
                p.time_scheme);

        calculer_parametres_derives(p);
        afficher_parametres(p);
        verifier_cfl(p);

        // ---- Vecteurs solution ----
        std::vector<double> U(p.N), U1(p.N);
        initialiser_U(U, p);
        ecrire_solution(U, 0, p);

        // ---- Vecteurs CSR de la matrice (remplis a chaque pas) ----
        std::vector<double> val;
        std::vector<int>    col_ind, row_ptr;

        // ---- Boucle en temps ----
        double t      = 0.0;
        int    Nstep  = static_cast<int>(std::ceil(p.Tf / p.dt));
        int    n_ecrit = 0;

        double dt_last = p.Tf - (Nstep - 1) * p.dt;
        if (dt_last < 1e-14 * p.dt) { --Nstep; dt_last = p.dt; }

        std::cout << "  Debut de la simulation : " << Nstep
                  << " pas de temps\n\n";

        const bool has_exact = (p.cas <= 2 || p.cas == 4);
        if (has_exact)
            std::cout << std::setw(8)  << "Pas"
                      << std::setw(14) << "t"
                      << std::setw(14) << "||e||_L2" << "\n";

        for (int n = 0; n < Nstep; ++n) {

            const double dt_cur = (n == Nstep - 1) ? dt_last : p.dt;
            Params p_cur = p;
            p_cur.dt     = dt_cur;

            // Construction de la matrice CSR (val, col_ind, row_ptr)
            construire_matrice(t, p_cur, val, col_ind, row_ptr);

            if (p.time_scheme == 1) {
                // ---- Explicite : U^{n+1} = A * U^n ----
                matvec(val, col_ind, row_ptr, U, U1);

            } else {
                // ---- Implicite : A * U^{n+1} = U^n ----
                U1 = U;   // point de depart pour BiCGSTAB
                const int iters = bicgstab(val, col_ind, row_ptr,
                                           U,   // membre droit b = U^n
                                           U1,  // solution x = U^{n+1}
                                           1e-10, 1000);
                if (n % 50 == 0)
                    std::cout << "  Pas " << std::setw(4) << n
                              << " : BiCGSTAB -> " << iters << " iter\n";
            }

            t += dt_cur;
            std::swap(U, U1);

            if ((n + 1) % WRITE_FREQ == 0 || n == Nstep - 1) {
                ++n_ecrit;
                ecrire_solution(U, n_ecrit, p_cur);
            }

            if (has_exact) {
                const bool dernier = (n == Nstep - 1);
                const bool periode = (n % (Nstep / 10 + 1) == 0);
                if (periode || dernier) {
                    const double err = erreur_L2(U, t, p);
                    std::cout << std::setw(8)  << n + 1
                              << std::setw(14) << std::scientific
                              << std::setprecision(6) << t
                              << std::setw(14) << err << "\n";
                }
            }
        }

        // ---- Bilan final ----
        std::cout << "\n  Simulation terminee.\n"
                  << "  t_final  = " << t << "\n"
                  << "  Fichiers = " << n_ecrit + 1 << "  (dans Sortie/)\n";
        if (has_exact)
            std::cout << "  ||e||_L2 = " << erreur_L2(U, t, p)
                      << "  (a t=" << t << ")\n";
        std::cout << "\n";

    } catch (const std::exception& e) {
        std::cerr << "\n[ERREUR] " << e.what() << "\n";
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
