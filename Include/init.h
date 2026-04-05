#ifndef INIT_H
#define INIT_H

// ============================================================
//  init.h  –  Parametres du probleme et lecture de para.dat
//  Projet CHP – Advection 2D par differences finies
// ============================================================

#include <string>

// Structure contenant tous les parametres du probleme
struct Params {
    int    cas;
    double xmin, xmax;
    double ymin, ymax;
    double Tf;
    int    Nx, Ny;
    double CFL;
    int    space_scheme; // 1=centre, 2=upwind
    int    time_scheme;  // 1=explicite, 2=implicite

    // Parametres derives
    double dx, dy, dt;
    double vx_cst, vy_cst;
    int    N;
};

// Lit le fichier filename et remplit les parametres individuellement.
// Lance std::runtime_error si le fichier est introuvable ou mal forme.
void Lecture(const std::string& filename,
             int&    cas,
             double& xmin, double& xmax,
             double& ymin, double& ymax,
             double& Tf,
             int&    Nx,   int&    Ny,
             double& CFL,
             int&    space_scheme,
             int&    time_scheme);

// Calcule les parametres derives (dx, dy, dt, N, vx_cst, vy_cst)
// a partir des parametres lus et les stocke dans p.
void calculer_parametres_derives(Params& p);

// Affiche un resume des parametres sur std::cout.
void afficher_parametres(const Params& p);

#endif // INIT_H
