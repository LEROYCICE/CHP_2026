// ============================================================
//  init.cpp  –  Lecture de para.dat et parametres derives
//  Projet CHP – Advection 2D par differences finies
// ============================================================

#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include "../Include/init.h"

// ----------------------------------------------------------
//  Lecture du fichier de parametres.
//
//  Format accepte : les valeurs peuvent etre sur une seule ligne
//  ou une par ligne, avec des commentaires (#) et des lignes
//  vides n'importe ou. Seules les valeurs non commentees sont lues,
//  dans l'ordre : cas xmin xmax ymin ymax Tf Nx Ny CFL
//                 space_scheme time_scheme
//
//  Exemples valides :
//    # cas
//    3
//    # xmin xmax
//    -1.0  1.0
//    ...
//  ou sur une ligne :
//    3 -1.0 1.0 -1.0 1.0 2.0 100 50 0.45 2 1
// ----------------------------------------------------------
void Lecture(const std::string& filename,
             int&    cas,
             double& xmin, double& xmax,
             double& ymin, double& ymax,
             double& Tf,
             int&    Nx,   int&    Ny,
             double& CFL,
             int&    space_scheme,
             int&    time_scheme)
{
    std::ifstream fichier(filename);
    if (!fichier.is_open())
        throw std::runtime_error("Lecture : impossible d'ouvrir '" + filename + "'");

    // On collecte tous les tokens (mots non commentes) dans un flux unique
    std::string contenu_utile;
    std::string ligne;

    while (std::getline(fichier, ligne)) {
        // Supprimer tout ce qui suit un '#' sur la ligne
        const std::size_t pos_comment = ligne.find('#');
        if (pos_comment != std::string::npos)
            ligne = ligne.substr(0, pos_comment);

        contenu_utile += ' ' + ligne;
    }

    // Lire les parametres depuis le flux agrege
    std::istringstream flux(contenu_utile);

    if (!(flux >> cas
               >> xmin >> xmax
               >> ymin >> ymax
               >> Tf
               >> Nx   >> Ny
               >> CFL
               >> space_scheme
               >> time_scheme))
        throw std::runtime_error(
            "Lecture : il manque des valeurs dans '" + filename + "'\n"
            "  Ordre attendu : cas xmin xmax ymin ymax Tf Nx Ny CFL "
            "space_scheme time_scheme");

    if (Nx <= 0 || Ny <= 0)
        throw std::runtime_error("Lecture : Nx et Ny doivent etre > 0");
    if (space_scheme < 1 || space_scheme > 2)
        throw std::runtime_error("Lecture : space_scheme doit etre 1 (centre) ou 2 (upwind)");
    if (time_scheme < 1 || time_scheme > 2)
        throw std::runtime_error("Lecture : time_scheme doit etre 1 (explicite) ou 2 (implicite)");
    if (cas < 1 || cas > 4)
        throw std::runtime_error("Lecture : cas doit etre 1, 2, 3 ou 4");
}

// ----------------------------------------------------------
//  Calcul des parametres derives apres lecture
// ----------------------------------------------------------
void calculer_parametres_derives(Params& p)
{
    p.dx = (p.xmax - p.xmin) / p.Nx;
    p.dy = (p.ymax - p.ymin) / p.Ny;
    p.N  = p.Nx * p.Ny;

    switch (p.cas) {
        case 1: p.vx_cst = 1.0; p.vy_cst = 0.0; break;
        case 2: p.vx_cst = 0.5; p.vy_cst = 0.5; break;
        case 3: p.vx_cst = 0.0; p.vy_cst = 0.0; break;
        case 4: p.vx_cst = 1.0; p.vy_cst = 1.0; break;
    }

    double vmax = 0.0;
    switch (p.cas) {
        case 1: vmax = 1.0;                                         break;
        case 2: vmax = std::sqrt(0.5*0.5 + 0.5*0.5);              break;
        case 3: vmax = std::sqrt(p.xmax*p.xmax + p.ymax*p.ymax);  break;
        case 4: vmax = std::sqrt(1.0*1.0 + 1.0*1.0);              break;
    }

    const double h = std::min(p.dx, p.dy);
    p.dt = (vmax > 1e-14) ? p.CFL * h / vmax : p.CFL * h;
}

// ----------------------------------------------------------
//  Affichage des parametres
// ----------------------------------------------------------
void afficher_parametres(const Params& p)
{
    std::cout << "==============================================\n"
              << "  Parametres du probleme\n"
              << "==============================================\n"
              << "  Cas test      : " << p.cas << "\n"
              << "  Domaine       : ["
                  << p.xmin << ", " << p.xmax << "] x ["
                  << p.ymin << ", " << p.ymax << "]\n"
              << "  Temps final   : Tf = " << p.Tf << "\n"
              << "  Grille        : Nx=" << p.Nx
                  << " x Ny=" << p.Ny
                  << "  (N=" << p.N << " inconnues)\n"
              << "  Pas d'espace  : dx=" << p.dx << "  dy=" << p.dy << "\n"
              << "  Pas de temps  : dt=" << p.dt << "\n"
              << "  CFL           : " << p.CFL << "\n"
              << "  Schema espace : "
                  << (p.space_scheme == 1 ? "Centre" : "Upwind") << "\n"
              << "  Schema temps  : "
                  << (p.time_scheme  == 1 ? "Explicite" : "Implicite") << "\n"
              << "==============================================\n";
}