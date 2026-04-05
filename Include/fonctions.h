#ifndef FONCTIONS_H
#define FONCTIONS_H

// ============================================================
//  fonctions.h  –  Champs de vitesse et solutions
//  Projet CHP – Advection 2D par differences finies
//
//  Cas 1 : Translation d'une gaussienne  v=(1,0)
//  Cas 2 : Translation d'un cylindre    v=(0.5,0.5)
//  Cas 3 : Rotation d'un cylindre       v=(-y,x)
//  Cas 4 : Sinus manufacture            v=(1,1)
// ============================================================

#include <vector>
#include "init.h"

double vx(double x, double y, double t, const Params& p);
double vy(double x, double y, double t, const Params& p);

double u0     (double x, double y,                const Params& p);
double u_exact(double x, double y, double t,      const Params& p);

void   initialiser_U(std::vector<double>& U,      const Params& p);
double erreur_L2    (const std::vector<double>& U, double t, const Params& p);

#endif // FONCTIONS_H
