// ============================================================
//  fonctions.cpp  –  Champs de vitesse, solution initiale
//                    et solution analytique
//  Projet CHP – Advection 2D par differences finies
//
//  Numerotation C (0-indexee) :
//    noeud (i, j) -> I = j * Nx + i
//    x_i = xmin + i * dx,  i = 0, ..., Nx-1
//    y_j = ymin + j * dy,  j = 0, ..., Ny-1
// ============================================================

#include <cmath>
#include <vector>
#include "../Include/fonctions.h"
#include "../Include/init.h"

// ----------------------------------------------------------
//  Champ de vitesse
// ----------------------------------------------------------

double vx(double x, double y, double t, const Params& p)
{
    (void)x; (void)t;
    switch (p.cas) {
        case 1:  return 1.0;
        case 2:  return 0.5;
        case 3:  return -y;        // v = (-y, x) : rotation solide
        case 4:  return 1.0;       // sinus manufacture, v=(1,1)
        default: return 0.0;
    }
}

double vy(double x, double y, double t, const Params& p)
{
    (void)y; (void)t;
    switch (p.cas) {
        case 1:  return 0.0;
        case 2:  return 0.5;
        case 3:  return  x;        // v = (-y, x) : rotation solide
        case 4:  return 1.0;       // sinus manufacture, v=(1,1)
        default: return 0.0;
    }
}

// ----------------------------------------------------------
//  Solution initiale u0(x, y)
// ----------------------------------------------------------

double u0(double x, double y, const Params& p)
{
    switch (p.cas) {
        case 1: {
            // Gaussienne centree en (0,0), sigma^2 = 0.0075
            const double sigma2 = 0.0075;
            return std::exp(-(x*x + y*y) / sigma2);
        }
        case 2:
        case 3: {
            // Cylindre : u=1 si ||(x,y)||_2 <= 0.4, 0 sinon
            return (x*x + y*y <= 0.4*0.4) ? 1.0 : 0.0;
        }
        case 4: {
            // Sinus manufacture : u0 = sin(pi*x)*sin(pi*y)
            // Solution exacte : u(x,y,t) = sin(pi*(x-t))*sin(pi*(y-t))
            // Domaine [-1,1]^2, v=(1,1), periodique de periode L=2
            return std::sin(M_PI * x) * std::sin(M_PI * y);
        }
        default:
            return 0.0;
    }
}

// ----------------------------------------------------------
//  Ramene val dans [lo, hi[ periodiquement
// ----------------------------------------------------------

double modulo_periodique(double val, double lo, double hi)
{
    const double L = hi - lo;
    val -= lo;
    val  = val - L * std::floor(val / L);
    return val + lo;
}

// ----------------------------------------------------------
//  Solution analytique u_exact(x, y, t)
// ----------------------------------------------------------

double u_exact(double x, double y, double t, const Params& p)
{
    switch (p.cas) {
        case 1: {
            const double x0 = modulo_periodique(x - 1.0*t, p.xmin, p.xmax);
            const double y0 = modulo_periodique(y,          p.ymin, p.ymax);
            return u0(x0, y0, p);
        }
        case 2: {
            const double x0 = modulo_periodique(x - 0.5*t, p.xmin, p.xmax);
            const double y0 = modulo_periodique(y - 0.5*t, p.ymin, p.ymax);
            return u0(x0, y0, p);
        }
        case 3: {
            // Rotation inverse de -t autour de l'origine
            const double x0 =  x*std::cos(t) + y*std::sin(t);
            const double y0 = -x*std::sin(t) + y*std::cos(t);
            return u0(x0, y0, p);
        }
        case 4: {
            // Formule exacte fermee
            return std::sin(M_PI*(x - t)) * std::sin(M_PI*(y - t));
        }
        default:
            return 0.0;
    }
}

// ----------------------------------------------------------
//  Initialisation du vecteur solution
// ----------------------------------------------------------

void initialiser_U(std::vector<double>& U, const Params& p)
{
    U.assign(p.N, 0.0);
    for (int j = 0; j < p.Ny; ++j) {
        const double y_j = p.ymin + j * p.dy;
        for (int i = 0; i < p.Nx; ++i) {
            const double x_i = p.xmin + i * p.dx;
            U[j * p.Nx + i]  = u0(x_i, y_j, p);
        }
    }
}

// ----------------------------------------------------------
//  Erreur L2 discrete
// ----------------------------------------------------------

double erreur_L2(const std::vector<double>& U, double t, const Params& p)
{
    double err = 0.0;
    for (int j = 0; j < p.Ny; ++j) {
        const double y_j = p.ymin + j * p.dy;
        for (int i = 0; i < p.Nx; ++i) {
            const double x_i  = p.xmin + i * p.dx;
            const double diff  = U[j*p.Nx + i] - u_exact(x_i, y_j, t, p);
            err += diff * diff;
        }
    }
    return std::sqrt(p.dx * p.dy * err);
}
