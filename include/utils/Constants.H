#ifndef _CONSTANTS_H_
#define _CONSTANTS_H_

#include <cmath>

//Fundamental electromagnetic constants
#define eps_0 (8.8541878188e-12)  
#define mu_0 (1.25663706127e-6)  
#define q_e -1.602176634e-19 // C
#define m_e 9.1093837139e-31 // kg

// Free parameters (must be defined elsewhere)
#define B_0 1e-2 // T
#define L_0 1e-3 // m
#define n_0 1e19 // m^-3

// Derived constants
#define c_0 (1/(std::sqrt(eps_0)*std::sqrt(mu_0)))
#define H_0 (B_0/mu_0)
#define E_0 (B_0/(std::sqrt(eps_0)*std::sqrt(mu_0)))
#define t_0 (L_0*std::sqrt(eps_0)*std::sqrt(mu_0))
#define T_0 (B_0*L_0/(std::sqrt(eps_0)*std::sqrt(mu_0)))
#define J_0 (B_0/(L_0*mu_0))
#define rho_0 (B_0*std::sqrt(eps_0)/(L_0*std::sqrt(mu_0)))
#define p_0 (std::pow(B_0, 2)*std::sqrt(eps_0)/(std::sqrt(mu_0)*n_0))
#define nu_0 (1/(L_0*std::sqrt(eps_0)*std::sqrt(mu_0)))
#define m_0 (std::pow(B_0, 2)*eps_0/n_0)
#define q_0 (B_0*std::sqrt(eps_0)/(L_0*std::sqrt(mu_0)*n_0))

// Dimensionless constants so the forms of the equations are familiar
#define Q_e (q_e/q_0)
#define M_e (m_e/m_0)
#define C_L 1
#define C_L2 1
#define C_h 1
#define Tau 10.

#define DENSITY_FACTOR .01 //use for testing at low cut-off frequency
#define B_FIELD_FACTOR 1 //use for testing at low cut-off frequency

#define eps0i 1. // inverse
#define mu0i 1.  // inverse

#define CFL 0.9

#define TINY 1e-25
#define BIGGY 1e25

#define SPACE_DIM 3

#define ZERO_PADDING 5

#define PRINT_WHEN_SAVING_DATA false
#endif