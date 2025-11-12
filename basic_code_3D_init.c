/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Relativistic jet propagation.

  This problem sets initial and boundary conditions for a jet
  propagating into a uniform medium with constant density and pressure.
  The computation can be carried using \c CYLINDRICAL or \c SPHERICAL
  coordinates.
  In the first case, the jet enters at the lower z-boundary with speed
  \f$ v_z = \beta\f$ while in the second case, a conical beam with a
  small aperture (\f$\theta = 5^\circ\f$) is injected with the same
  speed from the lower radial boundary.
  At the border of the nozzle, jet values are smoothly joined with
  ambient values using the Profile() function
  (in the current setting the profile is a sharp transition).
  The jet is pressure-matched so that the beam and ambient pressure
  coincide.
 
  The configuration is defined in terms of the following parameters:

  -# <tt>g_inputParam[BETA]</tt>:     the jet velocity;
  -# <tt>g_inputParam[RHO_IN]</tt>:   the jet density;
  -# <tt>g_inputParam[RHO_OUT]</tt>:  the ambient density.
  -# <tt>g_inputParam[PRESS_IN]</tt>: the jet pressure (also equal to ambient
                                     pressure)

  defined in \c pluto.ini.
  The \c TAUB equation of state is used.

  - Configurations #01 and #02 use \c CYLINDRICAL coordinates;
  - Configuration #03 employs \c SPHERICAL coordinates (see snapshot below)

  \image html rhd_jet.03.jpg "Density (log) for configuration #03 at t=200"

  \author A. Mignone (mignone@to.infn.it)
  \date   Feb 25, 2019
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

static double Profile(double x1, double x2, double x3, int nv);
static void GetJetValues (double *vjet, double x1, double x2, double x3);

/* ********************************************************************* */
void Init (double *v, double x1, double x2, double x3)
{
	
  const double eps = 1e-20;
  double R, r;
  R = sqrt(x1*x1 + x2*x2 + x3*x3);
  r = sqrt(x1*x1 + x2*x2);
  R = fmax(R, eps); r = fmax(r, eps);
  
  double lor;
  lor = 1 / sqrt(1 - pow(g_inputParam[BETA], 2.0));

  double z0, R0;
  z0 = 1.0;
  R0 = 1.0;
  
  if (atan2(r / x3) <= 0.2) {
    
    v[VX1] = sqrt(1 - 1 / pow(lor, 2.0)) * (x1 / R);
    v[VX2] = sqrt(1 - 1 / pow(lor, 2.0)) * (x2 / R);
    v[VX3] = sqrt(1 - 1 / pow(lor, 2.0)) * (x3 / R);

    v[RHO] = g_inputParam[RHO_IN] * pow(R / R0, -2.0);
    v[PRS] = g_inputParam[PRESS_IN] * pow(R / R0, -2.0*5.0/3.0);
  } else {
    v[VX1] = 0.0;
    v[VX2] = 0.0;
    v[VX3] = 0.0;

    v[RHO] = g_inputParam[RHO_OUT] * pow(x3 / z0, -0.5);
    v[PRS] = g_inputParam[PRESS_OUT] * pow(x3 / z0, -0.5);
  }

  #if NTRACER > 0
  v[TRC] = 0.0;
  #endif

}

/* ********************************************************************* */
void InitDomain (Data *d, Grid *grid)
/*! 
 * Assign initial condition by looping over the computational domain.
 * Called after the usual Init() function to assign initial conditions
 * on primitive variables.
 * Value assigned here will overwrite those prescribed during Init().
 *
 *
 *********************************************************************** */
{
}

/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/* 
 *
 *
 *********************************************************************** */
{

}
/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/*
 *
 *
 *********************************************************************** */
{
  #if GEOMETRY == CYLINDRICAL || GEOMETRY == CARTESIAN
  if (side == X3_BEG) {    int   nv, i, j, k;
    double vjet[NVAR], vout[NVAR];    X3_BEG_LOOP(k, j, i) {      double x1 = grid->x[IDIR][i];
      double x2 = grid->x[JDIR][j];
      double x3 = grid->x[KDIR][k];      NVAR_LOOP(nv) vout[nv] = d->Vc[nv][2*KBEG - k - 1][j][i];      double r  = sqrt(x1*x1 + x2*x2);
      double th = atan2(r, fmax(x3, 1e-20));
      const double th_jet = 0.2;      if (th <= th_jet) {
        GetJetValues(vjet, x1, x2, x3);
        NVAR_LOOP(nv){
          double prof = Profile(x1,x2,x3,nv);
          d->Vc[nv][k][j][i] = vout[nv] + (vjet[nv] - vout[nv]) * prof;
        }
      } else {
        NVAR_LOOP(nv) d->Vc[nv][k][j][i] = vout[nv];
      }
    }
  }
  #endif
}

/* ********************************************************************* */
void GetJetValues (double *vjet, double x1, double x2, double x3)
{

  double R, r;
  R = sqrt(x1*x1 + x2*x2 + x3*x3);
  r = sqrt(x1*x1 + x2*x2); 
  
  double lor;
  lor = 1 / sqrt(1 - pow(g_inputParam[BETA], 2.0));

  double z0, R0;
  z0 = 1.0;
  R0 = 1.0;

  vjet[VX1] = sqrt(1 - 1 / pow(lor, 2.0)) * (x1 / R);
  vjet[VX2] = sqrt(1 - 1 / pow(lor, 2.0)) * (x2 / R);
  vjet[VX3] = sqrt(1 - 1 / pow(lor, 2.0)) * (x3 / R);

  vjet[RHO] = g_inputParam[RHO_IN] * pow(R / R0, -2.0);
  vjet[PRS] = g_inputParam[PRESS_IN] * pow(R / R0, -2.0*5.0/3.0);

  #if NTRACER > 0
  vjet[TRC] = 1.0;
  #endif
  
}
 
/* ********************************************************************* */
double Profile(double x1, double x2, double x3, int nv)
/* 
 *
 *
 *********************************************************************** */
{
  const double eps = 1e-20;
  double r  = sqrt(x1*x1 + x2*x2);
  double th = atan2(r, fmax(x3, eps));  double theta_q, aq;  if (nv == RHO) {
    theta_q = 0.29;  aq = 10.0;
  } else if (nv == PRS) {
    theta_q = 0.29;  aq = 10.0;
  } else if (nv == VX1 || nv == VX2 || nv == VX3) {
    theta_q = 0.16;  aq =  8.0;
  } else {
    return 1.0;
  }  double arg = pow( th / fmax(theta_q, 1e-12), aq );
  return 1.0 / cosh(arg);
}
