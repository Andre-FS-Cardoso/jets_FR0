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
  double R, r;
  R = sqrt(x1*x1 + x2*x2 + x3*x3);
  r = sqrt(x1*x1 + x2*x2); 
  
  double lor;
  lor = 1 / sqrt(1 - pow(g_inputParam[BETA], 2.0));

  double vz, vr;
  vr = sqrt(1 - 1 / pow(lor, 2.0)) * (r / R);
  vz = sqrt(1 - 1 / pow(lor, 2.0)) * (x3 / R);
  
  if (atan(r / x3) <= 0.2) {
    
    v[VX1] = vr*(x1 / r);
    v[VX2] = vr*(x2 / r);
    v[VX3] = vz;

    v[RHO] = g_inputParam[RHO_IN] * pow(R / 1.0, -0.75);
    v[PRS] = g_inputParam[PRESS_IN] * pow(R / 1.0, -1.25);
  } else {
    v[VX1] = 0.0;
    v[VX2] = 0.0;
    v[VX3] = 0.0;

    v[RHO] = g_inputParam[RHO_OUT] * pow(x3 / 1.0, -0.5);
    v[PRS] = g_inputParam[PRESS_OUT] * pow(x3 / 1.0, -0.5);
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
  int   nv, i, j, k;
  double  r, z, vjet[NVAR], vout[NVAR];
  double x1, x2, x3;

  r = sqrt(x1*x1 + x2*x2); 

  #if GEOMETRY == CYLINDRICAL || GEOMETRY == CARTESIAN

  if (side == X2_BEG) {

    X2_BEG_LOOP(k, j, i) {

      x1 = grid->x[IDIR][i];
      x2 = grid->x[JDIR][j];
      x3 = grid->x[KDIR][k];

      NVAR_LOOP(nv) vout[nv] = d->Vc[nv][k][2 * JBEG - j - 1][i];

      if (atan(r / x3) <= 0.2) {
	GetJetValues(vjet, x1, x2, x3);
	NVAR_LOOP(nv) d->Vc[nv][k][j][i] = vout[nv] + (vjet[nv] - vout[nv])*Profile(x1,x2,x3,nv);
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

  double vz, vr;
  vr = sqrt(1 - 1 / pow(lor, 2.0)) * (r / R);
  vz = sqrt(1 - 1 / pow(lor, 2.0)) * (x3 / R);

  vjet[VX1] = vr*(x1 / r);
  vjet[VX2] = vr*(x2 / r);
  vjet[VX3] = vz;

  vjet[RHO] = g_inputParam[RHO_IN] * pow(R / 1.0, -0.75);
  vjet[PRS] = g_inputParam[PRESS_IN] * pow(R / 1.0, -1.25);

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
  double aq;
  double theta_q;
  
  double r;
  r = sqrt(x1*x1 + x2*x2);

  if (nv == RHO) {
    theta_q = 0.29;
    aq = 10.0;
  }
  if (nv == PRS) {
    theta_q = 0.29;
    aq = 10.0;
  }
  if (nv == VX1 || VX2) {
    theta_q = 0.16;
    aq = 8.0;
  }

  return 1.0/cosh(pow(r / (x3 * theta_q), aq));
}
