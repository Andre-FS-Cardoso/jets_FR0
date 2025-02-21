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

static double Profile(double r, int nv);
static void GetJetValues (double *vjet, double x1, double x2, double x3);

/* ********************************************************************* */
void Init (double *v, double x1, double x2, double x3)
{
  //double x1; // Em coordenadas cilíndricas, x1 é r.
  //double x2; // Em coordenadas cilíndricas, x2 é z.
  double R;
  R = sqrt(x1 * x1 + x2 * x2);

  double lor;
  lor = 1 / sqrt(1 - pow(g_inputParam[BETA], 2.0));
  if (x1 / x2 <= 0.2) {
    v[VX1] = sqrt(1 - 1 / pow(lor, 2.0)) * (x1 / R);
    v[VX2] = sqrt(1 - 1 / pow(lor, 2.0)) * (x2 / R);
    v[VX3] = 0.0;

    v[RHO] = g_inputParam[RHO_IN] * pow(R / 1.0, -2.0);
    v[PRS] = g_inputParam[PRESS_IN] * pow(R / 1.0, -2.0*5.0/3.0);
  } else {
    v[VX1] = 0.0;
    v[VX2] = 0.0;
    v[VX3] = 0.0;

    v[RHO] = g_inputParam[RHO_OUT] * pow(x2 / 1.0, -0.5);
    v[PRS] = g_inputParam[PRESS_OUT] * pow(x2 / 1.0, -0.5);
  }
  
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
void UserDefBoundary(const Data *d, RBox *box, int side, Grid *grid) 
{
  int nv, i, j, k;
  double r, z, vjet[NVAR], vout[NVAR];
  double x1, x2, x3;

  #if GEOMETRY == CYLINDRICAL

  // Tratando a fronteira inferior de Z (X2_BEG)
  if (side == X2_BEG) { // Limite inferior em Z

    

    X2_BEG_LOOP(k, j, i) {

      x1 = grid->x[IDIR][i]; // r
      x2 = grid->x[JDIR][j]; // z
      x3 = grid->x[KDIR][k];
      r = grid->x[IDIR][i];

      //R = sqrt(x1 * x1 + x2 * x2);

      //printf("PARA X2: r = %f, z = %f\n", r, z);

      if (x1 / x2 <= 0.2) { // Dentro do cone
        GetJetValues(vjet, x1, x2, x3);
        NVAR_LOOP(nv) d->Vc[nv][k][j][i] = vjet[nv];
      } else { // Fora do cone, aplica valores do ambiente
        NVAR_LOOP(nv) vout[nv] = d->Vc[nv][k][2 * JBEG - j - 1][i];
        vout[VX2] *= -1.0;
        NVAR_LOOP(nv) d->Vc[nv][k][j][i] = vout[nv] + (vjet[nv] - vout[nv])*Profile(r,nv);
      }
    }
  }

  #endif
}




/* ********************************************************************* */
void GetJetValues (double *vjet, double x1, double x2, double x3)
{
  //double x1; // Em coordenadas cilíndricas, x1 é r.
  //double x2; // Em coordenadas cilíndricas, x2 é z.
  double R;
  R = sqrt(x1 * x1 + x2 * x2);

  double lor;
  lor = 1 / sqrt(1 - pow(g_inputParam[BETA], 2.0));

  vjet[VX1] = sqrt(1 - 1 / pow(lor, 2.0)) * (x1 / R);
  vjet[VX2] = sqrt(1 - 1 / pow(lor, 2.0)) * (x2 / R);
  vjet[VX3] = 0.0;

  vjet[RHO] = g_inputParam[RHO_IN] * pow(R / 1.0, -2.0);
  vjet[PRS] = g_inputParam[PRESS_IN] * pow(R / 1.0, -2.0*5.0/3.0);

}

 
/* ********************************************************************* */
double Profile(double r, int nv)
/* 
 *
 *
 *********************************************************************** */
{
  int xn = 14;
  double r0 = 1.0;

  if (nv == RHO) r0 = 1.1;

  #if GEOMETRY == SPHERICAL
   r0 = 5.0/180.0*CONST_PI;
  #endif
  return 1.0/cosh(pow(r/r0,xn));
}
