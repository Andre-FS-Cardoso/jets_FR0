#define  PHYSICS                        RHD
#define  DIMENSIONS                     3
#define  GEOMETRY                       CARTESIAN
#define  BODY_FORCE                     NO
#define  COOLING                        NO
#define  RECONSTRUCTION                 LINEAR
#define  TIME_STEPPING                  RK2
#define  NTRACER                        1
#define  PARTICLES                      NO
#define  USER_DEF_PARAMETERS            5

/* -- physics dependent declarations -- */

#define  EOS                            TAUB
#define  ENTROPY_SWITCH                 NO
#define  RADIATION                      NO

/* -- user-defined parameters (labels) -- */

#define  BETA                           0
#define  RHO_IN                         1
#define  RHO_OUT                        2
#define  PRESS_IN                       3
#define  PRESS_OUT                      4

/* [Beg] user-defined constants (do not change this line) */

#define  SHOCK_FLATTENING               MULTID
#define  CHAR_LIMITING                  YES
#define  LIMITER                        MC_LIM
#define  EPS_PSHOCK_FLATTEN             2.5
#define  ARTIFICIAL_VISC                0.1
#define  UNIT_DENSITY                   CONST_mp
#define  UNIT_LENGTH                    0.5*CONST_pc
#define  UNIT_VELOCITY                  CONST_c

/* [End] user-defined constants (do not change this line) */
