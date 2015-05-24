#pragma OPENCL EXTENSION cl_khr_fp64: enable
#define SIPC_pi 3.14159265358979323846
#define SURFMAXN 4
#define LT 0.04

// materialProperties.h
#define MATBYTES 128

typedef struct scatProps {
    char   matname[MATBYTES] ;
    float  corlen      ;
    float  roughness   ;
    float  resistivity ;
    float  specular    ;
    float  diffuse     ;
    float  shinyness   ;
} scatProps ;

#define NMATERIALS 9
static constant scatProps materialProperties[NMATERIALS] = {
//   Name          corrLen      Roughness   Resistivity Specular    Diffuse     Shinyness
    {"MATERIAL",    100.0,      0.0,        0.0,        1.0,        0.0,        50.0        },
    {"ASPHALT",     0.5,        0.005,      1.0e18,     0.8,        0.2,        30.0        },
    {"BRICK",       0.1,        0.001,      1.0e18,     0.7,        0.3,        20.0        },
    {"CONCRETE",    0.2,        0.01,       120.0,      0.3,        0.7,        10.0        },
    {"METAL",       100.0,      0.0,        1.0e-8,     1.0,        0.0,        50.0        },
    {"ROOFING",     0.1,        0.1,        1.0e18,     0.6,        0.4,        40.0        },
    {"VEGETATION",  0.01,       0.1,        2000.0,     0.2,        0.8,        5.0         },
    {"WATER",       0.01,       0.1,        2.0e1,      1.0,        0.0,        50.0        },
    {"WOOD",        0.1,        0.001,      1.0e14,     0.6,        0.4,        10.0        }
} ;
// end of materialProperties.h


typedef struct { double x, y, z; } SPVector ;

typedef struct Triangle{
    int        id ;
    SPVector   AA ;
    SPVector   BB ;
    SPVector   CC ;
    SPVector   NN ;
    double     area ;
    double     globalToLocalMat[9];
    double     localToGlobalMat[9];
    int        matId;
} Triangle;

typedef struct Ray {
    SPVector org ;               // Origin of ray
    SPVector dir ;               // unit vector direction of ray
    double   pow ;               // power of ray at origin
    double   len ;               // length or distance travelled by ray to get to origin
    SPVector pol ;               // unit vector of direction of E field of ray
} Ray ;

typedef struct {
    float r;
    float i;
} SPCmplx;
/// polar structure
///
typedef struct {
    float pabs;
    float parg;
} SPCmplxPol;

typedef struct HitPoint {
    SPVector hit;       // Location of hitpoint in x,y,z
    int tri;            // index of triangle that this hit is on
} HitPoint ;

/// Make a cartesian complex number
///
#define CMPLX_F_MAKE(a,b,out){out.r = (float)a; out.i = (float)b;}
/// Add two complex numbers together
///
#define CMPLX_ADD(a,b,out) {out.r = a.r + b.r; out .i = a.i + b.i;}
/// Subtract two complex numbers
///
#define CMPLX_SUB(a,b,out) {out.r = a.r - b.r; out .i = a.i - b.i;}
/// Multiply two complex numbers together
///
#define CMPLX_MULT(a,b,out) {out.r = a.r * b.r - a.i * b.i; out.i = a.i * b.r + a.r * b.i;}
/// Complex scalar multiply. multiply complex number cmp by scalar x
///
#define CMPLX_SCMULT(x,cmp,out) {out.r = x * cmp.r; out.i = x * cmp.i;}
/// Perform complex division
///
#define CMPLX_DIV(a,b,out) {SPCmplxPol tmp1, tmp2, tmp3; CMPLX_CART2POL(a,tmp1); CMPLX_CART2POL(b,tmp2); CMPLX_PDIV(tmp1,tmp2,tmp3); CMPLX_POL2CART(tmp3,out);}
/// convert cartesian complex number to polar complex number
///
#define CMPLX_CART2POL(a,out) {out.pabs = sqrt (a.r * a.r + a.i * a.i); out.parg = atan2 (a.i,a.r);}
/// perform complex division of two polar complex numbers
///
#define CMPLX_PDIV(a,b,out) {out.pabs = a.pabs/b.pabs; out.parg = a.parg - b.parg;}
/// Convert polar complex number to cartesian complex number
///
#define CMPLX_POL2CART(a,out) {out.r = a.pabs * cos(a.parg); out.i = a.pabs * sin(a.parg);}

/// Find the dot product of two vectors
///
#define VECT_DOT(aVect,bVect) (aVect.x*bVect.x + aVect.y*bVect.y + aVect.z*bVect.z)
/// Create a vector
///
#define VECT_CREATE(a,b,c,outVect) {    \
    outVect.x = a;                      \
    outVect.y = b;                      \
    outVect.z = c;                      \
}
/// Find the equivalent unit vector v
///
#define VECT_NORM(aVect,outVect) {                              \
    double vect_unit_tmp = 1.0 / VECT_MOD(aVect);               \
    outVect.x = aVect.x*vect_unit_tmp;                          \
    outVect.y = aVect.y*vect_unit_tmp;                          \
    outVect.z = aVect.z*vect_unit_tmp;                          \
}

/// find a cross product of two vectors
///
#define VECT_CROSS(aVect,bVect,outVect) {                       \
    SPVector ___tmp ;                                           \
    ___tmp.x = aVect.y*bVect.z - aVect.z*bVect.y;               \
    ___tmp.y = aVect.z*bVect.x - aVect.x*bVect.z;               \
    ___tmp.z = aVect.x*bVect.y - aVect.y*bVect.x;               \
    outVect = ___tmp ;                                          \
}
/// Multiply a vector by a constant
///
#define VECT_SCMULT(inVect,inScal,outVect)  {   \
    outVect.x = (inVect.x)*inScal;              \
    outVect.y = (inVect.y)*inScal;              \
    outVect.z = (inVect.z)*inScal;              \
}
/// Find the modulus (or magnitdue) of a vector
///
#define VECT_MOD(aVect) (                                       \
    sqrt (aVect.x*aVect.x + aVect.y*aVect.y + aVect.z*aVect.z)  \
)
#define VECT_MAG(aVect)(                                        \
    sqrt (aVect.x*aVect.x + aVect.y*aVect.y + aVect.z*aVect.z)  \
)
/// Subtract two vectors
///
#define VECT_SUB(aVect,bVect,outVect) {                 \
    outVect.x = aVect.x - bVect.x;                      \
    outVect.y = aVect.y - bVect.y;                      \
    outVect.z = aVect.z - bVect.z;                      \
}

// Forward declarations here
//
int factorial(int n) ;
SPCmplx G_func4(double gamma);
SPCmplx G_func3(double gamma);
SPCmplx G_func2(double gamma);
SPCmplx G_func1(double gamma);
SPCmplx G_func0(double gamma);
void matmul(double *A,double *B, double *O, int Ax, int Ay,int Bx, int By);

static constant SPVector zz_hat = {0.0, 0.0, 1.0};
static constant double Z0 = 376.99111843077516; // Impedence of free space = 120 * PI

__kernel void POField(__global Triangle * tris, // input array of triangles
                      __global Ray *rays,       // input array of rays
                      int nrays,                // input number of incident rays to process
                      __global HitPoint *hitpoints,  // input array of hit points for each ray
                      SPVector RxPnt,           // input receiver location to calc field at
                      double k,                 // input k=2*PI/lambda
                      SPVector Vdir,            // input unit vector defining V pol at receiver
                      SPVector Hdir,            // input unit vector defining H pol at receiver
                      __global SPCmplx *EsVs,   // output array of scattered field strengths for V pol
                      __global SPCmplx *EsHs)   // output array of scattered field strengths for H pol
{

    int ind,n ;
    ind = get_global_id(0) ;

    if (ind >=0 && ind < nrays ) {

//        printf("[%d] Vdir : %f,%f,%f\n",ind,Vdir.x,Vdir.y,Vdir.z);
//        printf("[%d] Hdir : %f,%f,%f\n",ind,Hdir.x,Hdir.y,Hdir.z);
        
        // Some notes on variable terminology used in this function
        //
        // _il = incident ray in local coords
        // _ig = incident ray in global coords
        // _sl = scattered ray in local coords
        // _sg = scattered ray in global coords
        // uvw are direction cosines
        // E is Electric field
        // r is the magnitude of the vector defining the observation point in global coordinates (not the range)
        // J = surface current. J_l[3],J_g[3] are vectors in local and global coords
        // _par refers to parallel or V polarisation. _per refers to perpendicular polarisation
        // Jc_par, Jc_per are the complex surface current components in the parallel and perpendicular polarisation directions
        // E_par and E_per are the parallel and perpendicular E fields at the observation point in V/m
        //
        double Dp,Dq,D0;
        double Cp,Cq,C0;
        double u,v,w, x1,y1,z1,x2,y2,z2,x3,y3,z3,A;
        double r, Lt, magDp, magDq;
        double uvw_ig[3], uvw_il[3], sin_theta_il, cos_theta_il, tan_phi_il, sin_phi_il, cos_phi_il,E_ig[3], E_il[3] ;
        double Dq_minus_Dp, Eiphi_l, Eitheta_l, Rs, GamParr, GamPerp ;
        double J_l[3], J_g[3],J_par, J_per;
        double phs_ig, r_ig;
        double globalToLocalMat[9];
        double localToGlobalMat[9];
        
        for(int i=0; i<9; i++)globalToLocalMat[i] = tris[hitpoints[ind].tri].globalToLocalMat[i];
        for(int i=0; i<9; i++)localToGlobalMat[i] = tris[hitpoints[ind].tri].localToGlobalMat[i];
        
        int matId;
        
        SPVector uvw_sg, obsDir, rayDist, Raydir_l ;
        SPVector Eil, Eig, theta_l_hat, phi_l_hat, Jg;
        
        SPCmplx e_jDp, e_jDq, e_jD0, jDp, jD0, jDq, jDq_pow,jDp_min_jDq,jDp_a[SURFMAXN], jDq_a[SURFMAXN] ;
        SPCmplx C0_over_nplus1, jD_subs_pow, leftPart, G, Ic  ;
        SPCmplx Jc_par, Jc_per, Epar,Eper ;
        SPCmplx jkZ0_o_4PIr, e_jkr;
        SPCmplx sum, tmp, tmp1,t, t1, t2, t3,t4,t5, term1,term2,term3, braces ;
        
        // r is the magnitude of the vector defining the observation point in global coordinates (not the range)
        //
        VECT_SUB(RxPnt, hitpoints[ind].hit, obsDir);
        r = VECT_MAG(obsDir);
        VECT_NORM(obsDir, obsDir);
//        r  = VECT_MAG(RxPnt);
        VECT_SUB(hitpoints[ind].hit, rays[ind].org, rayDist);
        r_ig   = VECT_MAG(rayDist) + rays[ind].len;
        phs_ig = -(k * r_ig) - (SIPC_pi/2.0) ;
//        double p = phs_ig;
//        p = p / (2 * SIPC_pi) ;
//        int ip = (int)p;
//        p = (p - ip) * 2 * SIPC_pi;
//        printf("[%d] phs_ig : %f\n",ind,(180/SIPC_pi)*p);
        
        // Sort out the direction cosines for this ray, hitpoint and observation point
        // Assuming that the direction of the ray has been normalised then
        // the direction cosine is just the component of direction
        //
        uvw_sg.x = obsDir.x -rays[ind].dir.x ;
        uvw_sg.y = obsDir.y -rays[ind].dir.y ;
        uvw_sg.z = obsDir.z -rays[ind].dir.z ;
        
        // Calculate surface integral using Ludwig's integration algorithm [1], (modified for
        // triangular subregions by Pogorzelski [2] and then modified again for barycentric (simplex)
        // coordinates by Dos Santos [3].
        //
        u  = uvw_sg.x ;
        v  = uvw_sg.y ;
        w  = uvw_sg.z ;
        x1 = tris[hitpoints[ind].tri].AA.x;
        y1 = tris[hitpoints[ind].tri].AA.y ;
        z1 = tris[hitpoints[ind].tri].AA.z ;
        x2 = tris[hitpoints[ind].tri].BB.x ;
        y2 = tris[hitpoints[ind].tri].BB.y ;
        z2 = tris[hitpoints[ind].tri].BB.z ;
        x3 = tris[hitpoints[ind].tri].CC.x ;
        y3 = tris[hitpoints[ind].tri].CC.y ;
        z3 = tris[hitpoints[ind].tri].CC.z ;
        A  = tris[hitpoints[ind].tri].area ;
        
        Cp = Cq = 0;
        C0 = 1;
        Dp = k * ( ((x1 - x3) * u) + ((y1 - y3) * v) + ((z1 - z3) * w) ) ;
        Dq = k * ( ((x2 - x3) * u) + ((y2 - y3) * v) + ((z2 - z3) * w) ) ;
        D0 = k * ( (x3 * u) + (y3 * v) + (z3 * w) ) ;
        
        Lt      = LT ;          // Length of Taylor series region
        magDp   = fabs(Dp) ;
        magDq   = fabs(Dq) ;
        e_jD0.r = cosf(D0) ; e_jD0.i = sinf(D0) ;
        e_jDp.r = cosf(Dp) ; e_jDp.i = sinf(Dp) ;
        e_jDq.r = cosf(Dq) ; e_jDq.i = sinf(Dq) ;
        jDp.r   = 0;           jDp.i = Dp ;
        jD0.r   = 0;           jD0.i = D0 ;
        jDq.r   = 0;           jDq.i = Dq ;
        
        if ( magDp < Lt && magDq >= Lt ){  // Case 1
            sum.r = sum.i = 0;
            
            // Set up zeroeth term of sum
            //
            jD_subs_pow.r = 1.0;
            jD_subs_pow.i = 0.0 ;
            t1 = G_func0(-Dq) ;
            CMPLX_MULT(e_jDq, t1, t) ;
            C0_over_nplus1.r = -C0 ;
            C0_over_nplus1.i = 0 ;
            CMPLX_ADD(C0_over_nplus1, t, sum);
            
            // Now sum over n
            //
            // Loop unroll assume SURFMAXN=4, so that we dont have to keep testing
            //
            /*for (int n=1; n<=SURFMAXN; n++) {
                CMPLX_MULT(jD_subs_pow, jDp, t);
                jD_subs_pow = t ;
                CMPLX_SCMULT(1.0/(double)factorial(n), t, leftPart) ;
                t1 = G_func(n, -Dq) ;
                CMPLX_MULT(e_jDq, t1, t) ;
                C0_over_nplus1.r = -C0 / (n+1) ;
                C0_over_nplus1.i = 0.0 ;
                CMPLX_ADD(C0_over_nplus1, t, t1) ;
                CMPLX_MULT(leftPart, t1, t) ;
                CMPLX_ADD(sum, t, sum) ;
            }*/
            
            n=1;
            CMPLX_MULT(jD_subs_pow, jDp, t);
            jD_subs_pow = t ;
            CMPLX_SCMULT(1.0/(double)factorial(n), t, leftPart) ;
            t1 = G_func1(-Dq) ;
            CMPLX_MULT(e_jDq, t1, t) ;
            C0_over_nplus1.r = -C0 / (n+1) ;
            C0_over_nplus1.i = 0.0 ;
            CMPLX_ADD(C0_over_nplus1, t, t1) ;
            CMPLX_MULT(leftPart, t1, t) ;
            CMPLX_ADD(sum, t, sum) ;
            n=2;
            CMPLX_MULT(jD_subs_pow, jDp, t);
            jD_subs_pow = t ;
            CMPLX_SCMULT(1.0/(double)factorial(n), t, leftPart) ;
            t1 = G_func2(-Dq) ;
            CMPLX_MULT(e_jDq, t1, t) ;
            C0_over_nplus1.r = -C0 / (n+1) ;
            C0_over_nplus1.i = 0.0 ;
            CMPLX_ADD(C0_over_nplus1, t, t1) ;
            CMPLX_MULT(leftPart, t1, t) ;
            CMPLX_ADD(sum, t, sum) ;
            n=3;
            CMPLX_MULT(jD_subs_pow, jDp, t);
            jD_subs_pow = t ;
            CMPLX_SCMULT(1.0/(double)factorial(n), t, leftPart) ;
            t1 = G_func3(-Dq) ;
            CMPLX_MULT(e_jDq, t1, t) ;
            C0_over_nplus1.r = -C0 / (n+1) ;
            C0_over_nplus1.i = 0.0 ;
            CMPLX_ADD(C0_over_nplus1, t, t1) ;
            CMPLX_MULT(leftPart, t1, t) ;
            CMPLX_ADD(sum, t, sum) ;
            n=4;
            CMPLX_MULT(jD_subs_pow, jDp, t);
            jD_subs_pow = t ;
            CMPLX_SCMULT(1.0/(double)factorial(n), t, leftPart) ;
            t1 = G_func4(-Dq) ;
            CMPLX_MULT(e_jDq, t1, t) ;
            C0_over_nplus1.r = -C0 / (n+1) ;
            C0_over_nplus1.i = 0.0 ;
            CMPLX_ADD(C0_over_nplus1, t, t1) ;
            CMPLX_MULT(leftPart, t1, t) ;
            CMPLX_ADD(sum, t, sum) ;
            
            // Finished summation
            //
            CMPLX_SCMULT(2*A, e_jD0, t);
            CMPLX_DIV(t, jDq, t1);
            CMPLX_MULT(t1, sum, Ic);
            
        }else if ( magDp < Lt && magDq < Lt ){  // Case 2
            sum.r = sum.i = 0;
            jDp_a[0].r = 1 ; jDp_a[0].i = 0 ;
            jDq_a[0].r = 1 ; jDq_a[0].i = 0 ;
            
            for (int i=1; i<SURFMAXN; i++) {
                CMPLX_MULT(jDp_a[i-1], jDp, jDp_a[i]);
                CMPLX_MULT(jDq_a[i-1], jDq, jDq_a[i]);
            }
            for(int n=0; n<SURFMAXN; n++){
                for(int m=0; m<SURFMAXN; m++){
                    CMPLX_MULT(jDp_a[n],jDq_a[m],tmp);
                    CMPLX_SCMULT(1.f/((double)(factorial(m+n+2))), tmp, tmp1);
                    CMPLX_ADD(tmp1, sum, sum);
                }
            }
            
            // Finished summation
            //
            CMPLX_SCMULT(2*A, e_jD0, tmp);
            CMPLX_MULT(tmp, sum, Ic);
            
        } else if ( magDp >= Lt && magDq < Lt ){    // Case 3
            sum.r = sum.i = 0;
            jDq_pow.r = 1.0; jDq_pow.i = 0.0 ;
            G = G_func1(-Dp) ;
            sum = G ;
            
            // Loop unroll, assume SURFMAXN=4 so we dont have to keep testing
            //
            /*for(int n=1; n<SURFMAXN; n++){
                G = G_func4(n+1, -Dp) ;
                CMPLX_SCMULT(1.0/(double)(n+1), G, t);
                CMPLX_MULT(jDq_pow, jDq, t1);
                jDq_pow = t1 ;
                CMPLX_SCMULT(1/(double)factorial(n), jDq_pow, t1);
                CMPLX_MULT(t1, t, t2);
                CMPLX_ADD(sum, t2, sum);
            }*/
            n=1;
            G = G_func2(-Dp) ;
            CMPLX_SCMULT(1.0/(double)(n+1), G, t);
            CMPLX_MULT(jDq_pow, jDq, t1);
            jDq_pow = t1 ;
            CMPLX_SCMULT(1/(double)factorial(n), jDq_pow, t1);
            CMPLX_MULT(t1, t, t2);
            CMPLX_ADD(sum, t2, sum);
            n=2;
            G = G_func3(-Dp) ;
            CMPLX_SCMULT(1.0/(double)(n+1), G, t);
            CMPLX_MULT(jDq_pow, jDq, t1);
            jDq_pow = t1 ;
            CMPLX_SCMULT(1/(double)factorial(n), jDq_pow, t1);
            CMPLX_MULT(t1, t, t2);
            CMPLX_ADD(sum, t2, sum);
            n=3;
            G = G_func4(-Dp) ;
            CMPLX_SCMULT(1.0/(double)(n+1), G, t);
            CMPLX_MULT(jDq_pow, jDq, t1);
            jDq_pow = t1 ;
            CMPLX_SCMULT(1/(double)factorial(n), jDq_pow, t1);
            CMPLX_MULT(t1, t, t2);
            CMPLX_ADD(sum, t2, sum);
            
            // Finished summation
            //
            CMPLX_SCMULT(2*A, e_jD0, t);
            CMPLX_MULT(t, e_jDp, t1);
            CMPLX_MULT(t1, sum, Ic);
            
        } else if ( magDp >= Lt && magDq >= Lt && fabs(Dp - Dq) < Lt ){   // Case 4
            sum.r = sum.i = 0;
            // Set up zeroeth term of sum
            //
            CMPLX_SUB(jDp, jDq, jDp_min_jDq);
            jD_subs_pow.r = 1.0; jD_subs_pow.i = 0.0 ;
            // C0 = 1
            //
            t1 = G_func0(Dq) ;
            CMPLX_SCMULT(-1, t1, t1);
            CMPLX_ADD(t1, e_jDp, sum);
            
            // Now sum over n
            //
            // Loop unroll, assume SURFMAXN=4 so we dont have to keep testing
            //
            /*for (int n=1; n<=SURFMAXN; n++) {
                CMPLX_MULT(jD_subs_pow, jDp_min_jDq, t);
                jD_subs_pow = t ;
                CMPLX_SCMULT(1.0/(double)factorial(n), t, leftPart) ;
                t1 = G_func4(n, Dq) ;
                CMPLX_SCMULT(-1, t1, t2);
                CMPLX_SCMULT(1.0/(n+1.0), e_jDp, t3);
                CMPLX_ADD(t2, t3, t4);
                CMPLX_MULT(leftPart, t4, t5);
                CMPLX_ADD(t5, sum, sum) ;
            }*/
            n=1;
            CMPLX_MULT(jD_subs_pow, jDp_min_jDq, t);
            jD_subs_pow = t ;
            CMPLX_SCMULT(1.0/(double)factorial(n), t, leftPart) ;
            t1 = G_func1(Dq) ;
            CMPLX_SCMULT(-1, t1, t2);
            CMPLX_SCMULT(1.0/(n+1.0), e_jDp, t3);
            CMPLX_ADD(t2, t3, t4);
            CMPLX_MULT(leftPart, t4, t5);
            CMPLX_ADD(t5, sum, sum) ;
            n=2;
            CMPLX_MULT(jD_subs_pow, jDp_min_jDq, t);
            jD_subs_pow = t ;
            CMPLX_SCMULT(1.0/(double)factorial(n), t, leftPart) ;
            t1 = G_func2(Dq) ;
            CMPLX_SCMULT(-1, t1, t2);
            CMPLX_SCMULT(1.0/(n+1.0), e_jDp, t3);
            CMPLX_ADD(t2, t3, t4);
            CMPLX_MULT(leftPart, t4, t5);
            CMPLX_ADD(t5, sum, sum) ;
            n=3;
            CMPLX_MULT(jD_subs_pow, jDp_min_jDq, t);
            jD_subs_pow = t ;
            CMPLX_SCMULT(1.0/(double)factorial(n), t, leftPart) ;
            t1 = G_func3(Dq) ;
            CMPLX_SCMULT(-1, t1, t2);
            CMPLX_SCMULT(1.0/(n+1.0), e_jDp, t3);
            CMPLX_ADD(t2, t3, t4);
            CMPLX_MULT(leftPart, t4, t5);
            CMPLX_ADD(t5, sum, sum) ;
            n=4;
            CMPLX_MULT(jD_subs_pow, jDp_min_jDq, t);
            jD_subs_pow = t ;
            CMPLX_SCMULT(1.0/(double)factorial(n), t, leftPart) ;
            t1 = G_func4(Dq) ;
            CMPLX_SCMULT(-1, t1, t2);
            CMPLX_SCMULT(1.0/(n+1.0), e_jDp, t3);
            CMPLX_ADD(t2, t3, t4);
            CMPLX_MULT(leftPart, t4, t5);
            CMPLX_ADD(t5, sum, sum) ;
            // End of summation
            CMPLX_SCMULT(2*A, e_jD0, t);
            CMPLX_DIV(t, jDq, t1);
            CMPLX_MULT(t1, sum, Ic);
            
        } else {
            Dq_minus_Dp = Dq-Dp ;
            CMPLX_SCMULT((C0 / (Dp * Dq_minus_Dp)), e_jDp, term1);
            CMPLX_SCMULT((C0 / (Dq * Dq_minus_Dp)), e_jDq, term2);
            term3.r = C0 / (Dp * Dq) ; term3.i = 0;
            CMPLX_SUB(term1, term2, braces);
            CMPLX_SUB(braces, term3, braces) ;
            CMPLX_SCMULT(2*A, e_jD0, t) ;
            CMPLX_MULT(t, braces, Ic);
        }
        
        // At this point we have integrated the complex field over the triangular facet and have the
        // value stored in the complex number Ic
        //
        VECT_SCMULT(rays[ind].pol, sqrt(rays[ind].pow / (4 * SIPC_pi * r_ig * r_ig)), Eig) ;
//        printf("[%d] : pol : %f,%f,%f\n",ind,rays[ind].pol.x,rays[ind].pol.y,rays[ind].pol.z);
//        printf("[po %d] r_ig: %f input E: %e E corrected for r_ig: %e r to rx:%e total r: %e\n",
//               ind,r_ig,rays[ind].pow,sqrt(rays[ind].pow / (4 * SIPC_pi * r_ig * r_ig)),r,r+r_ig);
        
        // Assuming that the direction of the ray has been normalised then
        // the direction cosine is just the component of direction
        //
        uvw_ig[0] = -rays[ind].dir.x;
        uvw_ig[1] = -rays[ind].dir.y;
        uvw_ig[2] = -rays[ind].dir.z;
        
        matmul(globalToLocalMat, uvw_ig, uvw_il, 3, 3, 1, 3);
        
        sin_theta_il = sqrt(uvw_il[0]*uvw_il[0] +  uvw_il[1] * uvw_il[1] ) ;
        cos_theta_il = sqrt(1 - sin_theta_il*sin_theta_il) ;
        if(sin_theta_il <= 0.0001){
            tan_phi_il = 0.0 ;
            sin_phi_il = 1.0 ;
            cos_phi_il = 0.0 ;
        }else{
            tan_phi_il   = atan2(uvw_il[1] , uvw_il[0]) ;
            sin_phi_il   = uvw_il[1] / sin_theta_il ;
            cos_phi_il   = uvw_il[0] / sin_theta_il ;
        }
        
        // Also convert the direction of the E field into the coordinate system
        // of the triangle
        //
        VECT_CREATE(uvw_il[0], uvw_il[1], uvw_il[2], Raydir_l) ;
        E_ig[0] = Eig.x ;
        E_ig[1] = Eig.y ;
        E_ig[2] = Eig.z ;
        matmul(globalToLocalMat, E_ig, E_il, 3, 3, 1, 3) ;
        VECT_CREATE(E_il[0], E_il[1], E_il[2], Eil) ;
        if(fabs(VECT_DOT(zz_hat, Raydir_l)) >=0.99){
            VECT_CREATE(1, 0, 0, phi_l_hat) ;
        }else{
            VECT_CROSS (zz_hat, Raydir_l, phi_l_hat) ;
        }
        VECT_CROSS(Raydir_l, phi_l_hat, theta_l_hat) ;
        VECT_NORM(phi_l_hat, phi_l_hat) ;
        VECT_NORM(theta_l_hat, theta_l_hat) ;
        Eiphi_l   = VECT_DOT(Eil, phi_l_hat) ;
        Eitheta_l = VECT_DOT(Eil, theta_l_hat) ;
        
        // Calculate Gamma_parallel and Gamma_perpendicular
        //
        matId   = tris[hitpoints[ind].tri].matId ;
        Rs      = materialProperties[matId].resistivity ;
        GamParr = -1.0 * Z0 * cos_theta_il / (2*Rs + Z0*cos_theta_il) ;
        GamPerp = -1.0 * Z0 / ( 2.0*Rs*cos_theta_il + Z0);
        
        // Now calculate Jx_local and Jy_local and convert them to global
        // vectors.
        //
        J_l[0] = ((-1.0 * Eitheta_l * cos_phi_il * GamParr / Z0) + (Eiphi_l * sin_phi_il * GamPerp / Z0)) * cos_theta_il ;
        J_l[1] = ((-1.0 * Eitheta_l * sin_phi_il * GamParr / Z0) - (Eiphi_l * cos_phi_il * GamPerp / Z0)) * cos_theta_il ;
        J_l[2] = 0.0  ;
        matmul(localToGlobalMat, J_l, J_g, 3, 3, 1, 3) ;
        VECT_CREATE(J_g[0], J_g[1], J_g[2], Jg) ;
        
        // Find the component of the surface current in the parallel polarisation
        // direction
        //
        J_par = VECT_DOT(Vdir, Jg);
        J_per = VECT_DOT(Hdir, Jg);
        CMPLX_F_MAKE(J_par*cos(phs_ig), J_par*sin(phs_ig), Jc_par);
        CMPLX_F_MAKE(J_per*cos(phs_ig), J_per*sin(phs_ig), Jc_per);
        
        // Work out scaler (complex) component of E field
        //
        CMPLX_F_MAKE(0, -k*Z0/(4*SIPC_pi*r), jkZ0_o_4PIr) ;
        e_jkr.r = cos(-k*r);
        e_jkr.i = sin(-k*r);
        CMPLX_MULT(jkZ0_o_4PIr, e_jkr, tmp1);
        CMPLX_MULT(tmp1, Ic, tmp);
        CMPLX_MULT(Jc_par, tmp, Epar);
        CMPLX_MULT(Jc_per, tmp, Eper);
        
        EsVs[ind].r = Epar.r ; EsVs[ind].i = Epar.i ;
        EsHs[ind].r = Eper.r ; EsHs[ind].i = Eper.i ;
        
    }
    return ;
}

int factorial(int n)
{
    int c;
    int result = 1;
    
    for (c = 1; c <= n; c++)
        result = result * c;
    
    return result;
}

SPCmplx G_func0(double gamma){
    SPCmplx e_jgamma, jgamma, one;
    SPCmplx t1, ans ;
    e_jgamma.r = cos(gamma) ;
    e_jgamma.i = sin(gamma) ;
    jgamma.r   = 0 ;
    jgamma.i   = gamma ;
    one.r      = 1;
    one.i      = 0;
    CMPLX_SUB(e_jgamma, one, t1) ;
    CMPLX_DIV(t1, jgamma, ans) ;
    return ans ;
}

SPCmplx G_func1(double gamma){
    SPCmplx e_jgamma, jgamma, G, t1, ans;
    e_jgamma.r = cos(gamma) ;
    e_jgamma.i = sin(gamma) ;
    jgamma.r   = 0 ;
    jgamma.i   = gamma ;
    G = G_func0(gamma) ;
    CMPLX_SUB(e_jgamma, G, t1) ;
    CMPLX_DIV(t1, jgamma, ans) ;
    return ans;
}

SPCmplx G_func2(double gamma){
    SPCmplx e_jgamma, jgamma, G, t1, t2, ans;
    e_jgamma.r = cos(gamma) ;
    e_jgamma.i = sin(gamma) ;
    jgamma.r   = 0 ;
    jgamma.i   = gamma ;
    
    G = G_func1(gamma) ;
    CMPLX_SCMULT(2.0, G, t1);
    CMPLX_SUB(e_jgamma, t1, t2) ;
    CMPLX_DIV(t2, jgamma, ans) ;
    return ans;
}

SPCmplx G_func3(double gamma){
    SPCmplx e_jgamma, jgamma, G, t1, t2, ans;
    e_jgamma.r = cos(gamma) ;
    e_jgamma.i = sin(gamma) ;
    jgamma.r   = 0 ;
    jgamma.i   = gamma ;
    
    G = G_func2(gamma) ;
    CMPLX_SCMULT(3.0, G, t1);
    CMPLX_SUB(e_jgamma, t1, t2) ;
    CMPLX_DIV(t2, jgamma, ans) ;
    return ans;
}

SPCmplx G_func4(double gamma){
    SPCmplx e_jgamma, jgamma, G, t1, t2, ans;
    e_jgamma.r = cos(gamma) ;
    e_jgamma.i = sin(gamma) ;
    jgamma.r   = 0 ;
    jgamma.i   = gamma ;
    
    G = G_func3(gamma) ;
    CMPLX_SCMULT(4.0, G, t1);
    CMPLX_SUB(e_jgamma, t1, t2) ;
    CMPLX_DIV(t2, jgamma, ans) ;
    return ans;
}

// Matrix multiplication for a N size array
// Note : in order to remove the test for incompatable dimensions (for speed)
// it is essential that Ax = By
// Also the output matrix *O must have been allocated to the correct dimensions
//  double *A - matrix A
//  double *B - matrix B. By must equal Ax
//  double *O - output matrix. Dimensions must be Ay x Bx
//  int Ax  - Number of columns in A
//  int Ay  - Number of rows in A
//  int Bx  - number of columns in B
//  int By  - number of rows in B
//
void matmul(double *A, double *B, double *O, int Ax, int Ay,int Bx, int By)
{
    // Ax must be equal to By !!!
    //
    //    if (Ax != By) {
    //        printf("Error : Multiplication by incompatable matrices\n");
    //        exit(1);
    //    }
    
    int i,j,k;
    
    for (i=0;i<Ay;i++){
        for(j=0;j<Bx;j++){
            
            O[i*Bx+j] = 0;
            for(k=0;k<Ax;k++)
                O[i*Bx+j] += A[i*Ax+k] * B[k*Bx+j];
        }
    }
    return ;
}

