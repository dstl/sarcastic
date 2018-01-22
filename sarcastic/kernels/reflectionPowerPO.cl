

#include "materialProperties.h"

/// defines the structure for a vector
///
typedef union {
    struct { double x, y, z; };
    struct { double r, g, b; };
    struct { double cell[3]; };
    struct { double R, theta, phi; };
    struct { double lat, lon, alt; };
    
} SPVector ;

/// Create a vector
///
#define VECT_CREATE(a,b,c,outVect) {    \
    outVect.x = a;                      \
    outVect.y = b;                      \
    outVect.z = c;                      \
}

/// Multiply a vector by a constant
///
#define VECT_SCMULT(inVect,inScal,outVect)  {   \
    outVect.x = (inVect.x)*inScal;              \
    outVect.y = (inVect.y)*inScal;              \
    outVect.z = (inVect.z)*inScal;              \
}

/// Subtract two vectors
///
#define VECT_SUB(aVect,bVect,outVect) {                 \
    outVect.x = aVect.x - bVect.x;                      \
    outVect.y = aVect.y - bVect.y;                      \
    outVect.z = aVect.z - bVect.z;                      \
}

/// add two vectors
///
#define VECT_ADD(aVect,bVect,outVect) {                 \
    outVect.x = aVect.x + bVect.x;                      \
    outVect.y = aVect.y + bVect.y;                      \
    outVect.z = aVect.z + bVect.z;                      \
}

/// Find the modulus (or magnitdue) of a vector
///
#define VECT_MOD(aVect) (                                       \
    sqrt (aVect.x*aVect.x + aVect.y*aVect.y + aVect.z*aVect.z)  \
)
#define VECT_MAG(aVect)(                                        \
    sqrt (aVect.x*aVect.x + aVect.y*aVect.y + aVect.z*aVect.z)  \
)

/// Negate a vector
///
#define VECT_MINUS(inVect,outVect) {                            \
    VECT_CREATE(-inVect.x, -inVect.y, -inVect.z, outVect)       \
}

/// Find the dot product of two vectors
///
#define VECT_DOT(aVect,bVect) (aVect.x*bVect.x + aVect.y*bVect.y + aVect.z*bVect.z)

/// find a cross product of two vectors
///
#define VECT_CROSS(aVect,bVect,outVect) {                       \
    SPVector ___tmp ;                                           \
    ___tmp.x = aVect.y*bVect.z - aVect.z*bVect.y;               \
    ___tmp.y = aVect.z*bVect.x - aVect.x*bVect.z;               \
    ___tmp.z = aVect.x*bVect.y - aVect.y*bVect.x;               \
    outVect = ___tmp ;                                          \
}

/// Find the equivalent unit vector v
///
#define VECT_UNIT(aVect,outVect) {                              \
    double vect_unit_tmp = 1.0 / VECT_MOD(aVect);               \
    outVect.x = aVect.x*vect_unit_tmp;                          \
    outVect.y = aVect.y*vect_unit_tmp;                          \
    outVect.z = aVect.z*vect_unit_tmp;                          \
}
#define VECT_NORM(aVect,outVect) {                              \
    double vect_unit_tmp = 1.0 / VECT_MOD(aVect);               \
    outVect.x = aVect.x*vect_unit_tmp;                          \
    outVect.y = aVect.y*vect_unit_tmp;                          \
    outVect.z = aVect.z*vect_unit_tmp;                          \
}

typedef struct SPCmplx {
    float r;
    float i;
} cplxf;

typedef struct AABB {
    SPVector AA;
    SPVector BB;
} AABB;

typedef union {
    struct KdTreeLeaf {
        unsigned int flagDimAndOffset;
        // bits 0..1        : splitting dimension
        // bits 2..30       : offset bits
        // bits 31 (sign)   : flag whether node is a leaf
        float splitPosition;
        int Ropes[6];
        AABB aabb;
    } leaf;
    
    struct KdTreeBranch {
        unsigned int flagDimAndOffset;
        // bits 0..30       : offset to first child
        // bit 31 (sign)    : flag whether node is a leaf
        float splitPosition;
        int Ropes[6];
        AABB aabb;
    } branch;
    
} KdData ;
typedef struct Ray {
    SPVector org;    // Origin
    SPVector dir;    // Direction
    double   pow;    // Power for this ray
    double   len;    // Distance travelled to this ray's origin from transmission
    SPVector pol ;   // unit vector of direction of E field of ray
} Ray;

typedef struct rangeAndPower {
    double  range ;
    SPCmplx Es ;
} rangeAndPower ;

typedef struct Hit {
    double dist;
    int trinum;
    double u;
    double v;
} Hit;

typedef struct Triangle {
    int  triNum;    // Triangle ID
    double d;         // Constant of plane equation
    double nd_u;      // Normal.u / normal.k
    double nd_v;      // normal.v / normal.k
    int k;          // projection dimension
    double kbu;
    double kbv;
    double kbd;
    double kcu;
    double kcv;
    double kcd;
    int textureInd;
} Triangle;

typedef struct TriCoords {
    SPVector A ;      // Cartesian coordinates of triangle
    SPVector B ;
    SPVector Cc ;
} TriCoords ;

#define SURFMAXN 4
#define LT 0.04

int factorial(int n) ;
SPCmplx G_func(int n, double gamma) ;

__kernel void reflectPower(__global Triangle * Triangles,       // Array of Triangles. Each triangle references a texture
                           __global Texture  * Textures,        // Array of Textures containing specular, diffuse adn shinyness constants
                           __global Hit *hits,                  // Array of hit locations to x-ref with triangles (and then Textures) for material props
                           __global Ray           *rays,      // Rays arriving at triangle
                           const    int           nRays,      // Number of items in rays, Vrays
                           const    double        k,          // Wavenumber 2PI/lambda
                           const    SPVector      poldir,     // Direction of polarisation of receiver in global coords
                           __global Ray           *Vrays,     // Viewpoint unit vector. The vector direction to calculate power for (usually to receiver)
                           __global double        *ranges,    // Range to receiver for each shadow ray (precalculated in shadowRay generation)
                           __global rangeAndPower *rnp        // Output array of ray power at, and range to, reciever
                           )
{
    unsigned int modulo[5];
    int ind, k, ku, kv ;

    ind = get_global_id(0) ;
    
    if (ind >=0 && ind < nRays ) {
        
        modulo[0] = 0; modulo[1] = 1; modulo[2] = 2; modulo[3] = 0; modulo[4]=1;
        
        T   = Triangles[hits[ind].trinum];
        tex = globalMatProps[T.matInd] ;
        
        // Create surface normal for texture associated with this hitpoint
        //
        k  = T.k;
        ku = modulo[k+1];
        kv = modulo[k+2];
        
        N.cell[k]  = 1;
        N.cell[ku] = T.nd_u;
        N.cell[kv] = T.nd_v;
        VECT_NORM(N,N);

        
        // Sort out the direction cosines for this ray, hitpoint and observation point
        //
        SPVector uvw_sg ; // Direction cosines for scattering ray (HitPoint  to observation point)
        
        // Assuming that the direction of the ray has been normalised then
        // the direction cosine is just the component of direction
        //
        uvw_sg.x = obsDir.x -Vrays[ind].dir.x ;
        uvw_sg.y = obsDir.y -Vrays[ind].dir.y ;
        uvw_sg.z = obsDir.z -Vrays[ind].dir.z ;
        
        // Now need to perform two seperate calculations:
        // 1) The surface integral over triangle 'c' called Ic
        // 2) The surface current over the triangle in terms of Jx & Jy
        //
        // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        //
        // 1) Find the value of Ic
        //
        SPCmplx Ic = surfaceIntegral(k, tri, uvw_sg) ;
        
        // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        //
        // 2) Find the surface current over the triangle in terms of Jx & Jy. Multiply it
        // by the surfaceIntergral Ic and return the Electric field E (Cmplx amp & phase)
        // also return the polarisation direction of the Efield (in global coordinates)
        //
        SPCmplx Es ;
        rng            = ranges[ind] ;
        EField(k, rng, tri, ray, Ic, poldir, &Es) ;
        rnp[ind].Es    = Es ;
        rnp[ind].range = (rng + Vrays[ind].len) * 0.5 ;
        
        return ;
    }
}

SPCmplx surfaceIntegral (double k, triangle tri, SPVector uvw_sg){
    
    double Dp,Dq,D0;
    double Cp,Cq,C0;
    double u = uvw_sg.x ;
    double v = uvw_sg.y ;
    double w = uvw_sg.z ;
    double x1 = tri.AA.x ;
    double y1 = tri.AA.y ;
    double z1 = tri.AA.z ;
    double x2 = tri.BB.x ;
    double y2 = tri.BB.y ;
    double z2 = tri.BB.z ;
    double x3 = tri.CC.x ;
    double y3 = tri.CC.y ;
    double z3 = tri.CC.z ;
    double A  = tri.area ;
    
    Cp = Cq = 0;
    C0 = 1;
    Dp = k * ( ((x1 - x3) * u) + ((y1 - y3) * v) + ((z1 - z3) * w) ) ;
    Dq = k * ( ((x2 - x3) * u) + ((y2 - y3) * v) + ((z2 - z3) * w) ) ;
    D0 = k * ( (x3 * u) + (y3 * v) + (z3 * w) ) ;
    
    double Lt      = LT ;  // Length of Taylor series region
    double magDp   = fabs(Dp) ;
    double magDq   = fabs(Dq) ;
    
    SPCmplx e_jDp, e_jDq, e_jD0 ;
    SPCmplx jDp, jD0, jDq;
    SPCmplx Ic ;
    
    e_jD0.r = cosf(D0); e_jD0.i = sinf(D0) ;
    e_jDp.r = cosf(Dp); e_jDp.i = sinf(Dp) ;
    e_jDq.r = cosf(Dq); e_jDq.i = sinf(Dq) ;
    jDp.r = 0;          jDp.i = Dp ;
    jD0.r = 0;          jD0.i = D0 ;
    jDq.r = 0;          jDq.i = Dq ;
    
    if ( magDp < Lt && magDq >= Lt ){  // Case 1
        SPCmplx sum;
        sum.r = sum.i = 0;
        SPCmplx t, t1, C0_over_nplus1, jD_subs_pow, leftPart ;
        
        // Set up zeroeth term of sum
        //
        jD_subs_pow.r = 1.0;
        jD_subs_pow.i = 0.0 ;
        t1 = G_func(0, -Dq) ;
        CMPLX_MULT(e_jDq, t1, t) ;
        C0_over_nplus1.r = -C0 ;
        C0_over_nplus1.i = 0 ;
        CMPLX_ADD(C0_over_nplus1, t, sum);
        
        // Now sum over n
        //
        for (int n=1; n<=SURFMAXN; n++) {
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
        }
        
        // Finished summation
        //
        CMPLX_SCMULT(2*A, e_jD0, t);
        CMPLX_DIV(t, jDq, t1);
        CMPLX_MULT(t1, sum, Ic);
        
    }else if ( magDp < Lt && magDq < Lt ){  // Case 2
        
        SPCmplx jDp_a[SURFMAXN] ;
        SPCmplx jDq_a[SURFMAXN] ;
        SPCmplx sum; sum.r = sum.i = 0;
        SPCmplx tmp, tmp1;
        
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
        SPCmplx sum; sum.r = sum.i = 0;
        SPCmplx t, t1, t2, jDq_pow, G;
        
        jDq_pow.r = 1.0; jDq_pow.i = 0.0 ;
        G = G_func(1, -Dp) ;
        sum = G ;
        
        for(int n=1; n<SURFMAXN; n++){
            G = G_func(n+1, -Dp) ;
            CMPLX_SCMULT(1.0/(double)(n+1), G, t);
            CMPLX_MULT(jDq_pow, jDq, t1);
            jDq_pow = t1 ;
            CMPLX_SCMULT(1/(double)factorial(n), jDq_pow, t1);
            CMPLX_MULT(t1, t, t2);
            CMPLX_ADD(sum, t2, sum);
        }
        
        // Finished summation
        //
        CMPLX_SCMULT(2*A, e_jD0, t);
        CMPLX_MULT(t, e_jDp, t1);
        CMPLX_MULT(t1, sum, Ic);
        
    } else if ( magDp >= Lt && magDq >= Lt && fabs(Dp - Dq) < Lt ){   // Case 4
        SPCmplx t, t1,t2,t3,t4,t5, jDp_min_jDq, jD_subs_pow;
        SPCmplx sum; sum.r = sum.i = 0;
        SPCmplx leftPart;
        
        // Set up zeroeth term of sum
        //
        CMPLX_SUB(jDp, jDq, jDp_min_jDq);
        jD_subs_pow.r = 1.0; jD_subs_pow.i = 0.0 ;
        // C0 = 1
        //
        t1 = G_func(0, Dq) ;
        CMPLX_SCMULT(-1, t1, t1);
        CMPLX_ADD(t1, e_jDp, sum);
        
        // Now sum over n
        //
        for (int n=1; n<=SURFMAXN; n++) {
            CMPLX_MULT(jD_subs_pow, jDp_min_jDq, t);
            jD_subs_pow = t ;
            CMPLX_SCMULT(1.0/(double)factorial(n), t, leftPart) ;
            t1 = G_func(n, Dq) ;
            CMPLX_SCMULT(-1, t1, t2);
            CMPLX_SCMULT(1.0/(n+1.0), e_jDp, t3);
            CMPLX_ADD(t2, t3, t4);
            CMPLX_MULT(leftPart, t4, t5);
            CMPLX_ADD(t5, sum, sum) ;
        }
        
        // End of summation
        CMPLX_SCMULT(2*A, e_jD0, t);
        CMPLX_DIV(t, jDq, t1);
        CMPLX_MULT(t1, sum, Ic);
        
    } else {
        SPCmplx term1,term2,term3;
        double Dq_minus_Dp = Dq-Dp ;
        CMPLX_SCMULT((C0 / (Dp * Dq_minus_Dp)), e_jDp, term1);
        CMPLX_SCMULT((C0 / (Dq * Dq_minus_Dp)), e_jDq, term2);
        term3.r = C0 / (Dp * Dq) ; term3.i = 0;
        SPCmplx braces ;
        CMPLX_SUB(term1, term2,braces);
        CMPLX_SUB(braces, term3, braces) ;
        SPCmplx t ;
        CMPLX_SCMULT(2*A, e_jD0, t) ;
        CMPLX_MULT(t, braces, Ic);
    }
    
    return Ic ;
}

int factorial(int n)
{
    int c;
    int result = 1;
    
    for (c = 1; c <= n; c++)
        result = result * c;
    
    return result;
}

SPCmplx G_func(int n, double gamma){
    SPCmplx e_jgamma, jgamma, G, t1, t2, ans ;
    e_jgamma.r = cos(gamma) ;
    e_jgamma.i = sin(gamma) ;
    jgamma.r = 0 ; jgamma.i = gamma ;
    
    if (n==0) {
        SPCmplx one ;
        one.r = 1;
        one.i = 0;
        CMPLX_SUB(e_jgamma, one, t1) ;
        CMPLX_DIV(t1, jgamma, ans) ;
        
    }else{
        G = G_func(n-1, gamma) ;
        CMPLX_SCMULT(n, G, t1);
        CMPLX_SUB(e_jgamma, t1, t2) ;
        CMPLX_DIV(t2, jgamma, ans) ;
    }
    return ans ;
}

void EField(double k, double r, triangle tri, Ray ray, SPCmplx Ic, SPVector Es_parrdir, SPVector Es_perpdir, SPCmplx *Es_parr, SPCmplx *Es_perp ){
    
    SPVector Eig ;
    VECT_SCMULT(ray.pol, ray.pow, Eig) ;
    double uvw_ig[3], uvw_il[3]; // Direction cosines for incident Ray
    
    // Assuming that the direction of the ray has been normalised then
    // the direction cosine is just the component of direction
    //
    uvw_ig[0] = -ray.dir.x;
    uvw_ig[1] = -ray.dir.y;
    uvw_ig[2] = -ray.dir.z;
    
    matmul(tri.globalToLocalMat, uvw_ig, uvw_il, 3, 3, 1, 3);
    
    double sin_theta_il, cos_theta_il, tan_phi_il, sin_phi_il, cos_phi_il ;
    sin_theta_il = sqrt(uvw_il[0]*uvw_il[0] +  uvw_il[1] * uvw_il[1] ) ;
    cos_theta_il = sqrt(1 - sin_theta_il*sin_theta_il) ;
    if(sin_theta_il == 0){
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
    double E_ig[3], E_il[3] ;
    SPVector Eil, theta_l_hat, phi_l_hat, z_hat ;
    SPVector Raydir_l ;
    VECT_CREATE(uvw_il[0], uvw_il[1], uvw_il[2], Raydir_l) ;
    E_ig[0] = Eig.x ;
    E_ig[1] = Eig.y ;
    E_ig[2] = Eig.z ;
    matmul(tri.globalToLocalMat, E_ig, E_il, 3, 3, 1, 3) ;
    VECT_CREATE(E_il[0], E_il[1], E_il[2], Eil) ;
    VECT_CREATE(0, 0, 1, z_hat) ;
    if(fabs(VECT_DOT(z_hat, Raydir_l)) == 1.0){
        VECT_CREATE(1, 0, 0, phi_l_hat) ;
    }else{
        VECT_CROSS (z_hat, Raydir_l, phi_l_hat) ;
    }
    VECT_CROSS(Raydir_l, phi_l_hat, theta_l_hat) ;
    VECT_NORM(phi_l_hat, phi_l_hat) ;
    VECT_NORM(theta_l_hat, theta_l_hat) ;
    double Eiphi_l   = VECT_DOT(Eil, phi_l_hat) ;
    double Eitheta_l = VECT_DOT(Eil, theta_l_hat) ;
    
    // Calculate Gamma_parallel and Gamma_perpendicular
    //
    double Rs, GamParr, GamPerp ;
    Rs      = tri.Rs ;
    GamParr = -1.0 * Z0 * cos_theta_il / (2*Rs + Z0*cos_theta_il) ;
    GamPerp = -1.0 * Z0 / ( 2.0*Rs*cos_theta_il + Z0);
    
    // Now calculate Jx_local and Jy_local
    //
    double Jx_l, Jy_l ;
    Jx_l = ((-1.0 * Eitheta_l * cos_phi_il * GamParr / Z0) + (Eiphi_l * sin_phi_il * GamPerp / Z0)) * cos_theta_il ;
    Jy_l = ((-1.0 * Eitheta_l * sin_phi_il * GamParr / Z0) - (Eiphi_l * cos_phi_il * GamPerp / Z0)) * cos_theta_il ;
    
    // Convert the current vector to global coordinates
    //
    double J_l[3], J_g[3] ;
    J_l[0] = Jx_l ;
    J_l[1] = Jy_l ;
    J_l[2] = 0.0  ;
    matmul(tri.localToGlobalMat, J_l, J_g, 3, 3, 1, 3) ;
    SPVector Jg;
    VECT_CREATE(J_g[0], J_g[1], J_g[2], Jg) ;
    double A = VECT_DOT(Es_parrdir, Jg);
    double B = VECT_DOT(Es_perpdir, Jg);
    
    // Work out scaler (complex) component of E field
    //
    SPCmplx jkZ0_o_4PIr, tmp1, scaler, e_jkr;
    CMPLX_F_MAKE(0, -k*Z0/(4*SIPC_pi*r), jkZ0_o_4PIr) ;
    e_jkr.r = cos(-k*r);
    e_jkr.i = sin(-k*r);
    CMPLX_MULT(jkZ0_o_4PIr, e_jkr, tmp1);
    CMPLX_MULT(tmp1, Ic, scaler);
    SPCmplx par,per ;
    CMPLX_SCMULT(A, scaler, par) ;
    CMPLX_SCMULT(B, scaler, per) ;
    Es_parr->r = par.r ; Es_parr->i = par.i ;
    Es_perp->r = per.r ; Es_perp->i = per.i ;
    return ;
}


