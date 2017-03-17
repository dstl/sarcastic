#pragma OPENCL EXTENSION cl_khr_fp64: enable
#define SIPC_pi 3.14159265358979323846
#define SURFMAXN 4
#define LT 0.04

#define RAD2DEG(x)(180.0 * x / SIPC_pi)

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

/// Calculate complex magnitude
///
#define CMPLX_MAG(a) (sqrt(a.r * a.r + a.i * a.i))

/// Calculate the phase of a complex number
///
#define CMPLX_PHASE(a) (atan2(a.i, a.r))

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

// materialProperties.h
#define MATBYTES 128
typedef struct scatProps {
    char   matname[MATBYTES] ;
    float  corlen      ;
    float  roughness   ;
    float  Rs          ;
    float  Rm          ;
    float  specular    ;
    float  diffuse     ;
    float  shinyness   ;
} scatProps ;
#define NMATERIALS 9
static constant scatProps materialProperties[NMATERIALS] = {
//   Name          corrLen      Roughness   Rs          Rm      Specular    Diffuse     Shinyness
    {"MATERIAL",    100.0,      0.0,        0.0,        9e9,   1.0,        0.0,        50.0        },
    {"ASPHALT",     0.5,        0.005,      1.0e18,     9e9,   0.8,        0.2,        30.0        },
    {"BRICK",       0.1,        0.001,      1.0e18,     9e9,   0.7,        0.3,        20.0        },
    {"CONCRETE",    0.2,        0.01,       120.0,      9e9,   0.3,        0.7,        10.0        },
    {"METAL",       100.0,      0.0,        1.0e-8,     9e9,   1.0,        0.0,        50.0        },
    {"ROOFING",     0.1,        0.1,        1.0e18,     9e9,   0.6,        0.4,        40.0        },
    {"VEGETATION",  0.01,       0.1,        2000.0,     9e9,   0.2,        0.8,        5.0         },
    {"WATER",       0.01,       0.1,        2.0e1,      9e9,   1.0,        0.0,        50.0        },
    {"WOOD",        0.1,        0.001,      1.0e14,     9e9,   0.6,        0.4,        10.0        }
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

typedef struct {
    float pabs;
    float parg;
} SPCmplxPol;

typedef struct HitPoint {
    SPVector hit;       // Location of hitpoint in x,y,z
    int tri;            // index of triangle that this hit is on
} HitPoint ;

// Forward declarations here
//
int factorial(int n) ;
SPCmplx G_func4(double gamma);
SPCmplx G_func3(double gamma);
SPCmplx G_func2(double gamma);
SPCmplx G_func1(double gamma);
SPCmplx G_func0(double gamma);
SPVector translateVector(SPVector v, double *matrix) ;
void matmul(double *A,double *B, double *O, int Ax, int Ay,int Bx, int By) ;
void surfaceCurrents(SPVector Eig, SPVector Ri_hat, Triangle tri, SPVector Vdir, SPVector Hdir, double *Jpar, double *Jper) ;
SPCmplx ludwigIntegral(double Dp,double Dq,double D0, double A) ;

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
                      int firstBounce,          // if 1 then PO calcs use origin for field combination
                      __global SPCmplx *EsVs,   // output array of scattered field strengths for V pol
                      __global SPCmplx *EsHs)   // output array of scattered field strengths for H pol
{
    
    int ind ;
    double Dp,Dq,D0, u,v,w ;
    double x1,x2,x3,y1,y2,y3,z1,z2,z3,A ;
    double r_mag, ri_mag, cos_theta_i, sin_theta_i,cos_phi_i,sin_phi_i;
    double Rs, Rm;
    double Ei_phi_l, Hi_phi_l, Ei_theta_l, Hi_theta_l;
    double globalToLocalMat[9];
    double localToGlobalMat[9];
    double Jx,Jy,Mx,My, J_mag, M_mag, Jpolang, Mpolang, JnM_par, JnM_per, phs_ig, J_par,J_per, M_par,M_per ;
    SPVector Jl, Ml;
    SPCmplx Jc_par, Jc_per, minus_jkZ0_o_4PIr;
    SPVector Jg, Mg;
    SPCmplx Ic, t0_, t1_, t2_, e_minus_jkr, e_minus_jkri ;
    SPVector r, ri, r_hat, ri_hat, uvw_il, uvw_sl, g, h, Eig, Eil, Hil, phi_il_hat, theta_il_hat, hp ;
    Ray ray; Triangle tri;
    
    ind = get_global_id(0) ;

    if (ind >=0 && ind < nrays ) {

        SPVector origin;
        origin.x = (tris[hitpoints[ind].tri].AA.x + tris[hitpoints[ind].tri].BB.x + tris[hitpoints[ind].tri].CC.x) / 3;
        origin.y = (tris[hitpoints[ind].tri].AA.y + tris[hitpoints[ind].tri].BB.y + tris[hitpoints[ind].tri].CC.y) / 3;
        origin.z = (tris[hitpoints[ind].tri].AA.z + tris[hitpoints[ind].tri].BB.z + tris[hitpoints[ind].tri].CC.z) / 3;
        VECT_SUB(tris[hitpoints[ind].tri].AA,origin,tris[hitpoints[ind].tri].AA);
        VECT_SUB(tris[hitpoints[ind].tri].BB,origin,tris[hitpoints[ind].tri].BB);
        VECT_SUB(tris[hitpoints[ind].tri].CC,origin,tris[hitpoints[ind].tri].CC);

        hp  = hitpoints[ind].hit ;
        ray = rays[ind] ;
        tri = tris[hitpoints[ind].tri] ;
        Rs  = materialProperties[tri.matId].Rs ;
        Rm  = materialProperties[tri.matId].Rm ;
        for(int i=0; i<9; i++)globalToLocalMat[i] = tri.globalToLocalMat[i];
        for(int i=0; i<9; i++)localToGlobalMat[i] = tri.localToGlobalMat[i];
        
        VECT_SUB(RxPnt,hp,r);
        VECT_SUB(ray.org,hp,ri);
        
        SPVector illDir,obsDir;
        VECT_SUB(RxPnt,origin,obsDir);
        VECT_SUB(ray.org,origin,illDir);
        VECT_NORM(illDir,illDir);
        VECT_NORM(obsDir,obsDir);
        
        r_mag  = VECT_MAG(r);
        ri_mag = VECT_MAG(ri);
       
        VECT_SCMULT(r,  1.0/r_mag,  r_hat );
        VECT_SCMULT(ri, 1.0/ri_mag, ri_hat);
        
        // Direction cosine vectors for incident and scattering
        // in 'local' coordinates of triangle facet
        //
        uvw_il = translateVector(illDir, globalToLocalMat) ;
        uvw_sl = translateVector(obsDir,  globalToLocalMat) ;
        
        // Calculate the surface integral Ic
        //
        g.x = ri_hat.x ;
        g.y = ri_hat.y ;
        g.z = ri_hat.z ;

        h.x = r_hat.x ;
        h.y = r_hat.y ;
        h.z = r_hat.z ;
         
        u = h.x + g.x ;
        v = h.y + g.y ;
        w = h.z + g.z ;

        x1 = tri.AA.x ; y1 = tri.AA.y ; z1 = tri.AA.z ;
        x2 = tri.BB.x ; y2 = tri.BB.y ; z2 = tri.BB.z ;
        x3 = tri.CC.x ; y3 = tri.CC.y ; z3 = tri.CC.z ;
        A  = tri.area ;
        
        Dp = k * ( ((x1 - x3) * u) + ((y1 - y3) * v) + ((z1 - z3) * w) ) ;
        Dq = k * ( ((x2 - x3) * u) + ((y2 - y3) * v) + ((z2 - z3) * w) ) ;
        D0 = k * ( (x3 * u) + (y3 * v) + (z3 * w) ) ;
        
        Ic = ludwigIntegral(Dp, Dq, D0, A);

        // Find cos_theta_i, sin_phi_i and cos_phi_i
        //
        cos_theta_i = uvw_il.z ;
        sin_theta_i = sqrt(uvw_il.x*uvw_il.x +  uvw_il.y * uvw_il.y ) ;

        if( (cos_theta_i==0 && Rs == 0)  || (cos_theta_i == 0 && Rm == 0) ){
            cos_theta_i = 1.0e-8;
        }
        
        if(fabs(cos_theta_i) >= 0.9999){
            VECT_CREATE(1, 0, 0, phi_il_hat) ;
            VECT_CROSS(phi_il_hat, uvw_il, theta_il_hat) ;
            VECT_NORM(theta_il_hat, theta_il_hat) ;
            sin_phi_i = 1.0 ;
            cos_phi_i = 0.0 ;
        }else{
            VECT_CROSS (zz_hat, uvw_il, phi_il_hat) ;
            VECT_NORM(phi_il_hat, phi_il_hat) ;
            VECT_CROSS(phi_il_hat, uvw_il, theta_il_hat) ;
            VECT_NORM(theta_il_hat, theta_il_hat) ;
            sin_phi_i   = uvw_il.y / sin_theta_i ;
            cos_phi_i   = uvw_il.x / sin_theta_i ;
        }        

        // Find E and H fields in 'local' facet coordinate system
        //
        VECT_SCMULT(ray.pol, sqrt(ray.pow / (4 * SIPC_pi * (ri_mag+ray.len)*(ri_mag+ray.len))), Eig) ;

        Eil = translateVector(Eig, globalToLocalMat) ;
        
        // Eig is the Efield vector in the direction of the hitpoint.
        // Convert this to be in the direction of the origin. (assumes that waves are planar)
        //
        SPVector Eil_x_uvw_il,new_Eil_dir;
        double Eil_mag ;
        
        Eil_mag = VECT_MAG(Eil);
        VECT_CROSS(Eil, uvw_il, Eil_x_uvw_il);
        VECT_CROSS(uvw_il, Eil_x_uvw_il, new_Eil_dir);
        VECT_NORM(new_Eil_dir,new_Eil_dir);
        VECT_SCMULT(new_Eil_dir, Eil_mag, Eil);
        
        //////////////////////////////////////////

        VECT_CROSS(Eil, uvw_il, Hil);
        VECT_NORM(Hil, Hil) ;
        VECT_SCMULT(Hil, VECT_MAG(Eil)/Z0, Hil);
        
        Ei_phi_l   = VECT_DOT(Eil, phi_il_hat) ;
        Ei_theta_l = VECT_DOT(Eil, theta_il_hat);
        Hi_phi_l   = VECT_DOT(Hil, phi_il_hat) ;
        Hi_theta_l = VECT_DOT(Hil, theta_il_hat);
        
        // Calculate Surface currents Jx,Jy,Mx,My
        //
        Jx = 2 * cos_theta_i * (
                                (cos_phi_i * Ei_theta_l / ((Z0 * cos_theta_i) + (2 * Rs)))
                              - (sin_phi_i * Ei_phi_l   / (Z0 + (2 * Rs * cos_theta_i)))
                                ) ;
        Jy = 2 * cos_theta_i * (
                                (sin_phi_i * Ei_theta_l / ((Z0 * cos_theta_i) + (2 * Rs)))
                              + (cos_phi_i * Ei_phi_l   / (Z0 + (2 * Rs * cos_theta_i)))
                                ) ;
        Mx = 2 * Z0 * (
                        (cos_theta_i * cos_theta_i * cos_phi_i * Hi_theta_l / ( 1 + (2 * Z0 * Rm * cos_theta_i)) )
                      - (sin_phi_i * Hi_phi_l / (cos_theta_i + (2 * Z0 * Rm)))
                       ) ;
        My = 2 * Z0 * (
                        (-1 * cos_theta_i * cos_theta_i * sin_phi_i * Hi_theta_l / (1 + (2 * Z0 * Rm * cos_theta_i)) )
                      + (cos_phi_i * Hi_phi_l / (cos_theta_i + (2 * Z0 * Rm)))
                       ) ;
        VECT_CREATE(Jx, Jy, 0, Jl);
        VECT_CREATE(Mx, My, 0, Ml);
        
        // Translate surface currents into global coordinate system
        //
        Jg = translateVector(Jl, localToGlobalMat) ;
        Mg = translateVector(Ml, localToGlobalMat) ;

        J_mag = VECT_MAG(Jg);
        M_mag = VECT_MAG(Mg);
        J_par = fabs(VECT_DOT(Vdir, Jg));
        J_per = fabs(VECT_DOT(Hdir, Jg));
        M_par = fabs(VECT_DOT(Vdir, Mg));
        M_per = fabs(VECT_DOT(Hdir, Mg));
        
        Jpolang = atan2(J_par,J_per);
        Mpolang = atan2(M_par,M_per);
        J_par = J_mag*sin(Jpolang);
        J_per = J_mag*cos(Jpolang);
        M_par = M_mag*sin(Mpolang);
        M_per = M_mag*cos(Mpolang);
        
        // Sum current components together
        //
        JnM_par = J_par + M_par ;
        JnM_per = J_per + M_per ;
        
        // Calculate initial phase at origin of incident ray
        //
        phs_ig = -(k * ray.len) - (SIPC_pi/2.0) ;

        // Create complex current terms for V and H pol
        //
        CMPLX_F_MAKE(JnM_par*cos(phs_ig), JnM_par*sin(phs_ig), Jc_par);
        CMPLX_F_MAKE(JnM_per*cos(phs_ig), JnM_per*sin(phs_ig), Jc_per);
        
        // Create exponential (phase) terms
        //
        CMPLX_F_MAKE(0, -k*Z0/(4*SIPC_pi*r_mag), minus_jkZ0_o_4PIr) ;
        CMPLX_F_MAKE(cos(-k * r_mag),  sin(-k * r_mag) , e_minus_jkr ) ;
        CMPLX_F_MAKE(cos(-k * ri_mag), sin(-k * ri_mag), e_minus_jkri);

        // Multiply everything together
        //
        CMPLX_MULT(minus_jkZ0_o_4PIr, e_minus_jkr, t0_);
        CMPLX_MULT(t0_, e_minus_jkri, t1_);
        CMPLX_MULT(t1_, Ic, t2_);
        CMPLX_MULT(Jc_par, t2_, EsVs[ind]);
        CMPLX_MULT(Jc_per, t2_, EsHs[ind]);
        
    }
    return ;
}

SPCmplx ludwigIntegral(double Dp,double Dq,double D0, double A){
    
    int n;
    double Cp,Cq,C0;
    double Lt, magDp, magDq, Dq_minus_Dp;
    SPCmplx e_jDp, e_jDq, e_jD0, jDp, jD0, jDq, jDq_pow,jDp_min_jDq,jDp_a[SURFMAXN], jDq_a[SURFMAXN] ;
    SPCmplx sum, tmp, tmp1,t, t1, t2, t3,t4,t5, term1,term2,term3, braces ;
    SPCmplx C0_over_nplus1, jD_subs_pow, leftPart, G, Ic ;

    Cp = Cq = 0;
    C0 = 1;
    Lt      = LT ;          // Length of Taylor series region
    magDp   = fabs(Dp) ;
    magDq   = fabs(Dq) ;
    e_jD0.r = cos(D0) ; e_jD0.i = sin(D0) ;
    e_jDp.r = cos(Dp) ; e_jDp.i = sin(Dp) ;
    e_jDq.r = cos(Dq) ; e_jDq.i = sin(Dq) ;
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
        
    } else if ( magDp >= Lt && magDq >= Lt && fabs(Dq - Dp) < Lt ){   // Case 4

        sum.r = sum.i = 0;
        // Set up zeroeth term of sum
        //
        CMPLX_SUB(jDp, jDq, jDp_min_jDq);
        jD_subs_pow.r = 1.0; jD_subs_pow.i = 0.0 ;
        // C0 = 1
        //
        t1 = G_func0(Dq) ;
        CMPLX_SCMULT(-1, t1, t1);
        CMPLX_ADD(t1, e_jDq, sum);
        
        // Now sum over n
        //
        n=1;
        CMPLX_MULT(jD_subs_pow, jDp_min_jDq, t);
        jD_subs_pow = t ;
        CMPLX_SCMULT(1.0/(double)factorial(n), t, leftPart) ;
        t1 = G_func1(Dq) ;
        CMPLX_SCMULT(-1, t1, t2);
        CMPLX_SCMULT(1.0/(n+1.0), e_jDq, t3);
        CMPLX_ADD(t2, t3, t4);
        CMPLX_MULT(leftPart, t4, t5);
        CMPLX_ADD(t5, sum, sum) ;
        n=2;
        CMPLX_MULT(jD_subs_pow, jDp_min_jDq, t);
        jD_subs_pow = t ;
        CMPLX_SCMULT(1.0/(double)factorial(n), t, leftPart) ;
        t1 = G_func2(Dq) ;
        CMPLX_SCMULT(-1, t1, t2);
        CMPLX_SCMULT(1.0/(n+1.0), e_jDq, t3);
        CMPLX_ADD(t2, t3, t4);
        CMPLX_MULT(leftPart, t4, t5);
        CMPLX_ADD(t5, sum, sum) ;
        n=3;
        CMPLX_MULT(jD_subs_pow, jDp_min_jDq, t);
        jD_subs_pow = t ;
        CMPLX_SCMULT(1.0/(double)factorial(n), t, leftPart) ;
        t1 = G_func3(Dq) ;
        CMPLX_SCMULT(-1, t1, t2);
        CMPLX_SCMULT(1.0/(n+1.0), e_jDq, t3);
        CMPLX_ADD(t2, t3, t4);
        CMPLX_MULT(leftPart, t4, t5);
        CMPLX_ADD(t5, sum, sum) ;
        n=4;
        CMPLX_MULT(jD_subs_pow, jDp_min_jDq, t);
        jD_subs_pow = t ;
        CMPLX_SCMULT(1.0/(double)factorial(n), t, leftPart) ;
        t1 = G_func4(Dq) ;
        CMPLX_SCMULT(-1, t1, t2);
        CMPLX_SCMULT(1.0/(n+1.0), e_jDq, t3);
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
    
    return Ic;
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

void matmul(double *A, double *B, double *O, int Ax, int Ay,int Bx, int By)
{
    // Ax must be equal to By !!!
    //
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

SPVector translateVector(SPVector v, double *matrix){
    double vectAsArr[3], ansAsArr[3] ;
    SPVector ans;
    vectAsArr[0] = v.x;
    vectAsArr[1] = v.y;
    vectAsArr[2] = v.z;
    matmul(matrix, vectAsArr, ansAsArr, 3, 3, 1, 3) ;
    VECT_CREATE(ansAsArr[0], ansAsArr[1], ansAsArr[2], ans) ;
    return ans ;
}



