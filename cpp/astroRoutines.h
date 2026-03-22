// =========================================================================
// astroRoutines.h
// Shared astrodynamics and frame-transformation helpers for DA propagators.
// Included by all state propagation executables in cpp/.
// =========================================================================

#define _USE_MATH_DEFINES
#define WITH_ALGEBRAICMATRIX

#include <dace/dace.h>
#include <cmath>
#include <ctime>
#include <fstream>
#include <iomanip>
using namespace std;
using namespace DACE;

DA modulusLocal(DA numer, double denom)
{
    double numerf = numer.cons();
    double modf = remainder(numerf, denom);
    return numer - numerf + modf;
}

double modulusLocal(double numer, double denom)
{
    return remainder(numer, denom);
}

template<typename T> T true2eccAnomaly(const T theta, const T e)
{
    return 2.0 * atan2(sqrt(1. - e)*sin(theta / 2.), sqrt(1. + e) * cos(theta / 2.));
}

template<typename T> T ecc2trueAnomaly(const T E, const T e)
{
    return 2.0 * atan2(sqrt(1. + e)*sin(E / 2.), sqrt(1. - e) * cos(E / 2.));
}

template<typename T> T mean2eccAnomaly(const T M, const T e)
{
    T E = M;
    
    for (int i = 0; i < 20; i++) {
        E = M + e*sin(E);
    }
    return E;
}

template<typename T> T mean2trueAnomaly(const T M, const T e)
{
    T E = mean2eccAnomaly(M, e);
    
    return ecc2trueAnomaly(E, e);
}


template<typename T> T true2meanAnomaly(const T theta, const T e)
{
    T E = true2eccAnomaly(theta, e);
    
    return E - e*sin(E);
}


double min( double a, double b ) {
    
    double res = a;
    if(a>b){res = b;}
    return res;
}


double max( double a, double b ) {
    
    double res = a;
    if(b>a){res = b;}
    return res;
}


double normtmp( int N, vector<double> X ) {
    
    int I;
    double res = 0.0;
    
    for (I=0; I<N; I++) {
        if (X[I]<0) { X[I] = -X[I]; }
        res = max(res,X[I]);
    }
    
    return res;
}


template<typename T> AlgebraicVector<T> RK78(const int N, AlgebraicVector<T> Y0, const double X0, const double X1, AlgebraicVector<T>(*f)(const AlgebraicVector<T>&, double))
{
    
    double ERREST;
    double H0 = 0.001; double HS = 0.1; double H1 = 1000.0;
    double EPS = 1.e-12; double BS = 20*EPS;
    
    int I, J, K;
    
    // AlgebraicMatrix<T> Z(N,16);
    
    vector< vector<T> > Z;
    for( I = 0; I < N; I++ )
    {
        vector<T> tmp(16);
        Z.push_back(tmp);
    }
    
    AlgebraicVector<T> Y1(N);
    vector<double> Y1cons(N);
    
    double VIHMAX = 0.0, X, H;
    double RFNORM, HH0, HH1;
    
    double HSQR = 1.0/9.0;
    double A[13], C[13], D[13];
    double B[13][12];
    
    
    
    A[0] = 0.0; A[1] = 1.0/18.0; A[2] = 1.0/12.0; A[3] = 1.0/8.0; A[4] = 5.0/16.0; A[5] = 3.0/8.0;
    A[6] = 59.0/400.0; A[7] = 93.0/200.0; A[8] = 5490023248.0/9719169821.0; A[9] = 13.0/20.0; A[10] = 1201146811.0/1299019798.0; A[11] = 1.0;
    A[12] = 1.0;
    
    B[0][0] = 0.0; B[0][1] = 0.0; B[0][2] = 0.0; B[0][3] = 0.0; B[0][4] = 0.0;
    B[0][5] = 0.0; B[0][6] = 0.0; B[0][7] = 0.0; B[0][8] = 0.0; B[0][9] = 0.0;
    B[0][10] = 0.0; B[0][11] = 0.0;
    
    B[1][0] = 1.0/18.0; B[1][1] = 0.0; B[1][2] = 0.0; B[1][3] = 0.0; B[1][4] = 0.0;
    B[1][5] = 0.0; B[1][6] = 0.0; B[1][7] = 0.0; B[1][8] = 0.0; B[1][9] = 0.0;
    B[1][10] = 0.0; B[1][11] = 0.0;
    
    B[2][0] = 1.0/48.0; B[2][1] = 1.0/16.0; B[2][2] = 0.0; B[2][3] = 0.0; B[2][4] = 0.0;
    B[2][5] = 0.0; B[2][6] = 0.0; B[2][7] = 0.0; B[2][8] = 0.0; B[2][9] = 0.0;
    B[2][10] = 0.0; B[2][11] = 0.0;
    
    B[3][0] = 1.0/32.0; B[3][1] = 0.0; B[3][2] = 3.0/32.0; B[3][3] = 0.0; B[3][4] = 0.0;
    B[3][5] = 0.0; B[3][6] = 0.0; B[3][7] = 0.0; B[3][8] = 0.0; B[3][9] = 0.0;
    B[3][10] = 0.0; B[3][11] = 0.0;

    B[4][0] = 5.0/16.0; B[4][1] = 0.0; B[4][2] = -75.0/64.0; B[4][3] = 75.0/64.0; B[4][4] = 0.0;
    B[4][5] = 0.0; B[4][6] = 0.0; B[4][7] = 0.0; B[4][8] = 0.0; B[4][9] = 0.0;
    B[4][10] = 0.0; B[4][11] = 0.0;

    B[5][0] = 3.0/80.0; B[5][1] = 0.0; B[5][2] = 0.0; B[5][3] = 3.0/16.0; B[5][4] = 3.0/20.0;
    B[5][5] = 0.0; B[5][6] = 0.0; B[5][7] = 0.0; B[5][8] = 0.0; B[5][9] = 0.0;
    B[5][10] = 0.0; B[5][11] = 0.0;
    
    B[6][0] = 29443841.0/614563906.0; B[6][1] = 0.0; B[6][2] = 0.0; B[6][3] = 77736538.0/692538347.0; B[6][4] = -28693883.0/1125000000.0;
    B[6][5] = 23124283.0/1800000000.0; B[6][6] = 0.0; B[6][7] = 0.0; B[6][8] = 0.0; B[6][9] = 0.0;
    B[6][10] = 0.0; B[6][11] = 0.0;

    B[7][0] = 16016141.0/946692911.0; B[7][1] = 0.0; B[7][2] = 0.0; B[7][3] = 61564180.0/158732637.0; B[7][4] = 22789713.0/633445777.0;
    B[7][5] = 545815736.0/2771057229.0; B[7][6] = -180193667.0/1043307555.0; B[7][7] = 0.0; B[7][8] = 0.0; B[7][9] = 0.0;
    B[7][10] = 0.0; B[7][11] = 0.0;
    
    B[8][0] = 39632708.0/573591083.0; B[8][1] = 0.0; B[8][2] = 0.0; B[8][3] = -433636366.0/683701615.0; B[8][4] = -421739975.0/2616292301.0;
    B[8][5] = 100302831.0/723423059.0; B[8][6] = 790204164.0/839813087.0; B[8][7] = 800635310.0/3783071287.0; B[8][8] = 0.0; B[8][9] = 0.0;
    B[8][10] = 0.0; B[8][11] = 0.0;

    B[9][0] = 246121993.0/1340847787.0; B[9][1] = 0.0; B[9][2] = 0.0; B[9][3] = -37695042795.0/15268766246.0; B[9][4] = -309121744.0/1061227803.0;
    B[9][5] = -12992083.0/490766935.0; B[9][6] = 6005943493.0/2108947869.0; B[9][7] = 393006217.0/1396673457.0; B[9][8] = 123872331.0/1001029789.0; B[9][9] = 0.0;
    B[9][10] = 0.0; B[9][11] = 0.0;

    B[10][0] = -1028468189.0/846180014.0; B[10][1] = 0.0; B[10][2] = 0.0; B[10][3] = 8478235783.0/508512852.0; B[10][4] = 1311729495.0/1432422823.0;
    B[10][5] = -10304129995.0/1701304382.0; B[10][6] = -48777925059.0/3047939560.0; B[10][7] = 15336726248.0/1032824649.0; B[10][8] = -45442868181.0/3398467696.0; B[10][9] = 3065993473.0/597172653.0;
    B[10][10] = 0.0; B[10][11] = 0.0;

    B[11][0] = 185892177.0/718116043.0; B[11][1] = 0.0; B[11][2] = 0.0; B[11][3] = -3185094517.0/667107341.0; B[11][4] = -477755414.0/1098053517.0;
    B[11][5] = -703635378.0/230739211.0; B[11][6] = 5731566787.0/1027545527.0; B[11][7] = 5232866602.0/850066563.0; B[11][8] = -4093664535.0/808688257.0; B[11][9] = 3962137247.0/1805957418.0;
    B[11][10] = 65686358.0/487910083.0; B[11][11] = 0.0;

    B[12][0] = 403863854.0/491063109.0; B[12][1] = 0.0; B[12][2] = 0.0; B[12][3] = - 5068492393.0/434740067.0; B[12][4] = -411421997.0/543043805.0;
    B[12][5] = 652783627.0/914296604.0; B[12][6] = 11173962825.0/925320556.0; B[12][7] = -13158990841.0/6184727034.0; B[12][8] = 3936647629.0/1978049680.0; B[12][9] = -160528059.0/685178525.0;
    B[12][10] = 248638103.0/1413531060.0; B[12][11] = 0.0;

    C[0] = 14005451.0/335480064.0; C[1] = 0.0; C[2] = 0.0; C[3] = 0.0; C[4] = 0.0; C[5] = -59238493.0/1068277825.0;
    C[6] = 181606767.0/758867731.0; C[7] = 561292985.0/797845732.0; C[8] = -1041891430.0/1371343529.0; C[9] = 760417239.0/1151165299.0; C[10] = 118820643.0/751138087.0; C[11] = -528747749.0/2220607170.0;
    C[12] = 1.0/4.0;
  
    D[0] = 13451932.0/455176623.0; D[1] = 0.0; D[2] = 0.0; D[3] = 0.0; D[4] = 0.0; D[5] = -808719846.0/976000145.0;
    D[6] = 1757004468.0/5645159321.0; D[7] = 656045339.0/265891186.0; D[8] = -3867574721.0/1518517206.0; D[9] = 465885868.0/322736535.0; D[10] = 53011238.0/667516719.0; D[11] = 2.0/45.0;
    D[12] = 0.0;
    
    for( I = 0; I < N; I++ )
    {
        Z[I][0] = Y0[I];
        Z[I][1] = 0.0 ;
    }
    
    H = abs(HS) ; HH0 = abs(H0) ; HH1 = abs(H1) ;
    X = X0 ; RFNORM = 0.0 ; ERREST = 0.0 ;
    
    while(X != X1){
        
        // compute new stepsize
        if (RFNORM != 0) {H = H*min(4.0,exp(HSQR*log(EPS/RFNORM)));}
        if (abs(H)>abs(HH1)) { H = HH1; } else if ( abs(H)<abs(HH0)*0.99 ) {
            H = HH0;
            cout << "--- WARNING, MINIMUM STEPSIZE REACHED IN RK" << endl;
        }
        
        if ((X+H-X1)*H>0) { H = X1-X; }
        
        for (J = 0; J<13; J++) {
            
            for (I = 0; I<N; I++) {
                
                Y0[I] = 0.0 ; // EVALUATE RHS AT 13 POINTS
                
                for (K=0; K<J; K++) { Y0[I] = Y0[I] + Z[I][K+3]*B[J][K];}
                
                Y0[I] = H*Y0[I] + Z[I][0];
            }

            Y1 = f(Y0, X+H*A[J]);
            
            for (I = 0; I<N; I++) { Z[I][J+3] = Y1[I]; }
        }
        
        for (I = 0; I<N; I++) {
            
            Z[I][1] = 0.0 ; Z[I][2] = 0.0 ; // EXECUTE 7TH,8TH ORDER STEPS
            
            for (J = 0; J<13; J++) {
                Z[I][1] = Z[I][1] + Z[I][J+3]*D[J];
                Z[I][2] = Z[I][2] + Z[I][J+3]*C[J];
            }
            
            Y1[I] = (Z[I][2]-Z[I][1])*H;
            Z[I][2] = Z[I][2]*H+Z[I][0];
        }


        for (I = 0; I<N; I++) {Y1cons[I] = cons(Y1[I]);}
        
        RFNORM = normtmp(N,Y1cons); // ESTIMATE ERROR AND DECIDE ABOUT BACKSTEP
        
        if ((RFNORM>BS) && (abs(H/H0)>1.2)) {
            H = H/3.0;
            RFNORM = 0;
        }
        else {
            for (I = 0; I<N; I++) {Z[I][0] = Z[I][2];}
            X = X + H;
            VIHMAX = max(VIHMAX,H);
            ERREST = ERREST + RFNORM;
        }
    }
    
    for (I = 0; I<N; I++) {Y1[I] = Z[I][0];}
    
    return Y1;
    
}


template<typename T> T atmdensity(T h){
    
    const double p1 =   1.858e-11;
    const double p2 =  -4.097e-08;
    const double p3 =    4.33e-05;
    const double p4 =    -0.03884;
    const double p5 =      -15.91;
    
    T h2 = h*h;
    T h3 = h2*h;
    T h4 = h2*h2;
    T rho = exp(p1*h4 + p2*h3 + p3*h2 + p4*h + p5);
    
    return rho;
}

template <typename T> AlgebraicVector<T> kep2cart(const AlgebraicVector<T>& kep, const double mu = 398600.4415)
{
    /*member function to convert keplerian  classical element into Earth-Centred inertial reference frame element
     !< keplerian element kep = {a, e, i, RA, PA, TA}
     RA: rigth ascension of ascending node; PA: argument of periapsis; TA: true anomaly; i:orbital inclination
     !> return AlgebraicVector of ECI reference frame res = {x, y, z, dx, dy, dz}*/
    
    //const double mu = 398600.4415;
    
    T p = kep[0] * (1.0 - kep[1] * kep[1]);
    
    // position and velocity in perifocal refererence frame
    AlgebraicVector<T> rm(3), vm(3);
    rm[0] = p*cos(kep[5]) / (1.0 + kep[1] * cos(kep[5]));
    rm[1] = p*sin(kep[5]) / (1.0 + kep[1] * cos(kep[5]));
    rm[2] = 0.0;
    vm[0] = -1.0*sin(kep[5])*sqrt(mu / p);
    vm[1] = (kep[1] + cos(kep[5]))*sqrt(mu / p);
    vm[2] = 0.0;
    
    T cRA = cos(kep[3]);
    T sRA = sin(kep[3]);
    T cPA = cos(kep[4]);
    T sPA = sin(kep[4]);
    T ci = cos(kep[2]);
    T si = sin(kep[2]);
    
    T RR[3][3]; // rotational matrix from perifocal to eci reference frame
    RR[0][0] = cRA*cPA - sRA*ci*sPA;  RR[0][1] = -1.0*cRA*sPA - sRA*ci*cPA; RR[0][2] = sRA*si;
    RR[1][0] = sRA*cPA + cRA*ci*sPA;  RR[1][1] = -1.0*sRA*sPA + cRA*ci*cPA; RR[1][2] = -1.0*cRA*si;
    RR[2][0] = si*sPA;                RR[2][1] = si*cPA;                    RR[2][2] = ci;
    
    AlgebraicVector<T> rr(3), vv(3);
    for (unsigned int i = 0; i<3; i++){
        rr[i] = 0.0;
        vv[i] = 0.0;
        for (unsigned int j = 0; j<3; j++){
            rr[i] = rr[i] + RR[i][j] * rm[j];
            vv[i] = vv[i] + RR[i][j] * vm[j];
        }
    }
    
    AlgebraicVector<T> res(6);
    res[0] = rr[0];
    res[1] = rr[1];
    res[2] = rr[2];
    res[3] = vv[0];
    res[4] = vv[1];
    res[5] = vv[2];
    
    return res;
}


template<typename T> AlgebraicVector<T> cart2kep(const AlgebraicVector<T>& rv, const double mu = 398600.4418)
{
    //const double mu = 398600.4415;
    AlgebraicVector<T> kep(6);
    
    AlgebraicVector<T> rr(3), vv(3);
    for (int i = 0; i < 3; i++)
    {
        rr[i] = rv[i];
        vv[i] = rv[i + 3];
    }
    
    T r = rr.vnorm();
    T v = vv.vnorm();
    AlgebraicVector<T> h = cross(rr, vv);
    //cout << h << endl;
    
    kep[0] = mu / (2.0 * (mu / r - pow(v, 2) / 2.0));
    
    T h1sqr = pow(h[0], 2);
    T h2sqr = pow(h[1], 2);
    
    T RAAN;
    if (cons(h1sqr + h2sqr) == 0.0)
    {
        RAAN = 0.0;
    }
    else
    {
        T sinOMEGA = h[0] / sqrt(h1sqr + h2sqr);
        T cosOMEGA = -1.0*h[1] / sqrt(h1sqr + h2sqr);
        if (cons(cosOMEGA) >= 0.0)
        {
            if (cons(sinOMEGA) >= 0.0)
            {
                RAAN = asin(h[0] / sqrt(h1sqr + h2sqr));
            }
            else
            {
                RAAN = 2.0 * M_PI + asin(h[0] / sqrt(h1sqr + h2sqr));
            }
        }
        else
        {
            if (cons(sinOMEGA) >= 0.0)
            {
                RAAN = acos(-1.0*h[1] / sqrt(h1sqr + h2sqr));
            }
            else
            {
                RAAN = 2.0 * M_PI - acos(-1.0*h[1] / sqrt(h1sqr + h2sqr));
            }
        }
    }
    
    //RAAN = real(RAAN);
    
    AlgebraicVector<T> ee = 1.0 / mu*(cross(vv, h)) - rr / vnorm(rr);
    T e = vnorm(ee);
    
    T i = acos(h[2] / vnorm(h));
    
    kep[1] = e;
    kep[2] = i;
    kep[3] = RAAN;
    
    T omega;
    T theta;
    if (cons(e) <= 1.0e-8 && cons(i) < 1.0e-8)
    {
        e = 0.0;
        omega = atan2(rr[1], rr[0]);
        theta = 0.0;
        kep[4] = omega;
        kep[5] = theta;
        return kep;
    }
    
    if (cons(e) <= 1.0e-8 && cons(i) >= 1.0e-8)
    {
        omega = 0;
        AlgebraicVector<T> P(3), Q(3), W(3);
        P[0] = cos(omega)*cos(RAAN) - sin(omega)*sin(i)*sin(RAAN);
        P[1] = -1.0*sin(omega)*cos(RAAN) - cos(omega)*cos(i)*sin(RAAN);
        P[2] = sin(RAAN)*sin(i);
        Q[0] = cos(omega)*sin(RAAN) + sin(omega)*cos(i)*cos(RAAN);
        Q[1] = -1.0*sin(omega)*sin(RAAN) + cos(omega)*cos(i)*cos(RAAN);
        Q[2] = -1.0*cos(RAAN)*sin(i);
        W[0] = sin(omega)*sin(i);
        W[1] = cos(omega)*sin(i);
        W[2] = cos(i);
        AlgebraicVector<T> rrt = P*rr[0] + Q*rr[1] + W*rr[2];
        theta = atan2(rrt[1], rrt[0]);
        kep[4] = omega;
        kep[5] = theta;
        return kep;
    }
    
    T dotRxE = dot(rr, ee);
    T RxE = vnorm(rr)*vnorm(ee);
    if (abs(cons(dotRxE)) > abs(cons(RxE)) && abs(cons(dotRxE)) - abs(cons(RxE)) < abs(numeric_limits<double>::epsilon()*cons(dotRxE)))
    {
        dotRxE -= numeric_limits<double>::epsilon()*dotRxE;
    }
    theta = acos(dotRxE / RxE);
    
    if (cons(dot(rr, vv)) < 0.0)
    {
        theta = 2.0 * M_PI - theta;
    }
    
    if (cons(i) <= 1.0e-8 && cons(e) >= 1.0e-8)
    {
        i = 0.0;
        omega = atan2(ee[1], ee[0]);
        kep[4] = omega;
        kep[5] = theta;
        return kep;
    }
    
    T sino = rr[2] / r / sin(i);
    T coso = (rr[0] * cos(RAAN) + rr[1] * sin(RAAN)) / r;
    T argLat;
    if (cons(coso) >= 0.0)
    {
        if (cons(sino) >= 0.0)
        {
            argLat = asin(rr[2] / r / sin(i));
        }
        else
        {
            argLat = 2.0 * M_PI + asin(rr[2] / r / sin(i));
        }
    }
    else
    {
        if (cons(sino) >= 0.0)
        {
            argLat = acos((rr[0] * cos(RAAN) + rr[1] * sin(RAAN)) / r);
        }
        else
        {
            argLat = 2.0 * M_PI - acos((rr[0] * cos(RAAN) + rr[1] * sin(RAAN)) / r);
        }
    }
    //argLat = real(argLat);
    omega = argLat - theta;
    
    if (cons(omega) < 0.0)
    {
        omega = omega + 2.0 * M_PI;
    }
    //omega = real(omega);
    
    kep[4] = omega;
    kep[5] = theta;
    
    return kep;
}

template<typename T> AlgebraicVector<T> kep2delaunay(const AlgebraicVector<T>& kep, const double mu)
{
    // Keplerian orbital elements
    T a = kep[0];
    T e = kep[1];
    T i = kep[2];
    T RAAN = kep[3];
    T omega = kep[4];
    T M = kep[5];
    
    // Delaunay elements
    AlgebraicVector<T> delaunay(6);
    delaunay[0] = M; // l
    delaunay[1] = omega; // g
    delaunay[2] = RAAN; // h
    delaunay[3] = sqrt(mu*a); // L
    delaunay[4] = sqrt(1 - pow(e, 2)) * delaunay[3]; // G
    delaunay[5] = cos(i)*delaunay[4]; // H
    
    return delaunay;
}

template<typename T> AlgebraicVector<T> delaunay2kep(const AlgebraicVector<T>& delaunay, const double mu)
{
    // Delaunay elements
    T l = delaunay[0]; // l
    T g = delaunay[1]; // g
    T h = delaunay[2]; // h
    T L = delaunay[3]; // L
    T G = delaunay[4]; // G
    T H = delaunay[5]; // H
    
    // Keplerian orbital elements
    AlgebraicVector<T> kep(6);
    kep[0] = pow(L, 2) / mu; // a
    kep[1] = sqrt(1 - pow((G / L), 2)); // e
    kep[2] = acos(H / G); // i
    kep[3] = h; // RAAN
    kep[4] = g; // omega
    kep[5] = l; // M
    
    return kep;
}


template<typename T> AlgebraicVector<T> kep2hill(const AlgebraicVector<T>& kep, const double mu)

{
    AlgebraicVector<T> hill(6);
    T u, f, p;
    
    p = kep[0]*(1.0 - kep[1]*kep[1]);
    f = kep[5];
    
    hill[4] = sqrt(mu*p);
    hill[0] = p/(1.0 + kep[1]*cos(f));
    hill[1] = f + kep[4];
    hill[2] = kep[3];
    hill[3] = (hill[4]/p)*kep[1]*sin(f);
    hill[5] = hill[4]*cos(kep[2]);
    
    return hill;
}

template<typename T> AlgebraicVector<T> hill2cart(const AlgebraicVector<T>& hill, const double mu)
//  hill[] = {r, th, nu, R, Th, N}
//  cart[] = {x, y, z, vx, vy, vz}
{
    AlgebraicVector<T> u(3);
    AlgebraicVector<T> cart(6);
    
    T r, th, nu, R, Th, ci, si;
    int i;
    
    r  = hill[0];
    th = hill[1];
    nu = hill[2];
    R  = hill[3];
    Th = hill[4];
    ci = hill[5]/hill[4];
    si = sqrt(1.0 - ci*ci);
    
    u[0] = cos(th)*cos(nu) - ci*sin(th)*sin(nu);
    u[1] = cos(th)*sin(nu) + ci*sin(th)*cos(nu);
    u[2] = si*sin(th);
    
    for(i = 0; i < 3; ++i)
        
        cart[i] = r*u[i];
    cart[3] = (R*cos(th) - Th*sin(th)/r)*cos(nu) - (R*sin(th) + Th*cos(th)/r)*sin(nu)*ci;
    cart[4] = (R*cos(th) - Th*sin(th)/r)*sin(nu) + (R*sin(th) + Th*cos(th)/r)*cos(nu)*ci;
    cart[5] = (R*sin(th) + Th*cos(th)/r)*si;
    
    return cart;
}


template<typename T> AlgebraicVector<T> hill2kep(const AlgebraicVector<T>& hill, const double mu)
//  hill[] = {r, th, nu, R, Th, Nu}
//  cart[] = {x, y, z, vx, vy, vz}
{
    AlgebraicVector<T> kep(6);
    
    T r  = hill[0];
    T th = hill[1];
    T nu = hill[2];
    T R  = hill[3];
    T Th = hill[4];
    T Nu = hill[5];
    
    T i = acos(Nu/Th);
    T cs  =   (-1.0 + pow(Th,2)/(mu*r))*cos(th) + (R*Th*sin(th))/mu;
    T ss  =  -((R*Th*cos(th))/mu) + (-1.0 + pow(Th,2)/(mu*r))*sin(th);
    T e = sqrt(cs*cs+ss*ss);
    T p = Th*Th/mu;
    T costrue = 1.0/e*(p/r-1.0);
    T f = acos(costrue);
    
    if (cons(R)<0.0) {
        f = 2.0*M_PI-f;
    }
    T a = p/(1-e*e);
    
    kep[0] = a;
    kep[1] = e;
    kep[2] = i;
    kep[3] = nu;
    kep[4] = th-f;
    kep[5] = f;
    
    return kep;
    
}

template<typename T> AlgebraicVector<T> osculating2meanHill(AlgebraicVector<T> hillOsc, double mu, double J2, double rE)
{
    
    
    // Mean Delaunay elements
    const T r   = hillOsc[0]; // l
    const T th  = hillOsc[1]; // g
    const T nu  = hillOsc[2]; // h
    const T R   = hillOsc[3]; // L
    const T Th  = hillOsc[4]; // G
    const T Nu  = hillOsc[5]; // H
    
    T ci = Nu/Th;
    T si = sqrt(1.0-ci*ci);
    T cs  =   (-1.0 + pow(Th,2)/(mu*r))*cos(th) + (R*Th*sin(th))/mu;
    T ss  =  -((R*Th*cos(th))/mu) + (-1.0 + pow(Th,2)/(mu*r))*sin(th);
    T e = sqrt(cs*cs+ss*ss);
    T eta  = sqrt(1.0-e*e);
    
    T beta = 1.0/(1.0+eta);
    
    T p = Th*Th/mu;
    T costrue = 1/e*(p/r-1);
    
    T f = acos(costrue);
    
    if (cons(R)<0.0) {
        f = 2.0*M_PI-f;
    }
    
    T M = true2meanAnomaly(f,e);
    
    T phi  = f - M;
    
    const T rMean = r +  ((pow(rE,2)*beta*J2)/(2.*r) - (3*pow(rE,2)*beta*J2*pow(si,2))/(4.*r) +
                                 (pow(rE,2)*eta*J2*pow(mu,2)*r)/pow(Th,4) - (3*pow(rE,2)*eta*J2*pow(mu,2)*r*pow(si,2))/(2.*pow(Th,4)) +
                                 (pow(rE,2)*J2*mu)/(2.*pow(Th,2)) - (pow(rE,2)*beta*J2*mu)/(2.*pow(Th,2)) -
                                 (3.*pow(rE,2)*J2*mu*pow(si,2))/(4.*pow(Th,2)) + (3*pow(rE,2)*beta*J2*mu*pow(si,2))/(4.*pow(Th,2)) -
                                 (pow(rE,2)*J2*mu*pow(si,2)*cos(2*th))/(4.*pow(Th,2)));
    
    
    const T thMean = th + ((-3.*pow(rE,2)*J2*pow(mu,2)*phi)/pow(Th,4) + (15.*pow(rE,2)*J2*pow(mu,2)*phi*pow(si,2))/(4.*pow(Th,4)) -
                                  (5.*pow(rE,2)*J2*mu*R)/(2.*pow(Th,3)) - (pow(rE,2)*beta*J2*mu*R)/(2.*pow(Th,3)) +
                                  (3.*pow(rE,2)*J2*mu*R*pow(si,2))/pow(Th,3) + (3.*pow(rE,2)*beta*J2*mu*R*pow(si,2))/(4.*pow(Th,3)) -
                                  (pow(rE,2)*beta*J2*R)/(2.*r*Th) + (3.*pow(rE,2)*beta*J2*R*pow(si,2))/(4.*r*Th) +
                                  (-(pow(rE,2)*J2*mu*R)/(2.*pow(Th,3)) + (pow(rE,2)*J2*mu*R*pow(si,2))/pow(Th,3))*cos(2.*th) +
                                  (-(pow(rE,2)*J2*pow(mu,2))/(4.*pow(Th,4)) + (5.*pow(rE,2)*J2*pow(mu,2)*pow(si,2))/(8.*pow(Th,4)) +
                                   (pow(rE,2)*J2*mu)/(r*pow(Th,2)) - (3.*pow(rE,2)*J2*mu*pow(si,2))/(2.*r*pow(Th,2)))*sin(2.*th));
    
    const T nuMean = nu + ((3.*pow(rE,2)*ci*J2*pow(mu,2)*phi)/(2.*pow(Th,4)) + (3.*pow(rE,2)*ci*J2*mu*R)/(2.*pow(Th,3)) +
                                  (pow(rE,2)*ci*J2*mu*R*cos(2.*th))/(2.*pow(Th,3)) +
                                  ((pow(rE,2)*ci*J2*pow(mu,2))/(4.*pow(Th,4)) - (pow(rE,2)*ci*J2*mu)/(r*pow(Th,2)))*sin(2.*th));
    
    
    const T RMean = R  + (-(pow(rE,2)*beta*J2*R)/(2.*pow(r,2)) + (3.*pow(rE,2)*beta*J2*R*pow(si,2))/(4.*pow(r,2)) -
                                 (pow(rE,2)*eta*J2*pow(mu,2)*R)/(2.*pow(Th,4)) + (3.*pow(rE,2)*eta*J2*pow(mu,2)*R*pow(si,2))/(4.*pow(Th,4)) +
                                 (pow(rE,2)*J2*mu*pow(si,2)*sin(2.*th))/(2.*pow(r,2)*Th));
    
    
    const T ThMean = Th  + (((pow(rE,2)*J2*pow(mu,2)*pow(si,2))/(4.*pow(Th,3)) - (pow(rE,2)*J2*mu*pow(si,2))/(r*Th))*cos(2.*th) -
                                   (pow(rE,2)*J2*mu*R*pow(si,2)*sin(2.*th))/(2.*pow(Th,2)));
    
    const T NuMean = Nu +  0.;
    
    AlgebraicVector<T> hillMean(6);
    
    
    hillMean[0] = rMean;
    hillMean[1] = thMean;
    hillMean[2] = nuMean;
    hillMean[3] = RMean;
    hillMean[4] = ThMean;
    hillMean[5] = NuMean;
    
    return hillMean;
}


template<typename T> AlgebraicVector<T> mean2osculatingHill(AlgebraicVector<T> hillMean, double mu, double J2, double rE)
{
    
    // Mean Delaunay elements
    const T r  = hillMean[0]; // l
    const T th = hillMean[1]; // g
    const T nu = hillMean[2]; // h
    const T R  = hillMean[3]; // L
    const T Th = hillMean[4]; // G
    const T Nu = hillMean[5]; // H
    
    T ci = Nu/Th;
    T si = sqrt(1.0-ci*ci);
    T cs  =   (-1.0 + pow(Th,2)/(mu*r))*cos(th) + (R*Th*sin(th))/mu;
    T ss  =  -((R*Th*cos(th))/mu) + (-1.0 + pow(Th,2)/(mu*r))*sin(th);
    T e = sqrt(cs*cs+ss*ss);
    T eta  = sqrt(1.0-e*e);
    T beta = 1.0/(1.0+eta);
    
    T p = Th*Th/mu;
    T costrue = 1/e*(p/r-1);
    
    T f = acos(costrue);
    
    if (cons(R)<0.0) {
        f = 2.0*M_PI-f;
    }
    
    T M = true2meanAnomaly(f,e);
    
    T phi  = f - M;
    
    T rOsc = r -  ((pow(rE,2)*beta*J2)/(2.*r) - (3*pow(rE,2)*beta*J2*pow(si,2))/(4.*r) +
                          (pow(rE,2)*eta*J2*pow(mu,2)*r)/pow(Th,4) - (3*pow(rE,2)*eta*J2*pow(mu,2)*r*pow(si,2))/(2.*pow(Th,4)) +
                          (pow(rE,2)*J2*mu)/(2.*pow(Th,2)) - (pow(rE,2)*beta*J2*mu)/(2.*pow(Th,2)) -
                          (3.*pow(rE,2)*J2*mu*pow(si,2))/(4.*pow(Th,2)) + (3*pow(rE,2)*beta*J2*mu*pow(si,2))/(4.*pow(Th,2)) -
                          (pow(rE,2)*J2*mu*pow(si,2)*cos(2*th))/(4.*pow(Th,2)));
    
    
    T thOsc = th - ((-3.*pow(rE,2)*J2*pow(mu,2)*phi)/pow(Th,4) + (15.*pow(rE,2)*J2*pow(mu,2)*phi*pow(si,2))/(4.*pow(Th,4)) -
                           (5.*pow(rE,2)*J2*mu*R)/(2.*pow(Th,3)) - (pow(rE,2)*beta*J2*mu*R)/(2.*pow(Th,3)) +
                           (3.*pow(rE,2)*J2*mu*R*pow(si,2))/pow(Th,3) + (3.*pow(rE,2)*beta*J2*mu*R*pow(si,2))/(4.*pow(Th,3)) -
                           (pow(rE,2)*beta*J2*R)/(2.*r*Th) + (3.*pow(rE,2)*beta*J2*R*pow(si,2))/(4.*r*Th) +
                           (-(pow(rE,2)*J2*mu*R)/(2.*pow(Th,3)) + (pow(rE,2)*J2*mu*R*pow(si,2))/pow(Th,3))*cos(2.*th) +
                           (-(pow(rE,2)*J2*pow(mu,2))/(4.*pow(Th,4)) + (5.*pow(rE,2)*J2*pow(mu,2)*pow(si,2))/(8.*pow(Th,4)) +
                            (pow(rE,2)*J2*mu)/(r*pow(Th,2)) - (3.*pow(rE,2)*J2*mu*pow(si,2))/(2.*r*pow(Th,2)))*sin(2.*th));
    
    T nuOsc = nu - ((3.*pow(rE,2)*ci*J2*pow(mu,2)*phi)/(2.*pow(Th,4)) + (3.*pow(rE,2)*ci*J2*mu*R)/(2.*pow(Th,3)) +
                           (pow(rE,2)*ci*J2*mu*R*cos(2.*th))/(2.*pow(Th,3)) +
                           ((pow(rE,2)*ci*J2*pow(mu,2))/(4.*pow(Th,4)) - (pow(rE,2)*ci*J2*mu)/(r*pow(Th,2)))*sin(2.*th));
    
    
    T ROsc = R  - (-(pow(rE,2)*beta*J2*R)/(2.*pow(r,2)) + (3.*pow(rE,2)*beta*J2*R*pow(si,2))/(4.*pow(r,2)) -
                          (pow(rE,2)*eta*J2*pow(mu,2)*R)/(2.*pow(Th,4)) + (3.*pow(rE,2)*eta*J2*pow(mu,2)*R*pow(si,2))/(4.*pow(Th,4)) +
                          (pow(rE,2)*J2*mu*pow(si,2)*sin(2.*th))/(2.*pow(r,2)*Th));
    
    
    T ThOsc = Th  - (((pow(rE,2)*J2*pow(mu,2)*pow(si,2))/(4.*pow(Th,3)) - (pow(rE,2)*J2*mu*pow(si,2))/(r*Th))*cos(2.*th) -
                            (pow(rE,2)*J2*mu*R*pow(si,2)*sin(2.*th))/(2.*pow(Th,2)));
    
    T NuOsc = Nu +  0.;
    
    AlgebraicVector<T> hillOsc(6);
    
    hillOsc[0] = rOsc;
    hillOsc[1] = thOsc;
    hillOsc[2] = nuOsc;
    hillOsc[3] = ROsc;
    hillOsc[4] = ThOsc;
    hillOsc[5] = NuOsc;
    
    return hillOsc;
}

// Averaged J2
template<typename T> AlgebraicVector<T> averagedJ2rhs(AlgebraicVector<T> x, const double mu, const double J2, const double rE)
{
    // Delaunay elements
    T l = x[0]; // l
    T g = x[1]; // g
    T h = x[2]; // h
    T L = x[3]; // L
    T G = x[4]; // G
    T H = x[5]; // H
    
    T eta = G / L; // sqrt(1.0 - pow(e, 2));
    T ci = H / G; // cos(i);
    T si = sin(acos(ci));
    
    T dldt = pow(mu, 2) / pow(L, 3)
    + ((3.0 * J2*pow(rE, 2)*pow(mu, 4)) / (2.0*pow(L, 7)*pow(eta, 3)) -
              (9.0 * J2*pow(si, 2)*pow(rE, 2)*pow(mu, 4)) / (4.0*pow(L, 7)*pow(eta, 3)));
    
    T dgdt = ((3. * J2*pow(rE, 2)*pow(mu, 4)) / (2.*pow(L, 7)*pow(eta, 4)) -
                     (9. * J2*pow(si, 2)*pow(rE, 2)*pow(mu, 4)) / (4.*pow(L, 7)*pow(eta, 4)) +
                     (3. * pow(ci, 2)*J2*pow(rE, 2)*pow(mu, 4)) / (2.*G*pow(L, 6)*pow(eta, 3)));
    
    T dhdt = (-(3. * pow(ci, 2)*J2*pow(rE, 2)*pow(mu, 4)) / (2.*H*pow(L, 6)*pow(eta, 3)));
    
    AlgebraicVector<T> ff(6);
    ff[0] = dldt;
    ff[1] = dgdt;
    ff[2] = dhdt;
    ff[3] = 0;
    ff[4] = 0;
    ff[5] = 0;
    
    return ff;
    
}

template<typename T> AlgebraicVector<T> propJ2An(AlgebraicVector<T> xx0, T tof, double mu, double rE, double J2)
{
    
    AlgebraicVector<T> kep0(6), kep0Mean(6), kepfMean(6), del0Mean(6), delfMean, hill0(6), hill0Mean(6), hillfMean(6), hillf(6), kepf(6), xxf(6);
    
    
    kep0 = cart2kep(xx0, mu); //-> convert true to mean anomaly!
    
    // trasnform keplerian elements to Hill
    
    hill0 = kep2hill(kep0, mu);
    
    // from osculating to mean
    hill0Mean = osculating2meanHill(hill0, mu, J2, rE);
    
    kep0Mean = hill2kep(hill0Mean, mu);
    T meanAnomaly = true2meanAnomaly(kep0Mean[5], kep0Mean[1]);
    kep0Mean[5] = meanAnomaly;
    
    del0Mean = kep2delaunay(kep0Mean, mu);
    
    delfMean = averagedJ2rhs(del0Mean, mu, J2, rE);
    delfMean = delfMean*tof+del0Mean;
    
    kepfMean = delaunay2kep(delfMean, mu);
    
    T trueAnomaly = mean2trueAnomaly(kepfMean[5], kepfMean[1]);
    kepfMean[5] = trueAnomaly;
    
    hillfMean = kep2hill(kepfMean, mu);
    
    // transform mean to osculating
    hillf = mean2osculatingHill(hillfMean, mu, J2, rE);
    
    // transform keplerian elements to xxf
    xxf = hill2cart(hillf, mu);
    
    return xxf;
    
}

template<typename T> AlgebraicVector<T> propKepAn(AlgebraicVector<T> xx0, T t, double mu){
    
    int ord = DACE::DA::getMaxOrder();
    AlgebraicVector<T> rr0(3), vv0(3), xxf(6);
    for (int i=0; i<3; i++) {
        rr0[i] = xx0[i];
        vv0[i] = xx0[i+3];
    }
    AlgebraicVector<T> hh = cross(rr0,vv0);
    T h = vnorm(hh);
    T r0 = vnorm(rr0);
    T v0 = vnorm(vv0);
    
    T a = mu/(2*mu/r0 -v0*v0);
    T p = h*h/mu;
    T sigma0 = dot(rr0,vv0)/sqrt(mu);
    
    double tol = 1.0;
    int iter = 0;
    double scl = 1e-4;
    
    T F, Ft, G, Gt;
    
    if (cons(a)>0) {
        
        T MmM0 = t * sqrt(mu/a/a/a);
        T EmE0 = cons(MmM0);
        
        while ( tol>1e-13 || iter < ord + 1) {
            iter++;
            T fx0 = -(MmM0) + (EmE0) + (sigma0)/sqrt((a))*(1 - cos((EmE0))) - (1-(r0)/(a)) * sin((EmE0));
            T fxp = 1 + (sigma0)/sqrt((a)) * sin((EmE0)) - (1-(r0)/(a)) * cos((EmE0));
            tol = abs(cons(fx0/fxp));
            EmE0 = EmE0 - fx0/fxp;
        }
        
        T theta = 2*atan2((sqrt(a*p)*tan(EmE0/2)),( r0 + sigma0*sqrt(a)*tan(EmE0/2)));
        T r = p*r0 / (r0 + (p-r0)*cos(theta) - sqrt(p)*sigma0*sin(theta));
        
        //{compute the Lagrangian coefficients}
        F = 1 - a/r0 * (1 - cos(EmE0)) ;
        G = a*sigma0/sqrt(mu)*(1 - cos(EmE0)) + r0 * sqrt(a/mu) * sin(EmE0);
        Ft = - sqrt(mu*a)/(r*r0) * sin(EmE0);
        Gt = 1 - a/r * (1-cos(EmE0));
    }
    else {
        T NmN0 = t*sqrt(mu/(-a)/(-a)/(-a));
        T HmH0 = 0.0;
        
        while(tol>1e-14 || iter < ord + 1){
            iter ++;
            T fx0 = - (NmN0) - (HmH0) + (sigma0)/sqrt((-a)) * (-1 + cosh((HmH0))) + (1-(r0)/(a)) * sinh((HmH0)) ;
            T fxp = -1 + (sigma0)/sqrt((-a))*sinh((HmH0)) + (1-(r0)/(a))*cosh((HmH0));
            tol = abs(cons(fx0/fxp));
            HmH0 = HmH0 - fx0/fxp;
        }

        for( iter=0; iter<ord; iter++){
            T fx0 = - (NmN0) - HmH0 + (sigma0)/sqrt((-a)) * (-1 + cosh(HmH0)) + (1-(r0)/(a)) * sinh(HmH0) ;
            T fxp = -1 + (sigma0)/sqrt((-a))*sinh(HmH0) + (1-(r0)/(a))*cosh(HmH0);
            HmH0 = HmH0 - fx0/fxp;
        }
        
        //{DA expansion of HmH0 parameter}
        double Htemp, DE = 1.0;
        
        F = 1 - a/r0 * (1 - cosh(HmH0));
        G = a*sigma0/sqrt(mu)*(1 - cosh(HmH0)) + r0 * sqrt(-a/mu) * sinh(HmH0);
        AlgebraicVector<T> rv_temp(3);
        for (int i=0; i<3; i++) {
            rv_temp[i] = F * rr0[i] + G * vv0[i];
        }
        T r = vnorm(rv_temp);
        Ft = - sqrt(mu*(-a))/(r*r0) * sinh(HmH0);
        Gt = 1 - a/r*(1-cosh(HmH0));
    }
    
    for (int i=0 ; i<3; i++) {
        xxf[i] = F * rr0[i] + G * vv0[i];
        xxf[i+3] = Ft * rr0[i] + Gt * vv0[i];
    }
    return xxf;
}

template<typename T> AlgebraicVector<T> rhsp( const AlgebraicVector<T>& xx, double t )
{
    // physical parameters
    double mu = 398600.4418; //{[km^3/s^2]}
	double rE = 6378.137; //{km}
	const double J2 = 1.08262668e-3; //{-}
	const double J3 = -2.53648e-6; //{-}
	const double J4 = -1.6233e-6; //{-}
    
    // scaling quantitities
    const double Lsc = rE;
    const double Vsc = sqrt(mu/rE);
    const double Tsc = Lsc/Vsc;
    const double km = 1000.0; //{-}

    mu = mu/Lsc/Lsc/Lsc*Tsc*Tsc;
    rE = rE/Lsc;
    
    // spacecraft properties
	const double A = 10.0; //{-}
    const double CD = 2.2; //{-}
    const double m = 1000.0; //{-}
    
    AlgebraicVector<T> pos(3), res(6), vel(3);
    
    pos[0] = xx[0]; pos[1] = xx[1]; pos[2] = xx[2];
    vel[0] = xx[3]; vel[1] = xx[4]; vel[2] = xx[5];
    
    T r = pos.vnorm();
    T v = vel.vnorm();
    
    res[0] = xx[3];
    res[1] = xx[4];
    res[2] = xx[5];

    T x = pos[0];
	T y = pos[1];
	T z = pos[2];
	
	T mur3 = mu/r/r/r;
	T z2r2 = z/r*z/r;

	T J2rEr = 1.5*J2*rE/r*rE/r;
	T mur7Er3 = 5./2.*mur3*J3*rE/r*rE/r*rE/r/r;
	T mur7Er4 = 15./8.*mur3*J4*rE/r*rE/r*rE/r*rE/r;

	res[3] = -pos[0] * mur3 * (1. + J2rEr*(1. - 5.*z2r2));
	res[4] = -pos[1] * mur3 * (1. + J2rEr*(1. - 5.*z2r2));
	res[5] = -pos[2] * mur3 * (1. + J2rEr*(3. - 5.*z2r2));
    
    /*
    T h = (r-rE)*Lsc;
	T rho = atmdensity(h);
    T drag = -(1./2.*rho*CD*m/A*v*v)/km; // to compute it in km
    drag = drag/Lsc*Tsc*Tsc;             // into scaled
    */
    T drag = 0;
    
    res[3] = res[3] + (mur7Er3*x*z*(7.*z2r2-3.) + mur7Er4*x*(1.-14.*z2r2+21.*z2r2*z2r2) + drag*vel[0]/v);
	res[4] = res[4] + (mur7Er3*y*z*(7.*z2r2-3.) + mur7Er4*y*(1.-14.*z2r2+21.*z2r2*z2r2) + drag*vel[1]/v);
	res[5] = res[5] + (mur7Er3*r*r*(3./5. - 6.*z2r2+7.*z2r2*z2r2) +  mur7Er4*z*(5.-70./3.*z2r2+21.*z2r2*z2r2) + drag*vel[2]/v);
	
    return res;
    
}

template<typename T> AlgebraicVector<T> rhsm( const AlgebraicVector<T>& xx, double t )
{
    // physical parameters
    double mu = 398600.4418; //{[km^3/s^2]}
    double rE = 6378.137; //{km}
    const double J2 = 1.08262668e-3; //{-}
    const double J3 = -2.53648e-6; //{-}
    const double J4 = -1.6233e-6; //{-}
    
    // scaling quantitities
    const double Lsc = rE;
    const double Vsc = sqrt(mu/rE);
    const double Tsc = Lsc/Vsc;
    const double km = 1000.0; //{-}

    mu = mu/Lsc/Lsc/Lsc*Tsc*Tsc;
    rE = rE/Lsc;
    
    // spacecraft properties
    const double A = 10.0; //{-}
    const double CD = 2.2; //{-}
    const double m = 1000.0; //{-}
    
    AlgebraicVector<T> pos(3), res(6), vel(3);
    
    pos[0] = xx[0]; pos[1] = xx[1]; pos[2] = xx[2];
    vel[0] = xx[3]; vel[1] = xx[4]; vel[2] = xx[5];
    
    T r = pos.vnorm();
    T v = vel.vnorm();
    
    res[0] = xx[3];
    res[1] = xx[4];
    res[2] = xx[5];

    T x = pos[0];
    T y = pos[1];
    T z = pos[2];
    
    T mur3 = mu/r/r/r;
    T z2r2 = z/r*z/r;

    T J2rEr = 1.5*J2*rE/r*rE/r;
    T mur7Er3 = 5./2.*mur3*J3*rE/r*rE/r*rE/r/r;
    T mur7Er4 = 15./8.*mur3*J4*rE/r*rE/r*rE/r*rE/r;

    res[3] = -pos[0] * mur3 * (1. + J2rEr*(1. - 5.*z2r2));
    res[4] = -pos[1] * mur3 * (1. + J2rEr*(1. - 5.*z2r2));
    res[5] = -pos[2] * mur3 * (1. + J2rEr*(3. - 5.*z2r2));
    
   /*
    T h = (r-rE)*Lsc;
    T rho = atmdensity(h);
    T drag = -(1./2.*rho*CD*m/A*v*v)/km; // to compute it in km
    drag = drag/Lsc*Tsc*Tsc;             // into scaled
    */
    T drag = 0;
    
    res[3] = res[3] + (mur7Er3*x*z*(7.*z2r2-3.) + mur7Er4*x*(1.-14.*z2r2+21.*z2r2*z2r2) + drag*vel[0]/v);
    res[4] = res[4] + (mur7Er3*y*z*(7.*z2r2-3.) + mur7Er4*y*(1.-14.*z2r2+21.*z2r2*z2r2) + drag*vel[1]/v);
    res[5] = res[5] + (mur7Er3*r*r*(3./5. - 6.*z2r2+7.*z2r2*z2r2) +  mur7Er4*z*(5.-70./3.*z2r2+21.*z2r2*z2r2) + drag*vel[2]/v);
    
    return -res;
    
}


template<typename T> AlgebraicMatrix<T> rtn(const AlgebraicVector<T> xx)
{
    AlgebraicVector<T> rr(3), vv(3), tt(3), nn(3);
    AlgebraicMatrix<T> A(3,3);
    
    rr[0] = xx[0]; rr[1] = xx[1]; rr[2] = xx[2];
    vv[0] = xx[3]; vv[1] = xx[4]; vv[2] = xx[5];

    rr = rr/vnorm(rr);
    vv = vv/vnorm(vv);
    
    nn = cross(rr,vv);
    nn = nn/vnorm(nn);
    
    tt = cross(nn,rr);
    tt = tt/vnorm(tt);
    
    A.at(0,0) = rr[0]; A.at(0,1) = rr[1]; A.at(0,2) = rr[2];
    A.at(1,0) = tt[0]; A.at(1,1) = tt[1]; A.at(1,2) = tt[2];
    A.at(2,0) = nn[0]; A.at(2,1) = nn[1]; A.at(2,2) = nn[2];
    
    return A;
}


DA findTCA(const AlgebraicVector<DA> xx, const int nvar)
{
    AlgebraicVector<DA> rr(3), vv(3), MAPD(nvar), MAPI(nvar), dx(nvar);
    
    //relative velocities
    rr[0] = xx[0]; rr[1] = xx[1]; rr[2] = xx[2];
    vv[0] = xx[3]; vv[1] = xx[4]; vv[2] = xx[5];

    // we want rr to be orthogonal to vv, so dot(rr,vv) = 0
    DA rvdot = dot(rr,vv);
    double rvdotcons = cons(rvdot);
    
    // MAPD = (dot(rr,vv), dxx) = f(dxx, dt)
    MAPD[0] = rvdot-rvdotcons;
    for (int i = 1; i<nvar ; i++) {
        MAPD[i] = DA(i);}
    
    // MAPI = (dxx, dt) = f(dot(rr,vv), dxx)
    MAPI = MAPD.invert();
    
    // we need to evaluate the map in -cons(dot(rr,vv)), dxx
    dx[0] = - rvdotcons + 0*DA(1);
    for (int i = 1; i<nvar ; i++){
        dx[i] = DA(i);}
    MAPD = MAPI.eval(dx);
    
    // dt is the last row of the MAPD
    DA tca = MAPD[nvar-1];
    
    return tca;
}

AlgebraicVector<DA> EncounterPlane(AlgebraicVector<DA>VPrim,AlgebraicVector<DA>VSec, int axis){
    
    /*
    EncounterPlane.cpp -  Obtain the encounter plane of both satellites from
    the expression of their velocities

    DESCRIPTION:
     ---------------
    Obtain the encounter plane of both satellites from the expression of
    their velocities expressed in the same reference system

    INPUT:
    ---------------
    VPrim   Velocity of the maneouvrable satellite   [km/s]
    VSec    Velocity of the unmaneouvrable satellite [km/s]
    axis    which element of the trihedral to output (1 xi, 2 eta, 3 zeta)

    OUTPUT:
    ---------------
    xi,eta,zeta    Encounter plane unitary vectors  [-]

    REFERENCES:
    ---------------
    (none)
     
    CALLED FUNCTIONS:
    ---------------
    (none)

    AUTHORS:
    ---------------
    - Alvaro Martinez.
    - Carlos Belmonte.

    CHANGELOG:
    ---------------
    01/03/2020, Carlos Belmonte: First version and description.
    --------------------------------------------------------------------
     */

    AlgebraicVector<DA> relVel = VPrim-VSec;
    AlgebraicVector<DA> eta = relVel/relVel.vnorm();
    AlgebraicVector<DA> CrossVel = VSec.cross(VPrim);
    AlgebraicVector<DA> xi = CrossVel/CrossVel.vnorm();
    AlgebraicVector<DA> zeta = xi.cross(eta);
    zeta = zeta/vnorm(zeta);
    
    if(axis==1){
        return xi;
        
    } else if (axis==2){
        return eta;
        
    } else
    {
        return zeta;}
}


AlgebraicMatrix<DA> Bplane(AlgebraicVector<DA> vvs, AlgebraicVector<DA> vvd){
    
    /*
    EncounterPlane.cpp -  Obtain the encounter plane of both satellites from
    the expression of their velocities

    DESCRIPTION:
     ---------------
    Obtain the encounter plane of both satellites from the expression of
    their velocities expressed in the same reference system

    INPUT:
    ---------------
    vvs   Velocity of the maneouvrable satellite   [km/s]
    vvd    Velocity of the unmaneouvrable satellite [km/s]

    OUTPUT:
    ---------------
    Rotation matrix   [-]

    REFERENCES:
    ---------------
    (none)
     
    CALLED FUNCTIONS:
    ---------------
    (none)

    AUTHORS:
    ---------------
    - Alvaro Martinez.
    - Carlos Belmonte.

    CHANGELOG:
    ---------------
    01/03/2020, Carlos Belmonte: First version and description.
    --------------------------------------------------------------------
     */

    AlgebraicVector<DA> relVel = vvs-vvd;
    AlgebraicVector<DA> eta = relVel/relVel.vnorm();
    AlgebraicVector<DA> CrossVel = vvd.cross(vvs);
    AlgebraicVector<DA> xi = CrossVel/CrossVel.vnorm();
    AlgebraicVector<DA> zeta = xi.cross(eta);
    zeta = zeta/vnorm(zeta);

    AlgebraicMatrix<DA> toBplane(3,3);
    toBplane.at(0,0) = xi[0];     toBplane.at(0,1) = xi[1];     toBplane.at(0,2) = xi[2];
    toBplane.at(1,0) = eta[0];    toBplane.at(1,1) = eta[1];    toBplane.at(1,2) = eta[2];
    toBplane.at(2,0) = zeta[0];   toBplane.at(2,1) = zeta[1];   toBplane.at(2,2) = zeta[2];

    return toBplane;
}


// Rot Mat
AlgebraicMatrix<DA> RotMat(AlgebraicVector<DA>Vx1,AlgebraicVector<DA>Vy1,AlgebraicVector<DA>Vz1,AlgebraicVector<DA>Vx2,AlgebraicVector<DA>Vy2,AlgebraicVector<DA>Vz2){
    
    /*
    % RotMat.cpp -  Obtain the Covariance Matrix of both bodies involved
    % in the collision, expressed in the encounter plane reference frame
    %
    %% DESCRIPTION:
    %  ---------------
    %  It obtains the rotation matrix to change a vector or a matrix expressed
    %  in the system of reference 1 to express it in the system 2.
    %
    %% INPUT:
    %  ---------------
    %  Vx1,Vy1,Vz1 Unitary vectors of the system 1.
    %  Vx2,Vy2,Vz2 Unitary vectors of the system 2.
    %
    %% OUTPUT:
    %  ---------------
    %  R12     Rotation matrix to change from system of reference 1 to 2.
    %
    %% REFERENCES:
    %  ---------------
    %  (none)
    %
    %% CALLED FUNCTIONS:
    %  ---------------
    %  (none)
    %
    %% AUTHORS:
    %  ---------------
    %  - Alvaro Martinez.
    %  - Carlos Belmonte.
    %
    %% CHANGELOG:
    %  ---------------
    %  26/02/2020, Alvaro Martinez: First version and description.
    %
    %  -------------------------------------------------------------------- %
    */
    
    AlgebraicMatrix<DA> R12(3);
    R12.at(0,0) = Vx1.dot(Vx2);
    R12.at(0,1) = Vy1.dot(Vx2);
    R12.at(0,2) = Vz1.dot(Vx2);
    R12.at(1,0) = Vx1.dot(Vy2);
    R12.at(1,1) = Vy1.dot(Vy2);
    R12.at(1,2) = Vz1.dot(Vy2);
    R12.at(2,0) = Vx1.dot(Vz2);
    R12.at(2,1) = Vy1.dot(Vz2);
    R12.at(2,2) = Vz1.dot(Vz2);
    
    return R12;
}


int factorial(int n)
{
    // single line to find factorial
    return (n==1 || n==0) ? 1: n * factorial(n - 1);
}


DA CollProbChan(AlgebraicMatrix<DA> Cb, double r12, AlgebraicVector<DA> posRel){
    
    DA xi_exp = posRel[0];
    DA zeta_exp = posRel[1];
    
    DA sigma_xi = sqrt(Cb.at(0,0));
    DA sigma_zeta = sqrt(Cb.at(1,1));
    DA rho = Cb.at(0,1)/(sigma_xi*sigma_zeta);
    
    DA u = sqr(r12)/(sigma_xi*sigma_zeta*sqrt(1.0-sqr(rho)));
    DA v = (sqr((xi_exp/sigma_xi)) + sqr((zeta_exp/sigma_zeta)) - 2.0*rho*(xi_exp*zeta_exp/(sigma_xi*sigma_zeta)))/(1.0-sqr(rho));
    
    int n = 2;
    DA SecondLoop = 0;
    
    for (int m = 0; m <= n; m++){
        
        DA FirstLoop = 0;
        for (int k = 0; k <= m; k++){
            FirstLoop = FirstLoop + pow(u,k)/(pow(2.0,k) * factorial(k));
        }
        SecondLoop = SecondLoop + pow(v,m)/(pow(2.0,m) * factorial(m)) * (1.0 - exp(-u/2.0)*FirstLoop);
    }
    DA CollisionProbability = exp(-v/2.0)*SecondLoop;
    return CollisionProbability;

}

DA CollProbAlfano(AlgebraicMatrix<DA> Cb, double r12, AlgebraicVector<DA> posRel){

    DA xm, zm, sigmax, sigmaz, CollProb;
    AlgebraicVector<DA> xx(2,2);
    DA x = posRel[0];
    DA z = posRel[1];
    
    DA sigma_x = sqrt(Cb.at(0,0));
    DA sigma_z = sqrt(Cb.at(1,1));
    DA rho = Cb.at(0,1)/(sigma_x*sigma_z);
    DA theta = 1.0/2.0*atan(2.0*rho*sigma_x*sigma_z/(sigma_x*sigma_x-sigma_z*sigma_z));
    
    if (cons(sigma_z)>cons(sigma_x)){theta = theta + atan(1.0)*2.0;}
    AlgebraicMatrix<DA> R(2,2), C(2,2);
    
    R.at(0,0) = cos(theta); R.at(0,1) = sin(theta);
    R.at(1,0) = -sin(theta); R.at(1,1) = cos(theta);
    
    C = R*Cb*R.transpose();
    xx = R*posRel;
    xm = xx[0];
    zm = xx[1];

    sigmax = sqrt(C.at(0,0));
    sigmaz = sqrt(C.at(1,1));
    
   // int n = floor(cons(5.0*r12)/min(cons(sqrt(sigmaz)),cons(vnorm(posRel))));
    int n = 30;
    
    for (int i = 1; i<n+1; i++)
    {
        DA aux1 = ( zm + 2.0*r12/n*sqrt((n-i)*i) )/(sigmaz*sqrt(2.0));
        DA aux2 = (-zm + 2.0*r12/n*sqrt((n-i)*i) )/(sigmaz*sqrt(2.0));
        DA aux3 = -pow((r12*(2.0*i-n)/n + xm ),2 )/(2.0*sigmax*sigmax);
        CollProb = CollProb + (erf(aux1)+erf(aux2))*exp(aux3);

    }
    CollProb = r12*2.0/(sqrt(8.0*atan(1.0)*4.0)*sigmax*n)*CollProb;
    
    return CollProb;

}
