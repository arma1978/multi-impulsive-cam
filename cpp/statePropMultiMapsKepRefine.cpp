// =========================================================================
// statePropMultiMapsKepRefine.cpp
// Refinement-stage Kepler multi-map propagation around current maneuver guess.
// Used in iterative convex refinement loops.
// =========================================================================

#define _USE_MATH_DEFINES
#define WITH_ALGEBRAICMATRIX

#include <dace/dace.h>
#include <cmath>
#include <ctime>
#include <fstream>
#include <iomanip>
#include "astroRoutines.h"
using namespace std;
using namespace DACE;

int main( void )
{
    int i, j, k, order, nvar, maxImpNum;
    const double mu = 398600.4418; //{[km^3/s^2]}
    const double rE = 6378.137; //{km}
    const double J2 = 1.08262668e-3; //{-}
    
    // transform xx0 to keplerian elements
    double Lsc = rE;
    double Vsc = sqrt(mu/rE);
    double Tsc = Lsc/Vsc;
    double muSc = mu/Lsc/Lsc/Lsc*Tsc*Tsc;
    double PsR, PsT, PsN, PsRT, PsRN, PsTN;
    double PdR, PdT, PdN, PdRT, PdRN, PdTN;
    double disc;
    double r12;

    ifstream infile;
    infile.open("runtime/inputFullP.dat"); // debris initial state
    if (infile.is_open())
    {
        infile >> order;
        infile >> disc;
        infile >> r12;
        infile >> PsR;
        infile >> PsT;
        infile >> PsN;
        infile >> PsRT;
        infile >> PsRN;
        infile >> PsTN;
        infile >> PdR;
        infile >> PdT;
        infile >> PdN;
        infile >> PdRT;
        infile >> PdRN;
        infile >> PdTN;
        infile >> maxImpNum;
    }
    else
    {
        cout << "Input file not found!" << endl;
        return 1;
    }
    infile.close();
    
    // DA initialisation
    nvar = 7;
    DA::init( order, nvar );
    DA::setEps(1e-30);
    
    // DA vector initialisation
    AlgebraicVector<DA> xd0(6), xs0(6), xdfJ2(6), xsfJ2(6), rr(3), xx(6), vv(3);
    // dummy variables for reading
    AlgebraicVector<double> xdum(6);
    double tdum;
    
    // space debris state at t0
    infile.open("runtime/xd0Kep.dat"); // debris initial state
	if (infile.is_open())
	{
		infile >> xdum[0];
		infile >> xdum[1];
		infile >> xdum[2];
		infile >> xdum[3];
		infile >> xdum[4];
		infile >> xdum[5];
        infile >> tdum; //nominal time to TCA
	}
	else
	{
		cout << "Input file not found!" << endl;
		return 1;
	}
	infile.close();
    
    for (int i = 0; i < 6; i++)
    {
        xd0[i] = xdum[i] + 0*DA(i+1); //fake initialisation as DA
    }
    
    // scaling
    xd0[0] = xd0[0]/Lsc;
    xd0[1] = xd0[1]/Lsc;
    xd0[2] = xd0[2]/Lsc;
    xd0[3] = xd0[3]/Vsc;
    xd0[4] = xd0[4]/Vsc;
    xd0[5] = xd0[5]/Vsc;
    
    //initial state of spacecraft at t0
    DA t2c;
    DA dt = disc +0*DA(7); // five minute interval maybe read this from an input file
    infile.open("runtime/xxs.dat");
    infile >> tdum;
    int N  = min(floor(tdum/cons(dt)),maxImpNum) - 1; // removed one last stepto avoid maneuvering too close to collision!
    
    infile >> xdum[0];
    infile >> xdum[1];
    infile >> xdum[2];
    infile >> xdum[3];
    infile >> xdum[4];
    infile >> xdum[5];
    
    AlgebraicMatrix<double> dvv0(3,N+1);
    for (int i = 0; i < N+1; i++){
        for (int k = 0; k < 3; k++){
            infile >> dvv0.at(k,i);}
    }
    infile.close();
    ofstream outfile;
    outfile.open("runtime/mapsRefine.dat");
    outfile << setprecision(16);
    
    xs0[0] = (xdum[0] + DA(1))/Lsc;
    xs0[1] = (xdum[1] + DA(2))/Lsc;
    xs0[2] = (xdum[2] + DA(3))/Lsc;
    xs0[3] = (xdum[3] + dvv0.at(0,0) + DA(4))/Vsc;
    xs0[4] = (xdum[4] + dvv0.at(1,0) + DA(5))/Vsc;
    xs0[5] = (xdum[5] + dvv0.at(2,0) + DA(6))/Vsc;
    
    // N propagation
    for (int i = 0; i < N; i++){
        xsfJ2 = propKepAn(xs0, dt/Tsc, muSc);
        
        for (int j = 0; j<3 ; j++) {
            outfile << xsfJ2[j]*Lsc<<endl;
            xs0[j] = cons(xsfJ2[j]) + DA(j+1)/Lsc;}

        for (int j = 3; j<6 ; j++) {
            outfile << xsfJ2[j]*Vsc<<endl;
            xs0[j] = cons(xsfJ2[j]) + (dvv0.at(j-3,i+1) + DA(j+1))/Vsc;}
    }
    
    // missing interval to final time
   dt = tdum -N*dt;
    xsfJ2 = propKepAn(xs0, dt/Tsc, muSc);

   for (int j = 0; j<3 ; j++) {
       outfile << xsfJ2[j]*Lsc<<endl;
       xs0[j] = cons(xsfJ2[j]) + DA(j+1)/Lsc;}
   for (int j = 3; j<6 ; j++) {
       outfile << xsfJ2[j]*Vsc<<endl;
       xs0[j] = cons(xsfJ2[j]) + DA(j+1)/Vsc;}

    // final time espansion
    dt = (0.0+DA(7))/Tsc;
    DA tfc = tdum/Tsc+dt;
    
    xsfJ2  = propKepAn( xs0, dt, muSc); // propagated final state
    xdfJ2  = propKepAn( xd0, tfc, muSc); // propagated final state

    DA tca = findTCA(xsfJ2-xdfJ2, nvar);
    
    //evaluate the final states at tca (dependence is now dropped)
    AlgebraicVector<DA> dx(nvar);
    for (int i = 0; i < nvar-1; i++) {
        dx[i] = DA(i+1);}
    dx[nvar-1] = tca;
    xdfJ2 = xdfJ2.eval(dx);
    xsfJ2 = xsfJ2.eval(dx);
    tca   = tfc.eval(dx);

    // states at the encounter: spacececraft
    for (int j = 0; j<3 ; j++) {
        outfile << xsfJ2[j]*Lsc<<endl;
        xsfJ2[j] = xsfJ2[j]*Lsc;
    }
    for (int j = 3; j<6 ; j++) {
        outfile << xsfJ2[j]*Vsc<<endl;
        xsfJ2[j] = xsfJ2[j]*Vsc;
    }
    
    // states at the encounter: debris
    for (int j = 0; j<3 ; j++) {
        outfile << xdfJ2[j]*Lsc<<endl;
        xdfJ2[j] = xdfJ2[j]*Lsc;
    }
    for (int j = 3; j<6 ; j++) {
        outfile << xdfJ2[j]*Vsc<<endl;
        xdfJ2[j] = xdfJ2[j]*Vsc;
    }
    
    // time of closest approach
    outfile << tca*Tsc<<endl;
 
    xx = xsfJ2-xdfJ2; //relative position
    // from now on everything is in km km/s
    
    //relative velocities
    rr[0] = xx[0]; rr[1] = xx[1]; rr[2] = xx[2];
    vv[0] = xx[3]; vv[1] = xx[4]; vv[2] = xx[5];

   // cout << "value of TCA [s]" << endl << tca*Tsc << endl; // note this is valid for any variation of the final state
   // cout << "relative dustance [km]" << endl << vnorm(rr)*Lsc << endl;
    
    // collision probability section
    AlgebraicMatrix<DA> toRTNs(3,3), toRTNd(3,3), fromRTNs(3,3), fromRTNd(3,3), PdECI(3,3), PsECI(3,3);
    AlgebraicMatrix<double> Ps(3,3), Pd(3,3);
    // maybe read this from an input file
    Ps.at(0,0) = PsR;            Ps.at(0,1) = PsRT;            Ps.at(0,2) = PsRN;
    Ps.at(1,0) = PsRT;           Ps.at(1,1) = PsT;             Ps.at(1,2) = PsTN;
    Ps.at(2,0) = PsRN;           Ps.at(2,1) = PsTN;            Ps.at(2,2) = PsN;
    
    Pd.at(0,0) = PdR;            Pd.at(0,1) = PdRT;            Pd.at(0,2) = PdRN;
    Pd.at(1,0) = PdRT;           Pd.at(1,1) = PdT;             Pd.at(1,2) = PdTN;
    Pd.at(2,0) = PdRN;           Pd.at(2,1) = PdTN;            Pd.at(2,2) = PdN;
    
    toRTNs = rtn(xsfJ2);
    toRTNd = rtn(xdfJ2);
    fromRTNd = toRTNd.transpose(); // can be used to map back covariance defined in RTN
    fromRTNs = toRTNs.transpose();
    
    PsECI = fromRTNs*Ps*fromRTNs.transpose(); // covariance of spacecraft in ECI
    PdECI = fromRTNd*Pd*fromRTNd.transpose(); // covariance of debris in ECI

    //velocities
    AlgebraicVector<DA> vvsf(3), vvdf(3);
    vvsf[0] = xsfJ2[3]; vvsf[1] = xsfJ2[4]; vvsf[2] = xsfJ2[5];
    vvdf[0] = xdfJ2[3]; vvdf[1] = xdfJ2[4]; vvdf[2] = xdfJ2[5];

    // rotation matrix to Bplane
    AlgebraicMatrix<DA> toBplane = Bplane(vvsf, vvdf);

    AlgebraicMatrix<DA> Pb_3D = toBplane*(PsECI+PdECI)*toBplane.transpose();
    AlgebraicMatrix<DA> Pb_2D(2);
    Pb_2D.at(0,0) = Pb_3D.at(0,0); Pb_2D.at(0,1) = Pb_3D.at(0,2);
    Pb_2D.at(1,0) = Pb_3D.at(2,0); Pb_2D.at(1,1) = Pb_3D.at(2,2);
    
    AlgebraicVector<DA> rrb_3D = toBplane*rr;
    AlgebraicVector<DA> rrb_2D(2); rrb_2D[0]=rrb_3D[0]; rrb_2D[1] = rrb_3D[2];
    
    AlgebraicMatrix<DA> Nb_2D = Pb_2D.inv();
    AlgebraicMatrix<DA> Nb_3D = Pb_3D.inv();

    DA sqrMahalanobis2D = rrb_2D.dot(Nb_2D*rrb_2D);
    DA sqrMahalanobis3D = rrb_3D.dot(Nb_3D*rrb_3D);
    
    outfile << sqrMahalanobis2D <<endl;
        
    // states at the encounter: debris
    for (int j = 0; j<2 ; j++) {
        outfile << rrb_2D[j] <<endl;}
     
    AlgebraicVector<DA> Nb_2Drrb_2D(2);
    Nb_2Drrb_2D = Nb_2D*rrb_2D;
    for (int j = 0; j<2 ; j++) {
        outfile << Nb_2Drrb_2D[j] <<endl;}

   //elements of the covariance
    outfile << Pb_2D.at(0,0) <<endl;
    outfile << Pb_2D.at(0,1) <<endl;
    outfile << Pb_2D.at(1,0) <<endl;
    outfile << Pb_2D.at(1,1) <<endl;
    outfile.close();

}
