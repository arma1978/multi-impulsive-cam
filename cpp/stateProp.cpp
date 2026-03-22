// =========================================================================
// stateProp.cpp
// Single-step DA state propagation executable used by MATLAB workflows.
// Propagates chief/deputy states with Kepler/J2 dynamics from runtime inputs.
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
    int i, j;
    const double mu = 398600.4418; //{[km^3/s^2]}
    const double rE = 6378.137; //{km}
    const double J2 = 1.08262668e-3; //{-}
    
    // transform xx0 to keplerian elements
    double Lsc = rE;
    double Vsc = sqrt(mu/rE);
    double Tsc = Lsc/Vsc;
    double muSc = mu/Lsc/Lsc/Lsc*Tsc*Tsc;
    
    // DA initialisation
    int order, nvar;
    order = 1;
    nvar = 7;
    DA::init( order, nvar );
    DA::setEps(1e-30);
    
    // DA vector initialisation
    AlgebraicVector<DA> xd0(6), xs0(6), xdfJ2(6), xsfJ2(6),  xdfKep(6), xsfKep(6);
    // dummy variables for reading
    AlgebraicVector<double> xdum(6);
    double tdum;
    
    // space debris state at t0
    ifstream infile;
    infile.open("runtime/xd0.dat"); // debris initial state
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
    infile.open("runtime/xs0.dat");
	if (infile.is_open())
	{
		infile >> xdum[0];
		infile >> xdum[1];
		infile >> xdum[2];
		infile >> xdum[3];
		infile >> xdum[4];
		infile >> xdum[5];
        infile >> tdum; // nominal time to tca
	}
	else
	{
		cout << "Input file not found!" << endl;
		return 1;
	}
	infile.close();

    // spacecraft is DA
    for (int i = 0; i < 6; i++)
    {
        xs0[i] = xdum[i] + DA(i+1);
    }
    
    // scaling
    xs0[0] = xs0[0]/Lsc;
    xs0[1] = xs0[1]/Lsc;
    xs0[2] = xs0[2]/Lsc;
    xs0[3] = xs0[3]/Vsc;
    xs0[4] = xs0[4]/Vsc;
    xs0[5] = xs0[5]/Vsc;
    t2c  = (tdum+DA(7))/Tsc;

    // propagation using keplerian or J2 dynamics
    xdfJ2  = propJ2An( xd0, t2c, muSc, rE/Lsc, J2); // time dependence of final state on time
    xdfKep = propKepAn(xd0, t2c, muSc); // time dependence of final state on time
    
    xsfJ2  = propJ2An( xs0, t2c, muSc, rE/Lsc, J2); // propagated final state
    xsfKep = propKepAn(xs0, t2c, muSc); // propagated final state

    // check if the final condition (J2) is at TCA
    AlgebraicVector<DA> rr(3), vv(3), xx(6);
    xx = xsfJ2-xdfJ2;
    //relative velocities
    rr[0] = xx[0]; rr[1] = xx[1]; rr[2] = xx[2];
    vv[0] = xx[3]; vv[1] = xx[4]; vv[2] = xx[5];
    cout << "nominal value of dot product" << endl << cons(dot(rr,vv)) << endl;
    cout << "nominal value of relative distance [km]" << endl << cons(vnorm(rr))*Lsc << endl;

    //how to calculate tca
    DA tca = findTCA(xsfJ2-xdfJ2, nvar);
    
    //evaluate the final states at tca (dependence is now dropped)
    AlgebraicVector<DA> dx(nvar);
    for (int i = 0; i < nvar-1; i++)
    {dx[i] = DA(i+1);}
    dx[nvar-1] = tca;
    xdfJ2 = xdfJ2.eval(dx);
    xsfJ2 = xsfJ2.eval(dx);
    
    // check
    xx = xsfJ2-xdfJ2;
    //relative velocities
    rr[0] = xx[0]; rr[1] = xx[1]; rr[2] = xx[2];
    vv[0] = xx[3]; vv[1] = xx[4]; vv[2] = xx[5];

    cout << "value of dot product at TCA" << endl << dot(rr,vv) << endl; // note this is valide for any variation of the final state
    cout << "relative distance at TCA [km]" << endl << vnorm(rr)*Lsc << endl; // check if they are orthogonal
    cout << "TCA [s] from epoch0" << endl << t2c.eval(dx)*Tsc << endl; // check updated TCA
    
    // example of use of matrix -> useful for coll prob calculation
    AlgebraicMatrix<DA> toRTNs(3,3), toRTNd(3,3), fromRTNs(3,3), fromRTNd(3,3), PdECI(3,3), PsECI(3,3);
    AlgebraicMatrix<double> P(3,3);
    
    P.at(0,0) = 0.1/(Lsc*Lsc); P.at(0,1) = 0.0; P.at(0,2) = 0.0;
    P.at(1,0) = 0.0; P.at(1,1) = 1/(Lsc*Lsc); P.at(1,2) = 0.0;
    P.at(2,0) = 0.0; P.at(2,1) = 0.0; P.at(2,2) = 0.1/(Lsc*Lsc);
    
    toRTNs = rtn(xsfJ2);
    toRTNd = rtn(xdfJ2);
    fromRTNd = toRTNd.transpose(); // can be used to map back covariance defined in RTN
    fromRTNs = toRTNs.transpose();
    
    PsECI = fromRTNs*P*fromRTNs.transpose(); // covariance of spacecraft in ECI
    PdECI = fromRTNd*P*fromRTNd.transpose(); // covariance of debris in ECI

	// write the output
    ofstream outfile;
    outfile.open("runtime/xdfKep.dat");
    outfile << setprecision(16);
    for (i = 0; i<3 ; i++) {
        outfile << xdfKep[i]*Lsc<<endl;};
    for (i = 3; i<6 ; i++) {
        outfile << xdfKep[i]*Vsc<<endl;};
    outfile.close();
    
    outfile.open("runtime/xsfKep.dat");
	outfile << setprecision(16);
	for (i = 0; i<3 ; i++) {
	outfile << xsfKep[i]*Lsc<<endl;};
    for (i = 3; i<6 ; i++) {
        outfile << xsfKep[i]*Vsc<<endl;};
    outfile.close();

    outfile.open("runtime/xdfJ2.dat");
    outfile << setprecision(16);
    for (i = 0; i<3 ; i++) {
        outfile << xdfJ2[i]*Lsc<<endl;};
    for (i = 3; i<6 ; i++) {
        outfile << xdfJ2[i]*Vsc<<endl;};
    outfile.close();
    
    outfile.open("runtime/xsfJ2.dat");
    outfile << setprecision(16);
    for (i = 0; i<3 ; i++) {
    outfile << xsfJ2[i]*Lsc<<endl;};
    for (i = 3; i<6 ; i++) {
        outfile << xsfJ2[i]*Vsc<<endl;};
    outfile.close();
}
