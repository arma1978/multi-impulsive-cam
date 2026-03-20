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
    
    // DA vector initialisation
    AlgebraicVector<double> xd0(6), xs0(6), xdfJ2(6), xsfJ2(6),  xdfKep(6), xsfKep(6);
    // dummy variables for reading
    AlgebraicVector<double> xdum(6);
    double t2c, dt;
    
    // space debris state at t0
    ifstream infile;
    infile.open("runtime/xdTCA.dat"); // debris initial state
	if (infile.is_open())
	{
		infile >> xdum[0];
		infile >> xdum[1];
		infile >> xdum[2];
		infile >> xdum[3];
		infile >> xdum[4];
		infile >> xdum[5];
        infile >> t2c; //nominal time to TCA
        infile >> dt;
	}
	else
	{
		cout << "Input file not found!" << endl;
		return 1;
	}
	infile.close();
    
    // scaling
    xd0[0] = xdum[0]/Lsc;
    xd0[1] = xdum[1]/Lsc;
    xd0[2] = xdum[2]/Lsc;
    xd0[3] = xdum[3]/Vsc;
    xd0[4] = xdum[4]/Vsc;
    xd0[5] = xdum[5]/Vsc;
        
    //initial state of spacecraft at t0
    infile.open("runtime/xsTCA.dat");
	if (infile.is_open())
	{
		infile >> xdum[0];
		infile >> xdum[1];
		infile >> xdum[2];
		infile >> xdum[3];
		infile >> xdum[4];
		infile >> xdum[5];
        infile >> t2c; // nominal time to tca
	}
	else
	{
		cout << "Input file not found!" << endl;
		return 1;
	}
	infile.close();

    // scaling
    xs0[0] = xdum[0]/Lsc;
    xs0[1] = xdum[1]/Lsc;
    xs0[2] = xdum[2]/Lsc;
    xs0[3] = xdum[3]/Vsc;
    xs0[4] = xdum[4]/Vsc;
    xs0[5] = xdum[5]/Vsc;
    
    // propagation using keplerian or J2 dynamics
   // xdfJ2  = propJ2An(xd0, -t2c/Tsc, muSc, rE/Lsc, J2); // time dependence of state on time
    // integrate for time of flight
    xdfJ2 = RK78(6, xd0, 0.0, t2c/Tsc, rhsm);
    //xsfJ2  = propJ2An(xs0, -t2c/Tsc, muSc, rE/Lsc, J2); // propagated final state
    xsfJ2 = RK78(6, xs0, 0.0, t2c/Tsc, rhsm);
    xdfKep = propKepAn(xd0, -t2c/Tsc, muSc); // time dependence of final state on time
    xsfKep = propKepAn(xs0, -t2c/Tsc, muSc); // propagated final state

	// write the output
    ofstream outfile;
    outfile.open("runtime/xd0Kep.dat");
    outfile << setprecision(16);
    for (i = 0; i<3 ; i++) {
        outfile << xdfKep[i]*Lsc<<endl;};
    for (i = 3; i<6 ; i++) {
        outfile << xdfKep[i]*Vsc<<endl;};
    outfile << t2c <<endl;
    outfile.close();
    
    outfile.open("runtime/xs0Kep.dat");
	outfile << setprecision(16);
	for (i = 0; i<3 ; i++) {
	outfile << xsfKep[i]*Lsc<<endl;};
    for (i = 3; i<6 ; i++) {
        outfile << xsfKep[i]*Vsc<<endl;};
    outfile << t2c <<endl;
    outfile.close();

    outfile.open("runtime/xd0J2.dat");
    outfile << setprecision(16);
    for (i = 0; i<3 ; i++) {
        outfile << xdfJ2[i]*Lsc<<endl;};
    for (i = 3; i<6 ; i++) {
        outfile << xdfJ2[i]*Vsc<<endl;};
    outfile << t2c <<endl;
    outfile.close();
    
    outfile.open("runtime/xs0J2.dat");
    outfile << setprecision(16);
    for (i = 0; i<3 ; i++) {
    outfile << xsfJ2[i]*Lsc<<endl;};
    for (i = 3; i<6 ; i++) {
        outfile << xsfJ2[i]*Vsc<<endl;};
    outfile << t2c <<endl;
    outfile.close();
}
