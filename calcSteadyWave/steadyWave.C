/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "steadyWave.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

steadyWave::steadyWave(const dictionary& waveDict)
:
N_(waveDict.lookupOrDefault<int>("nModes",30)),
B_(N_,0.0),
E_(N_+1,0.0),
H_(readScalar(waveDict.lookup("height"))),
d_(readScalar(waveDict.lookup("depth"))),
uType_(waveDict.lookupOrDefault<int>("uType",2)),
nHsteps_(waveDict.lookupOrDefault<int>("nHeightSteps",20)),
pi_(Foam::acos(-1.0)),
g_(waveDict.lookupOrDefault<scalar>("g",9.81))
{

	if ( uType_ == 0 )
	{
		scalar lambda = readScalar(waveDict.lookup("length"));
		k_ = 2*pi_/lambda;
		omega_ = 0.0; //co-moving frame
		//u1or2 is not used for anything in this case
	}
	else if ( (uType_ == 1) || (uType_ == 2) )
	{
		//k_ is calculated by fenton()
		scalar tau = readScalar(waveDict.lookup("period"));
		omega_ = 2*pi_/tau;
		u1or2_ = waveDict.lookupOrDefault("u1or2",0.0);
	}
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool steadyWave::sameWave(const dictionary& waveDict)
{
	bool sameWave(true);
	
	scalar H(readScalar(waveDict.lookup("H")));
	scalar d(readScalar(waveDict.lookup("d")));
	int N(readScalar(waveDict.lookup("N")));
	int nHsteps(readScalar(waveDict.lookup("nHsteps")));
	int uType(readInt(waveDict.lookup("uType")));

	scalar prec(Foam::pow(10.0,1-scalar(IOstream::precision_)));
		
	sameWave = sameWave && (N == N_);
	sameWave = sameWave && (nHsteps == nHsteps_);
	sameWave = sameWave && (uType == uType_);
	sameWave = sameWave && (mag(H/d - H_/d_) < prec);

	scalar k(readScalar(waveDict.lookup("k")));
	scalar omega(readScalar(waveDict.lookup("omega")));
	scalar u1or2(readScalar(waveDict.lookup("u1or2")));		
		
	if ( uType_ == 0 )
	{
		sameWave = sameWave && (mag(k/d - k_/d_) < prec);
	}
	else if ( uType_ == 1 || uType_ == 2 )
	{
		sameWave = sameWave && (mag(omega - omega_) < prec);
		sameWave = sameWave && (mag(u1or2 - u1or2_) < prec);
	}
	else
	{
		Info << "uType has to be 0, 1 or 2!" << endl;		
	}
	
	return sameWave;
}

void steadyWave::fenton()
{

	scalar lort, g = mag(g_);

	if (uType_ == 0)
	{
		//k_ set to 2*pi_/lambda in init();
		scalar lambda = 2*pi_/k_;
		lort = lambda/d_; //Nondimensionalised wave length
		omega_ = 0.0;
		//u1or2_ not used in this case
	} else
	{
		scalar tau = 2*pi_/omega_;
		lort = tau * Foam::sqrt(g/d_); //Nondimensionalised wave period T*sqrt(g/d)
		k_ = Foam::pow(omega_,2.0)/ ( g * Foam::pow( Foam::tanh( Foam::pow( omega_*Foam::sqrt(d_/g) , 3.0/2.0 ) ) , 2.0/3.0) );
		u1or2_ /= Foam::sqrt(g*d_);
	}
	
//	Info << "Specifying initial guess for variables from 1st order theory (Fenton 1999, p. 14)" << endl; 
	scalar kd = k_*d_; //Mean depth normalized by wavenumber, k*d.
	
	scalarField keta(N_ + 1); //Nondimensionalised surface elevation above sea level, keta = k*eta
    forAll(keta,ii)
    {
		keta[ii] = kd + 0.5 * k_ * (H_/nHsteps_) * Foam::cos(ii*pi_/N_);
    }
	scalar u = Foam::sqrt(Foam::tanh(kd)); //Nondimensionalized mean horizontal velocity in co-moving frame, Ubar*sqrt(k/g).
	B_ = 0.0; //Fourier coefficients
	B_[0] = 0.5 * k_ * (H_/nHsteps_) / u;
	scalar q = kd * u; //Nondimensionalized volume flux in co-moving frame, q = Q*sqrt(k^3/g).
	scalar r = kd + 0.5 * Foam::pow(u,2.0); //Nondimensionalized constant in Bernoulli equation, R*k/g.
//	Info << "g = " << g << ", k_ guess = " << k_ << ", kd guess = " << kd << ", u = " << u << ", B[0] = " << B[0] << ", q = " << q << ", r = " << r << endl;

//	Info << "Putting guessed values into x vector" << endl;
	scalarField x(2*N_+5);
    forAll(keta,ii)
    {
		x[ii] = keta[ii];
    }
    forAll(B_,ii)
    {
		x[ii + N_ + 1] = B_[ii];
    }
	x[2*N_+1] = u; 
	x[2*N_+2] = kd;
	x[2*N_+3] = q;
	x[2*N_+4] = r;
	
	for (int ii=1; ii<=nHsteps_; ii++) 
	{
		findRoot(x, ii*H_/nHsteps_/d_, lort);
		Info << "Iteration " << ii << " of " << nHsteps_ << " with H/d = " << ii*H_/nHsteps_/d_ << " gave kd = " << x[2*N_+2] << endl;
	}

	//Converting x vector into named variables
    forAll(keta,ii)
    {
		keta[ii] = x[ii];
    }
    forAll(B_,ii)
    {
		B_[ii] = x[ii + N_ + 1];
    }
	u = x[2*N_+1];
	kd = x[2*N_+2];
	q = x[2*N_+3];
	r = x[2*N_+4];

	k_ = kd/d_;
	R_ = r*g/k_;

	//Calculating E_ using eq. 3.26 (with j running from 0 to N) and eq. 3.27 (with the RHS divided by N to correct the equation) in Fenton 1999.
	//With these corrections A_[j] = 2/N*E_{j+1} for j = [0,N-2], and A_[N-1] = E_{N}/N.
	scalarField etapp = keta/k_;
	etapp[0] *= 0.5;
	etapp[N_] *= 0.5;
	forAll(E_,jj)
	{
		forAll(etapp,mm)
		{
			E_[jj] += etapp[mm] * Foam::cos(jj*mm*pi_/N_);
		}
	}
	
	E_[0] *= 0.5;
	E_[N_] *= 0.5;
	E_ /= N_;

	uBar_ = u / Foam::sqrt(k_/g);
		
	Info << "Done calculating stream function wave" << endl;
}

void steadyWave::findRoot(scalarField& X, scalar h, scalar lort)
{
	scalar tol(1e-5); //HARDCODED TOLERANCE
	scalar err(2*tol);
	int nIter(0);
	while ( err > tol )
	{
		scalarField F0(fun(X,h,lort));
		scalarSquareMatrix M(-jac(X,h,lort));
		LUsolve(M, F0); //LUSolve puts dX into F0
		X += F0;
		err = max(mag(F0));
		nIter++;
		Info << "Root iteration: " <<  nIter << " Error = " << err << endl;
	}
}

scalarField steadyWave::fun(const scalarField& x, scalar h, scalar lort)
{

//	Info << "Converting from x vector into named variables" << endl;
	scalarField keta(N_+1);
    forAll(keta,ii)
    {
		keta[ii] = x[ii];
    }
	scalarField B(N_);
    forAll(B,ii)
    {
		int jj = ii + N_ + 1;
		B[ii] = x[jj];
    }
	scalar u = x[2*N_+1]; 
	scalar kd = x[2*N_+2];
	scalar q = x[2*N_+3];
	scalar r = x[2*N_+4]; 

	scalarField kinBC(N_+1); //Kinematic boundary condition, (3.7) in Fenton 1999
	scalarField dynBC(N_+1); //Dynamic boundary condition (3.8) in Fenton 1999
	scalar U, V, kx, ke;
	int j;
	scalar meanDepthEq = - (0.5/N_)*(keta[0] + keta[N_]) - kd; //Mean depth equation (3.9) in Fenton 1999
	forAll(kinBC,ii)
    {
		kx = ii*pi_/N_;
		ke = keta[ii];
		kinBC[ii] = - u * ke + q;
		dynBC[ii] = ke - r;
		U = - u;
		V = 0.0;
		forAll(B,jj)
		{
			j = jj + 1;
			kinBC[ii] += B[jj] * (Foam::sinh(j*ke) / Foam::cosh(j*kd)) * Foam::cos(j*kx);
			U += j * B[jj] * (Foam::cosh(j*ke) / Foam::cosh(j*kd)) * Foam::cos(j*kx);
			V += j * B[jj] * (Foam::sinh(j*ke) / Foam::cosh(j*kd)) * Foam::sin(j*kx);
		}
		dynBC[ii] += 0.5*Foam::pow(U,2.0) + 0.5*Foam::pow(V,2.0);
		meanDepthEq += keta[ii]/N_;
    }
	scalar waveHeightEq = keta[0] - keta[N_] - kd*h; //Wave height equation (3.10) in Fenton 1999
	
	scalar lastEq = 0;
	if (uType_ == 0) {
		scalar l = lort;
		lastEq = kd - 2.0*pi_/l; // Relative wavelength equation (3.11) in Fenton 1999
	}
	else {
		scalar t = lort;
		if (uType_ == 1) {
			lastEq = Foam::sqrt(kd)*u + kd*u1or2_ - 2.0*pi_/t; //Equation for mean current in laboratory frame (3.15) in Fenton 1999
		}
		else if (uType_ == 2) {
			lastEq = q/Foam::sqrt(kd) + kd*u1or2_ - 2.0*pi_/t; //Equation for mass transport velocity (3.16) in Fenton 1999
		}
	}	

//	Info << "Writing wave conditions into solution vector f" << endl;
	scalarField F(2*N_+5);
	forAll(kinBC,ii)
    {
		F[ii] = kinBC[ii];
		F[ii+N_+1] = dynBC[ii];
    }
	F[2*N_+2] = meanDepthEq;
	F[2*N_+3] = waveHeightEq;
	F[2*N_+4] = lastEq;

	return F;
}

scalarSquareMatrix steadyWave::jac(const scalarField& x, scalar h, scalar lort)
{
	scalar d(1e-10); //HARDCODED STEP SIZE
	scalarSquareMatrix J(x.size());
	forAll(x,ii)
	{
		scalarField dx(x.size(),0.0);
		dx[ii] = d;
		scalarField dF((fun(x+dx,h,lort) - fun(x,h,lort))/d);
		forAll(dF,jj)
		{
			J[jj][ii] = dF[jj];
		}
	}
	return J;
}

void steadyWave::writeDict(word outDictName)
{
	dictionary outDict(outDictName);

	outDict.add("H",H_,true);
	outDict.add("d",d_,true);
	outDict.add("N",N_,true);
	outDict.add("nHsteps",nHsteps_,true);
	outDict.add("uType",uType_,true);
	outDict.add("k",k_,true);
	outDict.add("omega",omega_,true);
	outDict.add("u1or2",u1or2_,true);
	outDict.add("R",R_,true);
	outDict.add("uBar",uBar_,true);
	outDict.add("B",B_,true);
	outDict.add("E",E_,true);
	outDict.add("lambda",2*Foam::acos(-1.0)/k_,true);
	if (omega_ > 0)
	{
		outDict.add("T",2*Foam::acos(-1.0)/omega_,true);
	}
	outDict.add("c",omega_/k_,true);
	outDict.write(OFstream(outDictName)(), false);

}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
