/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

Application
    cdugksFoam

Description
    cdugksFoam-Rykov diatomic molecule model based on Lianhua Zhu's and Qi Zhang's works

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvDVM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    //Add global valid option
    Foam::argList::addBoolOption("dvParallel", "Use discrete velocity domain decomposition\n");
    Foam::argList::addOption(
        "pd",
        "label",
        "num of phy domain"
    );
    //convergence monitor
    scalar TemperatureChange = 1.0;
    scalar rhoChange = 1.0;
    scalar Uchange = 1.0;

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "readTimeControlsExplicit.H"

    fvDVM dvm(rho, U, T, &argc, &argv, args);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    if(dvm.mpiReducer().rank() == 0)
        Info<< "\nStarting time loop\n" << endl;
    
    label It = 0;
    //while (runTime.run() && (TemperatureChange > convergeTol))
    while (runTime.run() && (Uchange > convergeTol))
    {
        #include "CourantNo.H" // calculate the Co num
        #include "readTimeControlsExplicit.H"
        #include "setDeltaTvar.H"

        runTime++;
        It++;

        if(dvm.mpiReducer().rank() == 0)
            Info<< "Time = " << runTime.timeName() << nl << endl;

        dvm.evolution();


        //write DF
        //dvm.writeDFonCells();
        runTime.write();
        if(dvm.mpiReducer().rank() == 0)
        {
            Info<< "Step =" << It << "  ExecutionTime = " << runTime.elapsedCpuTime() << " s"
                << "  ClockTime = " << runTime.elapsedClockTime() << " s"
                << nl << endl;
        }
        if(It%convergeCheckSteps == 0 && It >= convergeCheckSteps)
        {
                tmp<Foam::GeometricField<scalar, Foam::fvPatchField, Foam::volMesh> > 
                    deltaTem = mag(T-Told);
                tmp<Foam::GeometricField<scalar, Foam::fvPatchField, Foam::volMesh> > 
                    deltaRho = mag(rho-rhoOld);
                tmp<Foam::GeometricField<scalar, Foam::fvPatchField, Foam::volMesh> > 
                    deltaU = mag(U-Uold);
		// 保存上一步残差
		scalar TemperatureChangePrev = TemperatureChange;
		scalar rhoChangePrev         = rhoChange;
		scalar UchangePrev           = Uchange;

		// 计算当前步残差
		TemperatureChange = gSum(deltaTem()) / gSum(T);
		rhoChange         = gSum(deltaRho()) / gSum(rho);
		Uchange           = gSum(deltaU())   / gSum(mag(U)());

		// 取三者最大值作为综合残差
		//scalar residualPrev    = max(max(rhoChangePrev, UchangePrev), TemperatureChangePrev);
		//scalar residualCurrent = max(max(rhoChange,     Uchange),     TemperatureChange);

		// 取三者平方和开方作为综合残差（L2范数）
		scalar residualPrev    = Foam::sqrt(
                             rhoChangePrev         * rhoChangePrev
                           + UchangePrev           * UchangePrev
                           + TemperatureChangePrev * TemperatureChangePrev
                         );

		scalar residualCurrent = Foam::sqrt(
                             rhoChange         * rhoChange
                           + Uchange           * Uchange
                           + TemperatureChange * TemperatureChange
                         );

		dvm.adjustCFL(residualCurrent, residualPrev);

		Info << "Temperature changes = " << TemperatureChange << endl;
		Info << "Density     changes = " << rhoChange         << endl;
		Info << "Velocity    changes = " << Uchange << nl     << endl;
		//Info << "Current CFL = " << dvm.CFL() << nl  	      << endl;
                Told = T;
                rhoOld = rho;
                Uold = U;
        }
    }

    if(dvm.mpiReducer().rank() == 0)
        Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
