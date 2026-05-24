/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include <mpi.h>
#include "fvDVM.H"
#include "constants.H"
#include "fvm.H"
#include "calculatedMaxwellFvPatchField.H"
#include "symmetryModFvPatchField.H"
#include "pressureInFvPatchField.H"
#include "pressureOutFvPatchField.H"
#include "scalarIOList.H"
#include "fieldMPIreducer.H"

using namespace Foam::constant;
using namespace Foam::constant::mathematical;

#if FOAM_MAJOR <= 3
    #define BOUNDARY_FIELD_REF boundaryField()
#else
    #define BOUNDARY_FIELD_REF boundaryFieldRef()
#endif


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fvDVM, 0);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fvDVM::setDVgrid
(
    scalarField& weights,
    scalarField& Xis,
    scalar xiMin,
    scalar xiMax,
    label nXi
)
{
    // Read from file ./constant/Xis and ./constant/weights
    scalarIOList xiList
    (
        IOobject
        (
             "Xis",
             time_.caseConstant(),
             mesh_,
             IOobject::MUST_READ,
             IOobject::NO_WRITE
        )
    );

    scalarIOList weightList
    (
        IOobject
        (
             "weights",
             time_.caseConstant(),
             mesh_,
             IOobject::MUST_READ,
             IOobject::NO_WRITE
        )
    );

    for (label i = 0; i < nXi ; i++)
    {
        weights[i] = weightList[i];
        Xis[i] = xiList[i];
    }
}

void Foam::fvDVM::initialiseDV()
{

    scalarField weightsGlobal;
    vectorField XisGlobal;
    labelField  symmXtgID;
    labelField  symmYtgID;
    labelField  symmZtgID;


    if (vMesh_ == "no")
    {
        scalarField weights1D(nXiPerDim_);
        scalarField Xis(nXiPerDim_);

        //get discrete velocity points and weights
        setDVgrid
        (
            weights1D,
            Xis, 
            xiMin_.value(), 
            xiMax_.value(), 
            nXiPerDim_
        );

        if (mesh_.nSolutionD() == 3)    //3D(X & Y & Z)
        {
            nXiX_ = nXiY_ = nXiZ_ = nXiPerDim_;
            nXi_ = nXiX_*nXiY_*nXiZ_;

            weightsGlobal.setSize(nXi_);
            XisGlobal.setSize(nXi_);
            symmXtgID.setSize(nXi_);
            symmYtgID.setSize(nXi_);
            symmZtgID.setSize(nXi_);

            label i = 0;
            for (label iz = 0; iz < nXiZ_; iz++)
            {
                for (label iy = 0; iy < nXiY_; iy++)
                {
                    for (label ix = 0; ix < nXiZ_; ix++)
                    {
                        scalar weight = weights1D[iz]*weights1D[iy]*weights1D[ix];
                        vector xi(Xis[ix], Xis[iy], Xis[iz]);
                        weightsGlobal[i] = weight;
                        XisGlobal[i] = xi;
                        symmXtgID[i] = iz*nXiY_*nXiX_ + iy*nXiX_ + (nXiX_ - ix -1);
                        symmYtgID[i] = iz*nXiY_*nXiX_ + (nXiY_ - iy - 1)*nXiX_ + ix;
                        symmZtgID[i] = (nXiZ_ - iz -1)*nXiY_*nXiX_ + iy*nXiX_ + ix;
                        i++;
                    }
                }
            }
        }
        else
        {
            if (mesh_.nSolutionD() == 2)    //2D (X & Y)
            {
                nXiX_ = nXiY_ = nXiPerDim_;
                nXiZ_ = 1;
                nXi_ = nXiX_*nXiY_*nXiZ_;
                weightsGlobal.setSize(nXi_);
                XisGlobal.setSize(nXi_);
                symmXtgID.setSize(nXi_);
                symmYtgID.setSize(nXi_);
                symmZtgID.setSize(nXi_);
                label i = 0;
                for (label iy = 0; iy < nXiY_; iy++)
                {
                    for (label ix = 0; ix < nXiX_; ix++)
                    {
                        scalar weight = weights1D[iy]*weights1D[ix]*1;
                        vector xi(Xis[ix], Xis[iy], 0.0);
                        weightsGlobal[i] = weight;
                        XisGlobal[i] = xi;
                        symmXtgID[i] = iy*nXiX_ + (nXiX_ - ix -1);
                        symmYtgID[i] = (nXiY_ - iy - 1)*nXiX_ + ix;
                        symmZtgID[i] = 0;
                        i++;
                    }
                }
            }
            else    //1D (X)
            {
                nXiX_ = nXiPerDim_;
                nXiY_ = nXiZ_ = 1;
                nXi_ = nXiX_*nXiY_*nXiZ_;
                weightsGlobal.setSize(nXi_);
                XisGlobal.setSize(nXi_);
                symmXtgID.setSize(nXi_);
                symmYtgID.setSize(nXi_);
                symmZtgID.setSize(nXi_);
                label i = 0;
                for (label ix = 0; ix < nXiX_; ix++)
                {
                    scalar weight = weights1D[ix]*1*1;
                    vector xi(Xis[ix], 0.0, 0.0);
                    weightsGlobal[i] = weight;
                    XisGlobal[i] = xi;
                    symmXtgID[i] = (nXiX_ - ix -1);
                    symmYtgID[i] = 0;              
                    symmZtgID[i] = 0;              
                    i++;
                }
            }
        }
                
        if (mpiReducer_.rank() == 0)
        {
            Info << "fvDVM : Allocated " << XisGlobal.size()
                << " discrete velocities" << endl;
        }
        label nA = nXi_ / mpiReducer_.csize();//average
        label nB = nXi_ - nA * mpiReducer_.csize();//mod
        label nXiPart = nA + (label)(mpiReducer_.crank() < nB);//cyclic partition
        DV_.setSize(nXiPart);
        if (mpiReducer_.rank() %mpiReducer_.npd()==0) //first column
        {
            std::cout << "nrank    " << mpiReducer_.crank() <<std:: endl;
            std::cout << "nXisPart " << nXiPart << std::endl;
            //std::cout<<mpiReducer_.crank()<<" "<<mpiReducer_.rank()<<endl;
        }

        label chunk = 0;
        label gid = 0;
            
        forAll(DV_, i)
        {
            gid = chunk + mpiReducer_.crank();
            DV_.set
            (
                i,
                new discreteVelocity
                (
                    *this,
                    mesh_,
                    time_,
                    weightsGlobal[gid],
                    dimensionedVector("xi", dimLength / dimTime, XisGlobal[gid]),
                    i,
                    symmXtgID[gid],
                    symmYtgID[gid],
                    symmZtgID[gid]
                )
            );
            chunk += mpiReducer_.csize();          
        }
    }
    else
    {
        scalarField weights1D(nXiPerDim_);
	    scalarField Xis(nXiPerDim_);
        xiMax_.value()=0.0;
        Info << "Reading vMesh" << endl;
        //we have polyMesh and vMesh at the beginning
	    string mesh_name = args_.path() / word("constant") / word("polyMesh");//polyMesh dir
	    string vmesh_name = args_.path() / word("constant") / word("vMesh");//vMesh dir
	    string temp_name = args_.path() / word("constant") / word("tMesh");//temp name

        //1.change folder name, polyMesh-->tMesh,vMesh-->polyMesh
	    //as OF only recognizes certain path in the couse of reading mesh files
	
	    if(mpiReducer_.rank() == 0 || args_.optionFound("parallel")) 
        {
          std::rename(mesh_name.c_str(), temp_name.c_str());
          std::rename(vmesh_name.c_str(), mesh_name.c_str());
        }

        if (mpiReducer_.dvParallel())
            MPI_Barrier(MPI_COMM_WORLD);

        //2.read and construct vmesh 
        Foam::fvMesh vmesh
        (
            Foam::IOobject
            (
                Foam::fvMesh::defaultRegion,
                time_.timeName(),
                time_,
                Foam::IOobject::MUST_READ
            )
        );

        if (mpiReducer_.dvParallel())
            MPI_Barrier(MPI_COMM_WORLD);
        if(mpiReducer_.rank() == 0 || args_.optionFound("parallel")) 
        {
            std::rename(mesh_name.c_str(), vmesh_name.c_str());
            std::rename(temp_name.c_str(), mesh_name.c_str());
        }


        if (mesh_.nSolutionD() == 3)    //3D(X & Y & Z)
        {
            nXiX_ = nXiY_ = nXiZ_ = vmesh.C().size();
            nXi_ = vmesh.C().size();

            weightsGlobal.setSize(nXi_);
            XisGlobal.setSize(nXi_);
            symmXtgID.setSize(nXi_);
            symmYtgID.setSize(nXi_);
            symmZtgID.setSize(nXi_);

            label i;
            for (i = 0; i < nXi_; i++) 
            {
                vector xi(vmesh.C()[i].x(), vmesh.C()[i].y(), vmesh.C()[i].z());
                scalar weight(vmesh.V()[i]);
                weightsGlobal[i] = weight;
                XisGlobal[i] = xi;
                xiMax_.value()=max(xiMax_.value(),mag(xi));
            }
        }
        else
        {
            if (mesh_.nSolutionD() == 2)    //2D (X & Y)
            {
                nXiX_ = nXiY_ = vmesh.C().size();
                nXiZ_ = 1;
                nXi_ = vmesh.C().size();
                weightsGlobal.setSize(nXi_);
                XisGlobal.setSize(nXi_);
                symmXtgID.setSize(nXi_);
                symmYtgID.setSize(nXi_);
                symmZtgID.setSize(nXi_);
                label i;
                for (i = 0; i < nXi_; i++) 
                {
                    vector xi(vmesh.C()[i].x(), vmesh.C()[i].y(), 0.0);
                    scalar weight(vmesh.V()[i]);
                    weightsGlobal[i] = weight;
                    XisGlobal[i] = xi;
                    xiMax_.value()=max(xiMax_.value(),mag(xi));
                }
            }
            else    //1D (X)
            {
                nXiX_ = vmesh.C().size();
                nXiY_ = nXiZ_ = 1;
                nXi_ = vmesh.C().size();
                weightsGlobal.setSize(nXi_);
                XisGlobal.setSize(nXi_);
                symmXtgID.setSize(nXi_);
                symmYtgID.setSize(nXi_);
                symmZtgID.setSize(nXi_);
                label i;
                for (i = 0; i < nXi_; i++) 
                {
                    vector xi(vmesh.C()[i].x(), 0.0, 0.0);
                    scalar weight(vmesh.V()[i]);
                    weightsGlobal[i] = weight;
                    XisGlobal[i] = xi;
                    xiMax_.value()=max(xiMax_.value(),mag(xi));
                }
            }
        }


        if (mpiReducer_.rank() == 0)
        {
            Info << "fvDVM : Allocated " << XisGlobal.size()
                << " discrete velocities" << endl;
        }
        label nA = nXi_ / mpiReducer_.csize();//average
        label nB = nXi_ - nA * mpiReducer_.csize();//mod
        label nXiPart = nA + (label)(mpiReducer_.crank() < nB);//cyclic partition
        DV_.setSize(nXiPart);
        if (mpiReducer_.rank() %mpiReducer_.npd()==0) //first column
        {
            std::cout << "nrank    " << mpiReducer_.crank() <<std:: endl;
            std::cout << "nXisPart " << nXiPart << std::endl;
            //std::cout<<mpiReducer_.crank()<<" "<<mpiReducer_.rank()<<endl;
        }

        label chunk = 0;
        label gid = 0;
            
        forAll(DV_, i)
        {
            gid = chunk + mpiReducer_.crank();
            DV_.set
            (
                i,
                new discreteVelocity
                (
                    *this,
                    mesh_,
                    time_,
                    weightsGlobal[gid],
                    dimensionedVector("xi", dimLength / dimTime, XisGlobal[gid]),
                    i,
                    symmXtgID[gid],
                    symmYtgID[gid],
                    symmZtgID[gid]
                )
            );
            chunk += mpiReducer_.csize();          
        }
        
        //clear vmesh data structure to save memory
        vmesh.clearOut();
    }
    


}

void Foam::fvDVM::setCalculatedMaxwellRhoBC()
{
#if FOAM_MAJOR <= 3
    GeometricField<scalar, fvPatchField, volMesh>::GeometricBoundaryField& 
        rhoBCs = rhoVol_.boundaryField();
#else
    GeometricField<scalar, fvPatchField, volMesh>::Boundary& 
        rhoBCs = rhoVol_.BOUNDARY_FIELD_REF;
#endif
    forAll(rhoBCs, patchi)
    {
        if (rhoBCs[patchi].type() == "calculatedMaxwell")
        {

            const vectorField& SfPatch = mesh_.Sf().boundaryField()[patchi];
            calculatedMaxwellFvPatchField<scalar>& rhoPatch = 
                refCast<calculatedMaxwellFvPatchField<scalar> >(rhoBCs[patchi]);
            const vectorField& Upatch = Uvol_.boundaryField()[patchi];
            const scalarField& Tpatch = Tvol_.boundaryField()[patchi];

            forAll(rhoPatch, facei)
            {
                vector faceSf = SfPatch[facei];
                rhoPatch.inComingByRho()[facei] = 0; // set to zero
                forAll(DV_, dvi) // add one by one
                {
                    vector xi = DV_[dvi].xi().value();
                    scalar weight = DV_[dvi].weight();
                    if ( (xi & faceSf) < 0) //inComing
                    {
                        rhoPatch.inComingByRho()[facei] += 
                          - weight*(xi & faceSf)
                          *DV_[dvi].equilibriumMaxwellByRho
                          (
                              Upatch[facei], 
                              Tpatch[facei]
                          );
                    }
                }
            }

            if (mpiReducer_.dvParallel()&& mpiReducer_.npd() < mpiReducer_.nproc()) 
            {
				mpiReducer_.reduceField(rhoPatch.inComingByRho());
            }
        }

    }
}

void Foam::fvDVM::setSymmetryModRhoBC()
{
    //prepare the container (set size) to store all DF on the patchi
#if FOAM_MAJOR <= 3
    GeometricField<scalar, fvPatchField, volMesh>::GeometricBoundaryField& 
        rhoBCs = rhoVol_.boundaryField();
#else
    GeometricField<scalar, fvPatchField, volMesh>::Boundary& 
        rhoBCs = rhoVol_.BOUNDARY_FIELD_REF;
#endif
    forAll(rhoBCs, patchi)
    {
        label ps = rhoBCs[patchi].size();
        if (rhoBCs[patchi].type() == "symmetryMod")
        {
            symmetryModFvPatchField<scalar>& rhoPatch = 
                refCast<symmetryModFvPatchField<scalar> >(rhoBCs[patchi]);
            rhoPatch.dfContainer().setSize(ps*nXi_*2); //*2 means g and h
        }
    }
}

void Foam::fvDVM::updateGHbarPvol()
{
    forAll(DV_, DVid)
        DV_[DVid].updateGHvol();
}


void Foam::fvDVM::updateGHbarSurf()
{
    forAll(DV_, DVid)
        DV_[DVid].updateGHInterpolation();
}


void Foam::fvDVM::updateMaxwellWallRho()
{
#if FOAM_MAJOR <= 3
    GeometricField<scalar, fvPatchField, volMesh>::GeometricBoundaryField& 
        rhoBCs = rhoVol_.boundaryField();
#else
    GeometricField<scalar, fvPatchField, volMesh>::Boundary& 
        rhoBCs = rhoVol_.BOUNDARY_FIELD_REF;
#endif
    forAll(rhoBCs, patchi)
    {
        if (rhoBCs[patchi].type() == "calculatedMaxwell")
        {
            calculatedMaxwellFvPatchField<scalar>& rhoPatch = 
                refCast<calculatedMaxwellFvPatchField<scalar> >(rhoBCs[patchi]);
            if (mpiReducer_.dvParallel()&& mpiReducer_.npd() < mpiReducer_.nproc()) 
            {
				mpiReducer_.reduceField(rhoPatch.outGoing());
            }
        }
    }
    rhoVol_.correctBoundaryConditions();
}

void Foam::fvDVM::updateGHbarSurfMaxwellWallIn()
{
    forAll(DV_, DVid)
        DV_[DVid].updateGHbarSurfMaxwellWallIn();
}

void Foam::fvDVM::updateGHbarSurfSymmetryIn()
{
    //1. copy all DV's g/h to rho patch's dfContainer
    //2. MPI_Allgather the rho patch's dfContainer
    //if(args_.optionFound("dvParallel"))
    //{
    label rank  = mpiReducer_.rank();
    label nproc = mpiReducer_.nproc();
#if FOAM_MAJOR <= 3
    GeometricField<scalar, fvPatchField, volMesh>::GeometricBoundaryField& 
        rhoBCs = rhoVol_.boundaryField();
#else
    GeometricField<scalar, fvPatchField, volMesh>::Boundary& 
        rhoBCs = rhoVol_.BOUNDARY_FIELD_REF;
#endif
    forAll(rhoBCs, patchi)
    {
        label ps = rhoBCs[patchi].size();
        if (rhoBCs[patchi].type() == "symmetryMod")
        {
            symmetryModFvPatchField<scalar>& rhoPatch = 
                refCast<symmetryModFvPatchField<scalar> >(rhoBCs[patchi]);
            //compose the recvcout and displacement array
            labelField recvc(nproc);
            labelField displ(nproc);
            label chunck = nXi_/nproc;
            label left   = nXi_%nproc;
            forAll(recvc, i)
            {
                recvc[i] = 2*ps*(chunck + (i<left)) ;
                if(i<=left)
                    displ[i] = i*2*ps*(chunck + 1); // (i<=nXi_%nproc)
                else
                    displ[i] = 2*ps*(left*(chunck +1) + (i-left)*(chunck));
            }

            // check 12*28+15 dv's g
            label did = 1709;
            label pp  = did%nproc;
            label lid = did/nproc;
            //if(rank==pp)
            //{
                //Info << "processing by rank " << rank << endl;
                //Info << "12*28+15 outging g " << DV_[lid].gSurf()[0] << endl;
                //Info << "12*28+15 outging xi " << DV_[lid].xi() <<endl;
                //Info << "12*28+15 outging at boundary " << DV_[lid].gSurf().boundaryField()[patchi][0] << endl;
            //}
            // memcpy each dv's g/h to rho
            forAll(DV_, DVid)
            {
                //label shift = (nXi_ / nproc * rank + DVid)*2*ps;
                label shift = displ[rank] + DVid*2*ps;
                memcpy( (rhoPatch.dfContainer().data() + shift),
                        DV_[DVid].gSurf().boundaryField()[patchi].cdata(), ps*sizeof(scalar));
                memcpy( (rhoPatch.dfContainer().data() + shift + ps),
                        DV_[DVid].hSurf().boundaryField()[patchi].cdata(), ps*sizeof(scalar));
            }

            // check 
            //if(rank == pp)
                //Info << "dv gid 1709's g = " <<rhoPatch.dfContainer()[displ[pp]+lid*2*ps+32]<< endl;;


            //Allgather
            MPI_Allgatherv(
                //rhoPatch.dfContainer().data() + displ[rank],//2*ps*nXI_/nproc*rank, //send*
                MPI_IN_PLACE,
                2*ps*DV_.size(), //(how many DV i processed) * 2 * patch size
                MPI_DOUBLE,
                rhoPatch.dfContainer().data(),
                recvc.data(),
                displ.data(),
                MPI_DOUBLE,
                MPI_COMM_WORLD
                );
        }
    }
        forAll(DV_, DVid)
            DV_[DVid].updateGHbarSurfSymmetryIn();
}

void Foam::fvDVM::updateMacroSurf()
{
    // Init to zero before add one DV by one DV
    rhoSurf_ =  dimensionedScalar("0", rhoSurf_.dimensions(), 0);
    Usurf_ = dimensionedVector("0", Usurf_.dimensions(), vector(0, 0, 0));
    Tsurf_ = dimensionedScalar("0", Tsurf_.dimensions(), 0);
    qSurf_ = dimensionedVector("0", qSurf_.dimensions(), vector(0, 0, 0));
    stressSurf_ = dimensionedTensor
        (
            "0", 
            stressSurf_.dimensions(), 
            pTraits<tensor>::zero
        );

    surfaceVectorField rhoUsurf = rhoSurf_*Usurf_;
    surfaceScalarField rhoEsurf = rhoSurf_*magSqr(Usurf_);

    forAll(DV_, dvi)
    {
        discreteVelocity& dv = DV_[dvi];
        rhoSurf_  += dXiCellSize_*dv.weight()*dv.gSurf();
        rhoUsurf  += dXiCellSize_*dv.weight()*dv.gSurf()*dv.xi();
        rhoEsurf  += 0.5*dXiCellSize_*dv.weight()
           *(
                dv.gSurf()*magSqr(dv.xi()) 
              + dv.hSurf()
            );
    }

    if (mpiReducer_.dvParallel()&& mpiReducer_.npd() < mpiReducer_.nproc()) 
    {
		mpiReducer_.reduceField(rhoSurf_);
        mpiReducer_.reduceField(rhoUsurf);
        mpiReducer_.reduceField(rhoEsurf);
    }

    //- get Prim. from Consv.
    Usurf_ = rhoUsurf/rhoSurf_;

    Tsurf_ = (rhoEsurf - 0.5*rhoSurf_*magSqr(Usurf_))/((KInner_ + 3)/2.0*R_*rhoSurf_);

    updateTau(tauSurf_, Tsurf_, rhoSurf_);
    //- peculiar vel.

    surfaceVectorField c = Usurf_;

    //-get part heat flux 
    forAll(DV_, dvi)
    {
        discreteVelocity& dv = DV_[dvi];
        c = dv.xi() - Usurf_;
        qSurf_ += 0.5*dXiCellSize_*dv.weight()*c
            *(
                 magSqr(c)*dv.gSurf() 
               + dv.hSurf()
             );
        //- stressSurf is useless as we never update cell macro by macro flux 
        //- Comment out it as it is expansive
        //stressSurf_ += 
            //dXiCellSize_*dv.weight()*dv.gSurf()*c*c;
    }
    //- Get global heat flux, via MPI_Allreuce
    
    if (mpiReducer_.dvParallel()&& mpiReducer_.npd() < mpiReducer_.nproc())
    {
        mpiReducer_.reduceField(qSurf_);
    }

    //// correction for heat flux for implict solver from Rui zhang
    qSurf_ = tauSurf_/(tauSurf_ + 0.5*time_.deltaT()*Pr_)*qSurf_;
/*

#if FOAM_MAJOR <= 3
    GeometricField<scalar, fvPatchField, volMesh>::GeometricBoundaryField& 
        rhoBCs = rhoVol_.boundaryField();
#else
    GeometricField<scalar, fvPatchField, volMesh>::Boundary& 
        rhoBCs = rhoVol_.BOUNDARY_FIELD_REF;
#endif
    qWall_ = dimensionedVector("0", qWall_.dimensions(), vector(0, 0, 0));
    stressWall_ = dimensionedTensor
        (
            "0", 
            stressWall_.dimensions(), 
            pTraits<tensor>::zero
        );
    qfluxWall_ = dimensionedScalar("0", qfluxWall_.dimensions(), 0);
    presurefluxWall_ = dimensionedScalar("0", stressWall_.dimensions(), 0);
    stressfluxWall_ = dimensionedVector("0", stressWall_.dimensions(), vector(0, 0, 0));
    forAll(rhoBCs, patchi)
    {
        const vectorField n
                (
                        mesh_.Sf().boundaryField()[patchi]
                        /mesh_.magSf().boundaryField()[patchi]
                );        
        if (rhoBCs[patchi].type() == "calculatedMaxwell")
        {
            fvPatchField<vector>& qPatch = qWall_.BOUNDARY_FIELD_REF[patchi];
            fvPatchField<vector>& Upatch = Uvol_.BOUNDARY_FIELD_REF[patchi];
            fvPatchField<tensor>& stressPatch = stressWall_.BOUNDARY_FIELD_REF[patchi];
            fvPatchField<scalar>& qfluxPatch = qfluxWall_.BOUNDARY_FIELD_REF[patchi];
            fvPatchField<scalar>& presurefluxPatch = presurefluxWall_.BOUNDARY_FIELD_REF[patchi];
            fvPatchField<vector>& stressfluxPatch = stressfluxWall_.BOUNDARY_FIELD_REF[patchi];
            //- tau at surface use the tau at slip temperature as it is.
            fvsPatchField<scalar>&  tauPatch = tauSurf_.BOUNDARY_FIELD_REF[patchi];
            forAll(qPatch, facei)
            {
                forAll(DV_, dvi)
                {
                    scalar dXiCellSize = dXiCellSize_.value();
                    discreteVelocity& dv = DV_[dvi];
                    vector xi = dv.xi().value();
                    vector c = xi - Upatch[facei];
                    qPatch[facei] += 0.5*dXiCellSize*dv.weight()*c  //sometimes wall moves, then c != \xi
                        *(
                             magSqr(c)*dv.gSurf().boundaryField()[patchi][facei]
                           + dv.hSurf().boundaryField()[patchi][facei]
                         );
                    stressPatch[facei] += 
                        dXiCellSize*dv.weight()*dv.gSurf().boundaryField()[patchi][facei]*xi*xi;
                    qfluxPatch[facei] += 0.5*dXiCellSize*dv.weight()*(c&n[facei])
                                     *(
                                             magSqr(c)*dv.gSurf().boundaryField()[patchi][facei]
                                             + dv.hSurf().boundaryField()[patchi][facei]
                                     );
                    presurefluxPatch[facei] += dXiCellSize*dv.weight()*magSqr(c&n[facei])*dv.gSurf().boundaryField()[patchi][facei];
                    stressfluxPatch[facei] += dXiCellSize*(dv.weight()*(c&n[facei])* dv.gSurf().boundaryField()[patchi][facei]*c)^n[facei];
                }

            }
            if (mpiReducer_.dvParallel()&& mpiReducer_.npd() < mpiReducer_.nproc())
            {
                mpiReducer_.reduceField(qPatch);
                mpiReducer_.reduceField(stressPatch);
                mpiReducer_.reduceField(qfluxPatch);
                mpiReducer_.reduceField(presurefluxPatch);
                mpiReducer_.reduceField(stressfluxPatch);
            }
        }
    }*/

}


void Foam::fvDVM::updateGHsurf()
{
    forAll(DV_, DVid)
        DV_[DVid].updateGHsurf();

}

void Foam::fvDVM::updateWall()
{
#if FOAM_MAJOR <= 3
    GeometricField<scalar, fvPatchField, volMesh>::GeometricBoundaryField& 
        rhoBCs = rhoVol_.boundaryField();
#else
    GeometricField<scalar, fvPatchField, volMesh>::Boundary& 
        rhoBCs = rhoVol_.BOUNDARY_FIELD_REF;
#endif
    qWall_ = dimensionedVector("0", qWall_.dimensions(), vector(0, 0, 0));
    stressWall_ = dimensionedTensor
        (
            "0", 
            stressWall_.dimensions(), 
            pTraits<tensor>::zero
        );
    qfluxWall_ = dimensionedScalar("0", qfluxWall_.dimensions(), 0);
    presurefluxWall_ = dimensionedScalar("0", stressWall_.dimensions(), 0);
    stressfluxWall_ = dimensionedVector("0", stressWall_.dimensions(), vector(0, 0, 0));
    forAll(rhoBCs, patchi)
    {
        const vectorField n
                (
                        mesh_.Sf().boundaryField()[patchi]
                        /mesh_.magSf().boundaryField()[patchi]
                );        
        if (rhoBCs[patchi].type() == "calculatedMaxwell")
        {
            fvPatchField<vector>& qPatch = qWall_.BOUNDARY_FIELD_REF[patchi];
            fvPatchField<vector>& Upatch = Uvol_.BOUNDARY_FIELD_REF[patchi];
            fvPatchField<tensor>& stressPatch = stressWall_.BOUNDARY_FIELD_REF[patchi];
            fvPatchField<scalar>& qfluxPatch = qfluxWall_.BOUNDARY_FIELD_REF[patchi];
            fvPatchField<scalar>& presurefluxPatch = presurefluxWall_.BOUNDARY_FIELD_REF[patchi];
            fvPatchField<vector>& stressfluxPatch = stressfluxWall_.BOUNDARY_FIELD_REF[patchi];
            //- tau at surface use the tau at slip temperature as it is.
            fvsPatchField<scalar>&  tauPatch = tauSurf_.BOUNDARY_FIELD_REF[patchi];
            forAll(qPatch, facei)
            {
                forAll(DV_, dvi)
                {
                    scalar dXiCellSize = dXiCellSize_.value();
                    discreteVelocity& dv = DV_[dvi];
                    vector xi = dv.xi().value();
                    vector c = xi - Upatch[facei];
                    qPatch[facei] += 0.5*dXiCellSize*dv.weight()*c  //sometimes wall moves, then c != \xi
                        *(
                             magSqr(c)*dv.gSurf().boundaryField()[patchi][facei]
                           + dv.hSurf().boundaryField()[patchi][facei]
                         );
//                    stressPatch[facei] += 
//                        dXiCellSize*dv.weight()*dv.gSurf().boundaryField()[patchi][facei]*xi*xi;
                    stressPatch[facei] += 
                            dXiCellSize*dv.weight()*dv.gSurf().boundaryField()[patchi][facei]*c*c;
                    qfluxPatch[facei] += 0.5*dXiCellSize*dv.weight()*(c&n[facei])
                                     *(
                                             magSqr(c)*dv.gSurf().boundaryField()[patchi][facei]
                                             + dv.hSurf().boundaryField()[patchi][facei]
                                     );
                    presurefluxPatch[facei] += dXiCellSize*dv.weight()*magSqr(c&n[facei])*dv.gSurf().boundaryField()[patchi][facei];
                    stressfluxPatch[facei] += dXiCellSize*(dv.weight()*(c&n[facei])* dv.gSurf().boundaryField()[patchi][facei]*c)^n[facei];
                }

            }
            if (mpiReducer_.dvParallel()&& mpiReducer_.npd() < mpiReducer_.nproc())
            {
                mpiReducer_.reduceField(qPatch);
                mpiReducer_.reduceField(stressPatch);
                mpiReducer_.reduceField(qfluxPatch);
                mpiReducer_.reduceField(presurefluxPatch);
                mpiReducer_.reduceField(stressfluxPatch);
            }
        }
    }
}

void Foam::fvDVM::updateGHtildeVol()
{
    forAll(DV_, DVid)
        DV_[DVid].updateFlux();  /// delete "*dt/V[own]"
}

void Foam::fvDVM::updateMacroFlux()
{

    rhoFlux_ = dimensionedScalar("0", rhoFlux_.dimensions(), 0);
    rhoUFlux_ = dimensionedVector("0", rhoUFlux_.dimensions(), vector(0,0,0));
    rhoEFlux_ = dimensionedScalar("0", rhoEFlux_.dimensions(), 0);

    forAll(DV_, dvi)
    {
        discreteVelocity& dv = DV_[dvi];
        rhoFlux_ += dXiCellSize_*dv.weight()*dv.gFlux();
        rhoUFlux_ += dXiCellSize_*dv.weight()*dv.gFlux()*dv.xi();
        rhoEFlux_ += 0.5*dXiCellSize_*dv.weight()
            *(
                magSqr(dv.xi())*dv.gFlux() 
                + dv.hFlux()
             );

    }

    if (mpiReducer_.dvParallel() && mpiReducer_.npd() < mpiReducer_.nproc())
    {
        mpiReducer_.reduceField(rhoFlux_);
        mpiReducer_.reduceField(rhoUFlux_);
        mpiReducer_.reduceField(rhoEFlux_);
    }

}
/*
void Foam::fvDVM::PRSGSforwardsweep()
{
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();
    const vectorField Sf = mesh_.Sf();
    const vectorField C = mesh_.C();
    const scalarField V = mesh_.V();

    const scalar Dt = CFL_*time_.deltaTValue();
    scalar rij = 0;
    scalarField Di = V/Dt;
    label nei;

    scalar massFlux_= 0;
    vector momentFlux_= vector(0,0,0);
    scalar energyFlux_= 0;

    scalar newmassFlux_=0;
    vector newmomentFlux_= vector(0,0,0);
    scalar newenergyFlux_=0;

    #if FOAM_MAJOR <= 3
    GeometricField<scalar, fvPatchField, volMesh>::GeometricBoundaryField& 
        rhoBCs = rhoVol_.boundaryField();
    #else
    GeometricField<scalar, fvPatchField, volMesh>::Boundary& 
        rhoBCs = rhoVol_.BOUNDARY_FIELD_REF;
    #endif
    
    // boundary faces first
    for(label patchi = 0; patchi < rhoBCs.size(); ++patchi)
    {
        word type = rhoBCs[patchi].type();
        fvPatchField<scalar>& rhoSurfPatch = rhoVol_.BOUNDARY_FIELD_REF[patchi];
        const fvsPatchField<vector>& SfPatch = mesh_.Sf().boundaryField()[patchi];
        const fvsPatchField<vector>& CfPatch = mesh_.Cf().boundaryField()[patchi];
        const labelUList& faceCells = mesh_.boundary()[patchi].faceCells();

        if (type == "fixedValue" || type == "calculatedMaxwell" || type == "cyclic" || type == "processor")
        {
            //  check each boundary face in the patch
            forAll(rhoSurfPatch, facei)
            {
                vector normal = SfPatch[facei]/mag(SfPatch[facei]);
                rij = fabs(Usurf_[facei] & normal) + sqrt((KInner_ + 5.) / (KInner_ + 3.) * R_.value() * Tsurf_[facei]) 
                    + tauSurf_[facei] * R_.value() * Tsurf_[facei] / fabs((CfPatch[facei] - C[faceCells[facei]]) & normal);

                Di[faceCells[facei]] += 0.5 * rij * mag(SfPatch[facei]);
	        }
        }
	}
    
    for (label cellI = 0; cellI < mesh_.nCells(); ++cellI)
    {
        //Info << "Cell " << cellI << " has the following:";
        
        scalar Dj_rho = 0;
        vector Dj_rhoU = vector(0,0,0);
        scalar Dj_rhoE = 0;
        const labelList& faces = mesh_.cells()[cellI]; 
        forAll(faces, facei)
        {
            label faceIndex = faces[facei]; // 获取面索引

            if (faceIndex < mesh_.nInternalFaces())
            {
                vector normal = Sf[faceIndex]/mag(Sf[faceIndex]);

				if(owner[faceIndex] == cellI)
				{
					nei = neighbour[faceIndex]; // 获取邻居单元
				}
				else
				{
					nei = owner[faceIndex]; 
                    normal = -normal;	// 获取面的法向量
				}
                //Info << " Neighbour Cell " << neighbourCell;
                rij = fabs(Usurf_[faceIndex] & normal) + sqrt((KInner_ + 5.) / (KInner_ + 3.) * R_.value() * Tsurf_[faceIndex]) 
                        + 2.0 * tauSurf_[faceIndex] * R_.value() * Tsurf_[faceIndex] / fabs((C[nei] - C[cellI]) & normal);

                Di[cellI] += 0.5 * rij * mag(Sf[faceIndex]);

                scalar rho_ =  rhoVol_[nei];
                vector rhoU_ = rhoVol_[nei] * Uvol_[nei];
                scalar rhoE_ = 0.5 * rhoVol_[nei]*((3. + KInner_)* R_.value() * Tvol_[nei] + magSqr(Uvol_[nei]));

                scalar rho_new = rho_ + deltaRho_[nei];
                vector rhoU_new = rhoU_ + deltaRhoU_[nei];
                scalar rhoE_new = rhoE_ + deltaRhoE_[nei];

                EulerFlux(massFlux_, momentFlux_, energyFlux_, rho_, rhoU_, rhoE_, normal);
                EulerFlux(newmassFlux_, newmomentFlux_, newenergyFlux_, rho_new, rhoU_new, rhoE_new, normal);
                
                Dj_rho += 0.5*mag(Sf[faceIndex])*(newmassFlux_ - massFlux_ - rij * deltaRho_[nei]);
                Dj_rhoU += 0.5*mag(Sf[faceIndex])*(newmomentFlux_ - momentFlux_ - rij * deltaRhoU_[nei]);
                Dj_rhoE += 0.5*mag(Sf[faceIndex])*(newenergyFlux_ - energyFlux_- rij * deltaRhoE_[nei]);
            }
        }
        deltaRho_[cellI]  = (-rhoFlux_[cellI] - Dj_rho)/ Di[cellI];
        deltaRhoU_[cellI] = (-rhoUFlux_[cellI] - Dj_rhoU)/ Di[cellI];
        deltaRhoE_[cellI] = (-rhoEFlux_[cellI] - Dj_rhoE)/ Di[cellI];
    }

    deltaRho_.correctBoundaryConditions();
    deltaRhoU_.correctBoundaryConditions();
    deltaRhoE_.correctBoundaryConditions();
}*/

// Code Optimization


void Foam::fvDVM::PRSGSforwardsweep()
{
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();
    const vectorField Sf = mesh_.Sf();
    const vectorField C = mesh_.C();
    const scalarField V = mesh_.V();

    const scalar Dt = CFL_*time_.deltaTValue();
    scalar rij = 0;
    scalarField Di = V/Dt;
    label nei;

    const scalar sqrtConst = (KInner_ + 5.0) / (KInner_ + 3.0) * R_.value();

    #if FOAM_MAJOR <= 3
    GeometricField<scalar, fvPatchField, volMesh>::GeometricBoundaryField& 
        rhoBCs = rhoVol_.boundaryField();
    #else
    GeometricField<scalar, fvPatchField, volMesh>::Boundary& 
        rhoBCs = rhoVol_.BOUNDARY_FIELD_REF;
    #endif
    
    // boundary faces first
    for(label patchi = 0; patchi < rhoBCs.size(); ++patchi)
    {
        word type = rhoBCs[patchi].type();
        fvPatchField<scalar>& rhoSurfPatch = rhoVol_.BOUNDARY_FIELD_REF[patchi];
        const fvsPatchField<vector>& SfPatch = mesh_.Sf().boundaryField()[patchi];
        const fvsPatchField<vector>& CfPatch = mesh_.Cf().boundaryField()[patchi];
        const labelUList& faceCells = mesh_.boundary()[patchi].faceCells();

        if (type == "fixedValue" || type == "calculatedMaxwell" || type == "cyclic" || type == "processor" || type == "symmetryMod")
        {
            //  check each boundary face in the patch
            forAll(rhoSurfPatch, facei)
            {
                vector normal = SfPatch[facei]/mag(SfPatch[facei]);
                rij = fabs(Usurf_[facei] & normal) + sqrt(sqrtConst * Tsurf_[facei]) 
                    + tauSurf_[facei] * R_.value() * Tsurf_[facei] / fabs((CfPatch[facei] - C[faceCells[facei]]) & normal);

                Di[faceCells[facei]] += 0.5 * rij * mag(SfPatch[facei]);
                // if(type == "cyclic" || type == "processor")
                // {
                //     //需要补充更新相邻格子内的信息
                // }
	        }
        }
	}
    
    for (label cellI = 0; cellI < mesh_.nCells(); ++cellI)
    {
        //Info << "Cell " << cellI << " has the following:";
        
        scalar Dj_rho = 0;
        vector Dj_rhoU = vector(0,0,0);
        scalar Dj_rhoE = 0;
        const labelList& faces = mesh_.cells()[cellI]; 
        forAll(faces, facei)
        {
            label faceIndex = faces[facei]; // 获取面索引

            if (faceIndex < mesh_.nInternalFaces())
            {
                vector normal = Sf[faceIndex]/mag(Sf[faceIndex]);

				if(owner[faceIndex] == cellI)
				{
					nei = neighbour[faceIndex]; // 获取邻居单元
				}
				else
				{
					nei = owner[faceIndex]; 
                    normal = -normal;	// 获取面的法向量
				}
                //Info << " Neighbour Cell " << neighbourCell;
                rij = fabs(Usurf_[faceIndex] & normal) + sqrt(sqrtConst * Tsurf_[faceIndex]) 
                        + 2.0 * tauSurf_[faceIndex] * R_.value() * Tsurf_[faceIndex] / fabs((C[nei] - C[cellI]) & normal);

                Di[cellI] += 0.5 * rij * mag(Sf[faceIndex]);

                scalar rho_ =  rhoVol_[nei];
                vector rhoU_ = rhoVol_[nei] * Uvol_[nei];
                scalar rhoE_ = 0.5 * rhoVol_[nei]*((3. + KInner_)* R_.value() * Tvol_[nei] + magSqr(Uvol_[nei]));

                scalar rho_new = rho_ + deltaRho_[nei];
                vector rhoU_new = rhoU_ + deltaRhoU_[nei];
                scalar rhoE_new = rhoE_ + deltaRhoE_[nei];

                scalar massFlux_, energyFlux_, newmassFlux_, newenergyFlux_;
                vector momentFlux_, newmomentFlux_;

                EulerFlux(massFlux_, momentFlux_, energyFlux_, rho_, rhoU_, rhoE_, normal);
                EulerFlux(newmassFlux_, newmomentFlux_, newenergyFlux_, rho_new, rhoU_new, rhoE_new, normal);
                
                scalar halfSf = 0.5*mag(Sf[faceIndex]);
                Dj_rho += halfSf *(newmassFlux_ - massFlux_ - rij * deltaRho_[nei]);
                Dj_rhoU += halfSf *(newmomentFlux_ - momentFlux_ - rij * deltaRhoU_[nei]);
                Dj_rhoE += halfSf *(newenergyFlux_ - energyFlux_- rij * deltaRhoE_[nei]);
            }
        }
        deltaRho_[cellI]  = (-rhoFlux_[cellI] - Dj_rho)/ Di[cellI];
        deltaRhoU_[cellI] = (-rhoUFlux_[cellI] - Dj_rhoU)/ Di[cellI];
        deltaRhoE_[cellI] = (-rhoEFlux_[cellI] - Dj_rhoE)/ Di[cellI];
    }

    deltaRho_.correctBoundaryConditions();
    deltaRhoU_.correctBoundaryConditions();
    deltaRhoE_.correctBoundaryConditions();
}
/*
void Foam::fvDVM::PRSGSforwardsweep()
{
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();
    const vectorField& Sf = mesh_.Sf();
    const vectorField& C = mesh_.C();
    const scalarField& V = mesh_.V();

    const scalar Dt = CFL_*time_.deltaTValue();
    scalar rij = 0;
    scalarField Di = V/Dt;
    label nei;
    const scalar sqrtConst = (KInner_ + 5.0) / (KInner_ + 3.0) * R_.value();

    #if FOAM_MAJOR <= 3
    GeometricField<scalar, fvPatchField, volMesh>::GeometricBoundaryField& 
        rhoBCs = rhoVol_.boundaryField();
    #else
    GeometricField<scalar, fvPatchField, volMesh>::Boundary& 
        rhoBCs = rhoVol_.BOUNDARY_FIELD_REF;
    #endif

    // ---------------------------------------------------------------
    // 预先缓存 cyclic/processor patch 的邻居值
    // ---------------------------------------------------------------
    PtrList<scalarField> rhoNeiFields(rhoBCs.size());
    PtrList<vectorField> UNeiFields(rhoBCs.size());
    PtrList<scalarField> TNeiFields(rhoBCs.size());
    PtrList<scalarField> deltaRhoNeiFields(rhoBCs.size());
    PtrList<vectorField> deltaRhoUNeiFields(rhoBCs.size());
    PtrList<scalarField> deltaRhoENeiFields(rhoBCs.size());

    forAll(rhoBCs, patchi)
    {
        const word& ptype = rhoBCs[patchi].type();
        if (ptype == "cyclic" || ptype == "processor")
        {
            rhoNeiFields.set(patchi,
                new scalarField(rhoBCs[patchi].patchNeighbourField()));
            UNeiFields.set(patchi,
                new vectorField(Uvol_.boundaryField()[patchi].patchNeighbourField()));
            TNeiFields.set(patchi,
                new scalarField(Tvol_.boundaryField()[patchi].patchNeighbourField()));
            deltaRhoNeiFields.set(patchi,
                new scalarField(deltaRho_.boundaryField()[patchi].patchNeighbourField()));
            deltaRhoUNeiFields.set(patchi,
                new vectorField(deltaRhoU_.boundaryField()[patchi].patchNeighbourField()));
            deltaRhoENeiFields.set(patchi,
                new scalarField(deltaRhoE_.boundaryField()[patchi].patchNeighbourField()));
        }
    }

    // ---------------------------------------------------------------
    // boundary faces: 更新 Di（与原代码相同，顺序遍历）
    // ---------------------------------------------------------------
    for (label patchi = 0; patchi < rhoBCs.size(); ++patchi)
    {
        word type = rhoBCs[patchi].type();
        fvPatchField<scalar>& rhoSurfPatch = rhoVol_.BOUNDARY_FIELD_REF[patchi];
        const fvsPatchField<vector>& SfPatch = mesh_.Sf().boundaryField()[patchi];
        const fvsPatchField<vector>& CfPatch = mesh_.Cf().boundaryField()[patchi];
        const labelUList& faceCells = mesh_.boundary()[patchi].faceCells();

        if (type == "fixedValue" || type == "calculatedMaxwell"
         || type == "cyclic"     || type == "processor" || type == "symmetryMod")
        {
            forAll(rhoSurfPatch, facei)
            {
                vector normal = SfPatch[facei] / mag(SfPatch[facei]);
                rij = fabs(Usurf_[facei] & normal)
                    + sqrt(sqrtConst * Tsurf_[facei])
                    + tauSurf_[facei] * R_.value() * Tsurf_[facei]
                      / fabs((CfPatch[facei] - C[faceCells[facei]]) & normal);

                Di[faceCells[facei]] += 0.5 * rij * mag(SfPatch[facei]);
            }
        }
    }

    // ---------------------------------------------------------------
    // 主循环：forward sweep（从 0 到 nCells-1）
    // ---------------------------------------------------------------
    for (label cellI = 0; cellI < mesh_.nCells(); ++cellI)
    {
        scalar Dj_rho  = 0;
        vector Dj_rhoU = vector(0, 0, 0);
        scalar Dj_rhoE = 0;

        const labelList& faces = mesh_.cells()[cellI];
        forAll(faces, facei)
        {
            label faceIndex = faces[facei];

            if (faceIndex < mesh_.nInternalFaces())
            {
                // ---------- 内部面，逻辑不变 ----------
                vector normal = Sf[faceIndex] / mag(Sf[faceIndex]);

                if (owner[faceIndex] == cellI)
                {
                    nei = neighbour[faceIndex];
                }
                else
                {
                    nei = owner[faceIndex];
                    normal = -normal;
                }

                rij = fabs(Usurf_[faceIndex] & normal)
                    + sqrt(sqrtConst * Tsurf_[faceIndex])
                    + 2.0 * tauSurf_[faceIndex] * R_.value() * Tsurf_[faceIndex]
                      / fabs((C[nei] - C[cellI]) & normal);

                Di[cellI] += 0.5 * rij * mag(Sf[faceIndex]);

                scalar rho_  = rhoVol_[nei];
                vector rhoU_ = rhoVol_[nei] * Uvol_[nei];
                scalar rhoE_ = 0.5 * rhoVol_[nei] * (
                                   (3. + KInner_) * R_.value() * Tvol_[nei]
                                   + magSqr(Uvol_[nei]));

                scalar rho_new  = rho_  + deltaRho_[nei];
                vector rhoU_new = rhoU_ + deltaRhoU_[nei];
                scalar rhoE_new = rhoE_ + deltaRhoE_[nei];

                scalar massFlux_, energyFlux_, newmassFlux_, newenergyFlux_;
                vector momentFlux_, newmomentFlux_;

                EulerFlux(massFlux_,    momentFlux_,    energyFlux_,
                          rho_,    rhoU_,    rhoE_,    normal);
                EulerFlux(newmassFlux_, newmomentFlux_, newenergyFlux_,
                          rho_new, rhoU_new, rhoE_new, normal);

                scalar halfSf = 0.5 * mag(Sf[faceIndex]);
                Dj_rho  += halfSf * (newmassFlux_  - massFlux_  - rij * deltaRho_[nei]);
                Dj_rhoU += halfSf * (newmomentFlux_ - momentFlux_ - rij * deltaRhoU_[nei]);
                Dj_rhoE += halfSf * (newenergyFlux_ - energyFlux_ - rij * deltaRhoE_[nei]);
            }
            else
            {
                // ---------- 边界面：cyclic/processor 补充邻居贡献 ----------
                label patchi     = mesh_.boundaryMesh().whichPatch(faceIndex);
                label patchFacei = faceIndex - mesh_.boundaryMesh()[patchi].start();

                const word& ptype = rhoBCs[patchi].type();

                if (ptype == "cyclic" || ptype == "processor")
                {
                    const fvsPatchField<vector>& SfPatch =
                        mesh_.Sf().boundaryField()[patchi];
                    const fvsPatchField<vector>& CfPatch =
                        mesh_.Cf().boundaryField()[patchi];

                    vector normal = SfPatch[patchFacei] / mag(SfPatch[patchFacei]);
                    scalar halfSf = 0.5 * mag(SfPatch[patchFacei]);

                    scalar rho_  = rhoNeiFields[patchi][patchFacei];
                    vector U_nei = UNeiFields[patchi][patchFacei];
                    scalar T_nei = TNeiFields[patchi][patchFacei];

                    vector rhoU_ = rho_ * U_nei;
                    scalar rhoE_ = 0.5 * rho_ * (
                                       (3. + KInner_) * R_.value() * T_nei
                                       + magSqr(U_nei));

                    scalar deltaRho_nei  = deltaRhoNeiFields[patchi][patchFacei];
                    vector deltaRhoU_nei = deltaRhoUNeiFields[patchi][patchFacei];
                    scalar deltaRhoE_nei = deltaRhoENeiFields[patchi][patchFacei];

                    scalar rho_new  = rho_  + deltaRho_nei;
                    vector rhoU_new = rhoU_ + deltaRhoU_nei;
                    scalar rhoE_new = rhoE_ + deltaRhoE_nei;

                    scalar dist = fabs((CfPatch[patchFacei] - C[cellI]) & normal);
                    rij = fabs(Usurf_[faceIndex] & normal)
                        + sqrt(sqrtConst * Tsurf_[faceIndex])
                        + 2.0 * tauSurf_[faceIndex] * R_.value()
                          * Tsurf_[faceIndex] / dist;

                    scalar massFlux_, energyFlux_, newmassFlux_, newenergyFlux_;
                    vector momentFlux_, newmomentFlux_;

                    EulerFlux(massFlux_,    momentFlux_,    energyFlux_,
                              rho_,    rhoU_,    rhoE_,    normal);
                    EulerFlux(newmassFlux_, newmomentFlux_, newenergyFlux_,
                              rho_new, rhoU_new, rhoE_new, normal);

                    Dj_rho  += halfSf * (newmassFlux_  - massFlux_  - rij * deltaRho_nei);
                    Dj_rhoU += halfSf * (newmomentFlux_ - momentFlux_ - rij * deltaRhoU_nei);
                    Dj_rhoE += halfSf * (newenergyFlux_ - energyFlux_ - rij * deltaRhoE_nei);
                }
                // fixedValue/calculatedMaxwell: 不累加 Dj，Di 已在上方更新
            }
        }

        deltaRho_[cellI]  = (-rhoFlux_[cellI]  - Dj_rho)  / Di[cellI];
        deltaRhoU_[cellI] = (-rhoUFlux_[cellI] - Dj_rhoU) / Di[cellI];
        deltaRhoE_[cellI] = (-rhoEFlux_[cellI] - Dj_rhoE) / Di[cellI];
    }

    deltaRho_.correctBoundaryConditions();
    deltaRhoU_.correctBoundaryConditions();
    deltaRhoE_.correctBoundaryConditions();
}
void Foam::fvDVM::PRSGSbackwardsweep()
{
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();
    const vectorField& Sf = mesh_.Sf();
    const vectorField& C = mesh_.C();
    const scalarField& V = mesh_.V();

    const scalar Dt = CFL_*time_.deltaTValue();
    scalar rij = 0;
    scalarField Di = V/Dt;
    label nei;
    const scalar sqrtConst = (KInner_ + 5.0) / (KInner_ + 3.0) * R_.value();

    #if FOAM_MAJOR <= 3
    GeometricField<scalar, fvPatchField, volMesh>::GeometricBoundaryField& 
        rhoBCs = rhoVol_.boundaryField();
    #else
    GeometricField<scalar, fvPatchField, volMesh>::Boundary& 
        rhoBCs = rhoVol_.BOUNDARY_FIELD_REF;
    #endif

    // ---------------------------------------------------------------
    // 预先取好所有 patch 的 patchNeighbourField，避免在循环内重复调用
    // ---------------------------------------------------------------
    PtrList<scalarField> rhoNeiFields(rhoBCs.size());
    PtrList<vectorField> UNeiFields(rhoBCs.size());
    PtrList<scalarField> TNeiFields(rhoBCs.size());
    PtrList<scalarField> deltaRhoNeiFields(rhoBCs.size());
    PtrList<vectorField> deltaRhoUNeiFields(rhoBCs.size());
    PtrList<scalarField> deltaRhoENeiFields(rhoBCs.size());

    forAll(rhoBCs, patchi)
    {
        const word& ptype = rhoBCs[patchi].type();
        if (ptype == "cyclic" || ptype == "processor")
        {
            rhoNeiFields.set(patchi,
                new scalarField(rhoBCs[patchi].patchNeighbourField()));
            UNeiFields.set(patchi,
                new vectorField(Uvol_.boundaryField()[patchi].patchNeighbourField()));
            TNeiFields.set(patchi,
                new scalarField(Tvol_.boundaryField()[patchi].patchNeighbourField()));
            deltaRhoNeiFields.set(patchi,
                new scalarField(deltaRho_.boundaryField()[patchi].patchNeighbourField()));
            deltaRhoUNeiFields.set(patchi,
                new vectorField(deltaRhoU_.boundaryField()[patchi].patchNeighbourField()));
            deltaRhoENeiFields.set(patchi,
                new scalarField(deltaRhoE_.boundaryField()[patchi].patchNeighbourField()));
        }
    }

    // ---------------------------------------------------------------
    // 构建 boundary face -> (patchi, patchFacei) 的查找表
    // ---------------------------------------------------------------
    // OF6 中 mesh_.boundaryMesh() 提供了 whichPatch(faceIndex) 方法
    // faceIndex >= nInternalFaces 时使用

    // boundary faces first (Di 更新，保持原逻辑)
    for (label patchi = rhoBCs.size()-1; patchi >= 0; --patchi)
    {
        word type = rhoBCs[patchi].type();
        fvPatchField<scalar>& rhoSurfPatch = rhoVol_.BOUNDARY_FIELD_REF[patchi];
        const fvsPatchField<vector>& SfPatch = mesh_.Sf().boundaryField()[patchi];
        const fvsPatchField<vector>& CfPatch = mesh_.Cf().boundaryField()[patchi];
        const labelUList& faceCells = mesh_.boundary()[patchi].faceCells();

        if (type == "fixedValue" || type == "calculatedMaxwell" 
         || type == "cyclic"     || type == "processor" || type == "symmetryMod")
        {
            forAll(rhoSurfPatch, facei)
            {
                vector normal = SfPatch[facei] / mag(SfPatch[facei]);
                rij = fabs(Usurf_[facei] & normal) 
                    + sqrt(sqrtConst * Tsurf_[facei]) 
                    + tauSurf_[facei] * R_.value() * Tsurf_[facei] 
                      / fabs((CfPatch[facei] - C[faceCells[facei]]) & normal);

                Di[faceCells[facei]] += 0.5 * rij * mag(SfPatch[facei]);
            }
        }
    }

    // ---------------------------------------------------------------
    // 主循环：backward sweep
    // ---------------------------------------------------------------
    for (label cellI = mesh_.nCells() - 1; cellI >= 0; --cellI)
    {
        scalar Dj_rho  = 0;
        vector Dj_rhoU = vector(0, 0, 0);
        scalar Dj_rhoE = 0;

        const labelList& faces = mesh_.cells()[cellI];
        forAll(faces, facei)
        {
            label faceIndex = faces[facei];

            if (faceIndex < mesh_.nInternalFaces())
            {
                // ---------- 内部面，逻辑不变 ----------
                vector normal = Sf[faceIndex] / mag(Sf[faceIndex]);

                if (owner[faceIndex] == cellI)
                {
                    nei = neighbour[faceIndex];
                }
                else
                {
                    nei = owner[faceIndex];
                    normal = -normal;
                }

                rij = fabs(Usurf_[faceIndex] & normal) 
                    + sqrt(sqrtConst * Tsurf_[faceIndex]) 
                    + 2.0 * tauSurf_[faceIndex] * R_.value() * Tsurf_[faceIndex] 
                      / fabs((C[nei] - C[cellI]) & normal);

                Di[cellI] += 0.5 * rij * mag(Sf[faceIndex]);

                scalar rho_  = rhoVol_[nei];
                vector rhoU_ = rhoVol_[nei] * Uvol_[nei];
                scalar rhoE_ = 0.5 * rhoVol_[nei] * (
                                   (3. + KInner_) * R_.value() * Tvol_[nei] 
                                   + magSqr(Uvol_[nei]));

                scalar rho_new  = rho_  + deltaRho_[nei];
                vector rhoU_new = rhoU_ + deltaRhoU_[nei];
                scalar rhoE_new = rhoE_ + deltaRhoE_[nei];

                scalar massFlux_, energyFlux_, newmassFlux_, newenergyFlux_;
                vector momentFlux_, newmomentFlux_;

                EulerFlux(massFlux_,    momentFlux_,    energyFlux_,
                          rho_,    rhoU_,    rhoE_,    normal);
                EulerFlux(newmassFlux_, newmomentFlux_, newenergyFlux_,
                          rho_new, rhoU_new, rhoE_new, normal);

                scalar halfSf = 0.5 * mag(Sf[faceIndex]);
                Dj_rho  += halfSf * (newmassFlux_  - massFlux_  - rij * deltaRho_[nei]);
                Dj_rhoU += halfSf * (newmomentFlux_ - momentFlux_ - rij * deltaRhoU_[nei]);
                Dj_rhoE += halfSf * (newenergyFlux_ - energyFlux_ - rij * deltaRhoE_[nei]);
            }
            else
            {
                // ---------- 边界面 ----------
                // 找到该面属于哪个 patch 和 patch 内的局部索引
                label patchi    = mesh_.boundaryMesh().whichPatch(faceIndex);
                label patchFacei = faceIndex 
                                 - mesh_.boundaryMesh()[patchi].start();

                const word& ptype = rhoBCs[patchi].type();

                // 只有 cyclic/processor 才有邻居值需要处理
                if (ptype == "cyclic" || ptype == "processor")
                {
                    const fvsPatchField<vector>& SfPatch =
                        mesh_.Sf().boundaryField()[patchi];
                    const fvsPatchField<vector>& CfPatch =
                        mesh_.Cf().boundaryField()[patchi];

                    vector normal = SfPatch[patchFacei] / mag(SfPatch[patchFacei]);
                    scalar halfSf = 0.5 * mag(SfPatch[patchFacei]);

                    // 用预先缓存好的 patchNeighbourField 取邻居值
                    scalar rho_  = rhoNeiFields[patchi][patchFacei];
                    vector U_nei = UNeiFields[patchi][patchFacei];
                    scalar T_nei = TNeiFields[patchi][patchFacei];

                    vector rhoU_ = rho_ * U_nei;
                    scalar rhoE_ = 0.5 * rho_ * (
                                       (3. + KInner_) * R_.value() * T_nei
                                       + magSqr(U_nei));

                    scalar deltaRho_nei  = deltaRhoNeiFields[patchi][patchFacei];
                    vector deltaRhoU_nei = deltaRhoUNeiFields[patchi][patchFacei];
                    scalar deltaRhoE_nei = deltaRhoENeiFields[patchi][patchFacei];

                    scalar rho_new  = rho_  + deltaRho_nei;
                    vector rhoU_new = rhoU_ + deltaRhoU_nei;
                    scalar rhoE_new = rhoE_ + deltaRhoE_nei;

                    // rij：用 Cf 和 cellI 中心距离估算
                    scalar dist = fabs((CfPatch[patchFacei] - C[cellI]) & normal);
                    rij = fabs(Usurf_[faceIndex] & normal)   // 注意：这里需要你确认 Usurf_ 的索引方式
                        + sqrt(sqrtConst * Tsurf_[faceIndex])
                        + 2.0 * tauSurf_[faceIndex] * R_.value() 
                          * Tsurf_[faceIndex] / dist;

                    scalar massFlux_, energyFlux_, newmassFlux_, newenergyFlux_;
                    vector momentFlux_, newmomentFlux_;

                    EulerFlux(massFlux_,    momentFlux_,    energyFlux_,
                              rho_,    rhoU_,    rhoE_,    normal);
                    EulerFlux(newmassFlux_, newmomentFlux_, newenergyFlux_,
                              rho_new, rhoU_new, rhoE_new, normal);

                    Dj_rho  += halfSf * (newmassFlux_  - massFlux_  - rij * deltaRho_nei);
                    Dj_rhoU += halfSf * (newmomentFlux_ - momentFlux_ - rij * deltaRhoU_nei);
                    Dj_rhoE += halfSf * (newenergyFlux_ - energyFlux_ - rij * deltaRhoE_nei);
                }
                // fixedValue / calculatedMaxwell 等物理边界：
                // 不需要邻居 delta 贡献，Dj 不累加（Di 已在上方更新过）
            }
        }

        deltaRho_[cellI]  = (-rhoFlux_[cellI]  - Dj_rho)  / Di[cellI];
        deltaRhoU_[cellI] = (-rhoUFlux_[cellI] - Dj_rhoU) / Di[cellI];
        deltaRhoE_[cellI] = (-rhoEFlux_[cellI] - Dj_rhoE) / Di[cellI];
    }

    deltaRho_.correctBoundaryConditions();
    deltaRhoU_.correctBoundaryConditions();
    deltaRhoE_.correctBoundaryConditions();
}
*/

void Foam::fvDVM::PRSGSbackwardsweep()
{
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();
    const vectorField Sf = mesh_.Sf();
    const vectorField C = mesh_.C();
    const scalarField V = mesh_.V();

    const scalar Dt = CFL_*time_.deltaTValue();
    scalar rij = 0;
    scalarField Di = V/Dt;
    label nei;
    const scalar sqrtConst = (KInner_ + 5.0) / (KInner_ + 3.0) * R_.value();

    #if FOAM_MAJOR <= 3
    GeometricField<scalar, fvPatchField, volMesh>::GeometricBoundaryField& 
        rhoBCs = rhoVol_.boundaryField();
    #else
    GeometricField<scalar, fvPatchField, volMesh>::Boundary& 
        rhoBCs = rhoVol_.BOUNDARY_FIELD_REF;
    #endif
    
    // boundary faces first
    //forAll(rhoBCs, patchi)
    for(label patchi = rhoBCs.size()-1; patchi >=0; --patchi)
    {
        word type = rhoBCs[patchi].type();
        fvPatchField<scalar>& rhoSurfPatch = rhoVol_.BOUNDARY_FIELD_REF[patchi];
        const fvsPatchField<vector>& SfPatch = mesh_.Sf().boundaryField()[patchi];
        const fvsPatchField<vector>& CfPatch = mesh_.Cf().boundaryField()[patchi];
        const labelUList& faceCells = mesh_.boundary()[patchi].faceCells();

        if (type == "fixedValue" || type == "calculatedMaxwell" || type == "cyclic" || type == "processor" || type == "symmetryMod")
        {
            //  check each boundary face in the patch
            forAll(rhoSurfPatch, facei)
            {
                vector normal = SfPatch[facei]/mag(SfPatch[facei]);
                rij = fabs(Usurf_[facei] & normal) + sqrt(sqrtConst * Tsurf_[facei]) 
                    + tauSurf_[facei] * R_.value() * Tsurf_[facei] / fabs((CfPatch[facei] - C[faceCells[facei]]) & normal);

                Di[faceCells[facei]] += 0.5 * rij * mag(SfPatch[facei]);
                // if(type == "cyclic" || type == "processor")
                // {
                //     //需要补充更新相邻格子内的信息
                // }
	        }
        }
	}
    
    for (label cellI = mesh_.nCells() - 1; cellI >= 0; --cellI)
    {
        //Info << "Cell " << cellI << " has the following:";
        
        scalar Dj_rho = 0;
        vector Dj_rhoU = vector(0,0,0);
        scalar Dj_rhoE = 0;
        const labelList& faces = mesh_.cells()[cellI]; 
        forAll(faces, facei)
        {
            label faceIndex = faces[facei]; // 获取面索引

            if (faceIndex < mesh_.nInternalFaces())
            {
                vector normal = Sf[faceIndex]/mag(Sf[faceIndex]);

				if(owner[faceIndex] == cellI)
				{
					nei = neighbour[faceIndex]; // 获取邻居单元
				}
				else
				{
					nei = owner[faceIndex]; 
                    normal = -normal;	// 获取面的法向量
				}
                //Info << " Neighbour Cell " << neighbourCell;
                rij = fabs(Usurf_[faceIndex] & normal) + sqrt(sqrtConst * Tsurf_[faceIndex]) 
                        + 2.0 * tauSurf_[faceIndex] * R_.value() * Tsurf_[faceIndex] / fabs((C[nei] - C[cellI]) & normal);

                Di[cellI] += 0.5 * rij * mag(Sf[faceIndex]);

                scalar rho_ =  rhoVol_[nei];
                vector rhoU_ = rhoVol_[nei] * Uvol_[nei];
                scalar rhoE_ = 0.5 * rhoVol_[nei]*((3. + KInner_)* R_.value() * Tvol_[nei] + magSqr(Uvol_[nei]));

                scalar rho_new = rho_ + deltaRho_[nei];
                vector rhoU_new = rhoU_ + deltaRhoU_[nei];
                scalar rhoE_new = rhoE_ + deltaRhoE_[nei];

                scalar massFlux_, energyFlux_, newmassFlux_, newenergyFlux_;
                vector momentFlux_, newmomentFlux_;

                EulerFlux(massFlux_, momentFlux_, energyFlux_, rho_, rhoU_, rhoE_, normal);
                EulerFlux(newmassFlux_, newmomentFlux_, newenergyFlux_, rho_new, rhoU_new, rhoE_new, normal);
                
                scalar halfSf = 0.5*mag(Sf[faceIndex]);
                Dj_rho += halfSf *(newmassFlux_ - massFlux_ - rij * deltaRho_[nei]);
                Dj_rhoU += halfSf *(newmomentFlux_ - momentFlux_ - rij * deltaRhoU_[nei]);
                Dj_rhoE += halfSf *(newenergyFlux_ - energyFlux_- rij * deltaRhoE_[nei]);
            }
        }
        deltaRho_[cellI]  = (-rhoFlux_[cellI] - Dj_rho)/ Di[cellI];
        deltaRhoU_[cellI] = (-rhoUFlux_[cellI] - Dj_rhoU)/ Di[cellI];
        deltaRhoE_[cellI] = (-rhoEFlux_[cellI] - Dj_rhoE)/ Di[cellI];
    }

    deltaRho_.correctBoundaryConditions();
    deltaRhoU_.correctBoundaryConditions();
    deltaRhoE_.correctBoundaryConditions();
}
/*
void Foam::fvDVM::PRSGSbackwardsweep()
{
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();
    const vectorField Sf = mesh_.Sf();
    const vectorField C = mesh_.C();
    const scalarField V = mesh_.V();

    const scalar Dt = CFL_*time_.deltaTValue();
    scalar rij = 0;
    scalarField Di = V/Dt;
    label nei;

    scalar massFlux_= 0;
    vector momentFlux_= vector(0,0,0);
    scalar energyFlux_= 0;

    scalar newmassFlux_=0;
    vector newmomentFlux_= vector(0,0,0);
    scalar newenergyFlux_=0;

    #if FOAM_MAJOR <= 3
    GeometricField<scalar, fvPatchField, volMesh>::GeometricBoundaryField& 
        rhoBCs = rhoVol_.boundaryField();
    #else
    GeometricField<scalar, fvPatchField, volMesh>::Boundary& 
        rhoBCs = rhoVol_.BOUNDARY_FIELD_REF;
    #endif
    
    // boundary faces first
    //forAll(rhoBCs, patchi)
    for(label patchi = rhoBCs.size()-1; patchi >=0; --patchi)
    {
        word type = rhoBCs[patchi].type();
        fvPatchField<scalar>& rhoSurfPatch = rhoVol_.BOUNDARY_FIELD_REF[patchi];
        const fvsPatchField<vector>& SfPatch = mesh_.Sf().boundaryField()[patchi];
        const fvsPatchField<vector>& CfPatch = mesh_.Cf().boundaryField()[patchi];
        const labelUList& faceCells = mesh_.boundary()[patchi].faceCells();

        if (type == "fixedValue" || type == "calculatedMaxwell" || type == "cyclic" || type == "processor")
        {
            //  check each boundary face in the patch
            forAll(rhoSurfPatch, facei)
            {
                vector normal = SfPatch[facei]/mag(SfPatch[facei]);
                rij = fabs(Usurf_[facei] & normal) + sqrt((KInner_ + 5.) / (KInner_ + 3.) * R_.value() * Tsurf_[facei]) 
                    + tauSurf_[facei] * R_.value() * Tsurf_[facei] / fabs((CfPatch[facei] - C[faceCells[facei]]) & normal);

                Di[faceCells[facei]] += 0.5 * rij * mag(SfPatch[facei]);
	        }
        }
	}
    
    for (label cellI = mesh_.nCells() - 1; cellI >= 0; --cellI)
    {
        //Info << "Cell " << cellI << " has the following:";
        
        scalar Dj_rho = 0;
        vector Dj_rhoU = vector(0,0,0);
        scalar Dj_rhoE = 0;
        const labelList& faces = mesh_.cells()[cellI]; 
        forAll(faces, facei)
        {
            label faceIndex = faces[facei]; // 获取面索引

            if (faceIndex < mesh_.nInternalFaces())
            {
                vector normal = Sf[faceIndex]/mag(Sf[faceIndex]);

				if(owner[faceIndex] == cellI)
				{
					nei = neighbour[faceIndex]; // 获取邻居单元
				}
				else
				{
					nei = owner[faceIndex]; 
                    normal = -normal;	// 获取面的法向量
				}
                //Info << " Neighbour Cell " << neighbourCell;
                rij = fabs(Usurf_[faceIndex] & normal) + sqrt((KInner_ + 5.) / (KInner_ + 3.) * R_.value() * Tsurf_[faceIndex]) 
                        + 2.0 * tauSurf_[faceIndex] * R_.value() * Tsurf_[faceIndex] / fabs((C[nei] - C[cellI]) & normal);

                Di[cellI] += 0.5 * rij * mag(Sf[faceIndex]);

                scalar rho_ =  rhoVol_[nei];
                vector rhoU_ = rhoVol_[nei] * Uvol_[nei];
                scalar rhoE_ = 0.5 * rhoVol_[nei]*((3. + KInner_)* R_.value() * Tvol_[nei] + magSqr(Uvol_[nei]));

                scalar rho_new = rho_ + deltaRho_[nei];
                vector rhoU_new = rhoU_ + deltaRhoU_[nei];
                scalar rhoE_new = rhoE_ + deltaRhoE_[nei];

                EulerFlux(massFlux_, momentFlux_, energyFlux_, rho_, rhoU_, rhoE_, normal);
                EulerFlux(newmassFlux_, newmomentFlux_, newenergyFlux_, rho_new, rhoU_new, rhoE_new, normal);
                
                Dj_rho += 0.5*mag(Sf[faceIndex])*(newmassFlux_ - massFlux_ - rij * deltaRho_[nei]);
                Dj_rhoU += 0.5*mag(Sf[faceIndex])*(newmomentFlux_ - momentFlux_ - rij * deltaRhoU_[nei]);
                Dj_rhoE += 0.5*mag(Sf[faceIndex])*(newenergyFlux_ - energyFlux_- rij * deltaRhoE_[nei]);
            }
        }
        deltaRho_[cellI]  = (-rhoFlux_[cellI] - Dj_rho)/ Di[cellI];
        deltaRhoU_[cellI] = (-rhoUFlux_[cellI] - Dj_rhoU)/ Di[cellI];
        deltaRhoE_[cellI] = (-rhoEFlux_[cellI] - Dj_rhoE)/ Di[cellI];
    }

    deltaRho_.correctBoundaryConditions();
    deltaRhoU_.correctBoundaryConditions();
    deltaRhoE_.correctBoundaryConditions();
}
*/

void Foam::fvDVM::predicMacroVol()
{
    volVectorField rhoU_ = rhoVol_ * Uvol_;
    volScalarField rhoE_ = 0.5 * rhoVol_*((3. + KInner_)* R_ * Tvol_ + magSqr(Uvol_));

    rhoVol_ += deltaRho_;
    rhoU_ += deltaRhoU_;
    rhoE_ += deltaRhoE_;

    Uvol_ = rhoU_/rhoVol_;
    Tvol_ = (2.0 * rhoE_ - rhoVol_*magSqr(Uvol_))/(KInner_ + 3.) / rhoVol_/ R_;

    // rhoVol_ += deltaRho_;
    // volVectorField rhoU_ = (rhoVol_ * Uvol_) + deltaRhoU_;
    // volScalarField rhoE_ = (0.5 * rhoVol_*((3. + KInner_)* R_ * Tvol_ + magSqr(Uvol_))) + deltaRhoE_;
    // Uvol_ = rhoU_/rhoVol_;
    // Tvol_ = (2.0 * rhoE_ - rhoVol_*magSqr(Uvol_))/(KInner_ + 3.) / rhoVol_/ R_;

    //- Correct the macro field boundary conition
    //rhoVol_.correctBoundaryConditions();
    Uvol_.correctBoundaryConditions();
    Tvol_.correctBoundaryConditions();

    updateTau(tauVol_, Tvol_, rhoVol_);

    deltaRho_ = dimensionedScalar("0", rhoVol_.dimensions(), 0.0);
    deltaRhoU_ = dimensionedVector("0", deltaRhoU_.dimensions(), vector(0,0,0));
    deltaRhoE_ = dimensionedScalar("0", deltaRhoE_.dimensions(), 0.0);
}

void Foam::fvDVM::updateequilibriumShakhov()
{
    forAll(DV_, DVid)
        DV_[DVid].updateequilibriumShakhov();  
}

void Foam::fvDVM::forwardsweep()
{
    forAll(DV_, DVid)
        DV_[DVid].forwardsweep();  
}

void Foam::fvDVM::backwardsweep()
{
    forAll(DV_, DVid)
        DV_[DVid].backwardsweep();  
}

void Foam::fvDVM::updateGHVol()
{
    forAll(DV_, DVid)
        DV_[DVid].updateGHVol();  
}

void Foam::fvDVM::updateMacroVol()   
{
    volScalarField rhoold = rhoVol_;
    volVectorField rhoUold = rhoVol_ * Uvol_;
    volScalarField rhoEold = 0.5 * rhoVol_*((3. + KInner_)* R_ * Tvol_ + magSqr(Uvol_));
    

    volScalarField sumrhovol = rhoold;
    volVectorField sumrhoUvol = rhoUold;
    volScalarField sumrhoEvol = rhoEold;
    volVectorField sumq = qVol_;

    //     //- Old macros, used only if we update using macro fluxes.
    //     //Info << "update according to micro flux" << endl;
        
    sumrhovol = dimensionedScalar("0", rhoVol_.dimensions(), 0);
    sumrhoUvol = dimensionedVector("0", rhoUold.dimensions(), vector(0,0,0));
    sumrhoEvol = dimensionedScalar("0", rhoEold.dimensions(), 0);

    forAll(DV_, dvi)
    {
        discreteVelocity& dv = DV_[dvi];
        sumrhovol += dXiCellSize_*dv.weight()*(dv.gVol() - dv.gTildeVol());
        sumrhoUvol += dXiCellSize_*dv.weight()*(dv.gVol() - dv.gTildeVol())*dv.xi();
        sumrhoEvol += 0.5*dXiCellSize_*dv.weight()
            *(
                magSqr(dv.xi())*(dv.gVol() - dv.gTildeVol()) + (dv.hVol() - dv.hTildeVol())
             );
    }


    if (mpiReducer_.dvParallel() && mpiReducer_.npd() < mpiReducer_.nproc())
    {
        mpiReducer_.reduceField(sumrhovol);
        mpiReducer_.reduceField(sumrhoUvol);
        mpiReducer_.reduceField(sumrhoEvol);
    }
        
    rhoVol_ = rhoold + sumrhovol;

    Uvol_ = (rhoUold + sumrhoUvol) / rhoVol_;

    Tvol_ = (2.*(rhoEold + sumrhoEvol) - rhoVol_*magSqr(Uvol_))/((KInner_ + 3) * R_ * rhoVol_);

    // //- Correct the macro field boundary conition
    Uvol_.correctBoundaryConditions();
    Tvol_.correctBoundaryConditions();
    // //- Note for maxwell wall, the operation here update 
    // //- the boundary rho field but it's meaningless.
    updateTau(tauVol_, Tvol_, rhoVol_);


    //- peculiar vel.
    volVectorField c = Uvol_;
    //
    sumq = dimensionedVector("0", qVol_.dimensions(), vector(0, 0, 0));
    // stress_ = dimensionedTensor
    //     (
    //         "0",
    //         stress_.dimensions(),
    //         pTraits<tensor>::zero
    //     );
    forAll(DV_, dvi)
    {
        discreteVelocity& dv = DV_[dvi];
        c = dv.xi() - Uvol_;
        sumq += 0.5*dXiCellSize_*dv.weight()*c
            *(
                magSqr(c)*(dv.gVol() - dv.gTildeVol()) + (dv.hVol() - dv.hTildeVol())
             );
        // stress_ +=
        //     dXiCellSize_*dv.weight()*(dv.gVol() - dv.gTildeVol())*c*c;
    }

    // //- get global heat flux via MPI_Allreduce
    if (mpiReducer_.dvParallel()&& mpiReducer_.npd() < mpiReducer_.nproc())
    {
        mpiReducer_.reduceField(sumq);
        // mpiReducer_.reduceField(stress_);
    }
    qVol_ = sumq + qVol_ * (1 - Pr_);
}

void Foam::fvDVM::updatePressureInOutBC()
{
    // for pressureIn and pressureOut BC, the boundary value of Uvol(in/out) and Tvol(in/out) should be updated here!
    // boundary faces
#if FOAM_MAJOR <= 3
    GeometricField<scalar, fvPatchField, volMesh>::GeometricBoundaryField& 
        rhoBCs = rhoVol_.boundaryField();
#else
    GeometricField<scalar, fvPatchField, volMesh>::Boundary& 
        rhoBCs = rhoVol_.BOUNDARY_FIELD_REF;
#endif
    forAll(rhoBCs, patchi)
    {
        if (rhoBCs[patchi].type() == "pressureIn")
        {
            const fvsPatchField<vector>& SfPatch = mesh_.Sf().boundaryField()[patchi];
            const fvsPatchField<scalar>& magSfPatch = mesh_.magSf().boundaryField()[patchi];
            pressureInFvPatchField<scalar>& rhoPatch = 
                refCast<pressureInFvPatchField<scalar> >(rhoBCs[patchi]);
            fvPatchField<vector>& Upatch = Uvol_.BOUNDARY_FIELD_REF[patchi];
            const fvPatchField<scalar>& Tpatch = Tvol_.boundaryField()[patchi];
            const scalar pressureIn = rhoPatch.pressureIn();
            // now changed rho and U patch
            const labelUList& pOwner = mesh_.boundary()[patchi].faceCells();
            forAll(rhoPatch, facei)
            {
                const scalar  Tin = Tpatch[facei];
                // change density
                rhoPatch[facei] = pressureIn/R_.value()/Tin; // Accturally not changed at all :p

                // inner boundary cell data state data
                label own = pOwner[facei];
                vector Ui = Uvol_[own];
                scalar Ti = Tvol_[own];
                scalar rhoi = rhoVol_[own];
                scalar ai = sqrt(R_.value() * Ti * (KInner_ + 5)/(KInner_ + 3)); // sos

                // change normal velocity component based on the characteristics
                vector norm = SfPatch[facei]/magSfPatch[facei]; // boundary face normal vector
                scalar Un = Ui & norm; // normal component
                scalar UnIn = Un + (pressureIn - rhoi * R_.value() * Ti)/rhoi/ai; // change normal component
                Upatch[facei] = UnIn * norm + (Ui - Un * norm); // tangential component not changed.
            }
        }
        else if(rhoBCs[patchi].type() == "pressureOut")
        {
            const fvsPatchField<vector>& SfPatch = mesh_.Sf().boundaryField()[patchi];
            const fvsPatchField<scalar>& magSfPatch = mesh_.magSf().boundaryField()[patchi];
            pressureOutFvPatchField<scalar>& rhoPatch = 
                refCast<pressureOutFvPatchField<scalar> >(rhoBCs[patchi]);
            fvPatchField<vector>& Upatch = Uvol_.BOUNDARY_FIELD_REF[patchi];
            fvPatchField<scalar>& Tpatch = Tvol_.BOUNDARY_FIELD_REF[patchi];
            const scalar pressureOut = rhoPatch.pressureOut();
            // now changed rho and U patch
            const labelUList& pOwner = mesh_.boundary()[patchi].faceCells();
            forAll(rhoPatch, facei)
            {
                // inner cell data state data
                label own = pOwner[facei];
                vector Ui = Uvol_[own];
                scalar Ti = Tvol_[own];
                scalar rhoi = rhoVol_[own];
                scalar ai = sqrt(R_.value() * Ti * (KInner_ + 5)/(KInner_ + 3)); // sos

                // change outlet density
                rhoPatch[facei] = rhoi  +  (pressureOut - rhoi * R_.value() * Ti)/ai/ai; // Accturally not changed at all :p
                Tpatch[facei] = pressureOut/(R_.value() * rhoi);

                // change normal velocity component based on the characteristics
                vector norm = SfPatch[facei]/magSfPatch[facei]; // boundary face normal vector
                scalar Un = Ui & norm; // normal component
                scalar UnIn = Un + ( rhoi * R_.value() * Ti - pressureOut)/rhoi/ai; // change normal component
                Upatch[facei] = UnIn * norm + (Ui - Un * norm); // tangential component not changed.
            }
        }
    }
}

template<template<class> class PatchType, class GeoMesh>
void Foam::fvDVM::updateTau
(
 GeometricField<scalar, PatchType, GeoMesh>& tau,
 const GeometricField<scalar, PatchType, GeoMesh>& T, 
 const GeometricField<scalar, PatchType, GeoMesh>& rho
 )
{
    tau = muRef_*exp(omega_*log(T/Tref_))/rho/T/R_;
}


void Foam::fvDVM::writeDFonCell(label cellId)
{
    std::ostringstream convert;
    convert << cellId;
    scalarIOList df
    (
        IOobject
        (
             "DF"+convert.str(),
             "0",
             mesh_,
             IOobject::NO_READ,
             IOobject::AUTO_WRITE
        )
    );
    //set size of df
    df.setSize(nXi_);

    scalarList dfPart(DV_.size());
    //put cellId's DF to dfPart
    forAll(dfPart, dfi)
        dfPart[dfi] = DV_[dfi].gTildeVol()[cellId];

    label nproc = mpiReducer_.nproc();
    //gather
    //tmp list for recv
    scalarList dfRcv(nXi_);

    //Compose displc and recvc
    labelField recvc(nproc);
    labelField displ(nproc);
    label chunck = nXi_/nproc;
    label left   = nXi_%nproc;
    forAll(recvc, i)
    {
        recvc[i] = chunck + (i<left) ;
        if(i<=left)
            displ[i] = i*(chunck + 1); // (i<=nXi_%nproc)
        else
            displ[i] = left*(chunck +1) + (i-left)*(chunck);
    }
    MPI_Gatherv(dfPart.data(), dfPart.size(), MPI_DOUBLE,
            dfRcv.data(), recvc.data(), displ.data(),
            MPI_DOUBLE, 0, MPI_COMM_WORLD);
    //reposition
    if(mpiReducer_.rank() == 0)
    {
        forAll(df, i)
        {
            label p   = i%nproc;
            label ldi = i/nproc;
            df[i] = dfRcv[displ[p]+ldi];
        }
        df.write();
    }
}

void Foam::fvDVM::writeDFonCells()
{
    if (time_.outputTime())
        forAll(DFwriteCellList_, i)
            writeDFonCell(i);
}

void Foam::fvDVM::adjustCFL(scalar currentResidual, scalar previousResidual)
{
    const scalar cflMin = 10.0;
    const scalar cflMax = 10000.0;
    const scalar safetyFactor = 0.9;  // 防止震荡的安全系数
    const scalar alpha = 1.5;         // 灵敏度参数，范围 [1, 2]
    const scalar maxIncrease = 2.0;   // 单步最大增倍
    const scalar maxDecrease = 0.5;   // 单步最大减半

    if (previousResidual > SMALL && currentResidual > SMALL)
    {
        // 公式: CFL^{n+1} = CFL^n * (R^n / R^{n+1})^alpha
        // residualRatio = R^n / R^{n+1}，>1 表示残差在下降
        scalar residualRatio = previousResidual / currentResidual;

        scalar cflAdjustment = pow(residualRatio, alpha);

        // 应用安全系数，抑制过激调整
        cflAdjustment = 1.0 + safetyFactor * (cflAdjustment - 1.0);

        // 限制单步调整幅度
        cflAdjustment = max(maxDecrease, min(maxIncrease, cflAdjustment));

        CFL_ *= cflAdjustment;

        // 限制 CFL 在合法范围内
        CFL_ = max(cflMin, min(cflMax, CFL_));

        if (mpiReducer_.rank() == 0)
        {
            if (cflAdjustment > 1.1)
            {
                //Info << "Residual decreasing fast (R^n/R^{n+1}=" << residualRatio
                //     << "), increasing CFL to " << CFL_ << endl;
                Info << "Residual decreasing fast" << ", increasing CFL to " << CFL_ << endl;
            }
            else if (cflAdjustment < 0.9)
            {
                //Info << "Residual decreasing slow (R^n/R^{n+1}=" << residualRatio
                //     << "), reducing CFL to " << CFL_ << endl;
                Info << "Residual decreasing slow" << ", reducing CFL to " << CFL_ << endl;
            }
            else
            {
                //Info << "Residual stable (R^n/R^{n+1}=" << residualRatio
                //     << "), maintaining CFL at " << CFL_ << endl;
		Info << "Residual stable" << ", maintaining CFL at " << CFL_ << endl;
            }
        }
    }
    else
    {
        // 残差无效时保守减小 CFL
        CFL_ = max(cflMin, CFL_ * maxDecrease);

        if (mpiReducer_.rank() == 0)
        {
            Info << "Invalid residual detected, reducing CFL to " << CFL_ << endl;
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvDVM::fvDVM
(
 volScalarField& rho,
 volVectorField& U,
 volScalarField& T,
 int* argc,
 char*** argv,
 Foam::argList& args
 )
    :
    IOdictionary
    (
        IOobject
        (
        "DVMProperties",
        T.time().constant(),
        T.mesh(),
        IOobject::MUST_READ,
        IOobject::NO_WRITE
        )
    ),
    mesh_(rho.mesh()),
    time_(rho.time()),
    rhoVol_(rho),
    Uvol_(U),
    Tvol_(T),
    args_(args),
    fvDVMparas_(subOrEmptyDict("fvDVMparas")),
    gasProperties_(subOrEmptyDict("gasProperties")),
    nXiPerDim_(readLabel(fvDVMparas_.lookup("nDV"))),
    xiMax_(fvDVMparas_.lookup("xiMax")),
    xiMin_(fvDVMparas_.lookup("xiMin")),
    dXi_((xiMax_-xiMin_)/(nXiPerDim_ - 1)),
    dXiCellSize_
    (
        "dXiCellSize",
        pow(dimLength/dimTime, 3),
        scalar(1.0)
    ),
    vMesh_(fvDVMparas_.lookupOrDefault("vMesh", word("no"))),
    CFL_(fvDVMparas_.lookupOrDefault("CFL", 100)),
    //macroFlux_(fvDVMparas_.lookupOrDefault("macroFlux", word("no"))),
    //res_(fvDVMparas_.lookupOrDefault("res", 1.0e-12)),
    //checkSteps_(fvDVMparas_.lookupOrDefault("checkSteps", 100)),
    R_(gasProperties_.lookup("R")),
    omega_(readScalar(gasProperties_.lookup("omega"))),
    Tref_(gasProperties_.lookup("Tref")),
    muRef_(gasProperties_.lookup("muRef")),
    Pr_(readScalar(gasProperties_.lookup("Pr"))),
    KInner_((gasProperties_.lookupOrDefault("KInner", 0))),
    mpiReducer_(args, argc, argv), // args comes from setRootCase.H in dugksFoam.C;
    DV_(0),
    rhoFlux_
    (
        IOobject
        (
        "rhoFlux",
        mesh_.time().timeName(),
        mesh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar( "0", dimensionSet(1,-3,0,0,0,0,0), 0)
    ),
    rhoUFlux_
    (
        IOobject
        (
        "rhoUFlux",
        mesh_.time().timeName(),
        mesh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector( "0", dimensionSet(1,-2,-1,0,0,0,0), vector(0,0,0))
    ),
    rhoEFlux_
    (
        IOobject
        (
        "rhoEFlux",
        mesh_.time().timeName(),
        mesh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar( "0", dimensionSet(1,-1,-2,0,0,0,0), 0)
    ),
    deltaRho_
    (
        IOobject
        (
        "deltaRho",
        mesh_.time().timeName(),
        mesh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar( "0", rho.dimensions(), 0)
    ),
    deltaRhoU_
    (
        IOobject
        (
        "deltaRhoU",
        mesh_.time().timeName(),
        mesh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector( "0", dimensionSet(1,-2,-1,0,0,0,0), vector(0,0,0))
    ),
    deltaRhoE_
    (
        IOobject
        (
        "deltaRhoE",
        mesh_.time().timeName(),
        mesh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar( "0", dimensionSet(1,-1,-2,0,0,0,0), 0)
    ),
    rhoSurf_
    (
        IOobject
        (
        "rhoSurf",
        mesh_.time().timeName(),
        mesh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar( "0", rho.dimensions(), 0)
    ),
    Tsurf_
    (
        IOobject
        (
        "Tsurf",
        mesh_.time().timeName(),
        mesh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar( "0", T.dimensions(), 0)
    ),
    Usurf_
    (
        IOobject
        (
        "Usurf",
        mesh_.time().timeName(),
        mesh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector( "0", U.dimensions(), vector(0,0,0))
    ),
    qSurf_
    (
        IOobject
        (
        "qSurf",
        mesh_.time().timeName(),
        mesh_,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector( "0", dimMass/pow(dimTime,3), vector(0,0,0))
    ),
    stressSurf_
    (
        IOobject
        (
        "stressSurf",
        mesh_.time().timeName(),
        mesh_,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedTensor("0", dimensionSet(1,-1,-2,0,0,0,0), pTraits<tensor>::zero)
    ),
    qVol_
    (
        IOobject
        (
        "q",
        mesh_.time().timeName(),
        mesh_,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector( "0", dimMass/pow(dimTime,3), vector(0,0,0))
    ),
    stress_
    (
        IOobject
        (
        "stress",
        mesh_.time().timeName(),
        mesh_,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
        ),
        mesh_,
	    dimensionedTensor("0", dimensionSet(1,-1,-2,0,0,0,0), pTraits<tensor>::zero)
    ),
    tauVol_
    (
        IOobject
        (
        "tauVol",
        mesh_.time().timeName(),
        mesh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar( "0", dimTime, 0)
    ),
    tauSurf_
    (
        IOobject
        (
        "tauSurf",
        mesh_.time().timeName(),
        mesh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar( "0", dimTime, 0)
    ),
    qWall_
    (
        IOobject
        (
        "qWall",
        mesh_.time().timeName(),
        mesh_,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector( "0", dimMass/pow(dimTime,3), vector(0,0,0))
    ),
    stressWall_
    (
        IOobject
        (
        "stressWall",
        mesh_.time().timeName(),
        mesh_,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedTensor("0", dimensionSet(1,-1,-2,0,0,0,0), pTraits<tensor>::zero)
    ),
    qfluxWall_
    (
        IOobject
        (
        "qfluxWall",
        mesh_.time().timeName(),
        mesh_,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar( "0", dimMass/pow(dimTime,3), 0)
    ),
    presurefluxWall_
    (
        IOobject
        (
        "presureWall",
        mesh_.time().timeName(),
        mesh_,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar( "0", dimensionSet(1,-1,-2,0,0,0,0), 0)
    ),
    stressfluxWall_
    (
        IOobject
        (
        "stressfluxWall",
        mesh_.time().timeName(),
        mesh_,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector( "0", dimensionSet(1,-1,-2,0,0,0,0), vector(0,0,0))
    )

{
    DFwriteCellList_ = lookupOrDefault<labelList>("DFwriteCellList", labelList()),
    initialiseDV();
    setCalculatedMaxwellRhoBC();
    setSymmetryModRhoBC();
    // set initial rho in pressureIn/Out BC
    updatePressureInOutBC();
    updateTau(tauVol_, Tvol_, rhoVol_); //calculate the tau at cell when init
    Usurf_ = fvc::interpolate(Uvol_, "linear"); // for first time Dt calculation.
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvDVM::~fvDVM()
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fvDVM::evolution()
{
    // Info << "Begin evolution" << endl;
    updateGHbarPvol();
    // Info << "Done updateGHbarPvol " << endl;
    updateGHbarSurf();
    // Info << "Done updateGHbarSurf " << endl;
    // updateMaxwellWallRho();
    // //Info << "Done updateMaxwellWallRho " << endl;
    // updateGHbarSurfMaxwellWallIn();
    // //Info << "Done updateGHbarSurfMaxwellWallIn " << endl;
     updateGHbarSurfSymmetryIn();
    // //Info << "Done updateGHbarSurfSymmetryIn " << endl;
    updateMacroSurf();
    //Info << "Done updateMacroSurf " << endl;
    updateGHsurf();
    //Info << "Done updateGHsurf " << endl;
     updateMaxwellWallRho();
     //Info << "Done updateMaxwellWallRho " << endl;
     updateGHbarSurfMaxwellWallIn();
    // //Info << "Done updateGHbarSurfMaxwellWallIn " << endl;
	updateWall();
    updateGHtildeVol();
    //Info << "Done updateGHtildeVol " << endl;
    updateMacroFlux();
    //Info << "Done updateMacroFlux " << endl;
    for(label iter = 0; iter < 60; ++iter)
    {
        PRSGSforwardsweep();
        PRSGSbackwardsweep();
    }
    //Info << "Done updateSweep " << endl;
    predicMacroVol();

    updateequilibriumShakhov();

    for(label iter = 0; iter < 2; ++iter)
    {
        forwardsweep();
        backwardsweep();
    }

    updateGHVol();
    //Info << "Done updateMacroVol " << endl;
    updateMacroVol();
    //
    updatePressureInOutBC();
}


void Foam::fvDVM::getCoNum(scalar& maxCoNum, scalar& meanCoNum)
{
    scalar dt = time_.deltaTValue();
    scalarField UbyDx =
        mesh_.surfaceInterpolation::deltaCoeffs()
        *(mag(Usurf_) + sqrt(scalar(mesh_.nSolutionD()))*xiMax_);
    maxCoNum = gMax(UbyDx)*dt;
    meanCoNum = gSum(UbyDx)/UbyDx.size()*dt;
}



void Foam::fvDVM::EulerFlux
(
    scalar& massFlux,
    vector& momentumFlux,
    scalar& energyFlux,
    const scalar& rho,
    const vector& rhoU,
    const scalar& rhoE,
    const vector& n
)
{
    vector U = rhoU / rho;
    scalar p = (2.0 * rhoE - rho * magSqr(U)) / (KInner_ + 3.0);
    massFlux = rho * (U & n);
    momentumFlux = rhoU * (U & n) + p * n;
    energyFlux = (rhoE + p) * (U & n);
}



const fieldMPIreducer& Foam::fvDVM::mpiReducer() const
{
    return mpiReducer_;
}


// ************************************************************************* //
