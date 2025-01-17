/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Modified code Copyright (C) 2022 henning Scheufler
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

#include "isoAdvection.H"
#include "volFields.H"
#include "volPointInterpolation.H"
#include "fvcSurfaceIntegrate.H"
#include "fvcGrad.H"
#include "upwind.H"
#include "cellSet.H"
#include "meshTools.H"
#include "syncTools.H"
#include "profiling.H"


namespace Foam
{
    defineTypeNameAndDebug(isoAdvection, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::isoAdvection::isoAdvection
(
    geometricVoFMethod& geoVoF
)
:
    geoVoF_(geoVoF),
    alpha1_(geoVoF.alpha1()),
    alphaPhi1_(geoVoF.alphaPhi()),
    phi_(geoVoF.phi()),
        // General data
    mesh_(geoVoF.alpha1().mesh()),
    alpha1In_(geoVoF.alpha1().ref()),
    timeIndex_(-1),

    // Tolerances and solution controls
    nAlphaBounds_(geoVoF.modelDict().lookupOrDefault<label>("nAlphaBounds", 3)),
    isoFaceTol_(geoVoF.modelDict().lookupOrDefault<scalar>("isoFaceTol", 1e-8)),
    surfCellTol_(geoVoF.modelDict().lookupOrDefault<scalar>("surfCellTol", 1e-8)),
    writeIsoFacesToFile_(geoVoF.modelDict().lookupOrDefault("writeIsoFaces", false)),

    // Cell cutting data
    surfCells_(label(0.2*mesh_.nCells())),
    advectFace_(geoVoF.alpha1().mesh(), geoVoF.alpha1()),
    bsFaces_(label(0.2*mesh_.nBoundaryFaces())),
    bsx0_(bsFaces_.size()),
    bsn0_(bsFaces_.size()),
    bsUn0_(bsFaces_.size()),

    // Parallel run data
    procPatchLabels_(mesh_.boundary().size()),
    surfaceCellFacesOnProcPatches_(0)
{
   cutFaceAdvect::debug = debug;

    // Prepare lists used in parallel runs
    if (Pstream::parRun())
    {
        // Force calculation of required demand driven data (else parallel
        // communication may crash)
        mesh_.cellCentres();
        mesh_.cellVolumes();
        mesh_.faceCentres();
        mesh_.faceAreas();
        mesh_.magSf();
        mesh_.boundaryMesh().patchID();
        mesh_.cellPoints();
        mesh_.cellCells();
        mesh_.cells();

        // Get boundary mesh and resize the list for parallel comms
        setProcessorPatches();
    }
}


void Foam::isoAdvection::setProcessorPatches()
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();
    surfaceCellFacesOnProcPatches_.clear();
    surfaceCellFacesOnProcPatches_.resize(patches.size());

    // Append all processor patch labels to the list
    procPatchLabels_.clear();
    forAll(patches, patchi)
    {
        if
        (
            isA<processorPolyPatch>(patches[patchi])
            && patches[patchi].size() > 0
        )
        {
            procPatchLabels_.append(patchi);
        }
    }
}


void Foam::isoAdvection::extendMarkedCells
(
    bitSet& markedCell
) const
{
    // Mark faces using any marked cell
    bitSet markedFace(mesh_.nFaces());

    for (const label celli : markedCell)
    {
        markedFace.set(mesh_.cells()[celli]);  // set multiple faces
    }

    syncTools::syncFaceList(mesh_, markedFace, orEqOp<unsigned int>());

    // Update cells using any markedFace
    for (label facei = 0; facei < mesh_.nInternalFaces(); ++facei)
    {
        if (markedFace.test(facei))
        {
            markedCell.set(mesh_.faceOwner()[facei]);
            markedCell.set(mesh_.faceNeighbour()[facei]);
        }
    }
    for (label facei = mesh_.nInternalFaces(); facei < mesh_.nFaces(); ++facei)
    {
        if (markedFace.test(facei))
        {
            markedCell.set(mesh_.faceOwner()[facei]);
        }
    }
}


void Foam::isoAdvection::setDownwindFaces
(
    const label celli,
    DynamicLabelList& downwindFaces
) const
{
    DebugInFunction << endl;

    // Get necessary mesh data and cell information
    const labelList& own = mesh_.faceOwner();
    const cellList& cells = mesh_.cells();
    const cell& c = cells[celli];

    downwindFaces.clear();

    // Check all faces of the cell
    forAll(c, fi)
    {
        // Get face and corresponding flux
        const label facei = c[fi];
        const scalar phi = faceValue(phi_, facei);

        if (own[facei] == celli)
        {
            if (phi >= 0)
            {
                downwindFaces.append(facei);
            }
        }
        else if (phi < 0)
        {
            downwindFaces.append(facei);
        }
    }

    downwindFaces.shrink();
}


Foam::scalar Foam::isoAdvection::netFlux
(
    const surfaceScalarField& dVf,
    const label celli
) const
{
    scalar dV = 0;

    // Get face indices
    const cell& c = mesh_.cells()[celli];

    // Get mesh data
    const labelList& own = mesh_.faceOwner();

    forAll(c, fi)
    {
        const label facei = c[fi];
        const scalar dVff = faceValue(dVf, facei);

        if (own[facei] == celli)
        {
            dV += dVff;
        }
        else
        {
            dV -= dVff;
        }
    }

    return dV;
}


Foam::DynamicList<Foam::label>  Foam::isoAdvection::syncProcPatches
(
    surfaceScalarField& dVf,
    const surfaceScalarField& phi,
    bool returnSyncedFaces
)
{
    DynamicLabelList syncedFaces(0);
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    if (Pstream::parRun())
    {
        PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

        // Send
        forAll(procPatchLabels_, i)
        {
            const label patchi = procPatchLabels_[i];

            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>(patches[patchi]);

            UOPstream toNbr(procPatch.neighbProcNo(), pBufs);
            const scalarField& pFlux = dVf.boundaryField()[patchi];

            const List<label>& surfCellFacesOnProcPatch =
                surfaceCellFacesOnProcPatches_[patchi];

            const UIndirectList<scalar> dVfPatch
            (
                pFlux,
                surfCellFacesOnProcPatch
            );

            toNbr << surfCellFacesOnProcPatch << dVfPatch;
        }

        pBufs.finishedSends();


        // Receive and combine
        forAll(procPatchLabels_, patchLabeli)
        {
            const label patchi = procPatchLabels_[patchLabeli];

            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>(patches[patchi]);

            UIPstream fromNeighb(procPatch.neighbProcNo(), pBufs);
            List<label> faceIDs;
            List<scalar> nbrdVfs;

            fromNeighb >> faceIDs >> nbrdVfs;
            if (returnSyncedFaces)
            {
                List<label> syncedFaceI(faceIDs);
                for (label& faceI : syncedFaceI)
                {
                    faceI += procPatch.start();
                }
                syncedFaces.append(syncedFaceI);
            }

            if (debug)
            {
                Pout<< "Received at time = " << mesh_.time().value()
                    << ": surfCellFacesOnProcPatch = " << faceIDs << nl
                    << "Received at time = " << mesh_.time().value()
                    << ": dVfPatch = " << nbrdVfs << endl;
            }

            // Combine fluxes
            scalarField& localFlux = dVf.boundaryFieldRef()[patchi];

            forAll(faceIDs, i)
            {
                const label facei = faceIDs[i];
                localFlux[facei] = - nbrdVfs[i];
                if (debug && mag(localFlux[facei] + nbrdVfs[i]) > ROOTVSMALL)
                {
                    Pout<< "localFlux[facei] = " << localFlux[facei]
                        << " and nbrdVfs[i] = " << nbrdVfs[i]
                        << " for facei = " << facei << endl;
                }
            }
        }

        if (debug)
        {
            // Write out results for checking
            forAll(procPatchLabels_, patchLabeli)
            {
                const label patchi = procPatchLabels_[patchLabeli];
                const scalarField& localFlux = dVf.boundaryField()[patchi];
                Pout<< "time = " << mesh_.time().value() << ": localFlux = "
                    << localFlux << endl;
            }
        }

        // Reinitialising list used for minimal parallel communication
        forAll(surfaceCellFacesOnProcPatches_, patchi)
        {
            surfaceCellFacesOnProcPatches_[patchi].clear();
        }
    }

    return syncedFaces;
}


void Foam::isoAdvection::checkIfOnProcPatch(const label facei)
{
    if (!mesh_.isInternalFace(facei))
    {
        const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
        const label patchi = pbm.patchID()[facei - mesh_.nInternalFaces()];

        if (isA<processorPolyPatch>(pbm[patchi]) && pbm[patchi].size())
        {
            const label patchFacei = pbm[patchi].whichFace(facei);
            surfaceCellFacesOnProcPatches_[patchi].append(patchFacei);
        }
    }
}


void Foam::isoAdvection::applyBruteForceBounding()
{
    addProfilingInFunction(geometricVoF);
    bool alpha1Changed = false;

    scalar snapAlphaTol = geoVoF_.modelDict().lookupOrDefault("snapTol", 0);
    if (snapAlphaTol > 0)
    {
        alpha1_ =
            alpha1_
           *pos0(alpha1_ - snapAlphaTol)
           *neg0(alpha1_ - (1.0 - snapAlphaTol))
          + pos0(alpha1_ - (1.0 - snapAlphaTol));

        alpha1Changed = true;
    }

    if (geoVoF_.modelDict().lookupOrDefault("clip", true))
    {
        alpha1_ = min(scalar(1), max(scalar(0), alpha1_));
        alpha1Changed = true;
    }

    if (alpha1Changed)
    {
        alpha1_.correctBoundaryConditions();
    }
}

void Foam::isoAdvection::advect
(
    const volScalarField::Internal& Sp,
    const volScalarField::Internal& Su
)
{
    if (fv::localEulerDdt::enabled(mesh_))
    {
        const volScalarField& rDeltaT = fv::localEulerDdt::localRDeltaT(mesh_);
        advect(rDeltaT,Sp,Su);
    }
    else
    {
        scalar rDeltaT = 1/mesh_.time().deltaTValue();
        advect(rDeltaT,Sp,Su);
    }
}
// ************************************************************************* //
