/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2020 DLR
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

#include "ParkerYoung.H"
#include "fvc.H"
#include "leastSquareGrad.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ParkerYoung, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::ParkerYoung::gradSurf(const volScalarField& phi, const geometricVoF& surf)
{
    addProfilingInFunction(geometricVoF);
    leastSquareGrad<scalar> lsGrad("polyDegree1",mesh_.geometricD());

    zoneDistribute& exchangeFields = zoneDistribute::New(mesh_);

    exchangeFields.setUpCommforZone(surf.interfaceCell,true);

    Map<vector> mapCC
    (
        exchangeFields.getDatafromOtherProc(surf.interfaceCell, mesh_.C())
    );
    Map<scalar> mapPhi
    (
        exchangeFields.getDatafromOtherProc(surf.interfaceCell, phi)
    );

    DynamicField<vector> cellCentre(100);
    DynamicField<scalar> phiValues(100);

    const labelListList& stencil = exchangeFields.getStencil();

    forAll(surf.interfaceLabels, i)
    {
        const label celli = surf.interfaceLabels[i];

        cellCentre.clear();
        phiValues.clear();

        for (const label gblIdx : stencil[celli])
        {
            cellCentre.append
            (
                exchangeFields.getValue(mesh_.C(), mapCC, gblIdx)
            );
            phiValues.append
            (
                exchangeFields.getValue(phi, mapPhi, gblIdx)
            );
        }

        cellCentre -= mesh_.C()[celli];
        interfaceNormal_[i] = lsGrad.grad(cellCentre, phiValues);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ParkerYoung::ParkerYoung
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    interfaceNormal_(mesh.nCells()),
    isoFaceTol_(dict.lookupOrDefault<scalar>("isoFaceTol", 1e-8)),
    surfCellTol_(dict.lookupOrDefault<scalar>("surfCellTol", 1e-8)),
    sIterPLIC_(mesh_,surfCellTol_)
{
    
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::ParkerYoung::update(geometricVoF& newSurf,geometricVoF& oldSurf,Foam::timeState state)
{
    addProfilingInFunction(geometricVoF);
    if (state == timeState::oldState)
    {
        reconstruct(oldSurf);
    }
    else
    {
        reconstruct(newSurf);
    }

}

void Foam::ParkerYoung::reconstruct(geometricVoF& surf)
{
    addProfilingInFunction(geometricVoF);
    // const bool uptodate = alreadyReconstructed(forceUpdate);

    // if (uptodate && !forceUpdate)
    // {
    //     return;
    // }
    const volScalarField& alpha = surf.alpha();

    if (mesh_.topoChanging())
    {
        // Introduced resizing to cope with changing meshes
        if (surf.interfaceCell.size() != mesh_.nCells())
        {
            surf.interfaceCell.resize(mesh_.nCells());
        }
    }
    surf.interfaceCell = false;

    surf.interfaceLabels.clear();

    forAll(alpha, celli)
    {
        if (sIterPLIC_.isASurfaceCell(alpha[celli]))
        {
            surf.interfaceCell[celli] = true; // is set to false earlier
            surf.interfaceLabels.append(celli);
        }
    }
    interfaceNormal_.resize(surf.interfaceLabels.size());
    surf.centre = dimensionedVector("centre", dimLength, Zero);
    surf.normal  = dimensionedVector("normal", dimArea, Zero);

    gradSurf(alpha,surf);

    forAll(surf.interfaceLabels, i)
    {
        const label celli = surf.interfaceLabels[i];
        if (mag(interfaceNormal_[i]) == 0)
        {
            continue;
        }

        sIterPLIC_.vofCutCell
        (
            celli,
            alpha[celli],
            isoFaceTol_,
            100,
            interfaceNormal_[i]
        );

        if (sIterPLIC_.cellStatus() == 0)
        {
            surf.normal[celli] = sIterPLIC_.surfaceArea();
            surf.centre[celli] = sIterPLIC_.surfaceCentre();
            if (mag(surf.normal[celli]) == 0)
            {
                surf.normal[celli] = Zero;
                surf.centre[celli] = Zero;
            }
        }
        else
        {
            surf.normal[celli] = Zero;
            surf.centre[celli] = Zero;
        }
    }
}




// ************************************************************************* //
