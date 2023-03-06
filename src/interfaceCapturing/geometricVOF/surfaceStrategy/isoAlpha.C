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

#include "isoAlpha.H"
#include "profiling.H"

namespace Foam
{
    defineTypeNameAndDebug(isoAlpha, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::isoAlpha::isoAlpha
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    // Interpolation data
    ap_(mesh_.nPoints()),

    // Tolerances and solution controls
    isoFaceTol_(dict.lookupOrDefault<scalar>("isoFaceTol", 1e-8)),
    surfCellTol_(dict.lookupOrDefault<scalar>("surfCellTol", 1e-8)),
    sIterIso_(mesh_,ap_,surfCellTol_)
{

}

void Foam::isoAlpha::reconstruct(geometricVoF& surf)
{

    // const bool uptodate = alreadyReconstructed(forceUpdate);

    // if (uptodate && !forceUpdate)
    // {
    //     return;
    // }
    const volScalarField& alpha = surf.alpha();

    // Interpolating alpha1 cell centre values to mesh points (vertices)
    if (mesh_.topoChanging())
    {
        // Introduced resizing to cope with changing meshes
        if (ap_.size() != mesh_.nPoints())
        {
            ap_.resize(mesh_.nPoints());

        }
        if (surf.interfaceCell.size() != mesh_.nCells())
        {
            surf.interfaceCell.resize(mesh_.nCells());
        }
    }
    ap_ = volPointInterpolation::New(mesh_).interpolate(alpha);

    surf.interfaceLabels.clear();

    forAll(alpha,cellI)
    {
        if (sIterIso_.isASurfaceCell(alpha[cellI]))
        {
            surf.interfaceLabels.append(cellI);

            sIterIso_.vofCutCell
            (
                cellI,
                alpha[cellI],
                isoFaceTol_,
                100
            );

            if (sIterIso_.cellStatus() == 0)
            {
                surf.normal[cellI] = sIterIso_.surfaceArea();
                surf.centre[cellI] = sIterIso_.surfaceCentre();
                if (mag(surf.normal[cellI]) != 0)
                {
                    surf.interfaceCell[cellI] = true;
                }
                else
                {
                    surf.interfaceCell[cellI] = false;
                }
            }
            else
            {
                surf.normal[cellI] = Zero;
                surf.centre[cellI] = Zero;
                surf.interfaceCell[cellI] = false;
            }
         }
         else
         {
            surf.normal[cellI] = Zero;
            surf.centre[cellI] = Zero;
            surf.interfaceCell[cellI] = false;
         }
    }
}

void Foam::isoAlpha::update(geometricVoF& newSurf,geometricVoF& oldSurf,Foam::timeState state)
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
// ************************************************************************* //
