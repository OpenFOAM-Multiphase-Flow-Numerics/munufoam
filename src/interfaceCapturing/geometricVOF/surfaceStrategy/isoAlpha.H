/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 DHI
    Modified code Copyright (C) 2016-2019 OpenCFD Ltd.
    Modified code Copyright (C) 2018 Johan Roenby
    Modified code Copyright (C) 2019 DLR
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

Class
    Foam::isoAlpha

Description
   

SourceFiles
    isoAlpha.C

\*---------------------------------------------------------------------------*/

#ifndef isoAlpha_H
#define isoAlpha_H

#include "typeInfo.H"
#include "volFields.H"
#include "geometricVoF.H"
#include "interfaceCapturingMethod.H"
#include "volFields.H"

#include "surfaceIteratorIso.H"
#include "volPointInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class isoAlpha Declaration
\*---------------------------------------------------------------------------*/

class isoAlpha
{
    // Private Data

        //- Reference to mesh
        const fvMesh& mesh_;

        //- VOF field interpolated to mesh points
        scalarField ap_;

        // Switches and tolerances. Tolerances need to go into toleranceSwitches

        //- Tolerance for search of isoFace giving specified VOF value
        scalar isoFaceTol_;

        //- Tolerance for marking of surface cells:
        //  Those with surfCellTol_ < alpha1 < 1 - surfCellTol_
        scalar surfCellTol_;

        //- surfaceIterator finds the isovalue for specified VOF value
        surfaceIteratorIso sIterIso_;


    // Private Member Functions

        //- No copy construct
        isoAlpha(const isoAlpha&) = delete;

        //- No copy assignment
        void operator=(const isoAlpha&) = delete;

public:

    //- Runtime type information
    TypeName("isoAlpha");

    //- Constructors

        //- Construct given alpha, phi and velocity field. Note: phi should be
        // divergence free up to a sufficient tolerance
        isoAlpha
        (
            const fvMesh& mesh,
            const dictionary& dict
        );

    //- Destructor
    virtual ~isoAlpha() = default;

    void reconstruct(geometricVoF& surf);

    void update(geometricVoF& newSurf,geometricVoF& oldSurf,Foam::timeState state);


};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
