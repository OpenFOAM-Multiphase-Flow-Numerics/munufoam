/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
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

Class
    Foam::interface::ParkerYoung

Description
    Reconstructs an interface (centre and normal vector) consisting of planes
    to match the internal fluid distribution in cells. The interface normals
    are estimated by least square gradient scheme on the VoF Field (alpha).
    Also known as Parker Young method

    Original code supplied by Henning Scheufler, DLR (2019)

SourceFiles
    ParkerYoung.C

\*---------------------------------------------------------------------------*/

#ifndef ParkerYoung_H
#define ParkerYoung_H

#include "typeInfo.H"
#include "volFields.H"
#include "dimensionedScalar.H"
#include "autoPtr.H"
#include "surfaceIteratorPLIC.H"
#include "zoneDistribute.H"
#include "geometricVoF.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class ParkerYoung Declaration
\*---------------------------------------------------------------------------*/

class ParkerYoung
{
    // Private Data

        //- Reference to mesh
        const fvMesh& mesh_;

        //- Interpolation object from cell centres to points
        DynamicField<vector> interfaceNormal_;


    // Switches and tolerances. Tolerances need to go into toleranceSwitches

        //- Tolerance for search of isoFace giving specified VOF value
        scalar isoFaceTol_;

        //- Tolerance for marking of surface cells:
        //  Those with surfCellTol_ < alpha1 < 1 - surfCellTol_
        scalar surfCellTol_;

        //- SurfaceIterator finds the plane centre for specified VOF value
        surfaceIteratorPLIC sIterPLIC_;

        //- Compute gradient at the surfaces
        void gradSurf(const volScalarField& phi,const geometricVoF& surf);


    // Private Member Functions

        //- No copy construct
        ParkerYoung(const ParkerYoung&) = delete;

        //- No copy assignment
        void operator=(const ParkerYoung&) = delete;


public:

    //- Runtime type information
    TypeName("ParkerYoung");

    //- Construct from mesh and dict
    ParkerYoung
    (
        const fvMesh& mesh,
        const dictionary& dict
    );

    //- Destructor
    virtual ~ParkerYoung() = default;


    // Member Functions

    void reconstruct(geometricVoF& surf);

    void update(geometricVoF& newSurf,geometricVoF& oldSurf,Foam::timeState state);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


#endif

// ************************************************************************* //
