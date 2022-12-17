/*---------------------------------------------------------------------------*\
            Copyright (c) 2017-2019, German Aerospace Center (DLR)
-------------------------------------------------------------------------------
License
    This file is part of the VoFLibrary source code library, which is an
	unofficial extension to OpenFOAM.
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
    Foam::surfaceBase

Description
    estimates the curvature from the volume fraction field similiar to the
    mode of surfaceBase et al.

SourceFiles
    surfaceBase.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceBase_H
#define surfaceBase_H

#include "typeInfo.H"
#include "volFields.H"
#include "dimensionedScalar.H"
#include "autoPtr.H"
#include "DynamicField.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class surfaceBase Declaration
\*---------------------------------------------------------------------------*/

class surfaceBase
{

public:

    //- Runtime type information
    TypeName("surfaceBase");


    // Constructors

        //- Construct from dictionary
        surfaceBase
        (
            volScalarField& alpha1,
            const surfaceScalarField& phi,
            const volVectorField& U
        );


        //- Destructor
        virtual ~surfaceBase() = default;


    // public member

        //- Reference to the VoF Field
        volScalarField& alpha1_;

        //- Reference to the face fluxes
        const surfaceScalarField& phi_;

        //- Reference to the velocity field
        const volVectorField& U_;

        //- Interface area normals
        volVectorField normal_;

        //- Interface centres
        volVectorField centre_;

        //- Is interface cell?
        boolList interfaceCell_;

        //- List of cell labels that have an interface
        DynamicField<label> interfaceLabels_;

        //- Store timeindex/iteration to avoid multiple reconstruction
        mutable Pair<label> timeIndexAndIter_;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
