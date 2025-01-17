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
    Foam::geometricVoF

Description
    estimates the curvature from the volume fraction field similiar to the
    mode of geometricVoF et al.

SourceFiles
    geometricVoF.C

\*---------------------------------------------------------------------------*/

#ifndef geometricVoF_H
#define geometricVoF_H

#include "typeInfo.H"
#include "volFields.H"
#include "dimensionedScalar.H"
#include "autoPtr.H"
#include "surfaceBase.H"
#include "interfaceCapturingMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class geometricVoF Declaration
\*---------------------------------------------------------------------------*/

class geometricVoF
:
    public surfaceBase
{
    // private data members

        const volScalarField& alpha_;

public:

    //- Runtime type information
    TypeName("geometricVoF");


    // Constructors

        //- Construct from dictionary
        geometricVoF
        (
            volScalarField& alpha1,
            const surfaceScalarField& phi,
            const volVectorField& U,
            timeState state
        );


    //- Destructor
    virtual ~geometricVoF() = default;

    timeState state;

    
    const volScalarField& alpha() const
    {
        if (state == timeState::oldState)
        {
            return alpha_.oldTime();
        }
        return alpha_;
    }

    //- Interface area normals
    volVectorField normal;

    //- Interface centres
    volVectorField centre;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
