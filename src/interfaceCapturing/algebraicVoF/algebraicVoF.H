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
    Foam::algebraicVoF

Description
    estimates the curvature from the volume fraction field similiar to the
    mode of algebraicVoF et al.

SourceFiles
    algebraicVoF.C

\*---------------------------------------------------------------------------*/

#ifndef algebraicVoF_H
#define algebraicVoF_H

#include "typeInfo.H"
#include "volFields.H"
#include "dimensionedScalar.H"
#include "autoPtr.H"
#include "surfaceBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class algebraicVoF Declaration
\*---------------------------------------------------------------------------*/

class algebraicVoF
:
    public surfaceBase
{

private:


        // //- Disallow copy construct
        // algebraicVoF(const algebraicVoF&);

        // //- Disallow default bitwise assignment
        // void operator=(const algebraicVoF&);

public:

    //- Runtime type information
    TypeName("algebraicVoF");


    // Constructors

        //- Construct from dictionary
        algebraicVoF
        (
            volScalarField& alpha1,
            const surfaceScalarField& phi,
            const volVectorField& U
        );


    //- Destructor
    virtual ~algebraicVoF() = default;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
