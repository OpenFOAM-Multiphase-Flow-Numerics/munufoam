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
    Foam::standard

Description
   

SourceFiles
    standard.C

\*---------------------------------------------------------------------------*/

#ifndef standard_H
#define standard_H

#include "typeInfo.H"
#include "volFields.H"
#include "algebraicVoF.H"
#include "interfaceCapturingMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class standard Declaration
\*---------------------------------------------------------------------------*/

class standard
{

public:

    //- Runtime type information
    TypeName("standard");

    //- Constructors

        //- Construct given alpha, phi and velocity field. Note: phi should be
        // divergence free up to a sufficient tolerance
        standard
        (
            const fvMesh& mesh,
            const dictionary& dict
        );

    //- Destructor
    virtual ~standard() = default;


    void update(algebraicVoF& newSurf,algebraicVoF& oldSurf,Foam::timeState state);


};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
