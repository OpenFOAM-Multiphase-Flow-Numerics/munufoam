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
    Foam::implicitAdvection

Description
   

SourceFiles
    implicitAdvection.C

\*---------------------------------------------------------------------------*/

#ifndef implicitAdvection_H
#define implicitAdvection_H

#include "typeInfo.H"
#include "volFields.H"
#include "algebraicVoF.H"
#include "interfaceCapturingMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class implicitAdvection Declaration
\*---------------------------------------------------------------------------*/

class implicitAdvection
{
private:

    //- Reference to the VoF Field
    volScalarField& alpha1_;

    //- Reference to the face fluxes
    surfaceScalarField& alphaPhi_;

    //- Reference to the face fluxes
    const surfaceScalarField& phi_;

    //- Reference to the velocity field
    const volVectorField& U_;

    label nAlphaCorr_;

public:

    //- Runtime type information
    TypeName("implicitAdvection");

    //- Constructors

        //- Construct given alpha, phi and velocity field. Note: phi should be
        // divergence free up to a sufficient tolerance
        implicitAdvection
        (
            volScalarField& alpha1,
            surfaceScalarField& alphaPhi,
            const surfaceScalarField& phi,
            const volVectorField& U
        );


    //- Destructor
    virtual ~implicitAdvection() = default;

    void advect
    (
        algebraicVoF& algVoFNew,
        algebraicVoF& algVoFOld,
        const volScalarField::Internal& Sp,
        const volScalarField::Internal& Su
    );

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
