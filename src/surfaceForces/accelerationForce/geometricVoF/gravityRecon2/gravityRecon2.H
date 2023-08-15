/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
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
    Foam::gravityRecon2

Description

SourceFiles
    gravityRecon2.C
    newPhaseChangeModel.C

\*---------------------------------------------------------------------------*/

#ifndef gravityRecon2_H
#define gravityRecon2_H

#include "typeInfo.H"
#include "volFields.H"
#include "dimensionedScalar.H"
#include "autoPtr.H"
#include "uniformDimensionedFields.H"

#include "zoneDistribute.H"
#include "accelerationForceGeoVoFMethod.H"
#include "markInterfaceRegion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class gravityRecon2 Declaration
\*---------------------------------------------------------------------------*/

class gravityRecon2
:
    public accelerationForceGeoVoFMethod
{

private:


        const dictionary& gravityDict_;

        //- gravity vector
        dimensionedVector g_;

        //- Href Field
        uniformDimensionedScalarField hRef_;

    // Private Member Functions

        //- Disallow copy construct
        gravityRecon2(const gravityRecon2&);

        //- Disallow default bitwise assignment
        void operator=(const gravityRecon2&);

        vector closestDist(const point p, const vector n ,const vector cenrte);


        protected:

        //- Re-calculate the acceleration
        virtual void calculateAcc();

public:

    //- Runtime type information
    TypeName("gravityRecon2");


    // Constructors

        //- Construct from dictionary
        gravityRecon2
        (
            geometricVoFMethod& geoVoF,
            const dictionary& dict
        );

        gravityRecon2
        (
            interfaceCapturingMethod& geoVoF,
            const dictionary& dict
        );

    //- Destructor
    virtual ~gravityRecon2()
    {}

        //- return the pressure field depending
    virtual volScalarField& pressure
    (
        volScalarField& p,
        volScalarField& p_rgh
    );

    //- update static pressure field
    virtual void updatePressure
    (
        volScalarField& p,
        volScalarField& p_rgh,
        const volScalarField& rho
    );

    //- update dynamic pressure field
    virtual void updateRefPressure
    (
        volScalarField& p,
        volScalarField& p_rgh,
        const volScalarField& rho
    );

    virtual tmp<surfaceScalarField> accelerationForce();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
