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
    Foam::staticGravityPLIC

Description
    constant gravitional forces

SourceFiles
    staticGravityPLIC.C

\*---------------------------------------------------------------------------*/

#ifndef staticGravityPLIC_H
#define staticGravityPLIC_H

#include "typeInfo.H"
#include "volFields.H"
#include "dimensionedScalar.H"
#include "autoPtr.H"
#include "uniformDimensionedFields.H"
#include "geometricVoFMethod.H"

#include "accelerationForceGeoVoFMethod.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class staticGravityPLIC Declaration
\*---------------------------------------------------------------------------*/

class staticGravityPLIC
    : public accelerationForceGeoVoFMethod
{

        //- const access dict
        const dictionary& staticGravityPLICDict_;

        //- staticGravityPLIC vector
        dimensionedVector g_;

        //- Href Field
        uniformDimensionedScalarField hRef_;

    // Private Member Functions

        //- Disallow copy construct
        staticGravityPLIC(const staticGravityPLIC&);

        //- Disallow default bitwise assignment
        void operator=(const staticGravityPLIC&);


        protected:

        //- Re-calculate the interface curvature
        virtual void calculateAcc();


public:

    //- Runtime type information
    TypeName("staticGravityPLIC");


    // Constructors

        staticGravityPLIC
        (
            geometricVoFMethod& geoVoF,
            const dictionary& dict
        );

        staticGravityPLIC
        (
            interfaceCapturingMethod& geoVoF,
            const dictionary& dict
        );


    //- Destructor
    virtual ~staticGravityPLIC()
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
