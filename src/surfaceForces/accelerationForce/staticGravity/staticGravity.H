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
    Foam::staticGravity

Description
    constant gravitional forces

SourceFiles
    staticGravity.C

\*---------------------------------------------------------------------------*/

#ifndef staticGravity_H
#define staticGravity_H

#include "typeInfo.H"
#include "volFields.H"
#include "dimensionedScalar.H"
#include "autoPtr.H"
#include "uniformDimensionedFields.H"

#include "accelerationForceMethod.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class staticGravity Declaration
\*---------------------------------------------------------------------------*/

class staticGravity
    : public accelerationForceMethod
{

        //- const access dict
        const dictionary& staticGravityDict_;

        //- staticGravity vector
        dimensionedVector g_;

        //- Href Field
        uniformDimensionedScalarField hRef_;

    // Private Member Functions

        //- Disallow copy construct
        staticGravity(const staticGravity&);

        //- Disallow default bitwise assignment
        void operator=(const staticGravity&);


        protected:

            //- Re-calculate the interface curvature
            virtual void calculateAcc();


public:

    //- Runtime type information
    TypeName("staticGravity");


    // Constructors

        //- Construct from dictionary
        staticGravity
        (
            interfaceCapturingMethod& ICM,
            const dictionary& dict
        );


    //- Destructor
    virtual ~staticGravity()
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
