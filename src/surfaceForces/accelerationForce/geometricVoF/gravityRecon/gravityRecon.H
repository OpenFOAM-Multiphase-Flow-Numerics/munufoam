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
    Foam::gravityRecon

Description
    constant gravitional forces

SourceFiles
    gravityRecon.C

\*---------------------------------------------------------------------------*/

#ifndef gravityRecon_H
#define gravityRecon_H

#include "typeInfo.H"
#include "volFields.H"
#include "dimensionedScalar.H"
#include "autoPtr.H"
#include "uniformDimensionedFields.H"

#include "accelerationForceGeoVoFMethod.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class gravityRecon Declaration
\*---------------------------------------------------------------------------*/

class gravityRecon
    : public accelerationForceGeoVoFMethod
{

        //- const access dict
        const dictionary& gravityReconDict_;

        //- gravityRecon vector
        dimensionedVector g_;

        //- Href Field
        uniformDimensionedScalarField hRef_;

    // Private Member Functions

        //- Disallow copy construct
        gravityRecon(const gravityRecon&);

        //- Disallow default bitwise assignment
        void operator=(const gravityRecon&);


        protected:

        //- Re-calculate the interface curvature
        virtual void calculateAcc();


public:

    //- Runtime type information
    TypeName("gravityRecon");


    // Constructors

        //- Construct from dictionary
        gravityRecon
        (
            geometricVoFMethod& geoVoF,
            const dictionary& dict
        );

        gravityRecon
        (
            interfaceCapturingMethod& geoVoF,
            const dictionary& dict
        );



    //- Destructor
    virtual ~gravityRecon()
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
