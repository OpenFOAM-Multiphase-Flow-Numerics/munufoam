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
    Foam::accelerationForceGeoVoFMethod

Description
    abstract base class of the acceleration model
    handles dynamic and static formulation of the gravity/acceleration source terms

SourceFiles
    accelerationForceGeoVoFMethod.C

\*---------------------------------------------------------------------------*/

#ifndef accelerationForceGeoVoFMethod_H
#define accelerationForceGeoVoFMethod_H

#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "dimensionedScalar.H"
#include "autoPtr.H"
#include "Pair.H"

#include "IOdictionary.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "interfaceCapturingMethod.H"
#include "geometricVoFMethod.H"
#include "accelerationForceMethod.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class accelerationForceGeoVoFMethod Declaration
\*---------------------------------------------------------------------------*/

class accelerationForceGeoVoFMethod
:
    public accelerationForceMethod
{
protected:

    geometricVoFMethod& geoVoF_;


    geometricVoFMethod& checkCompatiablity(interfaceCapturingMethod& geoVoF);


    //- Re-calculate the acceleration
    virtual void calculateAcc();

private:

    // Private Member Functions

        //- Disallow copy construct
        accelerationForceGeoVoFMethod(const accelerationForceGeoVoFMethod&);

        //- Disallow default bitwise assignment
        void operator=(const accelerationForceGeoVoFMethod&);

public:

    //- Runtime type information
    TypeName("accelerationForceGeoVoFMethod");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            accelerationForceGeoVoFMethod,
            dictionary,
            (
                geometricVoFMethod& geoVoF,
                const dictionary& dict
            ),
            (geoVoF,dict)
        );

    // Selectors

        //- Return a reference to the selected phaseChange model
        static autoPtr<accelerationForceGeoVoFMethod> New
        (
            geometricVoFMethod& geoVoF,
            const dictionary& dict
        );

    // Constructors

        //- Construct from dictionary
        accelerationForceGeoVoFMethod
        (
            geometricVoFMethod& geoVoF,
            const dictionary& dict
        );

    //- Destructor
    virtual ~accelerationForceGeoVoFMethod()
    {}

    // Member Functions

        const surfaceScalarField& accf() const
        {
            return accf_;
        }

        const volScalarField& acc() const
        {
            return acc_;
        }

        //- return the pressure field depending
        virtual volScalarField& pressure
        (
            volScalarField& p,
            volScalarField& p_rgh
        ) = 0;

        //- update static pressure field
        virtual void updatePressure
        (
            volScalarField& p,
            volScalarField& p_rgh,
            const volScalarField& rho
        ) = 0;

        //- update dynamic pressure field
        virtual void updateRefPressure
        (
            volScalarField& p,
            volScalarField& p_rgh,
            const volScalarField& rho
        ) = 0;

        //- update the acceleration fields
        void correct()
        {
            calculateAcc();
        }

        //- return the accleration forces
        virtual tmp<surfaceScalarField> accelerationForce() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
