/*---------------------------------------------------------------------------*\
    Modified work | Copyright (c) 2017-2019, German Aerospace Center (DLR)
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
    Foam::geometricVoFMethod

Description
    Base class for advection schemes

SourceFiles
    geometricVoFMethod.C
    newgeometricVoFMethod.C

\*---------------------------------------------------------------------------*/

#ifndef geometricVoFMethod_H
#define geometricVoFMethod_H

#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "dimensionedScalar.H"
#include "autoPtr.H"
#include "IOdictionary.H"

#include "interfaceRepresentation.H"
#include "interfaceCapturingMethod.H"
#include "geometricVoF.H"




// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class geometricVoFMethod Declaration
\*---------------------------------------------------------------------------*/

class geometricVoFMethod
:
    public interfaceCapturingMethod
{

private:

        //- Disallow copy construct
        geometricVoFMethod(const geometricVoFMethod&);

        //- Disallow default bitwise assignment
        void operator=(const geometricVoFMethod&);

public:

    //- Runtime type information
    TypeName("geometricVoFMethod");

    //- Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        geometricVoFMethod,
        dictionary,
        (
            volScalarField& alpha1,
            const surfaceScalarField& phi,
            const volVectorField& U
        ),
        (alpha1, phi,U)
    );

    // Selectors

        //- Return a reference to the selected phaseChange model
        static autoPtr<geometricVoFMethod> New
        (
            volScalarField& alpha1,
            const surfaceScalarField& phi,
            const volVectorField& U
        );


    // Constructors

        //- Construct from dictionary
        geometricVoFMethod
        (
            volScalarField& alpha1,
            const surfaceScalarField& phi,
            const volVectorField& U
        );


    //- Destructor
    virtual ~geometricVoFMethod() = default;


    // Member Functions

        //- advection of the interface
        // virtual void advect() = 0;

        //- advection of the interface with source terms
        // virtual void advect
        // (
        //     const volScalarField::Internal& Sp,
        //     const volScalarField::Internal& Su
        // ) = 0;

        //- Reconstruct interface
        virtual const geometricVoF& surface(timeState state) = 0;

        // template < class Type >
        // tmp<surfaceScalarField> getRhoPhi
        // (
        //     const Type rho1,
        //     const Type rho2
        // ) const
        // {
        //     return tmp<surfaceScalarField>
        //     (
        //         new surfaceScalarField
        //         (
        //             "rhoPhi",
        //             (rho1 - rho2)*alphaPhi_ + rho2*phi_
        //         )
        //     );
        // }
        
        // const surfaceScalarField& alphaPhi() const
        // {
        //     return alphaPhi_;
        // }

        // virtual const interfaceRepresentation& surf() const  = override;

        // virtual interfaceRepresentation& surf() = override;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
