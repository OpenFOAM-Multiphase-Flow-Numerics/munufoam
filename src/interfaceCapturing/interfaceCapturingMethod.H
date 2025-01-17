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
    Foam::interfaceCapturingMethod

Description
    Base class for advection schemes

SourceFiles
    interfaceCapturingMethod.C
    newinterfaceCapturingMethod.C

\*---------------------------------------------------------------------------*/

#ifndef interfaceCapturingMethod_H
#define interfaceCapturingMethod_H

#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "dimensionedScalar.H"
#include "autoPtr.H"
#include "IOdictionary.H"

#include "surfaceBase.H"
#include "geometricZeroField.H"




// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

enum timeState
{
    oldState,
    newState
};

/*---------------------------------------------------------------------------*\
                Class interfaceCapturingMethod Declaration
\*---------------------------------------------------------------------------*/

class interfaceCapturingMethod
    : public IOdictionary
{

private:
        //-
        dictionary interfaceCapturingMethodCoeffs_;


protected:

        // reference to the VoF Field
        volScalarField& alpha1_;

        // reference to the face fluxes
        const surfaceScalarField& phi_;

        // reference to the velocity field
        const volVectorField& U_;

        //- Face volumetric transport
        surfaceScalarField alphaPhi_;

    // Private Member Functions

private:

        //- Disallow copy construct
        interfaceCapturingMethod(const interfaceCapturingMethod&);

        //- Disallow default bitwise assignment
        void operator=(const interfaceCapturingMethod&);


public:

    //-
    dictionary& modelDict();

    //-
    const dictionary& modelDict() const;


public:

    //- Runtime type information
    TypeName("interfaceCapturingMethod");

    //- Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        interfaceCapturingMethod,
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
        static autoPtr<interfaceCapturingMethod> New
        (
            volScalarField& alpha1,
            const surfaceScalarField& phi,
            const volVectorField& U
        );


    // Constructors

        //- Construct from dictionary
        interfaceCapturingMethod
        (
            volScalarField& alpha1,
            const surfaceScalarField& phi,
            const volVectorField& U
        );


    //- Destructor
    virtual ~interfaceCapturingMethod() = default;


    // Member Functions

        //- advection of the interface
        virtual void advect()
        {

            volScalarField::Internal Sp
            (
                IOobject
                (
                    "Sp",
                    alpha1_.mesh().time().timeName(),
                    alpha1_.mesh()
                ),
                alpha1_.mesh(),
                dimensionedScalar(dimless/dimTime, Zero)
            );

            volScalarField::Internal Su
            (
                IOobject
                (
                    "Su",
                    alpha1_.mesh().time().timeName(),
                    alpha1_.mesh()
                ),
                alpha1_.mesh(),
                dimensionedScalar(dimless/dimTime, Zero)
            );
            advect(Sp,Su);
        };

        //- advection of the interface with source terms
        virtual void advect
        (
            const volScalarField::Internal& Sp,
            const volScalarField::Internal& Su
        ) = 0;

        //- Reconstruct interface
        virtual const surfaceBase& surface(timeState state) = 0;

        template < class Type >
        tmp<surfaceScalarField> getRhoPhi
        (
            const Type rho1,
            const Type rho2
        ) const
        {
            return tmp<surfaceScalarField>
            (
                new surfaceScalarField
                (
                    "rhoPhi",
                    (rho1 - rho2)*alphaPhi_ + rho2*phi_
                )
            );
        }
        
        surfaceScalarField& alphaPhi()
        {
            return alphaPhi_;
        }

        const volVectorField& U() const
        {
            return U_;
        }

        const surfaceScalarField& phi() const
        {
            return phi_;
        }

        volScalarField& alpha1()
        {
            return alpha1_;
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
