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
    Foam::surfaceTensionForceModel

Description

SourceFiles
    surfaceTensionForceModel.C
    newsurfaceTensionForceModel.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceTensionForceModel_H
#define surfaceTensionForceModel_H

#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "dimensionedScalar.H"
#include "autoPtr.H"
#include "Pair.H"

#include "IOdictionary.H"
#include "volFields.H"
#include "surfaceFields.H"

#include "surfaceTensionModel.H"
#include "deltaFunctionModel.H"
#include "fvc.H"
#include "interfaceCapturingMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class surfaceTensionForceModel Declaration
\*---------------------------------------------------------------------------*/

class surfaceTensionForceModel
{

protected:

        //- reference to volume fraction field
        const volScalarField& alpha1_;

        //- reference to face velocity
        const surfaceScalarField& phi_;

        //- reference to velcoity
        const volVectorField& U_;

        autoPtr<surfaceTensionModel> surfTenModel_;
        
        //- surfaceTensionForce field
        surfaceScalarField surfaceTensionForce_;

private:

    // Private Member Functions

        static Ostream& printSurfaceForceModels
        (
            Ostream& os,
            const wordList& cmptNames,
            const wordList& surfaceForceNames
        );

        //- Disallow copy construct
        surfaceTensionForceModel(const surfaceTensionForceModel&);

        //- Disallow default bitwise assignment
        void operator=(const surfaceTensionForceModel&);

public:

    //- Runtime type information
    TypeName("surfaceTensionForceModel");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            surfaceTensionForceModel,
            dictionary,
            (
                const dictionary& dict,
                interfaceCapturingMethod& ICM
            ),
            (dict,ICM)
        );

    // Selectors

        //- Return a reference to the selected phaseChange model
        static autoPtr<surfaceTensionForceModel> New
        (
            const dictionary& dict,
            interfaceCapturingMethod& ICM
        );

    // Constructors

        //- Construct from dictionary
        surfaceTensionForceModel
        (
            const dictionary& dict,
            interfaceCapturingMethod& ICM
        );

    //- Destructor
    virtual ~surfaceTensionForceModel()
    {}

    // Member Functions


        virtual void correct() = 0;

        tmp<volScalarField> sigma() const
        {
            return surfTenModel_->sigma();
        }

        virtual const surfaceScalarField& surfaceTensionForce() const
        {
            return surfaceTensionForce_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
