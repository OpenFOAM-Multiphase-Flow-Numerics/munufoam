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
    Foam::algebraicVoFModel

Description
    Base class for advection schemes

SourceFiles
    algebraicVoFModel.C

\*---------------------------------------------------------------------------*/

#ifndef interfaceCapturingModel_H
#define interfaceCapturingModel_H

#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "algebraicVoFMethod.H"
#include "interfaceCapturingMethod.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "dimensionedScalar.H"
#include "autoPtr.H"
#include "IOdictionary.H"
#include "interfaceRepresentationModel.H"




// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class algebraicVoFModel Declaration
\*---------------------------------------------------------------------------*/
template<class Surface, class SurfaceStrategy, class AdvectionStrategy>
class algebraicVoFModel
: 
    public algebraicVoFMethod
{
    // Private Member Functions

    interfaceRepresentationModel<Surface,SurfaceStrategy> surface_;
    AdvectionStrategy advectionStrat_;

        //- Disallow copy construct
        algebraicVoFModel(const algebraicVoFModel&);

        //- Disallow default bitwise assignment
        void operator=(const algebraicVoFModel&);


public:

    //- Runtime type information
    TypeName("algebraicVoFModel");

    // Selectors

        // //- Return a reference to the selected phaseChange model
        // static autoPtr<algebraicVoFModel> New
        // (
        //     volScalarField& alpha1,
        //     const surfaceScalarField& phi,
        //     const volVectorField& U
        // );


    // Constructors
        algebraicVoFModel
        (
            volScalarField& alpha1,
            const surfaceScalarField& phi,
            const volVectorField& U
        );


    //- Destructor
    virtual ~algebraicVoFModel() = default;


    // Member Functions

        // //- advection of the interface
        void advect() override;

        virtual void advect
        (
            const volScalarField::Internal& Sp,
            const volScalarField::Internal& Su
        ) override;

        //- Reconstruct interface
        const Surface& surface(timeState state) override;


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#ifdef NoRepository
#   include "algebraicVoFModel.C"
#endif
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
