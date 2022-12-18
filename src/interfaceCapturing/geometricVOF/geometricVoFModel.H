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
    Foam::geometricVoFModel

Description
    Base class for advection schemes

SourceFiles
    geometricVoFModel.C

\*---------------------------------------------------------------------------*/

#ifndef interfaceCapturingModel_H
#define interfaceCapturingModel_H

#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "geometricVoFMethod.H"
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
                Class geometricVoFModel Declaration
\*---------------------------------------------------------------------------*/
template<class Surface, class SurfaceStrategy, class AdvectionStrategy>
class geometricVoFModel
: 
    public geometricVoFMethod
{
    // Private Member Functions

    interfaceRepresentationModel<Surface,SurfaceStrategy> surface_;
    AdvectionStrategy advectionStrat_;

        //- Disallow copy construct
        geometricVoFModel(const geometricVoFModel&);

        //- Disallow default bitwise assignment
        void operator=(const geometricVoFModel&);


public:

    //- Runtime type information
    TypeName("geometricVoFModel");

    // Selectors

        // //- Return a reference to the selected phaseChange model
        // static autoPtr<geometricVoFModel> New
        // (
        //     volScalarField& alpha1,
        //     const surfaceScalarField& phi,
        //     const volVectorField& U
        // );


    // Constructors
        geometricVoFModel
        (
            volScalarField& alpha1,
            const surfaceScalarField& phi,
            const volVectorField& U
        );


    //- Destructor
    virtual ~geometricVoFModel() = default;


    // Member Functions

        // //- advection of the interface
        void advect() override;

        //- Reconstruct interface
        const Surface& surface(timeState state) override;


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#ifdef NoRepository
#   include "geometricVoFModel.C"
#endif
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
