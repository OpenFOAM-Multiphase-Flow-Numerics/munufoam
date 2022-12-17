/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2020 DLR
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
    Foam::interfaceRepresentationModel

Description
    Original code supplied by Henning Scheufler, DLR (2019)

SourceFiles
    interfaceRepresentationModel.C

\*---------------------------------------------------------------------------*/

#ifndef interfaceRepresentationModel_H
#define interfaceRepresentationModel_H

#include "typeInfo.H"
#include "dimensionedScalar.H"
#include "volFields.H"
#include "DynamicField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class interfaceRepresentationModel Declaration
\*---------------------------------------------------------------------------*/

template<class Surface, class SurfaceStrategy>
class interfaceRepresentationModel
:
    public interfaceRepresentation
{
    // Private Data

        dictionary interfaceRepresentationModelCoeffs_;

        Surface oldSurface_;

        Surface newSurface_;

        SurfaceStrategy surfaceStrat_;

        //- No copy construct
        interfaceRepresentationModel(const interfaceRepresentationModel&) = delete;

        //- No copy assignment
        void operator=(const interfaceRepresentationModel&) = delete;


public:
    // Constructors

        TypeName("interfaceRepresentationModel");

        //- Construct from components
        interfaceRepresentationModel
        (
            volScalarField& alpha1,
            const surfaceScalarField& phi,
            const volVectorField& U,
            const dictionary& dict
        );


    //- Destructor
    virtual ~interfaceRepresentationModel() = default;


    // Member Functions

        //- Access to the model dictionary
        dictionary& modelDict()
        {
            return interfaceRepresentationModelCoeffs_;
        }

        //- Const access to the model dictionary
        const dictionary& modelDict() const
        {
            return interfaceRepresentationModelCoeffs_;
        }

        //- Reconstruct interface
        virtual const Surface& surface(timeState state) override;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#ifdef NoRepository
#   include "interfaceRepresentationModel.C"
#endif
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
