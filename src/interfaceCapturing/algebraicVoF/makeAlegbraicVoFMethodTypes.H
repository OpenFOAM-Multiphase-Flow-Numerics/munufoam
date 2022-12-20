/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#ifndef makeAlegbraicVoFMethodTypes_H
#define makeAlegbraicVoFMethodTypes_H

#include "addToRunTimeSelectionTable.H"
#include "interfaceCapturingMethod.H"
#include "dictionary.H"
#include "algebraicVoFModel.H"
#include "interfaceRepresentationModel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeAlegbraicVoFMethodTypes(ICM,Surface,SurfaceStrategy,AdvectionStrategy)    \
                                                                            \
    namespace Foam                                                          \
    {                                                                       \
        typedef Foam::interfaceRepresentationModel                          \
            <                                                               \
                Surface,                                                    \
                SurfaceStrategy                                             \
            >                                                               \
            Surface##_##SurfaceStrategy;                                    \
                                                                            \
                                                                            \
        typedef Foam::algebraicVoFModel                                     \
            <                                                               \
                Surface,                                                    \
                SurfaceStrategy,                                            \
                AdvectionStrategy                                           \
            >                                                               \
            Surface##_##SurfaceStrategy##_##AdvectionStrategy;              \
                                                                            \
        defineTemplateTypeNameAndDebug                                      \
        (                                                                   \
            Surface##_##SurfaceStrategy##_##AdvectionStrategy,              \
            0                                                               \
        );                                                                  \
                                                                            \
        addNamedToRunTimeSelectionTable                                     \
        (                                                                   \
            ICM,                                                            \
            Surface##_##SurfaceStrategy##_##AdvectionStrategy,              \
            dictionary,                                                     \
            Surface##_##SurfaceStrategy##_##AdvectionStrategy               \
        );                                                                  \
                                                                            \
        addNamedToRunTimeSelectionTable                                     \
        (                                                                   \
            interfaceCapturingMethod,                                       \
            Surface##_##SurfaceStrategy##_##AdvectionStrategy,              \
            dictionary,                                                     \
            Surface##_##SurfaceStrategy##_##AdvectionStrategy               \
        );                                                                  \
    }

#endif

// ************************************************************************* //

