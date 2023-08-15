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

#ifndef registerSurfaceTensionModel_H
#define registerSurfaceTensionModel_H

#include "addToRunTimeSelectionTable.H"
#include "surfaceTensionForceModel.H"
#include "dictionary.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#define registerBaseSurfaceTensionModel(surfaceTensionModel)                \
                                                                            \
                                                                            \
    defineTypeNameAndDebug(surfaceTensionModel, 0);                         \
                                                                            \
    addToRunTimeSelectionTable                                              \
    (                                                                       \
        surfaceTensionForceModel,                                          \
        surfaceTensionModel,                                                \
        dictionary                                                          \
    );                                                                      \
                                                                            \
    surfaceTensionForceModel::registerModel register##surfaceTensionModel  \
    (                                                                       \
        #surfaceTensionModel,                                               \
        "any"                                                               \
    );


#define registerDerivedSurfaceTensionModel(modelType,surfaceTensionModel,ICTYPE)\
                                                                            \
                                                                            \
    defineTypeNameAndDebug(surfaceTensionModel, 0);                         \
                                                                            \
    addToRunTimeSelectionTable                                              \
    (                                                                       \
        surfaceTensionForceModel,                                          \
        surfaceTensionModel,                                                \
        dictionary                                                          \
    );                                                                      \
                                                                            \
    addToRunTimeSelectionTable                                              \
    (                                                                       \
        modelType,                                                          \
        surfaceTensionModel,                                                \
        dictionary                                                          \
    );                                                                      \
                                                                            \
    surfaceTensionForceModel::registerModel register##surfaceTensionModel  \
    (                                                                       \
        #surfaceTensionModel,                                               \
        #ICTYPE                                                             \
    );


#endif

// ************************************************************************* //

