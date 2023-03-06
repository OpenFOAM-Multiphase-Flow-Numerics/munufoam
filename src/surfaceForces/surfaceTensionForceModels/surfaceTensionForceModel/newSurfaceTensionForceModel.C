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

\*---------------------------------------------------------------------------*/

#include "surfaceTensionForceModel.H"
#include "interfaceCapturingMethod.H"
#include "geometricVoFMethod.H"
#include "stringOps.H"



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


Foam::autoPtr<Foam::surfaceTensionForceModel>
Foam::surfaceTensionForceModel::New
(
    const dictionary& dict,
    interfaceCapturingMethod& ICM
)
{

    word surfaceTensionForceModelTypeName
    (
        dict.get<word>("surfaceTensionForceModel")
    );

    Info<< "Selecting surfaceTension model "
        << surfaceTensionForceModelTypeName << endl;

 
    // allowing to correctly specify the type without looking up the interfaceCapturingMethod
    auto ctorPtr = dictionaryConstructorTable(surfaceTensionForceModelTypeName);

    // // get the typeName of the interfaceCapturingMethod
    // const interfaceCapturingMethod& ICM =
    //     alpha1.mesh().lookupObject<interfaceCapturingMethod>(interfaceCapturingMethod::typeName);

    if (!ctorPtr)
    {
        // concrete type failed this interfaceCapturingMethod ICM with surfaceTensionForceModel
        word combinedType = ICM.type() + "_" + surfaceTensionForceModelTypeName;
        ctorPtr = dictionaryConstructorTable(combinedType);
    }

    if (!ctorPtr)
    {
        // check ICM category with surfaceTensionForceModel
        const word ICMCategory = stringOps::split<word>(ICM.type(), '_').str(0);
        word combinedType = ICMCategory + "_" + surfaceTensionForceModelTypeName;
        ctorPtr = dictionaryConstructorTable(combinedType);
    }

    if (!ctorPtr)
    {
        // check with base class interfaceCapturingMethod with surfaceTensionForceModel
        word combinedType = "interfaceCapturingMethod_" + surfaceTensionForceModelTypeName;
        ctorPtr = dictionaryConstructorTable(combinedType);
    }

    if (!ctorPtr)
    {
        FatalErrorInLookup
        (
            "surfaceTensionForceModel",
            surfaceTensionForceModelTypeName,
            *dictionaryConstructorTablePtr_
        );

        printSurfaceForceModels
        (
            FatalIOError,
            wordList
            ({
                "interfaceCapturingMethod",
                "interfaceRepresentation",
                "surfaceTensionForceModel"
            }),
            dictionaryConstructorTablePtr_->sortedToc()
        ) << exit(FatalIOError);

    }

    return autoPtr<surfaceTensionForceModel>
    (
        ctorPtr
        (
            dict,
            ICM
        )
    );

}


// ************************************************************************* //
