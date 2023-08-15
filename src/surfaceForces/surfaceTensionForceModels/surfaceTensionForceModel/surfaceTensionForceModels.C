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

    if (!ctorPtr)
    {
        surfaceTensionForceModel::registerModel::printCompatibilityTable
        (
            FatalIOError,
            wordList
            ({
                "accelerationModel",
                "interfaceCapturingMethod"
            }),
            *surfaceTensionForceModel::compatibilityTable_
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
