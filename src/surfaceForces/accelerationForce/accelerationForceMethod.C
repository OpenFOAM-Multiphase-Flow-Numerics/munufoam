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

#include "accelerationForceMethod.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(accelerationForceMethod, 0);
    defineRunTimeSelectionTable(accelerationForceMethod, dictionary);
    std::unique_ptr<HashTable<word>> accelerationForceMethod::compatibilityTable_(new HashTable<word>());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Foam::accelerationForceMethod::registerModel
// (
//     word modelName,
//     word compatibleICM
// )
// :
//     modelName_(modelName),
//     compatibleICM_(compatibleICM),
// {

// }


Foam::accelerationForceMethod::accelerationForceMethod
(
    interfaceCapturingMethod& ICM,
    const dictionary& dict
)
:
    accf_
    (
        IOobject
        (
            "ghf",
            ICM.alpha1().mesh().time().timeName(),
            ICM.alpha1().mesh()
        ),
        ICM.alpha1().mesh(),
        dimensionedScalar("accf_", dimAcceleration*dimLength, 0.0)
    ),
    acc_
    (
        IOobject
        (
            "gh",
            ICM.alpha1().mesh().time().timeName(),
            ICM.alpha1().mesh()
        ),
        ICM.alpha1().mesh(),
        dimensionedScalar("acc_", dimAcceleration*dimLength, 0.0)
    )
{

}


// * * * * * * * * * * * * * * Public Access Member Functions  * * * * * * * //


void Foam::accelerationForceMethod::calculateAcc()
{
    notImplemented("bool Foam::accelerationForceMethod::calculateAcc()");;
}



// ************************************************************************* //
