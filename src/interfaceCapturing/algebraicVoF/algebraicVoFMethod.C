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

\*---------------------------------------------------------------------------*/

#include "algebraicVoFMethod.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(algebraicVoFMethod, 0);
    defineRunTimeSelectionTable(algebraicVoFMethod, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::algebraicVoFMethod::algebraicVoFMethod
(
    volScalarField& alpha1,
    const surfaceScalarField& phi,
    const volVectorField& U
)
:
  interfaceCapturingMethod(alpha1,phi,U)
{


}



Foam::autoPtr<Foam::algebraicVoFMethod>
Foam::algebraicVoFMethod::New
(
    volScalarField& alpha1,
    const surfaceScalarField& phi,
    const volVectorField& U
)
{
    IOdictionary fvSolutionDict
    (
        IOobject
        (
            "fvSolution",
            alpha1.time().system(),
            alpha1.db(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE,
            false
        )
    );
    const dictionary& alphaDict = fvSolutionDict.subDict("solvers").subDict(alpha1.name());

    word algebraicVoFMethodTypeName = 
    (
        alphaDict.get<word>("interfaceType") + "_"
      + alphaDict.get<word>("interfaceRepresentation")  + "_"
      + alphaDict.get<word>("interfaceCapturingScheme")
    );

    Info<< "Selecting interfaceCapturingScheme: "
        << algebraicVoFMethodTypeName << endl;

    auto* ctorPtr = dictionaryConstructorTable(algebraicVoFMethodTypeName);

    if (!ctorPtr)
    {
        FatalErrorInLookup
        (
            "algebraicVoFMethodScheme",
            algebraicVoFMethodTypeName,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalError);
    }

    return autoPtr<algebraicVoFMethod>(ctorPtr( alpha1, phi,U));
}



// ************************************************************************* //
