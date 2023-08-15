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

#include "accelerationForceGeoVoFMethod.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(accelerationForceGeoVoFMethod, 0);
    defineRunTimeSelectionTable(accelerationForceGeoVoFMethod, dictionary);
}

// * * * * * * * * * * * * * * Protected Access Member Functions   * * * * * //

Foam::geometricVoFMethod& Foam::accelerationForceGeoVoFMethod::checkCompatiablity(interfaceCapturingMethod& geoVoF)
{
    try
    {
        geometricVoFMethod& ICMPtr = dynamic_cast<geometricVoFMethod&>(geoVoF);
        return ICMPtr;
    }
    catch(std::bad_cast const&)
    {

        accelerationForceMethod::registerModel::printCompatibilityTable
        (
            FatalIOError,
            wordList
            ({
                "accelerationModel",
                "interfaceCapturingMethod"
            }),
            *accelerationForceMethod::compatibilityTable_
        ) << exit(FatalIOError);
    }

    return dynamic_cast<geometricVoFMethod&>(geoVoF);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::accelerationForceGeoVoFMethod::accelerationForceGeoVoFMethod
(
    geometricVoFMethod& geoVoF,
    const dictionary& dict
)
:
  accelerationForceMethod(geoVoF,dict),
  geoVoF_(geoVoF)
{

}


// * * * * * * * * * * * * * * Public Access Member Functions  * * * * * * * //


void Foam::accelerationForceGeoVoFMethod::calculateAcc()
{
    notImplemented("bool Foam::accelerationForceGeoVoFMethod::calculateAcc()");;
}



// ************************************************************************* //
