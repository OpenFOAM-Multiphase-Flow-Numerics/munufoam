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

#include "surfaceBase.H"

#include "alphaContactAngleTwoPhaseFvPatchScalarField.H"
#include "mathematicalConstants.H"
#include "surfaceInterpolate.H"
#include "fvcDiv.H"
#include "fvcGrad.H"

namespace Foam
{
    defineTypeNameAndDebug(surfaceBase, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::surfaceBase::surfaceBase
(
    volScalarField& alpha1,
    const surfaceScalarField& phi,
    const volVectorField& U
)
:
    alpha1_(alpha1),
    phi_(phi),
    U_(U),
    normal_
    (
        IOobject
        (
            IOobject::groupName("interfaceNormal", alpha1.group()),
            alpha1_.mesh().time().timeName(),
            alpha1_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        //     dict.getOrDefault("writeFields",false)
        //   ? IOobject::AUTO_WRITE
        //   : IOobject::NO_WRITE
        ),
        alpha1_.mesh(),
        dimensionedVector(dimArea, Zero)
    ),
    centre_
    (
        IOobject
        (
            IOobject::groupName("interfaceCentre", alpha1.group()),
            alpha1_.mesh().time().timeName(),
            alpha1_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        //     dict.getOrDefault("writeFields",false)
        //   ? IOobject::AUTO_WRITE
        //   : IOobject::NO_WRITE
        ),
        alpha1_.mesh(),
        dimensionedVector(dimLength, Zero)
    ),
    interfaceCell_(alpha1_.mesh().nCells(), false),
    interfaceLabels_(0.2*alpha1_.mesh().nCells()),
    timeIndexAndIter_(0, 0)
{

}


// * * * * * * * * * * * * * * Public Access Member Functions  * * * * * * * * * * * * * * //



// ************************************************************************* //
