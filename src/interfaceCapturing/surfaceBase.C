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
    alpha1(alpha1),
    phi(phi),
    U(U),
    interfaceCell(alpha1.mesh().nCells(), false),
    interfaceLabels(0.2*alpha1.mesh().nCells()),
    timeIndexAndIter(0, 0)
{

}


// * * * * * * * * * * * * * * Public Access Member Functions  * * * * * * * * * * * * * * //



// ************************************************************************* //
