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

#include "staticGravity.H"
#include "gravityMeshObject.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(staticGravity, 0);
    addToRunTimeSelectionTable(accelerationForceMethod,staticGravity, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::staticGravity::staticGravity
(
    interfaceCapturingMethod& ICM,
    const dictionary& dict
)
:
    accelerationForceMethod
    (
        ICM,
        dict
    ),
    staticGravityDict_(dict),
    g_
    (
        "gravity",
        dimAcceleration,
        vector(0,0,0)
    ),
    hRef_
    (
        IOobject
        (
            "hRef",
            acc_.mesh().time().constant(),
            acc_.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        dimensionedScalar(dimLength, Zero)
    )
{
    calculateAcc();
}


// * * * * * * * * * * * * * * Public Access Member Functions  * * * * * * * //

void Foam::staticGravity::calculateAcc()
{
    // read only if mesh changed would be clever
    const fvMesh& mesh = acc_.mesh();
    const uniformDimensionedVectorField& g = meshObjects::gravity::New(mesh.time());
    g_.value() = g.value();

    dimensionedScalar ghRef
    (
        mag(g_.value()) > SMALL
      ? g_ & (cmptMag(g_.value())/mag(g_.value()))*hRef_
      : dimensionedScalar("ghRef", g_.dimensions()*dimLength, 0)
    );

    acc_ = (g_ & mesh.C()) - ghRef;
    accf_ = (g_ & mesh.Cf()) - ghRef;
}


Foam::volScalarField& Foam::staticGravity::pressure
(
    volScalarField& p,
    volScalarField& p_rgh
)
{
    return p;
}


void Foam::staticGravity::updatePressure
(
    volScalarField& p,
    volScalarField& p_rgh,
    const volScalarField& rho
)
{
    p_rgh = p - rho*acc_;
}


void Foam::staticGravity::updateRefPressure
(
    volScalarField& p,
    volScalarField& p_rgh,
    const volScalarField& rho
)
{
    p == p_rgh + rho*acc_;
}


Foam::tmp<Foam::surfaceScalarField> Foam::staticGravity::accelerationForce()
{
    const fvMesh& mesh = acc_.mesh();
    const volScalarField& rho = mesh.lookupObject<volScalarField>("rho");
    return fvc::interpolate(rho)*(g_ & mesh.Sf())/mesh.magSf();
}

// ************************************************************************* //
