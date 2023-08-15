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

#include "staticGravityPLIC.H"
#include "gravityMeshObject.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"
#include "interfaceCapturingMethod.H"
#include "cutFacePLIC.H"
#include "registerAccelerationModels.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    registerDerivedAccelerationModel(accelerationForceGeoVoFMethod,staticGravityPLIC,geometricVoFMethod);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::staticGravityPLIC::staticGravityPLIC
(
    geometricVoFMethod& geoVoF,
    const dictionary& dict
)
:
    accelerationForceGeoVoFMethod
    (
        geoVoF,
        dict
    ),
    staticGravityPLICDict_(dict),
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
            geoVoF.alpha1().mesh().time().constant(),
            geoVoF.alpha1().mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        dimensionedScalar(dimLength, Zero)
    )
{
    calculateAcc();
}


Foam::staticGravityPLIC::staticGravityPLIC
(
    interfaceCapturingMethod& geoVoF,
    const dictionary& dict
)
:
    staticGravityPLIC
    (
        checkCompatiablity(geoVoF),
        dict
    )
{

}

// * * * * * * * * * * * * * * Public Access Member Functions  * * * * * * * //

void Foam::staticGravityPLIC::calculateAcc()
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


Foam::volScalarField& Foam::staticGravityPLIC::pressure
(
    volScalarField& p,
    volScalarField& p_rgh
)
{
    return p;
}


void Foam::staticGravityPLIC::updatePressure
(
    volScalarField& p,
    volScalarField& p_rgh,
    const volScalarField& rho
)
{
    p_rgh = p - rho*acc_;
}


void Foam::staticGravityPLIC::updateRefPressure
(
    volScalarField& p,
    volScalarField& p_rgh,
    const volScalarField& rho
)
{
    // p_rgh = p - rho*acc_;
    p == p_rgh + rho*acc_;
}


Foam::tmp<Foam::surfaceScalarField> Foam::staticGravityPLIC::accelerationForce()
{
    const fvMesh& mesh = acc_.mesh();
    const volScalarField& rho = mesh.lookupObject<volScalarField>("rho");

    const geometricVoF& surf = geoVoF_.surface(timeState::newState);
    const boolList& interfaceCell = surf.interfaceCell;
    const volVectorField& centre = surf.centre;
    forAll(acc_,celli)
    {
        if (!interfaceCell[celli])
        {
            acc_[celli] = (g_.value() & mesh.C()[celli]); // - ghRef;
        }
        else
        {
            acc_[celli] = (g_.value() & centre[celli]); // - ghRef;
        }

    }
    surfaceScalarField rhof(fvc::interpolate(rho));

    // Mark faces using any marked cell
    bitSet markedFace(mesh.nFaces());

    for (const label celli : surf.interfaceLabels)
    {
        markedFace.set(mesh.cells()[celli]);  // set multiple faces
    }

    // syncTools::syncFaceList(mesh, markedFace, orEqOp<unsigned int>());
    cutFacePLIC cutFace(mesh);

    const labelList& owner = mesh.faceOwner();
    const labelList& neighbour = mesh.faceNeighbour();

    forAll(rhof,faceI)
    {
        
        if (markedFace.found(faceI))
        {
            label celli = -1;
            bool ownSurf = interfaceCell[owner[faceI]];
            bool neiSurf = interfaceCell[neighbour[faceI]];
            if (ownSurf)
            {
                celli = owner[faceI];
            }
            if (neiSurf)
            {
                celli = neighbour[faceI];
            }
            vector n = surf.normal[celli];
            if(mag(n) != 0)
            {
                n /= mag(n);
                label cutStatus = cutFace.calcSubFace
                (
                    faceI,
                    n,
                    surf.centre[celli]
                );
            }
            scalar alphaf = mag(cutFace.subFaceArea())/mesh.magSf()[faceI];
            if (!ownSurf)
            {
                alphaf = surf.alpha()[owner[faceI]];
            }
            if (!neiSurf)
            {
                alphaf = surf.alpha()[neighbour[faceI]];
            }

            // Info << "celli " << celli << endl;
            // Info << "alphaf " << alphaf << endl;
            // Info << "rhof " << rhof[faceI] << " rhof " << 1000*alphaf + (1-alphaf)*1  << endl;
            rhof[faceI] = 1000*alphaf + (1-alphaf)*1;
        }
    }
    return rhof*fvc::snGrad(acc_);
}

// ************************************************************************* //
