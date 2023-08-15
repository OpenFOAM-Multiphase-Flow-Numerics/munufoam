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

#include "gravityRecon.H"
#include "gravityMeshObject.H"
#include "fvc.H"

#include "accelerationForceGeoVoFMethod.H"
#include "geometricVoFMethod.H"
#include "cutFacePLIC.H"
#include "registerAccelerationModels.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    registerDerivedAccelerationModel(accelerationForceGeoVoFMethod,gravityRecon,geometricVoFMethod);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::gravityRecon::gravityRecon
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
    gravityReconDict_(dict),
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


Foam::gravityRecon::gravityRecon
(
    interfaceCapturingMethod& geoVoF,
    const dictionary& dict
)
:
    gravityRecon
    (
        checkCompatiablity(geoVoF),
        dict
    )
{

}


// * * * * * * * * * * * * * * Public Access Member Functions  * * * * * * * //

void Foam::gravityRecon::calculateAcc()
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


Foam::volScalarField& Foam::gravityRecon::pressure
(
    volScalarField& p,
    volScalarField& p_rgh
)
{
    return p_rgh;
}


void Foam::gravityRecon::updatePressure
(
    volScalarField& p,
    volScalarField& p_rgh,
    const volScalarField& rho
)
{
    p == p_rgh + rho*acc_;
}


void Foam::gravityRecon::updateRefPressure
(
    volScalarField& p,
    volScalarField& p_rgh,
    const volScalarField& rho
)
{
    p_rgh = p - rho*acc_;
}


Foam::tmp<Foam::surfaceScalarField> Foam::gravityRecon::accelerationForce()
{
    const fvMesh& mesh = acc_.mesh();
    const geometricVoF& surf = geoVoF_.surface(timeState::newState);
    const boolList& interfaceCell = surf.interfaceCell;
    const volVectorField& centre = surf.centre;
    const volScalarField& rho = mesh.lookupObject<volScalarField>("rho");

    accf_ = (g_ & mesh.Cf());

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

    forAll(accf_,faceI)
    {
        
        if (markedFace.found(faceI))
        {
            bool ownSurf = interfaceCell[owner[faceI]];
            bool neiSurf = interfaceCell[neighbour[faceI]];
            vector faceCentre = Zero;
            if (ownSurf)
            {
                faceCentre = centre[owner[faceI]];
            }
            if (neiSurf)
            {
                faceCentre = centre[neighbour[faceI]];
            }

            if (ownSurf & neiSurf)
            {
                faceCentre = 0.5*centre[owner[faceI]] + 0.5*centre[neighbour[faceI]];
            }
            accf_[faceI] = g_.value() & faceCentre;
        }
    }
    
    // return -accf_*fvc::snGrad(rho);
    return -dimensionedScalar("rhoJump",dimDensity,999)*accf_*fvc::snGrad(surf.alpha());
}

// ************************************************************************* //
