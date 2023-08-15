/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

#include "gravityRecon2.H"
#include "addToRunTimeSelectionTable.H"


#include "gravityMeshObject.H"
#include "plane.H"
#include "fvc.H"

#include "interfaceCapturingMethod.H"
#include "geometricVoFMethod.H"
#include "cutFacePLIC.H"
#include "registerAccelerationModels.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    registerDerivedAccelerationModel(accelerationForceGeoVoFMethod,gravityRecon2,geometricVoFMethod);
}

Foam::vector Foam::gravityRecon2::closestDist(const point p, const vector n ,const vector centre)
{
    vector normal = n/mag(n);
    return p - normal*((p - centre) & normal);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::gravityRecon2::gravityRecon2
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
    gravityDict_(dict),
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
    // calculateAcc();
}

Foam::gravityRecon2::gravityRecon2
(
    interfaceCapturingMethod& geoVoF,
    const dictionary& dict
)
:
    gravityRecon2
    (
        checkCompatiablity(geoVoF),
        dict
    )
{

}

// * * * * * * * * * * * * * * Public Access Member Functions  * * * * * * * * * * * * * * //


Foam::volScalarField& Foam::gravityRecon2::pressure
(
    volScalarField& p,
    volScalarField& p_rgh
)
{
    return p_rgh;
}


void Foam::gravityRecon2::updatePressure
(
    volScalarField& p,
    volScalarField& p_rgh,
    const volScalarField& rho
)
{
    p == p_rgh + rho*acc_;
}


void Foam::gravityRecon2::updateRefPressure
(
    volScalarField& p,
    volScalarField& p_rgh,
    const volScalarField& rho
)
{
    p_rgh = p - rho*acc_;
}


void Foam::gravityRecon2::calculateAcc()
{
    const fvMesh& mesh = acc_.mesh();
    zoneDistribute& exchangeFields = zoneDistribute::New(mesh);
    const uniformDimensionedVectorField& g = meshObjects::gravity::New(mesh.time());
    g_.value() = g.value();
    dimensionedScalar hRef("hRef",dimLength, gravityDict_.lookupOrDefault("hRef",0));


    dimensionedScalar ghRef
    (
        mag(g_.value()) > SMALL
      ? g_ & (cmptMag(g_.value())/mag(g_.value()))*hRef
      : dimensionedScalar("ghRef", g_.dimensions()*dimLength, 0)
    );

    acc_ = (g_ & mesh.C()) - ghRef;
    accf_ = ((g_ & mesh.Cf()) - ghRef);

    const geometricVoF& surf = geoVoF_.surface(timeState::newState);

    const volVectorField& faceCentre = surf.centre;
    const volVectorField& faceNormal = surf.normal;

    boolList nextToInterface(mesh.nCells(),false);
    markInterfaceRegion markIF(mesh);

    markIF.markCellsNearSurf(surf.interfaceCell,1,nextToInterface);

    exchangeFields.setUpCommforZone(nextToInterface,true);

    Map<vector > mapCentres =
        exchangeFields.getDatafromOtherProc(nextToInterface,faceCentre);
    Map<vector > mapNormal =
        exchangeFields.getDatafromOtherProc(nextToInterface,faceNormal);

    const labelListList& stencil = exchangeFields.getStencil();

    forAll(surf.interfaceCell,celli)
    {
        if(mag(faceNormal[celli]) != 0)
        {
            vector closeP = closestDist(mesh.C()[celli],-faceNormal[celli],faceCentre[celli]);
            acc_[celli] = closeP & g_.value();
        }
        else if(nextToInterface[celli])
        {
            // the to the interface
            vector averageSurfP = vector::zero;
            scalar avgWeight = 0;
            const point p = mesh.C()[celli];

            forAll(stencil[celli],i)
            {
                const label& gblIdx = stencil[celli][i];
                vector n = -exchangeFields.getValue(faceNormal,mapNormal,gblIdx);
                if (mag(n) != 0)
                {
                    n /= mag(n);
                    vector c =
                        exchangeFields.getValue(faceCentre,mapCentres,gblIdx);
                    vector distanceToIntSeg = (c - p);
                    vector closeP = closestDist(p,n,c);
                    scalar weight = 0;

                    if (mag(distanceToIntSeg) != 0)
                    {
                        distanceToIntSeg /= mag(distanceToIntSeg);
                        weight = sqr(mag(distanceToIntSeg & n));
                    }
                    else // exactly on the center
                    {
                        weight = 1;
                    }
                    averageSurfP += closeP * weight;
                    avgWeight += weight;
                }
            }

            if (avgWeight != 0)
            {
                averageSurfP /= avgWeight;
                acc_[celli] = averageSurfP & g_.value();
            }

        }
        else
        {
            // do nothing
        }

    }

    acc_.correctBoundaryConditions();
    accf_ = fvc::interpolate(acc_);

}


Foam::tmp<Foam::surfaceScalarField> Foam::gravityRecon2::accelerationForce()
{
    const fvMesh& mesh = acc_.mesh();
    const geometricVoF& surf = geoVoF_.surface(timeState::newState);
    const boolList& interfaceCell = surf.interfaceCell;
    const volVectorField& centre = surf.centre;
    const volScalarField& rho = mesh.lookupObject<volScalarField>("rho");
    // acc_ = (g_ & mesh.C());
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
    return -(fvc::snGrad(acc_*rho) - rhof*fvc::snGrad(acc_));
}


// ************************************************************************* //
