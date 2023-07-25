/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2020 DLR
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

#include "plicRDF.H"
#include "interpolationCellPoint.H"
#include "fvc.H"
#include "leastSquareGrad.H"
// #include "addToRunTimeSelectionTable.H"
#include "alphaContactAngleTwoPhaseFvPatchScalarField.H"
#include "profiling.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

    defineTypeNameAndDebug(plicRDF, 0);

}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::plicRDF::interpolateNormal(const geometricVoF& surf)
{
    addProfilingInFunction(geometricVoF);
    scalar dt = mesh_.time().deltaTValue();
    zoneDistribute& exchangeFields = zoneDistribute::New(mesh_);
    const volScalarField& alpha = surf.alpha();

    leastSquareGrad<scalar> lsGrad("polyDegree1",mesh_.geometricD());

    exchangeFields.setUpCommforZone(surf.interfaceCell,false);

    Map<vector> mapCentre
    (
        exchangeFields.getDatafromOtherProc(surf.interfaceCell, surf.centre)
    );
    Map<vector> mapNormal
    (
        exchangeFields.getDatafromOtherProc(surf.interfaceCell, surf.normal)
    );

    Map<vector> mapCC
    (
        exchangeFields.getDatafromOtherProc(surf.interfaceCell, mesh_.C())
    );
    Map<scalar> mapAlpha
    (
        exchangeFields.getDatafromOtherProc(surf.interfaceCell, alpha)
    );

    DynamicField<vector > cellCentre(100);
    DynamicField<scalar > alphaValues(100);

    DynamicList<vector> foundNormals(30);

    const labelListList& stencil = exchangeFields.getStencil();

    forAll(surf.interfaceLabels, i)
    {
        const label celli = surf.interfaceLabels[i];
        vector estimatedNormal = vector::zero;
        scalar weight = 0;
        foundNormals.clear();
        forAll(stencil[celli], i)
        {
            const label& gblIdx = stencil[celli][i];
            vector n =
                exchangeFields.getValue(surf.normal, mapNormal, gblIdx);
            point p = mesh_.C()[celli]-U_[celli]*dt;
            if (mag(n) != 0)
            {
                n /= mag(n);
                vector centre =
                    exchangeFields.getValue(surf.centre, mapCentre, gblIdx);
                vector distanceToIntSeg = (tensor::I- n*n) & (p - centre);
                estimatedNormal += n /max(mag(distanceToIntSeg), SMALL);
                weight += 1/max(mag(distanceToIntSeg), SMALL);
                foundNormals.append(n);
            }
        }

        if (weight != 0 && mag(estimatedNormal) != 0)
        {
            estimatedNormal /= weight;
            estimatedNormal /= mag(estimatedNormal);
        }

        bool tooCoarse = false;

        if (foundNormals.size() > 1 && mag(estimatedNormal) != 0)
        {
            forAll(foundNormals, i)
            {
                // all have the length of 1
                // to coarse if normal angle is bigger than 10 deg
                if ((estimatedNormal & foundNormals[i]) <= 0.98)
                {
                    tooCoarse = true;
                }
            }
        }
        else
        {
            tooCoarse = true;
        }

        // if a normal was found and the interface is fine enough
        // smallDist is always smallDist
        if (mag(estimatedNormal) != 0 && !tooCoarse)
        {
            interfaceNormal_[i] = estimatedNormal;
        }
        else
        {
            cellCentre.clear();
            alphaValues.clear();

            forAll(stencil[celli],i)
            {
                const label& gblIdx = stencil[celli][i];
                cellCentre.append
                (
                    exchangeFields.getValue(mesh_.C(), mapCC, gblIdx)
                );
                alphaValues.append
                (
                    exchangeFields.getValue(alpha, mapAlpha, gblIdx)
                );
            }
            cellCentre -= mesh_.C()[celli];
            interfaceNormal_[i] = lsGrad.grad(cellCentre, alphaValues);
        }

    }
}

void Foam::plicRDF::gradSurf(const volScalarField& phi,const geometricVoF& surf)
{
    addProfilingInFunction(geometricVoF);
    leastSquareGrad<scalar> lsGrad("polyDegree1", mesh_.geometricD());
    zoneDistribute& exchangeFields = zoneDistribute::New(mesh_);

    exchangeFields.setUpCommforZone(surf.interfaceCell, false);

    Map<vector> mapCC
    (
        exchangeFields.getDatafromOtherProc(surf.interfaceCell, mesh_.C())
    );
    Map<scalar> mapPhi
    (
        exchangeFields.getDatafromOtherProc(surf.interfaceCell, phi)
    );

    DynamicField<vector> cellCentre(100);
    DynamicField<scalar> phiValues(100);

    const labelListList& stencil = exchangeFields.getStencil();

    forAll(surf.interfaceLabels, i)
    {
        const label celli = surf.interfaceLabels[i];

        cellCentre.clear();
        phiValues.clear();

        for (const label gblIdx : stencil[celli])
        {
            cellCentre.append
            (
                exchangeFields.getValue(mesh_.C(), mapCC, gblIdx)
            );
            phiValues.append
            (
                exchangeFields.getValue(phi, mapPhi, gblIdx)
            );
        }

        cellCentre -= mesh_.C()[celli];
        interfaceNormal_[i] = lsGrad.grad(cellCentre, phiValues);
    }
}


void Foam::plicRDF::setInitNormals(bool interpolate,geometricVoF& surf)
{
    addProfilingInFunction(geometricVoF);
    zoneDistribute& exchangeFields = zoneDistribute::New(mesh_);


    const volScalarField& alpha = surf.alpha();

    forAll(alpha, celli)
    {
        if (sIterPLIC_.isASurfaceCell(alpha[celli]))
        {
            surf.interfaceCell[celli] = true; // is set to false earlier
            surf.interfaceLabels.append(celli);
        }
    }
    interfaceNormal_.setSize(surf.interfaceLabels.size());

    RDF_.markCellsNearSurf(surf.interfaceCell, 1);
    const boolList& nextToInterface = RDF_.nextToInterface();
    exchangeFields.updateStencil(nextToInterface);

    if (interpolate)
    {
        interpolateNormal(surf);
    }
    else
    {
        gradSurf(alpha,surf);
    }
}


void Foam::plicRDF::calcResidual
(
    List<normalRes>& normalResidual,
    const geometricVoF& surf
)
{
    addProfilingInFunction(geometricVoF);
    zoneDistribute& exchangeFields = zoneDistribute::New(mesh_);
    exchangeFields.setUpCommforZone(surf.interfaceCell,false);

    Map<vector> mapNormal
    (
        exchangeFields.getDatafromOtherProc(surf.interfaceCell, surf.normal)
    );

    const labelListList& stencil = exchangeFields.getStencil();

    forAll(surf.interfaceLabels, i)
    {
        const label celli = surf.interfaceLabels[i];
        if (mag(surf.normal[celli]) == 0 || mag(interfaceNormal_[i]) == 0)
        {
            normalResidual[i].celli = celli;
            normalResidual[i].normalResidual = 0;
            normalResidual[i].avgAngle = 0;
            continue;
        }

        scalar avgDiffNormal = 0;
        scalar maxDiffNormal = GREAT;
        scalar weight= 0;
        const vector cellNormal = surf.normal[celli]/mag(surf.normal[celli]);
        forAll(stencil[celli],j)
        {
            const label gblIdx = stencil[celli][j];
            vector normal =
                exchangeFields.getValue(surf.normal, mapNormal, gblIdx);

            if (mag(normal) != 0 && j != 0)
            {
                vector n = normal/mag(normal);
                scalar cosAngle = max(min((cellNormal & n), 1), -1);
                avgDiffNormal += acos(cosAngle) * mag(normal);
                weight += mag(normal);
                if (cosAngle < maxDiffNormal)
                {
                    maxDiffNormal = cosAngle;
                }
            }
        }

        if (weight != 0)
        {
            avgDiffNormal /= weight;
        }
        else
        {
            avgDiffNormal = 0;
        }

        vector newCellNormal = normalised(interfaceNormal_[i]);

        scalar normalRes = (1 - (cellNormal & newCellNormal));
        normalResidual[i].celli = celli;
        normalResidual[i].normalResidual = normalRes;
        normalResidual[i].avgAngle = avgDiffNormal;
    }
}

void Foam::plicRDF::centreAndNormalBC(geometricVoF& surf)
{
    addProfilingInFunction(geometricVoF);
    scalar convertToRad = Foam::constant::mathematical::pi/180.0;

    // check if face is cut
    cutFacePLIC cutFace(mesh_);

    const volScalarField::Boundary& abf = surf.alpha().boundaryField();
    volVectorField::Boundary& cbf = surf.centre.boundaryFieldRef();
    volVectorField::Boundary& nbf = surf.normal.boundaryFieldRef();

    const fvBoundaryMesh& boundary = mesh_.boundary();

    // we need a surfaceVectorField to compute theta
    surfaceVectorField normalf(fvc::interpolate(surf.normal));

    forAll(boundary, patchi)
    {
        if (isA<alphaContactAngleTwoPhaseFvPatchScalarField>(abf[patchi]))
        {
            forAll(normalf.boundaryFieldRef()[patchi],i)
            {
                const label celli = boundary[patchi].faceCells()[i];
                vector n = surf.normal[celli];
                if(mag(n) != 0)
                {
                    n /= mag(n);
                    normalf.boundaryFieldRef()[patchi][i] = n;
                }
            }
        }
    }

    forAll(boundary, patchi)
    {
        if (isA<alphaContactAngleTwoPhaseFvPatchScalarField>(abf[patchi]))
        {
            alphaContactAngleTwoPhaseFvPatchScalarField& acap =
                const_cast<alphaContactAngleTwoPhaseFvPatchScalarField&>
                (
                    refCast<const alphaContactAngleTwoPhaseFvPatchScalarField>
                    (
                        abf[patchi]
                    )
                );

            fvsPatchVectorField& nHatp = normalf.boundaryFieldRef()[patchi];
            const scalarField theta
            (
                convertToRad*acap.theta(U_.boundaryField()[patchi], nHatp)
            );

            const vectorField nf
            (
                boundary[patchi].nf()
            );

            // Reset nHatp to correspond to the contact angle
            forAll(nbf[patchi],i)
            {
                const label celli = boundary[patchi].faceCells()[i];
                const label faceI = boundary[patchi].start() + i;
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

                    if(cutStatus == 0)
                    {
                        //const point cutEdgeCentre = average(cutFace.surfacePoints());

                        // project Normal on the face
                        vector projN = (tensor::I - nf[i]*nf[i]) & n;

                        // normalise
                        projN /= mag(projN) + SMALL;

                        vector nTheta = sin(theta[i])*nf[i] - cos(theta[i])*projN;
                        vector nHat =  cos(theta[i])*nf[i] + sin(theta[i])*projN;

                        cbf[patchi][i] = surf.centre[celli] + 2*nTheta/boundary[patchi].deltaCoeffs()[i]; // should point outside of the domain
                        nbf[patchi][i] = nHat*mag(surf.normal[celli]);

                    }

                }
                else
                {
                    cbf[patchi][i] = vector::zero;
                    nbf[patchi][i] = vector::zero;
                }
            }

            acap.evaluate();
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::plicRDF::plicRDF
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),

    interfaceNormal_(0.2*mesh_.nCells()),

    isoFaceTol_(dict.getOrDefault<scalar>("isoFaceTol", 1e-8)),
    surfCellTol_(dict.getOrDefault<scalar>("surfCellTol", 1e-8)),
    tol_(dict.getOrDefault<scalar>("tol" , 1e-6)),
    relTol_(dict.getOrDefault<scalar>("relTol" , 0.1)),
    iteration_(dict.getOrDefault<label>("iterations" , 5)),
    interpolateNormal_(dict.getOrDefault<bool>("interpolateNormal", true)),
    RDF_(reconstructedDistanceFunction::New(mesh)),
    sIterPLIC_(mesh_,surfCellTol_),
    U_(mesh_.lookupObject<volVectorField>(dict.getOrDefault<word>("UName", "U")))
{

}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


void Foam::plicRDF::update(geometricVoF& newSurf,geometricVoF& oldSurf,Foam::timeState state)
{
    addProfilingInFunction(geometricVoF);
    if (state == timeState::oldState)
    {
        reconstruct(oldSurf);
    }
    else
    {
        reconstruct(newSurf);
    }

}


void Foam::plicRDF::reconstruct(geometricVoF& surf)
{
    addProfilingInFunction(geometricVoF);
    zoneDistribute& exchangeFields = zoneDistribute::New(mesh_);

    const volScalarField& alpha = surf.alpha();

    if (mesh_.topoChanging())
    {
        // Introduced resizing to cope with changing meshes
        if (surf.interfaceCell.size() != mesh_.nCells())
        {
            surf.interfaceCell.resize(mesh_.nCells());
        }
    }
    surf.interfaceCell = false;
    surf.interfaceLabels.clear();

    // Sets surf.interfaceCell and interfaceNormal
    setInitNormals(interpolateNormal_,surf);

    surf.centre = dimensionedVector("centre", dimLength, Zero);
    surf.normal = dimensionedVector("normal", dimArea, Zero);

    // nextToInterface is update on setInitNormals
    const boolList& nextToInterface_ = RDF_.nextToInterface();

    bitSet tooCoarse(mesh_.nCells(),false);

    for (int iter=0; iter<iteration_; ++iter)
    {
        forAll(surf.interfaceLabels, i)
        {
            const label celli = surf.interfaceLabels[i];
            if (mag(interfaceNormal_[i]) == 0 || tooCoarse.test(celli))
            {
                continue;
            }
            sIterPLIC_.vofCutCell
            (
                celli,
                alpha[celli],
                isoFaceTol_,
                100,
                interfaceNormal_[i]
            );

            if (sIterPLIC_.cellStatus() == 0)
            {

                surf.normal[celli] = sIterPLIC_.surfaceArea();
                surf.centre[celli] = sIterPLIC_.surfaceCentre();
                if (mag(surf.normal[celli]) == 0)
                {
                    surf.normal[celli] = Zero;
                    surf.centre[celli] = Zero;
                }
            }
            else
            {
                surf.normal[celli] = Zero;
                surf.centre[celli] = Zero;
            }
        }

        surf.normal.correctBoundaryConditions();
        surf.centre.correctBoundaryConditions();
        List<normalRes> normalResidual(surf.interfaceLabels.size());

        surfaceVectorField::Boundary nHatb(mesh_.Sf().boundaryField());
        nHatb *= 1/(mesh_.magSf().boundaryField());

        {
            centreAndNormalBC(surf);
            RDF_.constructRDF
            (
                nextToInterface_,
                surf.centre,
                surf.normal,
                exchangeFields,
                false
            );
            // RDF_.updateContactAngle(alpha, U_, nHatb);
            gradSurf(RDF_,surf);
            calcResidual(normalResidual,surf);
        }

        label resCounter = 0;
        scalar avgRes = 0;
        scalar avgNormRes = 0;

        forAll(normalResidual,i)
        {

            const label celli = normalResidual[i].celli;
            const scalar normalRes= normalResidual[i].normalResidual;
            const scalar avgA = normalResidual[i].avgAngle;

            if (avgA > 0.26 && iter > 0) // 15 deg
            {
                tooCoarse.set(celli);
            }
            else
            {
                avgRes += normalRes;
                scalar normRes = 0;
                scalar discreteError = 0.01*sqr(avgA);
                if (discreteError != 0)
                {
                    normRes= normalRes/max(discreteError, tol_);
                }
                else
                {
                    normRes= normalRes/tol_;
                }
                avgNormRes += normRes;
                resCounter++;

            }
        }

        reduce(avgRes,sumOp<scalar>());
        reduce(avgNormRes,sumOp<scalar>());
        reduce(resCounter,sumOp<label>());

        if (resCounter == 0) // avoid division  by zero and leave loop
        {
            resCounter = 1;
            avgRes = 0;
            avgNormRes = 0;
        }


        if (iter == 0)
        {
            DebugInfo
                << "initial residual absolute = "
                << avgRes/resCounter << nl
                << "initial residual normalized = "
                << avgNormRes/resCounter << nl;
        }

        if
        (
            (
                (avgNormRes/resCounter < relTol_ || avgRes/resCounter < tol_)
             && (iter >= 1 )
            )
         || iter + 1  == iteration_
        )
        {
            DebugInfo
                << "iterations = " << iter << nl
                << "final residual absolute = "
                << avgRes/resCounter << nl
                << "final residual normalized = " << avgNormRes/resCounter
                << endl;

            break;
        }
    }
}



// ************************************************************************* //
