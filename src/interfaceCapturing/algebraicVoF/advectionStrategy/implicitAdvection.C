/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Modified code Copyright (C) 2022 henning Scheufler
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

#include "implicitAdvection.H"
#include "fvCFD.H"
#include "CMULES.H"

namespace Foam
{
    defineTypeNameAndDebug(implicitAdvection, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::implicitAdvection::implicitAdvection
(
    volScalarField& alpha1,
    surfaceScalarField& alphaPhi,
    const surfaceScalarField& phi,
    const volVectorField& U
)
:
    alpha1_(alpha1),
    alphaPhi_(alphaPhi),
    phi_(phi),
    U_(U),
    nAlphaCorr_(1)
{
    const dictionary& alphaControls = alpha1_.mesh().solverDict(alpha1_.name());
    nAlphaCorr_= alphaControls.get<label>("nAlphaCorr");
}



void Foam::implicitAdvection::advect
(
    algebraicVoF& algVoFNew,
    algebraicVoF& algVoFOld,
    const volScalarField::Internal& Sp,
    const volScalarField::Internal& Su
)
{
    word alphaScheme("div(phi,alpha)");

    const fvMesh& mesh = alpha1_.mesh();

    for (int aCorr=0; aCorr<nAlphaCorr_; aCorr++)
    {

        fvScalarMatrix alpha1Eqn
        (
            // (
            //     LTS_
            //     ? fv::localEulerDdtScheme<scalar>(mesh).fvmDdt(alpha1_)
            //     : fv::EulerDdtScheme<scalar>(mesh).fvmDdt(alpha1_)
            // )
            fv::EulerDdtScheme<scalar>(mesh).fvmDdt(alpha1_)
            + fvm::div(phi_,alpha1_,alphaScheme)

            ==
            Su + fvm::Sp(Sp , alpha1_)
        );


        alpha1Eqn.solve();
    }

    Info<< "Phase-1 volume fraction = "
        << alpha1_.weightedAverage(mesh.Vsc()).value()
        << "  Min(" << alpha1_.name() << ") = " << min(alpha1_).value()
        << "  Max(" << alpha1_.name() << ") = " << max(alpha1_).value()
        << endl;
}
// ************************************************************************* //
