#include "fvCFD.H"
//#include "fvMesh.H"
#include "interfaceCapturingMethod.H"
#include "geometricVoFMethod.H"
// #include "isoAdvection.H"
// #include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
// #include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
// #include "turbulentTransportModel.H"
// #include "pimpleControl.H"
// #include "fvOptions.H"
// #include "CorrectPhi.H"
// #include "fvcSmooth.H"
// #include "surfaceForces.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{


    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    advector->advect();

    Info << "toc ICM" << interfaceCapturingMethod::dictionaryConstructorTablePtr_->sortedToc() << endl;
    Info << "toc geometricVoF" << geometricVoFMethod::dictionaryConstructorTablePtr_->sortedToc() << endl;

    mesh.lookupObject<interfaceCapturingMethod>("interfaceCapturingMethod");
    mesh.lookupObject<geometricVoF>("interfaceCapturingMethod");



    // #include "initCorrectPhi.H"
    // #include "createUfIfPresent.H"

    // const bool overwrite = args.found("overwrite");

    // turbulence->validate();

    // #include "CourantNo.H"
    // #include "setInitialDeltaT.H"

    // // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    // Info<< "\nStarting time loop\n" << endl;

    // while (runTime.run())
    // {
    //     #include "readDyMControls.H"
    //     #include "CourantNo.H"
    //     #include "alphaCourantNo.H"
    //     #include "setDeltaT.H"

    //     ++runTime;

    //     // helps with the initialization of AMR
    //     if(overwrite) 
    //     {
    //         runTime.setTime(runTime.value() - runTime.deltaTValue(), 1);
    //         runTime.writeAndEnd();
    //     }

    //     Info<< "Time = " << runTime.timeName() << nl << endl;

    //     // --- Pressure-velocity PIMPLE corrector loop
    //     while (pimple.loop())
    //     {
    //         if (pimple.firstIter() || moveMeshOuterCorrectors)
    //         {
    //             advector->surf().reconstruct();

    //             mesh.update();

    //             if (mesh.changing())
    //             {
    //                 // gets recompute by surfaces forces
    //                 // gh = (g & mesh.C()) - ghRef;
    //                 // ghf = (g & mesh.Cf()) - ghRef;
    //                 advector->surf().mapAlphaField();
    //                 alpha2 = 1.0 - alpha1;
    //                 alpha2.correctBoundaryConditions();
    //                 rho == alpha1*rho1 + alpha2*rho2;
    //                 rho.correctBoundaryConditions();
    //                 rho.oldTime() = rho;
    //                 alpha2.oldTime() = alpha2;

    //                 MRF.update();

    //                 if (correctPhi)
    //                 {
    //                     // Calculate absolute flux
    //                     // from the mapped surface velocity
    //                     phi = mesh.Sf() & Uf();

    //                     #include "correctPhi.H"

    //                     // Make the flux relative to the mesh motion
    //                     fvc::makeRelative(phi, U);

    //                     mixture.correct();
    //                 }

    //                 if (checkMeshCourantNo)
    //                 {
    //                     #include "meshCourantNo.H"
    //                 }
    //             }
    //         }

    //         if(overwrite)
    //         {
    //             continue;
    //         }

    //         #include "alphaControls.H"
    //         #include "alphaEqnSubCycle.H"

    //         mixture.correct();

    //         surfForces.correct();

    //         if (pimple.frozenFlow())
    //         {
    //             continue;
    //         }

    //         #include "UEqn.H"

    //         // --- Pressure corrector loop
    //         while (pimple.correct())
    //         {
    //             #include "pEqn.H"
    //         }

    //         if (pimple.turbCorr())
    //         {
    //             turbulence->correct();
    //         }
    //     }

    //     runTime.write();

    //     runTime.printExecutionTime(Info);
    // }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //