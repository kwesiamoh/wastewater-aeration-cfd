/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM
   \\    /   O peration     | Version: v2412
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Application
    aerationBioFoam

Description
    Steady frozen-flow biological scalar transport solver for aeration-tank
    BOD5 and COD calculations using time-averaged hydrodynamic and DO fields.
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    Info<< "\nCreating frozen-flow reduced ASM1 fields\n" << endl;

    #include "createFields.H"


    // Advance from the stored CFD mean-field time (300)
    // to the biological-model output time (301).

    ++runTime;

    Info<< "\nReduced ASM1 output time = "
        << runTime.timeName() << nl << endl;


    // --------------------------------------------------------------------- //
    // Reduced aerobic ASM1 calculation
    // --------------------------------------------------------------------- //

    #include "ASM1Eqns.H"


    // --------------------------------------------------------------------- //
    // Positivity protection
    // --------------------------------------------------------------------- //

    SI.max(zeroConcentration);
    SS.max(zeroConcentration);
    XI.max(zeroConcentration);
    XS.max(zeroConcentration);
    XBH.max(zeroConcentration);
    XP.max(zeroConcentration);
    SO.max(zeroConcentration);


    SI.correctBoundaryConditions();
    SS.correctBoundaryConditions();
    XI.correctBoundaryConditions();
    XS.correctBoundaryConditions();
    XBH.correctBoundaryConditions();
    XP.correctBoundaryConditions();
    SO.correctBoundaryConditions();


    // --------------------------------------------------------------------- //
    // Final diagnostic fields
    // --------------------------------------------------------------------- //

    oxygenFactor =
        SO/(KOH + SO);

    substrateCOD =
        SS + XS;

    solubleCOD =
        SI + SS;

    mixedLiquorCOD =
        SI + SS + XI + XS + XBH + XP;


    oxygenFactor.correctBoundaryConditions();
    substrateCOD.correctBoundaryConditions();
    solubleCOD.correctBoundaryConditions();
    mixedLiquorCOD.correctBoundaryConditions();


    Info<< "\nWriting reduced ASM1 fields at time "
        << runTime.timeName() << nl << endl;

    runTime.write();


    Info<< "\nFinal reduced ASM1 ranges:" << nl

        << "    SI  = "
        << gMin(SI.primitiveField()) << " to "
        << gMax(SI.primitiveField()) << " kg COD/m3" << nl

        << "    SS  = "
        << gMin(SS.primitiveField()) << " to "
        << gMax(SS.primitiveField()) << " kg COD/m3" << nl

        << "    XI  = "
        << gMin(XI.primitiveField()) << " to "
        << gMax(XI.primitiveField()) << " kg COD/m3" << nl

        << "    XS  = "
        << gMin(XS.primitiveField()) << " to "
        << gMax(XS.primitiveField()) << " kg COD/m3" << nl

        << "    XBH = "
        << gMin(XBH.primitiveField()) << " to "
        << gMax(XBH.primitiveField()) << " kg COD/m3" << nl

        << "    XP  = "
        << gMin(XP.primitiveField()) << " to "
        << gMax(XP.primitiveField()) << " kg COD/m3" << nl

        << "    SO  = "
        << gMin(SO.primitiveField()) << " to "
        << gMax(SO.primitiveField()) << " kg O2/m3" << nl

        << "    substrateCOD = "
        << gMin(substrateCOD.primitiveField()) << " to "
        << gMax(substrateCOD.primitiveField()) << " kg COD/m3" << nl

        << "    solubleCOD = "
        << gMin(solubleCOD.primitiveField()) << " to "
        << gMax(solubleCOD.primitiveField()) << " kg COD/m3"
        << nl << endl;


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
