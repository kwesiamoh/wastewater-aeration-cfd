/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM
   \\    /   O peration     |
    \\  /    A nd           | Reduced biological reactive transport solver
     \\/     M anipulation  |
-------------------------------------------------------------------------------

Application
    bioReactiveTransportFoam

Description
    One-way coupled reactive scalar solver for activated-sludge
    aeration-tank simulations.

    Frozen CFD fields provide:
        - liquid-phase transport
        - turbulent scalar mixing
        - oxygen saturation
        - gas-liquid oxygen-transfer coefficient

    Transported biological fields:
        Ss  readily biodegradable soluble COD
        Xs  slowly biodegradable COD
        DO  dissolved oxygen

    Heterotrophic biomass XBH is prescribed as an operating condition.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "simpleControl.H"


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Reduced activated-sludge COD and dissolved-oxygen "
        "reactive transport solver."
    );

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    simpleControl simple(mesh);

    #include "createFields.H"


    Info<< "\nStarting biological reactive transport\n" << endl;

    #include "CourantNo.H"


    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;


        while (simple.correctNonOrthogonal())
        {
            // -------------------------------------------------------------
            // Slowly biodegradable COD hydrolysis
            // -------------------------------------------------------------

            DOPos = max(DO, zeroConc);
            XsPos = max(Xs, zeroConc);

            DOMonod =
                DOPos/(KOH + DOPos);

            XsRatio =
                XsPos/XBH;

            kHyd =
                kh
               *DOMonod
               /(KX + XsRatio);


            fvScalarMatrix XsEqn
            (
                fvm::ddt(alphaWaterMean, Xs)
              + fvm::div(phi, Xs)
			  - fvm::Sp(fvc::div(phi), Xs)
              - fvm::laplacian(alphaDeffXs, Xs)
              + fvm::Sp(alphaWaterMean*kHyd, Xs)
            );

            XsEqn.relax();
            XsEqn.solve();


            XsPos = max(Xs, zeroConc);

            rhoHyd =
                kHyd*XsPos;


            // -------------------------------------------------------------
            // Readily biodegradable COD utilisation
            // -------------------------------------------------------------

            SsPos = max(Ss, zeroConc);
            DOPos = max(DO, zeroConc);

            DOMonod =
                DOPos/(KOH + DOPos);


            kSsUptake =
                (muH/YH)
               *XBH
               *DOMonod
               /(Ks + SsPos);


            fvScalarMatrix SsEqn
            (
                fvm::ddt(alphaWaterMean, Ss)
              + fvm::div(phi, Ss)
			  - fvm::Sp(fvc::div(phi), Ss)
              - fvm::laplacian(alphaDeffSs, Ss)
              + fvm::Sp(alphaWaterMean*kSsUptake, Ss)
             ==
                alphaWaterMean*rhoHyd
            );

            SsEqn.relax();
            SsEqn.solve();


            // -------------------------------------------------------------
            // Dissolved oxygen
            // -------------------------------------------------------------

            SsPos = max(Ss, zeroConc);
            DOPos = max(DO, zeroConc);

            SsMonod =
                SsPos/(Ks + SsPos);


            kO2Uptake =
                oxygenPerGrowth
               *muH
               *XBH
               *SsMonod
               /(KOH + DOPos);


            fvScalarMatrix DOEqn
            (
                fvm::ddt(alphaWaterMean, DO)
              + fvm::div(phi, DO)
			  - fvm::Sp(fvc::div(phi), DO)
              - fvm::laplacian(alphaDeffDO, DO)
              + fvm::Sp(oxygenTransferCoeffPostMean, DO)
              + fvm::Sp(alphaWaterMean*kO2Uptake, DO)
             ==
                oxygenTransferCoeffPostMean*DOsatLocalMean
            );

            DOEqn.relax();
            DOEqn.solve();


            // -------------------------------------------------------------
            // Diagnostic biological rates
            // -------------------------------------------------------------

            SsPos = max(Ss, zeroConc);
            DOPos = max(DO, zeroConc);

            SsMonod =
                SsPos/(Ks + SsPos);

            DOMonod =
                DOPos/(KOH + DOPos);


            rhoH =
                muH
               *SsMonod
               *DOMonod
               *XBH;

            SsUptakeRate =
                rhoH/YH;

            oxygenUptakeRate =
                oxygenPerGrowth*rhoH;
        }


        Info<< "    Ss min/max = "
            << gMin(Ss) << " "
            << gMax(Ss) << nl;

        Info<< "    Xs min/max = "
            << gMin(Xs) << " "
            << gMax(Xs) << nl;

        Info<< "    DO min/max = "
            << gMin(DO) << " "
            << gMax(DO) << nl << endl;


        runTime.write();
    }


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
