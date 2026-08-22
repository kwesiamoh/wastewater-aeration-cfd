// ************************************************************************* //
// Reduced aerobic ASM1 equations
//
// Scope:
//     - prescribed heterotrophic biomass XBH
//     - aerobic heterotrophic growth
//     - aerobic hydrolysis of XS -> SS
//     - dissolved oxygen consumption by heterotrophic growth
//     - bubble-mediated oxygen transfer
//
// Not included:
//     - anoxic growth / denitrification
//     - nitrification
//     - dynamic biomass inventory / RAS
//     - biomass decay products
//
// This is therefore an ASM1-derived aerobic carbon subset,
// not a complete ASM1 plant model.
// ************************************************************************* //


Info<< "Solving reduced aerobic ASM1 transport\n" << endl;


// ------------------------------------------------------------------------- //
// Stoichiometric oxygen demand associated with aerobic heterotrophic growth
//
// ASM1 process 1:
//
//     rho1 = muH
//            SS/(KS + SS)
//            SO/(KOH + SO)
//            XBH
//
// SS consumption:
//
//     -(1/YH) rho1
//
// oxygen consumption:
//
//     -((1 - YH)/YH) rho1
// ------------------------------------------------------------------------- //

dimensionedScalar oxygenPerGrowth
(
    "oxygenPerGrowth",
    dimless,
    (1.0 - YH.value())/YH.value()
);


// ------------------------------------------------------------------------- //
// Diagnostic process-rate fields
// ------------------------------------------------------------------------- //

volScalarField rhoGrowth
(
    IOobject
    (
        "rhoGrowth",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "zero",
        dimDensity/dimTime,
        0.0
    )
);


volScalarField rhoHydrolysis
(
    IOobject
    (
        "rhoHydrolysis",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "zero",
        dimDensity/dimTime,
        0.0
    )
);


volScalarField oxygenUptakeH
(
    IOobject
    (
        "oxygenUptakeH",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "zero",
        dimDensity/dimTime,
        0.0
    )
);


// ------------------------------------------------------------------------- //
// Nonlinear coefficient fields
// ------------------------------------------------------------------------- //

volScalarField substrateFactor
(
    IOobject
    (
        "substrateFactor",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    SS/(KS + SS)
);


volScalarField hydrolysisCoeff
(
    IOobject
    (
        "hydrolysisCoeff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    kh
   *oxygenFactor
   *XBH
   /(KX*XBH + XS)
);


volScalarField growthSSCoeff
(
    IOobject
    (
        "growthSSCoeff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    (muH/YH)
   *oxygenFactor
   *XBH
   /(KS + SS)
);


volScalarField oxygenGrowthCoeff
(
    IOobject
    (
        "oxygenGrowthCoeff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    oxygenPerGrowth
   *muH
   *substrateFactor
   *XBH
   /(KOH + SO)
);


// ------------------------------------------------------------------------- //
// Picard iteration
//
// The reaction rates are nonlinear in SS, XS and SO.
// Each transport equation is linearised using the current iterate.
//
// Under-relaxation is applied to the equation matrices.
// ------------------------------------------------------------------------- //

const label nBioCorrectors = 500;
const scalar bioTolerance = 1e-6;
const scalar bioRelax = 1.0;

bool bioConverged = false;


for (label corr = 0; corr < nBioCorrectors; ++corr)
{
    scalarField SSold(SS.primitiveField());
    scalarField XSold(XS.primitiveField());
    scalarField SOold(SO.primitiveField());


    // --------------------------------------------------------------------- //
    // Update oxygen limitation
    // --------------------------------------------------------------------- //

    oxygenFactor =
        SO/(KOH + SO);

    oxygenFactor.correctBoundaryConditions();


    // --------------------------------------------------------------------- //
    // ASM1 aerobic hydrolysis
    //
    // Reduced aerobic form of process rho7:
    //
    //     rho7 =
    //         kh
    //         [SO/(KOH + SO)]
    //         [(XS/XBH)/(KX + XS/XBH)]
    //         XBH
    //
    // which can be written:
    //
    //     rho7 =
    //         kh
    //         fO
    //         XBH
    //         XS/(KX*XBH + XS)
    //
    // Hydrolysis converts XS -> SS on a COD basis.
    // --------------------------------------------------------------------- //

    hydrolysisCoeff =
        kh
       *oxygenFactor
       *XBH
       /(KX*XBH + XS);


    fvScalarMatrix XSEqn
    (
        fvm::div(alphaPhiBio, XS)
      - fvm::laplacian(alphaBio*DeffParticulate, XS)
      + fvm::Sp(alphaBio*hydrolysisCoeff, XS)
    );

    XSEqn.relax(bioRelax);
    XSEqn.solve();

    XS.max(zeroConcentration);
    XS.correctBoundaryConditions();


    // Recalculate the hydrolysis rate from the updated XS field

    rhoHydrolysis =
        kh
       *oxygenFactor
       *XBH
       *XS
       /(KX*XBH + XS);

    rhoHydrolysis.max
    (
        dimensionedScalar
        (
            "zeroRate",
            dimDensity/dimTime,
            0.0
        )
    );


    // --------------------------------------------------------------------- //
    // Readily biodegradable substrate SS
    //
    // ASM1 aerobic heterotrophic growth:
    //
    //     rho1 =
    //         muH
    //         SS/(KS + SS)
    //         SO/(KOH + SO)
    //         XBH
    //
    // SS stoichiometry:
    //
    //     -(1/YH) rho1
    //
    // Hydrolysis contributes +rho7 to SS.
    //
    // The growth sink is linearised as:
    //
    //     [muH/YH * fO * XBH/(KS + SS_old)] SS
    // --------------------------------------------------------------------- //

    growthSSCoeff =
        (muH/YH)
       *oxygenFactor
       *XBH
       /(KS + SS);


    fvScalarMatrix SSEqn
    (
        fvm::div(alphaPhiBio, SS)
      - fvm::laplacian(alphaBio*DeffSoluble, SS)
      + fvm::Sp(alphaBio*growthSSCoeff, SS)
     ==
        alphaBio*rhoHydrolysis
    );

    SSEqn.relax(bioRelax);
    SSEqn.solve();

    SS.max(zeroConcentration);
    SS.correctBoundaryConditions();


    // --------------------------------------------------------------------- //
    // Update growth kinetics using the new SS field
    // --------------------------------------------------------------------- //

    substrateFactor =
        SS/(KS + SS);

    substrateFactor.correctBoundaryConditions();


    rhoGrowth =
        muH
       *substrateFactor
       *oxygenFactor
       *XBH;


    rhoGrowth.max
    (
        dimensionedScalar
        (
            "zeroRate",
            dimDensity/dimTime,
            0.0
        )
    );


    // --------------------------------------------------------------------- //
    // Dissolved oxygen
    //
    // Bubble oxygen source:
    //
    //     KLa_eff (SOsat - SO)
    //
    // where oxygenTransferCoeffBio already contains the liquid-phase
    // weighting inherited from the original aeration CFD.
    //
    // Heterotrophic oxygen demand:
    //
    //     ((1-YH)/YH) rho1
    //
    // The Monod oxygen term is linearised in SO.
    // --------------------------------------------------------------------- //

    oxygenGrowthCoeff =
        oxygenPerGrowth
       *muH
       *substrateFactor
       *XBH
       /(KOH + SO);


    fvScalarMatrix SOEqn
    (
        fvm::div(alphaPhiBio, SO)
      - fvm::laplacian(alphaBio*DeffSoluble, SO)

      + fvm::Sp
        (
            oxygenTransferCoeffBio,
            SO
        )

      + fvm::Sp
        (
            alphaBio*oxygenGrowthCoeff,
            SO
        )

     ==
        oxygenTransferCoeffBio*DOsatLocalMean
    );

    SOEqn.relax(bioRelax);
    SOEqn.solve();

    SO.max(zeroConcentration);
    SO.correctBoundaryConditions();


    // --------------------------------------------------------------------- //
    // Update final nonlinear rates for this corrector
    // --------------------------------------------------------------------- //

    oxygenFactor =
        SO/(KOH + SO);

    oxygenFactor.correctBoundaryConditions();


    substrateFactor =
        SS/(KS + SS);


    rhoGrowth =
        muH
       *substrateFactor
       *oxygenFactor
       *XBH;


    rhoHydrolysis =
        kh
       *oxygenFactor
       *XBH
       *XS
       /(KX*XBH + XS);


    oxygenUptakeH =
        oxygenPerGrowth*rhoGrowth;


    // --------------------------------------------------------------------- //
    // Nonlinear convergence
    // --------------------------------------------------------------------- //

    const scalar dSS =
        gMax
        (
            mag
            (
                SS.primitiveField()
              - SSold
            )
        );


    const scalar dXS =
        gMax
        (
            mag
            (
                XS.primitiveField()
              - XSold
            )
        );


    const scalar dSO =
        gMax
        (
            mag
            (
                SO.primitiveField()
              - SOold
            )
        );


    const scalar SSscale =
        max(gMax(SSold), scalar(1e-12));

    const scalar XSscale =
        max(gMax(XSold), scalar(1e-12));

    const scalar SOscale =
        max(gMax(SOold), scalar(1e-12));


    const scalar relSS = dSS/SSscale;
    const scalar relXS = dXS/XSscale;
    const scalar relSO = dSO/SOscale;


    const scalar maxRelativeChange =
        max
        (
            relSS,
            max(relXS, relSO)
        );


    Info<< "ASM1 corrector " << corr + 1
        << ": relChange SS=" << relSS
        << " XS=" << relXS
        << " SO=" << relSO
        << " max=" << maxRelativeChange
        << nl;


    if (maxRelativeChange < bioTolerance)
    {
        Info<< "Reduced ASM1 nonlinear iteration converged after "
            << corr + 1 << " correctors."
            << nl << endl;

        bioConverged = true;
        break;
    }
}


if (!bioConverged)
{
    WarningInFunction
        << "Reduced ASM1 nonlinear iteration reached "
        << nBioCorrectors
        << " correctors without meeting tolerance "
        << bioTolerance
        << endl;
}


// ------------------------------------------------------------------------- //
// Final process rates and diagnostic concentration fields
// ------------------------------------------------------------------------- //

oxygenFactor =
    SO/(KOH + SO);


substrateFactor =
    SS/(KS + SS);


rhoGrowth =
    muH
   *substrateFactor
   *oxygenFactor
   *XBH;


rhoHydrolysis =
    kh
   *oxygenFactor
   *XBH
   *XS
   /(KX*XBH + XS);


oxygenUptakeH =
    oxygenPerGrowth*rhoGrowth;


// Wastewater carbon diagnostics

substrateCOD =
    SS + XS;

solubleCOD =
    SI + SS;

mixedLiquorCOD =
    SI + SS + XI + XS + XBH + XP;


substrateCOD.correctBoundaryConditions();
solubleCOD.correctBoundaryConditions();
mixedLiquorCOD.correctBoundaryConditions();


Info<< nl
    << "Reduced ASM1 final field ranges:" << nl

    << "    SS min/max = "
    << gMin(SS.primitiveField()) << " "
    << gMax(SS.primitiveField()) << nl

    << "    XS min/max = "
    << gMin(XS.primitiveField()) << " "
    << gMax(XS.primitiveField()) << nl

    << "    SO min/max = "
    << gMin(SO.primitiveField()) << " "
    << gMax(SO.primitiveField()) << nl

    << "    substrateCOD min/max = "
    << gMin(substrateCOD.primitiveField()) << " "
    << gMax(substrateCOD.primitiveField()) << nl

    << "    solubleCOD min/max = "
    << gMin(solubleCOD.primitiveField()) << " "
    << gMax(solubleCOD.primitiveField())
    << nl << endl;


// ************************************************************************* //
