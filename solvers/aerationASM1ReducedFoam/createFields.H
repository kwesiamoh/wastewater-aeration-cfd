// ************************************************************************* //
// Frozen hydrodynamics + reduced aerobic ASM1 fields
// ************************************************************************* //


// ------------------------------------------------------------------------- //
// Frozen CFD fields
// ------------------------------------------------------------------------- //

Info<< "Reading frozen alphaUwaterMean\n" << endl;

volVectorField alphaUwaterMean
(
    IOobject
    (
        "alphaUwaterMean",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);


Info<< "Reading frozen alpha.waterMean\n" << endl;

volScalarField alphaWaterMean
(
    IOobject
    (
        "alpha.waterMean",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);


Info<< "Reading frozen nut.waterMean\n" << endl;

volScalarField nutWaterMean
(
    IOobject
    (
        "nut.waterMean",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);


// DOMean is retained as a reference field from the completed
// aeration-only CFD simulation.  It is NOT the DO field solved by ASM1.

Info<< "Reading frozen DOMean reference field\n" << endl;

volScalarField DOMean
(
    IOobject
    (
        "DOMean",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);


Info<< "Reading frozen DOsatLocalMean\n" << endl;

volScalarField DOsatLocalMean
(
    IOobject
    (
        "DOsatLocalMean",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


Info<< "Reading frozen oxygenTransferCoeffPostMean\n" << endl;

volScalarField oxygenTransferCoeffMean
(
    IOobject
    (
        "oxygenTransferCoeffPostMean",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);


// ------------------------------------------------------------------------- //
// Reduced ASM1 properties
//
// All concentration-based ASM1 parameters are expressed consistently
// with the OpenFOAM concentration fields in kg/m3.
//
// Numerically:
//     1 kg/m3 = 1000 mg/L
//
// Rate constants are supplied in SI units (1/s) in asm1Properties.
// ------------------------------------------------------------------------- //

Info<< "Reading constant/asm1Properties\n" << endl;

IOdictionary asm1Properties
(
    IOobject
    (
        "asm1Properties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);


// Effective molecular diffusivity for soluble ASM1 components

dimensionedScalar DmSoluble
(
    "DmSoluble",
    dimViscosity,
    asm1Properties
);


// Effective molecular diffusivity for particulate ASM1 components

dimensionedScalar DmParticulate
(
    "DmParticulate",
    dimViscosity,
    asm1Properties
);


// Turbulent Schmidt number

dimensionedScalar ScT
(
    "ScT",
    dimless,
    asm1Properties
);


// Heterotrophic maximum specific growth rate

dimensionedScalar muH
(
    "muH",
    dimless/dimTime,
    asm1Properties
);


// Heterotrophic decay coefficient

dimensionedScalar bH
(
    "bH",
    dimless/dimTime,
    asm1Properties
);


// Hydrolysis rate coefficient

dimensionedScalar kh
(
    "kh",
    dimless/dimTime,
    asm1Properties
);


// Half-saturation coefficient for readily biodegradable substrate

dimensionedScalar KS
(
    "KS",
    dimDensity,
    asm1Properties
);


// Oxygen half-saturation coefficient for heterotrophs

dimensionedScalar KOH
(
    "KOH",
    dimDensity,
    asm1Properties
);


// Hydrolysis saturation coefficient

dimensionedScalar KX
(
    "KX",
    dimless,
    asm1Properties
);


// Heterotrophic yield coefficient

dimensionedScalar YH
(
    "YH",
    dimless,
    asm1Properties
);


// Fraction of biomass decay forming particulate inert products

dimensionedScalar fP
(
    "fP",
    dimless,
    asm1Properties
);


// ------------------------------------------------------------------------- //
// Numerical concentration regularisation
// ------------------------------------------------------------------------- //

dimensionedScalar zeroConcentration
(
    "zeroConcentration",
    dimDensity,
    0.0
);


dimensionedScalar smallConcentration
(
    "smallConcentration",
    dimDensity,
    SMALL
);


// ------------------------------------------------------------------------- //
// Liquid-fraction regularisation
//
// alphaBio is used only to keep the transport matrix well-conditioned
// in cells where the time-averaged water fraction approaches zero.
//
// Biological interpretation/visualisation should still be restricted
// to the liquid region using alpha.waterMean.
// ------------------------------------------------------------------------- //

dimensionedScalar alphaWaterMin
(
    "alphaWaterMin",
    dimless,
    1e-6
);


volScalarField alphaBio
(
    IOobject
    (
        "alphaBio",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    max(alphaWaterMean, alphaWaterMin)
);


// ------------------------------------------------------------------------- //
// Tank patch IDs
// ------------------------------------------------------------------------- //

const label inletPatch =
    mesh.boundaryMesh().findPatchID("inlet");

const label outletPatch =
    mesh.boundaryMesh().findPatchID("outlet");

const label atmospherePatch =
    mesh.boundaryMesh().findPatchID("atmosphere");

const label sideWallsPatch =
    mesh.boundaryMesh().findPatchID("sideWalls");

const label floorPatch =
    mesh.boundaryMesh().findPatchID("floor");

const label diffusersPatch =
    mesh.boundaryMesh().findPatchID("diffusers");


if
(
    inletPatch < 0
 || outletPatch < 0
 || atmospherePatch < 0
 || sideWallsPatch < 0
 || floorPatch < 0
 || diffusersPatch < 0
)
{
    FatalErrorInFunction
        << "One or more required tank boundary patches were not found."
        << nl
        << "Required patches: inlet outlet atmosphere "
        << "sideWalls floor diffusers"
        << exit(FatalError);
}


// ------------------------------------------------------------------------- //
// Bubble-mediated oxygen-transfer mask
//
// The original dispersed-bubble oxygen-transfer model gives an
// unphysical interfacial-area contribution in cells directly adjacent
// to the resolved atmosphere/free-surface patch.
//
// Probe results showed:
//     alpha.waterMean ~ 0.652 at z = 5.0775 m
//     oxygen-transfer coefficient ~ 0.195 1/s
//
// whereas one layer below:
//     alpha.waterMean ~ 0.991
//     oxygen-transfer coefficient ~ 0.00487 1/s.
//
// Therefore only cells directly adjacent to "atmosphere" are excluded
// from diffuser-bubble oxygen transfer.
//
// Direct free-surface reaeration is not modelled here.
// ------------------------------------------------------------------------- //

Info<< "Constructing bubble oxygen-transfer mask\n" << endl;


volScalarField bubbleOTMask
(
    IOobject
    (
        "bubbleOTMask",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "one",
        dimless,
        1.0
    )
);


const labelUList& atmosphereCells =
    mesh.boundary()[atmospherePatch].faceCells();


forAll(atmosphereCells, facei)
{
    bubbleOTMask[atmosphereCells[facei]] = 0.0;
}


bubbleOTMask.boundaryFieldRef()[atmospherePatch] = scalar(0);
bubbleOTMask.correctBoundaryConditions();

// Make mask boundary values consistent with neighbouring internal cells.
// The atmosphere remains explicitly excluded.

forAll(bubbleOTMask.boundaryField(), patchi)
{
    if (patchi == atmospherePatch)
    {
        bubbleOTMask.boundaryFieldRef()[patchi] = scalar(0);
    }
    else
    {
        bubbleOTMask.boundaryFieldRef()[patchi] =
            bubbleOTMask.boundaryField()[patchi].patchInternalField();
    }
}


// Corrected bubble-mediated oxygen-transfer coefficient

volScalarField oxygenTransferCoeffBio
(
    IOobject
    (
        "oxygenTransferCoeffBio",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    oxygenTransferCoeffMean*bubbleOTMask
);


oxygenTransferCoeffBio.boundaryFieldRef()[atmospherePatch] = scalar(0);
oxygenTransferCoeffBio.correctBoundaryConditions();


// Make derived-field boundary values consistent with neighbouring
// internal cells. These patch values are diagnostics only; the
// oxygen-transfer source is volumetric.

forAll(oxygenTransferCoeffBio.boundaryField(), patchi)
{
    if (patchi == atmospherePatch)
    {
        oxygenTransferCoeffBio.boundaryFieldRef()[patchi] = scalar(0);
    }
    else
    {
        oxygenTransferCoeffBio.boundaryFieldRef()[patchi] =
            oxygenTransferCoeffBio.boundaryField()[patchi].patchInternalField();
    }
}


// ------------------------------------------------------------------------- //
// Conservative frozen water-phase flux
//
// Start from:
//
//     alphaUwaterMean = mean(alpha.water * U.water)
//
// A direct interpolation to the faces is not exactly conservative and,
// in the validation run, produced an unphysical atmosphere flux.
//
// The procedure is therefore:
//
// 1. reconstruct a raw face flux,
// 2. enforce the physical biological-transport boundary fluxes,
// 3. rescale the outlet profile to exactly match inlet flow,
// 4. perform a Poisson projection to obtain a divergence-free field.
// ------------------------------------------------------------------------- //

Info<< "Constructing raw frozen water-phase flux\n" << endl;


surfaceScalarField alphaPhiRaw
(
    IOobject
    (
        "alphaPhiRaw",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fvc::flux(alphaUwaterMean)
);


// No wastewater advection through impermeable/non-hydraulic boundaries

alphaPhiRaw.boundaryFieldRef()[atmospherePatch] = scalar(0);
alphaPhiRaw.boundaryFieldRef()[sideWallsPatch]   = scalar(0);
alphaPhiRaw.boundaryFieldRef()[floorPatch]       = scalar(0);
alphaPhiRaw.boundaryFieldRef()[diffusersPatch]   = scalar(0);


// Preserve inlet flow and balance outlet exactly

const scalar inletFlux =
    gSum(alphaPhiRaw.boundaryField()[inletPatch]);

const scalar outletFluxRaw =
    gSum(alphaPhiRaw.boundaryField()[outletPatch]);

const scalar targetOutletFlux =
    -inletFlux;


if (mag(outletFluxRaw) < SMALL)
{
    FatalErrorInFunction
        << "Raw outlet flux is zero; cannot balance frozen flow."
        << exit(FatalError);
}


alphaPhiRaw.boundaryFieldRef()[outletPatch] *=
    targetOutletFlux/outletFluxRaw;


Info<< "Conditioned frozen-flow boundary fluxes:" << nl
    << "    inlet      = "
    << gSum(alphaPhiRaw.boundaryField()[inletPatch]) << nl

    << "    outlet     = "
    << gSum(alphaPhiRaw.boundaryField()[outletPatch]) << nl

    << "    atmosphere = "
    << gSum(alphaPhiRaw.boundaryField()[atmospherePatch]) << nl

    << "    sideWalls  = "
    << gSum(alphaPhiRaw.boundaryField()[sideWallsPatch]) << nl

    << "    floor      = "
    << gSum(alphaPhiRaw.boundaryField()[floorPatch]) << nl

    << "    diffusers  = "
    << gSum(alphaPhiRaw.boundaryField()[diffusersPatch])
    << nl << endl;


// ------------------------------------------------------------------------- //
// Conservative flux projection
//
// Solve:
//
//     laplacian(psiFlux) = div(alphaPhiRaw)
//
// and construct:
//
//     alphaPhiBio = alphaPhiRaw - flux(psiFlux)
//
// giving:
//
//     div(alphaPhiBio) ~= 0
//
// zeroGradient on psiFlux preserves the imposed boundary fluxes.
// ------------------------------------------------------------------------- //

volScalarField psiFlux
(
    IOobject
    (
        "psiFlux",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "zero",
        dimensionSet(0, 2, -1, 0, 0, 0, 0),
        0.0
    ),
    zeroGradientFvPatchScalarField::typeName
);


mesh.setFluxRequired(psiFlux.name());


fvScalarMatrix psiFluxEqn
(
    fvm::laplacian(psiFlux)
 ==
    fvc::div(alphaPhiRaw)
);


psiFluxEqn.setReference(0, 0.0);
psiFluxEqn.solve();


surfaceScalarField alphaPhiBio
(
    IOobject
    (
        "alphaPhiBio",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    alphaPhiRaw - psiFluxEqn.flux()
);


volScalarField divAlphaPhiBio
(
    IOobject
    (
        "divAlphaPhiBio",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fvc::div(alphaPhiBio)
);


Info<< "Projected alphaPhiBio divergence min/max = "
    << gMin(divAlphaPhiBio.primitiveField()) << " "
    << gMax(divAlphaPhiBio.primitiveField())
    << nl << endl;


// ------------------------------------------------------------------------- //
// Reduced aerobic ASM1 state variables
//
// SI   = soluble inert organic matter
// SS   = readily biodegradable substrate
// XI   = particulate inert organic matter
// XS   = slowly biodegradable substrate
// XBH  = active heterotrophic biomass
// XP   = particulate products arising from biomass decay
// SO   = dissolved oxygen
//
// Carbon variables are expressed as kg COD/m3.
// SO is expressed as kg O2/m3.
// ------------------------------------------------------------------------- //

Info<< "Reading ASM1 field SI\n" << endl;

volScalarField SI
(
    IOobject
    (
        "SI",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


Info<< "Reading ASM1 field SS\n" << endl;

volScalarField SS
(
    IOobject
    (
        "SS",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


Info<< "Reading ASM1 field XI\n" << endl;

volScalarField XI
(
    IOobject
    (
        "XI",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


Info<< "Reading ASM1 field XS\n" << endl;

volScalarField XS
(
    IOobject
    (
        "XS",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


Info<< "Reading ASM1 field XBH\n" << endl;

volScalarField XBH
(
    IOobject
    (
        "XBH",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


Info<< "Reading ASM1 field XP\n" << endl;

volScalarField XP
(
    IOobject
    (
        "XP",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


Info<< "Reading ASM1 dissolved oxygen field SO\n" << endl;

volScalarField SO
(
    IOobject
    (
        "SO",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


// ------------------------------------------------------------------------- //
// ASM1 aerobic oxygen limitation
//
//     fO = SO/(KOH + SO)
//
// This field is recalculated in the reaction loop as SO changes.
// The initial field is created here for output and diagnostics.
// ------------------------------------------------------------------------- //

volScalarField oxygenFactor
(
    IOobject
    (
        "oxygenFactor",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    SO/(KOH + max(SO, zeroConcentration))
);


// ------------------------------------------------------------------------- //
// Effective transport diffusivities
//
// Soluble components:
//     SI, SS, SO
//
// Particulate/suspended components:
//     XI, XS, XBH, XP
//
// Turbulent dispersion dominates molecular diffusion in this aeration tank.
// ------------------------------------------------------------------------- //

volScalarField DeffSoluble
(
    IOobject
    (
        "DeffSoluble",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    DmSoluble + nutWaterMean/ScT
);


volScalarField DeffParticulate
(
    IOobject
    (
        "DeffParticulate",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    DmParticulate + nutWaterMean/ScT
);


// ------------------------------------------------------------------------- //
// COD diagnostics
//
// mixedLiquorCOD includes biomass and particulate sludge.  It must not be
// interpreted as final effluent COD because the CFD domain does not include
// the secondary clarifier.
//
// substrateCOD represents the biodegradable substrate pool relevant to
// carbon removal.
//
// solubleCOD represents soluble wastewater COD in the reduced model.
// ------------------------------------------------------------------------- //

volScalarField mixedLiquorCOD
(
    IOobject
    (
        "mixedLiquorCOD",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    SI + SS + XI + XS + XBH + XP
);


volScalarField substrateCOD
(
    IOobject
    (
        "substrateCOD",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    SS + XS
);


volScalarField solubleCOD
(
    IOobject
    (
        "solubleCOD",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    SI + SS
);


// ************************************************************************* //