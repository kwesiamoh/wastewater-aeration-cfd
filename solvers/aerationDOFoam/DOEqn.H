// ************************************************************************* //
// Real-time dissolved oxygen transport
//
// alpha1     = instantaneous alpha.air
// alpha2     = instantaneous alpha.water
// U1         = instantaneous U.air
// U2         = instantaneous U.water
// alphaPhi2  = instantaneous conservative water-phase flux
// ************************************************************************* //

Info<< "Solving dissolved oxygen transport" << endl;


// ------------------------------------------------------------------------- //
// Effective oxygen diffusivity in the liquid
// ------------------------------------------------------------------------- //

tmp<volScalarField> tDeffDO
(
    DmDO + nutWater/ScT
);


// ------------------------------------------------------------------------- //
// Instantaneous gas-liquid slip velocity
// ------------------------------------------------------------------------- //

volScalarField slipSpeed
(
    IOobject
    (
        "slipSpeed",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mag(U1 - U2)
);


// Protect division in cells approaching pure gas
volScalarField alphaWaterSafe
(
    IOobject
    (
        "alphaWaterSafe",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    max(alpha2, alphaWaterMin)
);


// ------------------------------------------------------------------------- //
// Interfacial area per LIQUID volume
//
//              6 alpha_g
//     a = ---------------------
//          d_b (1 - alpha_g)
//
// Since alpha2 = 1 - alpha1:
//
//              6 alpha1
//     a = ----------------
//          d_b alpha2
// ------------------------------------------------------------------------- //

volScalarField interfacialAreaLiquid
(
    IOobject
    (
        "interfacialAreaLiquid",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    scalar(6)*alpha1/(dBubble*alphaWaterSafe)
);


// ------------------------------------------------------------------------- //
// Higbie liquid-side mass-transfer coefficient
//
//                     D_L |Ur|
//     kL = 2 sqrt( ---------------- )
//                     pi d_b
// ------------------------------------------------------------------------- //

volScalarField kL
(
    IOobject
    (
        "kL",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    scalar(2)
   *sqrt
    (
        DmDO*slipSpeed
       /(scalar(3.14159265358979323846)*dBubble)
    )
);


// Wastewater-corrected local volumetric mass-transfer coefficient
volScalarField kLa
(
    IOobject
    (
        "kLa",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    alphaFactor*kL*interfacialAreaLiquid
);


// ------------------------------------------------------------------------- //
// Local saturation concentration.
//
// First-order pressure correction from the reference saturation condition.
// p is the instantaneous absolute pressure field from twoPhaseEulerFoam.
// ------------------------------------------------------------------------- //

volScalarField DOsatLocal
(
    IOobject
    (
        "DOsatLocal",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
   	betaFactor*DOsatRef*(p/pRefDO)
);


// Conservative equation is written per mixture volume.
// kLa above is defined per liquid volume, hence multiply by alpha2.
// ------------------------------------------------------------------------- //

volScalarField oxygenTransferCoeff
(
    IOobject
    (
        "oxygenTransferCoeff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    alpha2*kLa
);


// ------------------------------------------------------------------------- //
// Dissolved oxygen:
//
// d(alpha_w DO)/dt
// + div(alphaPhi_w DO)
// - div(alpha_w Deff grad(DO))
// = alpha_w kLa (DOsat - DO)
// ------------------------------------------------------------------------- //

fvScalarMatrix DOEqn
(
    fvm::ddt(alpha2, DO)
  + fvm::div(alphaPhi2, DO)
  - fvm::laplacian(alpha2*tDeffDO(), DO)
  + fvm::Sp(oxygenTransferCoeff, DO)
 ==
    oxygenTransferCoeff*DOsatLocal
);

DOEqn.relax();
DOEqn.solve();

DO.max(zeroConcentration);


Info<< "DO min/max = "
    << gMin(DO.primitiveField()) << " "
    << gMax(DO.primitiveField()) << nl;

Info<< "kLa min/max = "
    << gMin(kLa.primitiveField()) << " "
    << gMax(kLa.primitiveField()) << nl << endl;
