Info<< "Reading biodegradable soluble COD field Ss\n" << endl;

volScalarField Ss
(
    IOobject
    (
        "Ss",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


Info<< "Reading slowly biodegradable COD field Xs\n" << endl;

volScalarField Xs
(
    IOobject
    (
        "Xs",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


Info<< "Reading dissolved oxygen field DO\n" << endl;

volScalarField DO
(
    IOobject
    (
        "DO",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


// Frozen CFD fields are stored in biological case time 0

Info<< "Reading frozen mean liquid transport field alphaUwaterMean\n" << endl;

volVectorField alphaUwaterMean
(
    IOobject
    (
        "alphaUwaterMean",
        "0",
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);


Info<< "Reading frozen mean liquid fraction alpha.waterMean\n" << endl;

volScalarField alphaWaterMean
(
    IOobject
    (
        "alpha.waterMean",
        "0",
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);


Info<< "Reading frozen turbulent viscosity nut.waterMean\n" << endl;

volScalarField nutWaterMean
(
    IOobject
    (
        "nut.waterMean",
        "0",
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);


Info<< "Reading frozen oxygen saturation field DOsatLocalMean\n" << endl;

volScalarField DOsatLocalMean
(
    IOobject
    (
        "DOsatLocalMean",
        "0",
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);


Info<< "Reading frozen oxygen-transfer coefficient field\n" << endl;

volScalarField oxygenTransferCoeffPostMean
(
    IOobject
    (
        "oxygenTransferCoeffPostMean",
        "0",
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);


Info<< "Reading biologicalProperties\n" << endl;

IOdictionary biologicalProperties
(
    IOobject
    (
        "biologicalProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);


// Molecular diffusivities and turbulent Schmidt number

dimensionedScalar DmSs
(
    "DmSs",
    dimViscosity,
    biologicalProperties
);

dimensionedScalar DmXs
(
    "DmXs",
    dimViscosity,
    biologicalProperties
);

dimensionedScalar DmDO
(
    "DmDO",
    dimViscosity,
    biologicalProperties
);

dimensionedScalar ScT
(
    "ScT",
    dimless,
    biologicalProperties
);


// Heterotrophic biological parameters

dimensionedScalar muH
(
    "muH",
    dimless/dimTime,
    biologicalProperties
);

dimensionedScalar Ks
(
    "Ks",
    dimDensity,
    biologicalProperties
);

dimensionedScalar KOH
(
    "KOH",
    dimDensity,
    biologicalProperties
);

dimensionedScalar YH
(
    "YH",
    dimless,
    biologicalProperties
);

dimensionedScalar XBH
(
    "XBH",
    dimDensity,
    biologicalProperties
);

dimensionedScalar kh
(
    "kh",
    dimless/dimTime,
    biologicalProperties
);

dimensionedScalar KX
(
    "KX",
    dimless,
    biologicalProperties
);


// Frozen liquid-phase superficial flux:
// alphaUwaterMean = mean(alpha.water * U.water)

Info<< "Reading frozen mean water-phase face flux alphaPhi.waterMean\n" << endl;

surfaceScalarField phi
(
    IOobject
    (
        "alphaPhi.waterMean",
        "0",
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);


// Effective scalar diffusivities

volScalarField alphaDeffSs
(
    IOobject
    (
        "alphaDeffSs",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    alphaWaterMean*(DmSs + nutWaterMean/ScT)
);

volScalarField alphaDeffXs
(
    IOobject
    (
        "alphaDeffXs",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    alphaWaterMean*(DmXs + nutWaterMean/ScT)
);

volScalarField alphaDeffDO
(
    IOobject
    (
        "alphaDeffDO",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    alphaWaterMean*(DmDO + nutWaterMean/ScT)
);


// Zero concentration used for bounded kinetic evaluation

const dimensionedScalar zeroConc
(
    "zeroConc",
    dimDensity,
    0.0
);


volScalarField SsPos
(
    IOobject
    (
        "SsPos",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    max(Ss, zeroConc)
);

volScalarField XsPos
(
    IOobject
    (
        "XsPos",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    max(Xs, zeroConc)
);

volScalarField DOPos
(
    IOobject
    (
        "DOPos",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    max(DO, zeroConc)
);


// Monod limitation fields

volScalarField SsMonod
(
    IOobject
    (
        "SsMonod",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    SsPos/(Ks + SsPos)
);

volScalarField DOMonod
(
    IOobject
    (
        "DOMonod",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    DOPos/(KOH + DOPos)
);

volScalarField XsRatio
(
    IOobject
    (
        "XsRatio",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    XsPos/XBH
);


// Semi-implicit kinetic coefficients [1/s]

volScalarField kHyd
(
    IOobject
    (
        "kHyd",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    kh*DOMonod/(KX + XsRatio)
);


volScalarField kSsUptake
(
    IOobject
    (
        "kSsUptake",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    (muH/YH)*XBH*DOMonod/(Ks + SsPos)
);


const scalar oxygenPerGrowth =
    (1.0 - YH.value())/YH.value();


volScalarField kO2Uptake
(
    IOobject
    (
        "kO2Uptake",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    oxygenPerGrowth
   *muH
   *XBH
   *SsMonod
   /(KOH + DOPos)
);


// Diagnostic biological reaction fields

volScalarField rhoHyd
(
    IOobject
    (
        "rhoHyd",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    kHyd*XsPos
);


volScalarField rhoH
(
    IOobject
    (
        "rhoH",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    muH*SsMonod*DOMonod*XBH
);


volScalarField SsUptakeRate
(
    IOobject
    (
        "SsUptakeRate",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    rhoH/YH
);


volScalarField oxygenUptakeRate
(
    IOobject
    (
        "oxygenUptakeRate",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    oxygenPerGrowth*rhoH
);
