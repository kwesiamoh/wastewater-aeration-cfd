// ************************************************************************* //
// Dissolved oxygen transport
// ************************************************************************* //

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


// Instantaneous water turbulent viscosity.
// This field is updated by the multiphase turbulence model.
const volScalarField& nutWater =
    mesh.lookupObject<volScalarField>("nut.water");


Info<< "Reading constant/DOProperties\n" << endl;

IOdictionary DOProperties
(
    IOobject
    (
        "DOProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);


// Molecular diffusivity of oxygen in water
dimensionedScalar DmDO
(
    "DmDO",
    dimViscosity,
    DOProperties
);


// Turbulent Schmidt number
dimensionedScalar ScT
(
    "ScT",
    dimless,
    DOProperties
);


// Oxygen saturation concentration at reference pressure
dimensionedScalar DOsatRef
(
    "DOsatRef",
    dimDensity,
    DOProperties
);


// Reference pressure for saturation correction
dimensionedScalar pRefDO
(
    "pRefDO",
    dimPressure,
    DOProperties
);


// Wastewater oxygen-transfer correction factor
dimensionedScalar alphaFactor
(
    "alphaFactor",
    dimless,
    DOProperties
);

dimensionedScalar betaFactor
(
    "betaFactor",
    dimless,
    DOProperties
);

// Effective bubble diameter
dimensionedScalar dBubble
(
    "dBubble",
    dimLength,
    DOProperties
);


dimensionedScalar alphaWaterMin
(
    "alphaWaterMin",
    dimless,
    1e-6
);


dimensionedScalar zeroConcentration
(
    "zeroConcentration",
    dimDensity,
    0
);
