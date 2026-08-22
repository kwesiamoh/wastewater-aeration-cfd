// ************************************************************************* //
// COD transport
// ************************************************************************* //

Info<< "Solving biodegradable COD transport" << endl;

tmp<volScalarField> tDeffCOD
(
    DmCOD + nutWaterMean/ScT
);

fvScalarMatrix CODbioEqn
(
    fvm::div(alphaPhiBio, CODbio)
  - fvm::laplacian(alphaBio*tDeffCOD(), CODbio)
  + fvm::Sp(alphaBio*kCOD*oxygenFactor, CODbio)
);

CODbioEqn.solve();


Info<< "Solving inert COD transport" << endl;

fvScalarMatrix CODinertEqn
(
    fvm::div(alphaPhiBio, CODinert)
  - fvm::laplacian(alphaBio*tDeffCOD(), CODinert)
);

CODinertEqn.solve();

Info<< "COD equations solved" << nl << endl;
