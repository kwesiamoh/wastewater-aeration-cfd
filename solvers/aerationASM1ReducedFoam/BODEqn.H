// ************************************************************************* //
// Steady BOD5-equivalent transport
// ************************************************************************* //

Info<< "Solving BOD5 transport" << endl;

tmp<volScalarField> tDeffBOD
(
    DmBOD + nutWaterMean/ScT
);

fvScalarMatrix BOD5Eqn
(
    fvm::div(alphaPhiBio, BOD5)
  - fvm::laplacian(alphaBio*tDeffBOD(), BOD5)
  + fvm::Sp(alphaBio*kBOD*oxygenFactor, BOD5)
);

BOD5Eqn.solve();

Info<< "BOD5 equation solved" << nl << endl;
