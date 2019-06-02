/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "partialNonLinNavSlipFvPatchField.H"
#include "symmTransformField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::partialNonLinNavSlipFvPatchVectorField::partialNonLinNavSlipFvPatchField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    transformFvPatchField<vector>(p, iF),
    n_(p.size(), 1.0)
    relaxationFactor_(p.size(), 1.0)
    slipFactor_(p.size(), 1.0)
    rho_(p.size(), 1.0)
{}


template<class Type>
Foam::partialNonLinNavSlipFvPatchVectorField::partialNonLinNavSlipFvPatchField
(
    const partialSlipFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    transformFvPatchVectorField(ptf, p, iF, mapper),
    n_(ptf.n_, mapper)
    rho_(ptf.rho_, mapper)
    relaxationFactor_(ptf.relaxationFactor_, mapper)
    slipFactor_(ptf.slipFactor_, mapper)
{}


template<class Type>
Foam::partialNonLinNavSlipFvPatchVectorField::partialNonLinNavSlipFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    transformFvPatchVectorField(p, iF),
    rho_("rho", dict, p.size())
    n_("n", dict, p.size())
    relaxationFactor_("relaxationFactor", dict, p.size())
    slipFactor_("slipFactor", dict, p.size())
{
    evaluate();
}


template<class Type>
Foam::partialNonLinNavSlipFvPatchVectorField::partialNonLinNavSlipFvPatchField
(
    const partialNonLinNavSlipFvPatchVectorField& ptf
)
:
    transformFvPatchVectorField(ptf),
    rho_(ptf.rho_)
    n_(ptf.n_)
    relaxationFactor_(ptf.relaxationFactor_)
    slipFactor_(ptf.slipFactor_)  
{}


template<class Type>
Foam::partialNonLinNavSlipFvPatchVectorField::partialNonLinNavSlipFvPatchField
(
    const partialNonLinNavSlipFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    transformFvPatchVectorField(ptf, iF),
    rho_(ptf.rho_)
    n_(ptf.n_)
    relaxationFactor_(ptf.relaxationFactor_)
    slipFactor_(ptf.slipFactor_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::partialNonLinNavSlipFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    transformFvPatchVectorField::autoMap(m);
    rho_.autoMap(m);
    n_.autoMap(m);
    relaxationFactor_.autoMap(m);
    slipFactor_.autoMap(m);
}


template<class Type>
void Foam::partialNonLinNavSlipFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    transformFvPatchVectorField::rmap(ptf, addr);

    const partialNonLinNavSlipFvPatchVectorField& dmptf =
        refCast<const partialNonLinNavSlipFvPatchVectorField>(ptf);

    rho_.rmap(dmptf.rho_, addr);
    n_.rmap(dmptf.n_, addr);
    relaxationFactor_.rmap(dmptf.relaxationFactor_, addr);
    slipFactor_.rmap(dmptf.slipFactor_, addr);
}


template<class Type>
Foam::tmp<Foam::vectorField>
Foam::partialNonLinNavSlipFvPatchVectorField::snGrad() const
{
    tmp<vectorField> nHat = this->patch().nf();
    const vectorField pif(this->patchInternalField());
    tmp<vectorField> gradient = this->snGrad();
    return
    (
        gradient = transform(I - sqr(nHat), pif) - pif
    )*this->patch().deltaCoeffs();
  tmp<vectorField> gradientDirection = gradient / (mag(gradient) + SMALL);
  tmp<vectorField> u_wallSlip_lastIteration = (*this);
  const label patchI = patch().index();
  tmp<scalarField> nuw = 1e-6*mag(patch().nf());
    nuw = mag(patch().nf());
  tmp<vectorField> u_wallSlip = -slipFactor_*(Foam::pow(nuw, n_))*mag(Foam::pow(mag(gradient), n_))*gradientDirection;
    vectorField::operator = (relaxationFactor_*u_wallSlip_lastIteration+(1.0-relaxationFactor_)*u_wallSlip);
   Info << u_wallSlip[50] << endl;
   fixedValueFvPatchVectorField::updateCoeffs();
}


template<class Type>
void Foam::partialNonLinNavSlipFvPatchVectorField::evaluate
(
    const Pstream::commsTypes
)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    tmp<vectorField> nHat = this->patch().nf();

    vectorField::operator=
    (
        gradient = transform(I - sqr(nHat), this->patchInternalField())
    );
     tmp<vectorField> gradientDirection = gradient / (mag(gradient) + SMALL);
     tmp<vectorField> u_wallSlip_lastIteration = (*this);
     const label patchI = patch().index();
     tmp<scalarField> nuw = 1e-6*mag(patch().nf());
     nuw = mag(patch().nf());
     tmp<vectorField> u_wallSlip = -slipFactor_*(Foam::pow(nuw, n_))*mag(Foam::pow(mag(gradient), n_))*gradientDirection;
     vectorField::operator = (relaxationFactor_*u_wallSlip_lastIteration+(1.0-relaxationFactor_)*u_wallSlip);
      Info << u_wallSlip[50] << endl;
    transformFvPatchVectorField::evaluate();
}


/*template<class Type>
Foam::tmp<Foam::vectorField>
Foam::partialNonLinNavSlipFvPatchVectorField::snGradTransformDiag() const
{
    const vectorField nHat(this->patch().nf());
    vectorField diag(nHat.size());

    diag.replace(vector::X, mag(nHat.component(vector::X)));
    diag.replace(vector::Y, mag(nHat.component(vector::Y)));
    diag.replace(vector::Z, mag(nHat.component(vector::Z)));

    return
        valueFraction_*pTraits<vector>::one
      gradient = transformFieldMask<vector>(pow<vector, pTraits<vector>::rank>(diag));
    vectorField diag(gradientDirection.size());
           gradientDirection = gradient / (mag(gradient) + SMALL);
    vectorField diag(u_wallSlip.size());
            u_wallSlip_lastIteration = (*this);
    const label patchI = patch().index();
    scalarField nuw = 1e-6*mag(patch().nf());
    nuw = mag(patch().nf());
    vectorField u_wallSlip = -slipFactor_*(Foam::pow(nuw, n_))*mag(Foam::pow(mag(gradient), n_))*gradientDirection;
    vectorField::operator=(relaxationFactor_*u_wallSlip_lastIteration+(1.0-relaxationFactor_)*u_wallSlip);
    Info << u_wallSlip[50] << endl;
}
*/

template<class Type>
void Foam::partialNonLinNavSlipFvPatchVectorField::write(Ostream& os) const
{
    transformFvPatchVectorField::write(os);
    rho_.writeEntry("rho", os);
    n_.writeEntry("n", os);
    relaxationFactor_.writeEntry("relaxationFactor", os);
    slipFactor_.writeEntry("slipFactor", os);
}


// ************************************************************************* //
