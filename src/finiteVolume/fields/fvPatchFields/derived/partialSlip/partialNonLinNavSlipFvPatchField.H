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

Class
    Foam::partialNonLinNavSlipFvPatchField

Group
    grpWallBoundaryConditions grpGenericBoundaryConditions

Description
    This boundary condition provides a partial slip condition.  The amount of
    slip is controlled by a user-supplied field.

See also
    Foam::transformFvPatchField

SourceFiles
    partial NonLinNavSlipFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef partialNonLinNavSlipFvPatchField_H
#define partialNonLinNavSlipFvPatchField_H

#include "transformFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class partialSlipFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class partialNonLinNavSlipFvPatchField
:
    public transformFvPatchVectorField
{
    // Private data

        //-  value used for boundary condition
        scalarField rho_;
        scalarField n_;
        scalarField relaxationFactor_;
        scalarField slipFactor_;


public:

    //- Runtime type information
    TypeName("partialNonLinNavSlip");


    // Constructors

        //- Construct from patch and internal field
        partialNonLinNavSlipFvPatchField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        partialNonLinNavSlipFvPatchField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given partialSlipFvPatchField onto a new patch
        partialNonLinNavSlipFvPatchField
        (
            const partialSlipFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        partialNonLinNavSlipFvPatchField
        (
            const partialNonLinNavSlipFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new partialNonLinNavSlipFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        partialNonLinNavSlipFvPatchField
        (
            const partialSlipFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new partialSlipFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Attributes

            //- Return false: this patch field is not altered by assignment
            virtual bool assignable() const
            {
                return false;
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchVectorField&,
                const labelList&
            );


        // Return defining fields

            virtual scalarField& rho()
            {
                return rho_;
            }
            virtual scalarField& n()
            {
                return n_;
            }
            virtual scalarField& relaxationFactor()
            {
                return relaxationFactor_;
            }
            virtual scalarField& slipFactor()
            {
                return slipFactor_;
            }

            virtual const scalarField& rho() const
            {
                return rho_;
            }
            virtual const scalarField& n() const
            {
                return n_;
            }
            virtual const scalarField& relaxationFactor() const
            {
                return relaxationFactor_;
            }
            virtual const scalarField& slipFactor() const
            {
                return slipFactor_;
            }


        // Evaluation functions

            //- Return gradient at boundary
            virtual tmp<vectorField> snGrad() const;

            //- Evaluate the patch field
            virtual void evaluate
            (
                const Pstream::commsTypes commsType=Pstream::blocking
            );

            //- Return face-gradient transform diagonal
            virtual tmp<vectorField> snGradTransformDiag() const;


        //- Write
        virtual void write(Ostream&) const;


    // Member operators

        virtual void operator=(const UList<T>&) {}

        virtual void operator=(const fvPatchField<vector>&) {}
        virtual void operator+=(const fvPatchField<vector>&) {}
        virtual void operator-=(const fvPatchField<vector>&) {}
        virtual void operator*=(const fvPatchField<scalar>&) {}
        virtual void operator/=(const fvPatchField<scalar>&) {}

        virtual void operator+=(const Field<vector>&) {}
        virtual void operator-=(const Field<vector>&) {}

        virtual void operator*=(const Field<scalar>&) {}
        virtual void operator/=(const Field<scalar>&) {}

        virtual void operator=(const vector&) {}
        virtual void operator+=(const vector&) {}
        virtual void operator-=(const vector&) {}
        virtual void operator*=(const scalar) {}
        virtual void operator/=(const scalar) {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "partialNonLinNavSlipFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
