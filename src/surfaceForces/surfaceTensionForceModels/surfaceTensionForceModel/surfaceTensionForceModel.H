/*---------------------------------------------------------------------------*\
            Copyright (c) 2017-2019, German Aerospace Center (DLR)
-------------------------------------------------------------------------------
License
    This file is part of the VoFLibrary source code library, which is an
	unofficial extension to OpenFOAM.
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
    Foam::surfaceTensionForceModel

Description

SourceFiles
    surfaceTensionForceModel.C
    newsurfaceTensionForceModel.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceTensionForceModel_H
#define surfaceTensionForceModel_H

#include "typeInfo.H"
#include "wordIOList.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "dimensionedScalar.H"
#include "autoPtr.H"
#include "Pair.H"

#include "IOdictionary.H"
#include "volFields.H"
#include "surfaceFields.H"

#include "surfaceTensionModel.H"
#include "fvc.H"
#include "interfaceCapturingMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class surfaceTensionForceModel Declaration
\*---------------------------------------------------------------------------*/

class surfaceTensionForceModel
{

protected:

        //- reference to volume fraction field
        const volScalarField& alpha1_;

        //- reference to face velocity
        const surfaceScalarField& phi_;

        //- reference to velcoity
        const volVectorField& U_;

        autoPtr<surfaceTensionModel> surfTenModel_;
        
        //- surfaceTensionForce field
        surfaceScalarField surfaceTensionForce_;

private:

    // Private Member Functions

        static Ostream& printSurfaceForceModels
        (
            Ostream& os,
            const wordList& cmptNames,
            const wordList& surfaceForceNames
        );

        //- Disallow copy construct
        surfaceTensionForceModel(const surfaceTensionForceModel&);

        //- Disallow default bitwise assignment
        void operator=(const surfaceTensionForceModel&);

public:

    static std::unique_ptr<HashTable<word>> compatibilityTable_;

    class registerModel
    {
        public:
            registerModel(word modelName,word compatibleICM)
            {
                surfaceTensionForceModel::compatibilityTable_->insert(modelName,compatibleICM);
            }

            static Ostream& printCompatibilityTable
            (
                Ostream& os,
                const wordList& cmptNames,
                const HashTable<word>& modelNames
            )
            {
                   // Build a table of constituent parts by split name into constituent parts
                    // - remove incompatible entries from the list
                    // - note: row-0 contains the names of constituent parts (ie, the header)

                    DynamicList<wordList> outputTbl;
                    outputTbl.resize(modelNames.size()+1);

                    label rowi = 0;

                    // Header
                    outputTbl[rowi] = cmptNames;
                    if (!outputTbl[rowi].empty())
                    {
                        ++rowi;
                    }

                    forAllConstIters(*surfaceTensionForceModel::compatibilityTable_, iter)
                    {
                        outputTbl[rowi] = wordList({iter.key(),iter()});
                        if (!outputTbl[rowi].empty())
                        {
                            ++rowi;
                        }
                    }

                    if (rowi > 1)
                    {
                        os  << "The selected accelerationForceModel does not excist or is not compatible with the interfaceCapturing method: " << nl
                            << "Please select one of the following cominations:" << nl
                            << endl;
                        outputTbl.resize(rowi);
                        Foam::printTable(outputTbl, os);
                    }

                    return os; 
            }

    };

    //- Runtime type information
    TypeName("surfaceTensionForceModel");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            surfaceTensionForceModel,
            dictionary,
            (
                const dictionary& dict,
                interfaceCapturingMethod& ICM
            ),
            (dict,ICM)
        );

    // Selectors

        //- Return a reference to the selected phaseChange model
        static autoPtr<surfaceTensionForceModel> New
        (
            const dictionary& dict,
            interfaceCapturingMethod& ICM
        );

    // Constructors

        //- Construct from dictionary
        surfaceTensionForceModel
        (
            const dictionary& dict,
            interfaceCapturingMethod& ICM
        );

    //- Destructor
    virtual ~surfaceTensionForceModel()
    {}

    // Member Functions


        virtual void correct() = 0;

        tmp<volScalarField> sigma() const
        {
            return surfTenModel_->sigma();
        }

        virtual const surfaceScalarField& surfaceTensionForce() const
        {
            return surfaceTensionForce_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
