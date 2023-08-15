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
    Foam::accelerationForceMethod

Description
    abstract base class of the acceleration model
    handles dynamic and static formulation of the gravity/acceleration source terms

SourceFiles
    accelerationForceMethod.C
    accelerationForceMethods.C

\*---------------------------------------------------------------------------*/

#ifndef accelerationForceMethod_H
#define accelerationForceMethod_H

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
#include "interfaceCapturingMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class accelerationForceMethod Declaration
\*---------------------------------------------------------------------------*/

class accelerationForceMethod
{
protected:



        //- ghf
        surfaceScalarField accf_; 

        //- gh
        volScalarField acc_;

        //- Re-calculate the acceleration
        virtual void calculateAcc();

private:


    // Private Member Functions

        //- Disallow copy construct
        accelerationForceMethod(const accelerationForceMethod&);

        //- Disallow default bitwise assignment
        void operator=(const accelerationForceMethod&);

public:

    static std::unique_ptr<HashTable<word>> compatibilityTable_;

    class registerModel
    {
        public:
            registerModel(word modelName,word compatibleICM)
            {
                accelerationForceMethod::compatibilityTable_->insert(modelName,compatibleICM);
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

                    forAllConstIters(*accelerationForceMethod::compatibilityTable_, iter)
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
    TypeName("accelerationForceMethod");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            accelerationForceMethod,
            dictionary,
            (
                interfaceCapturingMethod& ICM,
                const dictionary& dict
            ),
            (ICM,dict)
        );

    // Selectors

        //- Return a reference to the selected phaseChange model
        static autoPtr<accelerationForceMethod> New
        (
            interfaceCapturingMethod& ICM,
            const dictionary& dict
        );

    // Constructors

        //- Construct from dictionary
        accelerationForceMethod
        (
            interfaceCapturingMethod& ICM,
            const dictionary& dict
        );

    //- Destructor
    virtual ~accelerationForceMethod()
    {}

    // Member Functions

        // static void registerModel(word ModelName,word compatibleICM)
        // {
        //     compatibilityTable_->insert(ModelName,compatibleICM);
        // }

        const surfaceScalarField& accf() const
        {
            return accf_;
        }

        const volScalarField& acc() const
        {
            return acc_;
        }

        //- return the pressure field depending
        virtual volScalarField& pressure
        (
            volScalarField& p,
            volScalarField& p_rgh
        ) = 0;

        //- update static pressure field
        virtual void updatePressure
        (
            volScalarField& p,
            volScalarField& p_rgh,
            const volScalarField& rho
        ) = 0;

        //- update dynamic pressure field
        virtual void updateRefPressure
        (
            volScalarField& p,
            volScalarField& p_rgh,
            const volScalarField& rho
        ) = 0;

        //- update the acceleration fields
        void correct()
        {
            calculateAcc();
        }

        //- return the accleration forces
        virtual tmp<surfaceScalarField> accelerationForce() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
