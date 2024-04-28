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
    Foam::geoVoFSurface

Description

Author
    Henning Scheufler, DLR, all rights reserved.

SourceFiles
    geoVoFSurface.C
    geoVoFSurfaceTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef geoVoFSurface_H
#define geoVoFSurface_H

#include "PackedBoolList.H"
#include "volFields.H"
#include "slicedVolFields.H"
#include "MeshedSurface.H"
#include "MeshedSurfacesFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;
class plane;
class treeBoundBox;
class triSurface;

/*---------------------------------------------------------------------------*\
                       Class geoVoFSurface Declaration
\*---------------------------------------------------------------------------*/

class geoVoFSurface
:
    public meshedSurface
{
    // Private typedefs for convenience
    typedef meshedSurface MeshStorage;

    // Private data

        //- Reference to mesh
        const fvMesh& mesh_;

        //- For every triangle the original cell in mesh
        labelList meshCells_;


public:


    //- Runtime type information
    TypeName("geoVoFSurface");


    // Constructors

        //- Construct from lookup
        geoVoFSurface(const fvMesh& mesh);


    // Member Functions

        //- For every face, the original cell in mesh
        const labelList& meshCells() const
        {
            return meshCells_;
        }

        //- Interpolates cCoords,pCoords. Uses the references to the original
        //  fields used to create the iso surface.
        template<class Type>
        tmp<Field<Type> > interpolate
        (
            const GeometricField<Type, fvPatchField, volMesh>& cCoords,
            const Field<Type>& pCoords
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "geoVoFSurfaceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
