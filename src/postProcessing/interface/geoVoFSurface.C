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

\*---------------------------------------------------------------------------*/

#include "geoVoFSurface.H"
#include "geometricVoFMethod.H"
#include "cutCellPLIC.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(geoVoFSurface, 0);

}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::geoVoFSurface::geoVoFSurface
(
    const fvMesh& mesh
)
:
    MeshStorage(),
    mesh_(mesh)
{

    // the method are registered under the interfaceCaputuringMethod name
    geometricVoFMethod& geoVoF = mesh_.lookupObjectRef<geometricVoFMethod>(interfaceCapturingMethod::typeName);
    const geometricVoF& surf = geoVoF.surface(timeState::newState);


    cutCellPLIC cellCut(mesh_);

    const volVectorField& normal = surf.normal;
    const volVectorField& centre = surf.centre;

    const boolList& geoVoFSurfaceCells = surf.interfaceCell;

    DynamicList< List<point> > facePts;
    DynamicList< label > geoVoFSurfaceCellAdressing(0.1*mesh_.nCells());

    forAll(geoVoFSurfaceCells,cellI)
    {
        if(geoVoFSurfaceCells[cellI])
        {
            if(mag(normal[cellI]) != 0)
            {
                vector n = -normal[cellI]/mag(normal[cellI]);

                scalar cutVal = (centre[cellI]-mesh_.C()[cellI]) & n;

                cellCut.calcSubCell(cellI,cutVal,n);
                const auto& fPoints = cellCut.facePoints();
                if (fPoints.size() >= 3)
                {
                    facePts.append(fPoints);
                    geoVoFSurfaceCellAdressing.append(cellI);
                }
            }
        }

    }

    meshCells_.setSize(geoVoFSurfaceCellAdressing.size());

    forAll(meshCells_,i)
    {
        meshCells_[i] = geoVoFSurfaceCellAdressing[i];
    }

    // Transfer to mesh storage
    {
        faceList faces(facePts.size());


        label nPoints = 0;
        forAll(facePts,i)
        {
            face f(facePts[i].size());
            forAll(f,fi)
            {
                f[fi] = nPoints + fi;
            }
            faces[i] = f;

            nPoints += facePts[i].size();
        }
        pointField points(nPoints);

        nPoints = 0; // reuse
        forAll(facePts,i)
        {
            forAll(facePts[i],fi)
            {
                points[nPoints] = facePts[i][fi];
                nPoints++;
            }

        }

        // geoVoFSurface has no zones
        surfZoneList zones(0);

        MeshStorage updated(std::move(points), std::move(faces), surfZoneList());

        this->MeshStorage::transfer(updated);
    }
}


// ************************************************************************* //
