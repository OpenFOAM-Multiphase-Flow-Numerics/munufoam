Installation
============

You can install numufoam by following these steps:

Clone the numufoam repository:

   .. code-block:: bash

      git clone https://github.com/OpenFOAM-Multiphase-Flow-Numerics/numufoam.git

Navigate to the numufoam directory:

   .. code-block:: bash

      cd numufoam

numufoam uses Cmake to build, thus the standard Cmake procedure should work, however, we recommend using one of the provided Cmake presets detailed below `below <Building with Cmake Presets>`_. From a build directory, you can execute:

   .. code-block:: bash

        ./Allwmake
        ./get-gmsh.sh # will install gmsh version 4.9.3 as gmshv493

