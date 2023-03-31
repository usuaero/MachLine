# MachLine Studies
Since MachLine is first and foremost a research code, we have found it desirable to include a collection of the studies we have conducted on MachLine. These studies cover a wide range of aspects from paneling sensitivity and grid convergence to accuracy compared to wind tunnel data.

Each study consists of a set of meshes, one or more input files, and one or more Python scripts used to run MachLine automatically, parse the data, and process results. For each, the Python scripts should be run from the MachLine root directory. Do not worry about copying the MachLine executable; the Python scripts are set up to automatically find the executable in the MachLine root directory.

The results of these studies are reproduced in various of the technical works documenting the development of MachLine.

## For Developers
The studies should be organized by the type of geometry/flow being organized. So, for example, rather than lumping a bunch of grid convergence studies for various geometries together, keep the sphere convergence study with the other sphere studies, and the double-wedge wing convergence study with the comparison against shock-expansion theory. This should help with keeping storage requirements down, since the mesh files are probably the largest out of all the files being stored.

At the top of any given study script, there should be a global variable called ```RERUN_MACHLINE```. This should allow the user to specify whether MachLine is to be rerun or the data should just be pulled from previously-generated data and report files.

To keep things consistent, try to refer to your mesh densities with the tags "ultra-coarse", "coarse", "medium", "fine", and "ultra-fine".

Common chunks of Python code can be put in scripts in the main studies/ directory and then imported into individual study scripts.