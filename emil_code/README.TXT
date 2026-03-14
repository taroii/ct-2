To install w/o GPU tv grad. descent routines:

(1) install packages: numpy, matplotlib, ipython

------

(2) General notes:
* Check the first few lines of all python programs, because
  that's where all the system dependent stuff is.

------

(3) compile fortran files:
Fortran subroutines are supposed to be compiled with gnu's
gfortran command which is a fortran 95 compiler (see help for
f2py on how to specify gnu95 if it is not your default).
If you only have g77, then you have to go through all the fortran
files and replace the floor(x) command with, say,
int(x + a large number)- the same large number.

If fortran 95 is your default here are the compilation commands:
f2py -c ray_proj3D.f -m ray_proj3D
f2py -c tomo3DfortranAux.f -m tomo3Dfort
f2py -c totalvar.f -m totalvar
f2py -c view_proj3D.f -m view_proj3D
f2py -c view_proj3D_tomo.f -m view_proj3D_tomo

------

(4) tests:
Not much yet.
You can type:
ipython -pylab tomo3D.py

that will have some projection tests, but more thorough testing
will follow


##########################################################################
ADDITIONAL INSTRUCTIONS FOR THE GPU ROUTINES WILL GO HERE
##########################################################################






EXAMPLES FOLDER:

These are TV reconstruction programs. Some use the
very flexible linogram3D data object which can model just about
any situation, and some use the sinogram3D data object which
allows arbitrary source trajectories and detector position and
orientations, but it is designed for a flat panel detector.

Running each program produces a folder with results in the data folder.

gridRandom10Views_eps_0p01.py:
run with:
ipython gridRandom10Views_eps_0p01.py
OR
python gridRandom10Views_eps_0p01.py

- Uses linogram data object
- 10 x-ray source views randomly selected, full sampling on a flat-panel
  detector.
- Phantom is an array of dots.
- data tolerance set to 0.01, the program aims to find the min. TV for this
  tolerance.


gridRandom10VRays_eps_0p01.py:
run with:
ipython gridRandom10VRays_eps_0p01.py
OR
python gridRandom10VRays_eps_0p01.py

- Uses linogram data object
- Projection ray selection is totally random, but there are the same
  number of total rays as previous example
- Phantom is an array of dots.
- data tolerance set to 0.01, the program aims to find the min. TV for this
  tolerance.


jawPhantomCCB.py:
run with:
ipython jawPhantomCCB.py
OR
python jawPhantomCCB.py

- Uses sinogram3D data object
- X-ray source scans 180 degrees on a circular orbit
- Phantom is the FORBILD jaw phantom
- data tolerance set to whatever the data error is at iteration 500;
  the program aims to find the min. TV for this tolerance.





