README
===================


This is a somewhat "lighter" version of KORAL, a two-temperature, radiative GRMHD code written by Aleksander Sadowski

To download the development gpu branch:

.. code-block:: bash

  git clone -b gpu https://github.com/achael/koral_lite.git

To install on Adroit (without MPI) using the default test problem make sure the following are in your .bashrc

.. code-block:: bash
  module load openmpi/gcc/1.10.2/64
  module load hdf5/gcc/openmpi-1.10.2/1.10.0
  module load fftw/gcc/openmpi-1.10.2/3.3.4
  module load gsl/2.6

Then compile:
.. code-block:: bash

  ./mser.sh

To run the test problem:

.. code-block:: bash

  mkdir dumps
  ./ko
