README
===================


This is a somewhat "lighter" version of KORAL, a two-temperature, radiative GRMHD code written by Aleksander Sadowski

To download the development gpu branch:
.. code-block:: bash

  git clone -b gpu https://github.com/achael/koral_lite.git

To install on Adroit (without MPI) using the default test problem:
.. code-block:: bash

  ./mser.sh

To run the test problem:
.. code-block:: bash

  mkdir dumps
  ./ko
