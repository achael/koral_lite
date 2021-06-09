README
===================


This is a somewhat "lighter" version of KORAL, a two-temperature, radiative GRMHD code written by Aleksander Sadowski

To download the development gpu branch:

.. code-block:: bash

  git clone -b gpu https://github.com/achael/koral_lite.git

To install on **Adroit** make sure the following are in your .bashrc

.. code-block:: bash

  module load openmpi/gcc/1.10.2/64
  module load hdf5/gcc/openmpi-1.10.2/1.10.0
  module load fftw/gcc/openmpi-1.10.2/3.3.4
  module load gsl/2.6
  module load cudatoolkit/11.0
  
To install on **Della** make sure the following are in your .bashrc

.. code-block:: bash

  module load openmpi/gcc/4.1.0
  module load hdf5/gcc/openmpi-4.1.0/1.10.6
  module load fftw/gcc/openmpi-4.1.0/3.3.9
  module load cudatoolkit/11.1
  module load gsl/2.6

On Adroit/Della, we also define the following commands in .bashrc to set up interactive jobs 

.. code-block:: bash

  alias sinter='salloc -p gpu --reservation=hackathon --nodes=1 --ntasks=1 --mem=32G --time=04:00:00 /bin/bash'
  alias sinter='galloc -p gpu --reservation=hackathon --nodes=1 --ntasks=1 --mem=32G --time=04:00:00 --gres=gpu:1 /bin/bash'

and then start a gpu interactive session with: 
  
.. code-block:: bash

  source ~/.bashrc
  ginter

Then compile KORAL:

.. code-block:: bash

  ./mser.sh

Make sure you create a dumps directory, otherwise you will get a segmentation fault. The test problem will start at t=0 every time (reading dumps0000.dat if available). Empty this directory if you change the problem size. 

.. code-block:: bash

  mkdir dumps

To run the test problem on the cpu:

.. code-block:: bash

  ./ko
  

To run the test problem on the gpu:

.. code-block:: bash

  srun -n 1 ./ko_gpu

