.. _hicTransform:

hicTransform
============

.. argparse::
   :ref: hicexplorer.hicTransform.parse_arguments
   :prog: hicTransform


Background
===========

hicTransform transforms a given input matrix to a matrix with the defined method applied on.

- obs_exp
- obs_exp_norm
- obs_exp_lieberman
- obs_exp_non_zero
- pearson
- covariance

All expected values are computed per genomic distances. 

Usage
-----

.. code:: bash

    $ hicTransform -m matrix.cool --method obs_exp -o obs_exp.cool

For all images data from Rao 2014 was used.

Observed / Expected
-------------------

All values, including non-zero values, are used to compute the expected values per genomic distance. 

.. math::

    m_{i,j} = \frac{}
    expected\ value\ distance\ j = \frac{ \sum diagonal(j) }{|diagonal(j)|}

.. image:: ../../images/obs_exp.png

Observed / Expected norm
------------------------

The expected matrix is computed in the same way Homer software computes it with the option '-norm' set,
to conserve this reference, HiCExplorer names this expected matrix computation 'norm'.

sum(diagonal(i-j)) * sum(row(j)) * sum(row(i)) / sum(matrix)

.. image:: ../../images/obs_exp_norm.png

Observed / Expected lieberman
-----------------------------

The expected matrix is computed in the way Lieberman-Aiden used it in the 2009 publication, it is quite similar 
to obs/exp matrix computation.

sum(diagonal(j) / (length of chromosome - j))

.. image:: ../../images/obs_exp_lieberman.png

Observed / Expected non zero
----------------------------

Only non-zero values are used to compute the expected values per genomic distance, i.e. only non-zero values are taken into account
 for the denominator. 

sum(diagonal(j) / number of non-zero elements in diagonal(j)

.. image:: ../../images/obs_exp_norm.png

Pearson correlation matrix
--------------------------

Pearson_i,j = C_i,j / sqrt(C_i,i * C_j,j) and C is the covariance matrix


.. image:: ../../images/pearson.png


.. image:: ../../images/obs_exp_pearson.png

The first image shows the Pearson correlation on the original interaction matrix, the second one shows 
the Person correlation matrix on an observed/expected matrix. A consecutive computation like this is used in 
the A/B compartment computation.


Covariance matrix
-----------------

Cov_i,j = E[M_i, M_j] - my_i * my_j where M is the input matrix and my the mean.
