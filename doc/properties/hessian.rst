Analytic GFN1-xTB Hessian and vibrational analysis
==================================================

``tblite`` provides an analytic Cartesian Hessian for non-periodic GFN1-xTB.
The implementation differentiates the analytic GFN1-xTB gradient and solves the
first-order self-consistent response equations; production Hessians therefore do
not use finite differences.

Supported scope
---------------

The analytic Hessian supports:

* closed-shell and collinear open-shell GFN1-xTB,
* zero or finite electronic temperature,
* the GFN1 halogen correction,
* ALPB and GBSA solvation, including the supported CDS contributions,
* harmonic vibrational analysis with translational and rotational projection.

The following are deliberately outside the supported scope:

* GFN2-xTB and other Hamiltonians,
* periodic systems,
* external electric fields,
* ddX solvation models.

Unsupported combinations fail explicitly rather than silently omitting a
second-derivative contribution.

Command line
------------

Evaluate the Cartesian Hessian with::

   tblite run --method gfn1 --hess coord.xyz

Run the Hessian followed by harmonic vibrational analysis with::

   tblite run --method gfn1 --freq coord.xyz

The Hessian is reported in Hartree/Bohr\ :sup:`2`. Frequencies are reported in
cm\ :sup:`-1`; imaginary modes use negative frequencies. The vibrational
analysis projects five rigid-body modes for linear molecules and six for
nonlinear molecules.

Python
------

The Python interface returns the Hessian and vibrational information in the
usual result object::

   from tblite.interface import Calculator

   calc = Calculator("GFN1-xTB", numbers, positions)
   res = calc.hessian()
   hessian = res.get("hessian")
   frequencies = res.get("frequencies")
   modes = res.get("normal-modes")
   n_imaginary = res.get("n-imaginary")

``hessian`` has shape ``(3*nat, 3*nat)``. Cartesian normal modes are returned as
one displacement vector per vibrational mode.

Numerical considerations
------------------------

The electronic response is defined only for a differentiable electronic state.
An exact occupied/virtual frontier degeneracy at zero electronic temperature can
make the state-specific response non-unique. In that situation the Hessian
routine reports an error; a small finite electronic temperature can be used when
a state-averaged response is appropriate.

The implementation validates the coupled response equations internally and
symmetrizes only round-off-level antisymmetric noise. Tests compare analytic
Hessian columns against central differences of the analytic gradient and also
check finiteness, Hessian symmetry, translational invariance, and rotational
covariance.

Vibrational classification
--------------------------

By default, a vibrational frequency below ``-20 cm^-1`` is counted as imaginary.
The count can be used to distinguish a local minimum (zero significant imaginary
modes) from a first-order saddle point (one significant imaginary mode), but this
classification should be interpreted together with the optimized geometry and
the corresponding normal mode.
