# Soft solid with a liquid drop

In this code, we will let a viscous or viscoelastic liquid drop rest on a soft solid film until it reaches an equilibrium state. The gravity in this case should be in the -x direction only.
First run this code and then proceed with the code: softsliding.c.

The equilibrium state will depend on (\alpha, Bo, Ec)... So, it needs to be run for every case that we are interested in for this project.

## changelog Sep 24, 2024 v4.0 ($Ec \gg 1$, use with $De \gg \sqrt{Ec}$ for elastic solids)

#axi
Here, we combine v2.1 and v3.0. -- this is the prefereed version for viscoelastic films at large Ec, De. In the limit $De \to \infty$, this code should give an elastic film response.
- must ensure that $De \gg \sqrt{Ec}$

## changelog Sep 22, 2024 v2.1 (arbitary De)

#axi
Here, we oblitrate v3.0 and then add viscoelasticity. Hopefully, with large enough De, the system would mimic purely elastic behavior. -- this version is prefered for viscoelastic films at arbitary De.

## changelog Sep 22, 2024 v3.0 (very slow!)

#axi
This code supports very large Ec. This is done by renormalizing the equations such that the elastic modulus is always 1. (it is very slow!)

## changelog Aug 23, 2024 v2.0 ($Ec \leq 10$ works well)

This code is in axi. The log-conform-elastic.h is already compatible with axi.  -- this version is prefered for purely elastic films with Ec < 1 (maybe even 10).

## changelog Aug 10, 2024 v1.1

* Decreased domain size. The drop is now at the center of the domain. Run this version to get the equilibrium shape of the drop+film, export the interface as .dat file that could be directly/indirectly read into the softsliding.c code. ## TODO: Jnandeep....

## changelog Aug 5, 2024 v1.0

* This code uses the reduced gravity formulation described here: [https://www.annualreviews.org/content/journals/10.1146/annurev-fluid-122316-045034](https://www.annualreviews.org/content/journals/10.1146/annurev-fluid-122316-045034)