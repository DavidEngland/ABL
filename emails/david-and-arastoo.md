David and Arastoo,
I am trying to make our results better using only curvature as we think about putting a proposal together for NSF/DOE.  The attached plot shows results for D - 2000* Ricurvature rather than the D=0.3+300*Ricurvature we had in our submitted paper.
Part of the disagreement is near the surface. The temperatures in the coarser grid are colder than the 2-m grid.  We already use similarity theory to get the fluxes - u* and theta*u*. Below is a page from England McNider 1995.
I was thinking that because are grids are so large it might be better to use similarity to specify temperature and wind at the first level. However, since we use similarity to get u* and theta*but we use the first level temperature and wind. If we now use similarity to get temperature or wind we have an inconsistency since we are trying to solve for two variables u* and u(1) with only one equation and theta* and T(1)  with only one equation. Do you think I can solve for the two variables iteratively. This is what we used to do (before David's paper)  to get L from u* then u* depends on new L.
Equations  50, 51a and 51b from E&Mc 1995 are give below. Would need to first solve for u* and theta* 51 a and 51b. Then invert these and solve for theta(1) and V(1). Then recompute Ri(1/2) and solve for new u* and theta* until they converge.
It does take us backwards in that we lose the efficiency of David's approach.
Alternatively we could make grids smaller near the surface to improve the coarse grid prediction. This may be more efficient than the iteration approach.
Any thoughts on this before I start. It may be easier and more computationally eff in icient to just to add a few more layers near the surface for the coarse grids. May try adding grids first.
Dick
