The Inverse Logarithmic Mean: $\frac{\ln(z/z_0)}{z-z_0}$
While we often focus on heights ($z_g$, $z_L$), the expression $\frac{\ln(z/z_0)}{z-z_0}$ is the Layer-Averaged Gradient ($\langle \partial U / \partial z \rangle$).
1. Physical Definition
In a neutral surface layer, the wind profile is $U(z) = \frac{u_*}{\kappa} \ln\left(\frac{z}{z_0}\right)$.
The gradient (shear) is $\frac{\partial U}{\partial z} = \frac{u_*}{\kappa} \cdot \frac{1}{z}$.
The Bulk Shear measured across a grid layer from $z_0$ to $z$ is:

$$\frac{\Delta U}{\Delta z} = \frac{U(z) - U(z_0)}{z - z_0} = \frac{u_*}{\kappa} \left[ \frac{\ln(z/z_0)}{z - z_0} \right]$$
2. Why it matters for ML and Bias
This expression is the exact scaling factor between the friction velocity ($u_*$) and the observed wind difference.
The Jensen Link: Because $1/z$ is convex, the average value of $1/z$ over the interval $[z_0, z]$ is always greater than the value of $1/z$ at the midpoint:
$$\frac{\ln(z/z_0)}{z-z_0} > \frac{1}{(z + z_0)/2}$$
The Numerical Correction: If an NWP model assumes the gradient is simply the wind difference divided by the distance, it is implicitly using this inverse logarithmic mean.
The Curvature Connection: If you multiply this by the geometric mean height $z_g$, you get:
$$\text{Scaled Bias Term} = z_g \cdot \frac{\ln(z/z_0)}{z-z_0}$$

This term often appears in the derivation of the Bias Ratio ($B$), acting as the geometric correction for how the log-layer "bends" over the grid height.
3. Summary for Students
When you see $\frac{\ln(z/z_0)}{z-z_0}$, don't just see a formula; see the Average Steepness of the wind profile across that layer.
It is the "True Bulk Gradient."
It is mathematically equal to $1/z_L$.
It is the reason why coarse-grid models ($large \Delta z$) naturally struggle to "see" the intense shear happening right at the surface.
