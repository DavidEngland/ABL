Biazar:  
In general, I don't think stretching vertical grid is a good idea. I think we need a vertical grid structure with finer grid spacing where our fundamental target variables have high spatial gradient (such as potential temperature, ozone, ...), and relaxing it in other area. But, since we use fixed grid structure, I've been doing this manually in WRF for predefined estimated elevations (such as surface, PBL, and tropopause). I came up with some formulation to do this automatically by just giving it the estimated heights, number of levels, surface and top pressure (or height), but I can't find it now. Doing this, I was able to capture a tropopause folding event with exceptional agreement with lidar observation. I think non-uniform vertical grid spacing is another way of better capturing vertical mixing in high gradient regions.

My response:
What we really want is a vertical distribution of levels that follows the regions where the key variables have strong structure, for example potential temperature, ozone, shear, or moisture. In that sense, I agree that the better principle is not “stretch the grid,” but “concentrate resolution where the physics demands it.”

Mathematically, I think of it this way: if the local truncation or representation error scales like some function of curvature, then the grid spacing should be smallest where quantities such as

$$
\left| \frac{\partial \theta}{\partial z} \right|,\quad
\left| \frac{\partial^2 \theta}{\partial z^2} \right|,\quad
\left| \frac{\partial U}{\partial z} \right|,\quad
\left| \frac{\partial q}{\partial z} \right|
$$

are large, rather than being prescribed by a single monotone stretching law from the surface upward.

So I would now separate the ideas this way:

1. A stretched grid is useful for problems like shallow stable layers or LLJs, where the dominant gradients are often near the surface.
2. A targeted nonuniform grid is more general, because it can place levels around estimated heights such as the surface inversion, PBL top, jet nose, or tropopause.
3. In that sense, your WRF approach is more physically adaptive, even if the grid itself is fixed during a given run.

I also think your tropopause-folding example is important, because it shows that the benefit is not just numerical elegance. Better vertical placement can materially improve the representation of mixing and layer structure when sharp gradients are present aloft, not just near the ground.

What I am converging toward is a framework where the grid is designed from estimated feature heights and perhaps feature widths. Conceptually, that means using a level-density function $w(z)$ and distributing levels so that more of them fall where $w(z)$ is large. For example,

$$
w(z) = w_0 + \sum_i a_i \exp\!\left[-\frac{(z-h_i)^2}{2\sigma_i^2}\right]
$$

where $h_i$ are estimated important heights, $\sigma_i$ are their vertical influence scales, and $a_i$ controls how strongly we refine there. A simple stretched grid is then just a limiting special case where the weighting is biased primarily toward the surface.

So my current view is:

- you are right that global stretching is too restrictive as a general strategy
- the better idea is feature-targeted nonuniform spacing
- for LLJ-focused or shallow SBL problems, surface-biased stretching can still be a useful first approximation
- ultimately, the vertical grid should reflect expected gradient structure, not just distance from the surface