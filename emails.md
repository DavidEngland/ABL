From:  Dick McNider
David,

Glad you sent this email. BTW how much of your golfing picture was real. I assume that is you real face, although your beard was immaculately trimmed.

The equations you present show you must have continued to think about this . I was impressed.

I am attaching our new manuscript. Also,  response to reviewers. Especially Reviewer #3 where I gave you most all the credit for the mathematics in our 1995 paper.

The manuscript is basically a problem I have been working on for nearly forty five years. When I told Martha that she just rolled her eyes.

I have a question for you. In the paper even after trying to make our correction dependent on curvature in Ri we don't do that well  near the surface (see figure 11).  I think part of this is that with coarse grid we don't do a good job of getting the curvature near the surface.  See figure 10.

We mention that perhaps using similarity near the surface (below first model level) to calculate curvature in Ri might be a better path.

I was going to take a shot at this but since you are the expert thought I would ask you first!

My question to you is - can we use the similarity forms to find curvature in Ri in Ri. It seems that if you define Ri =( g/theta)( dtheta /dz)/ (( dv/dz)*(dv/dz)) that we should be able to use similarity to get profiles of V and theta. Then take second derivative of Ri.

What do you think. Also, to make it easier we approximate the stability function as an exponential.

From:  Arastoo P. Biazar

Thank you, David. It appears that you have put in quite a bit of work on this problem. I just glanced at your notes, and if I understand them correctly, you have addressed a couple of issues that we have been struggling with. First, are you suggesting that near-surface stability should be treated differently from the elevated stratification? Second, what would be the implications of using the geometric mean height in numerical models (different vertical grid spacing)?

Best,
Arastoo

From:  Me (David E. England, PhD)

I also made some significant progress on a 45 year old problem that I found as an undergrad at UNA.  Found out that if lower barriers by half that clustering will occur without any affinity.  Gave rise to a "hyperbolic strip" and hyperbolic weights.  Now have a very well paying job when I can get the work and a few other projects.

Not long ago I looked at MOST, sea-ice model and even quasi-hydrostatic correction as well as some other areas.  Been working on what I call the "Hasse-Stirling Framework" c.f.:  <https://github.com/DavidEngland/hasse-stirling/blob/main/Hasse-Stirling-Operator-Table.md>

Where interpret means, is just like finite differences, it is an estimate.  Our case, like a log transform, average in log space, but if transform back...Should look like a line plotted in log coordinates.  Just showing where the 2 point sqrt(z1*z2) comes from.  Took a while to recall Ri_1/2.  Other questions require discussion, no short, quick or easy way about it.  Quantum computing should allow for more and better.

Attaching some preliminary work, but I need to still get up to speed on some areas.  Having trouble viewing some articles.  Attached work is just for \phi(\zeta)=(1-\beta*\zeta)^-\alpha profiles.

BTW, most work is run through an AI LLM that have been helping to train.  Providing the right prompts is the key.

The hat is real in my FB golf photo.  Have a Harris tweed, but not with elbow patches.  Buddy Ricky Crawford, a Scottish soccer player at UAH visited the isles of Harris and Lewis.  They still weave it at home.
Cheers,
Dave

From:  Dick

Thanks David. Yes. Our results depend on the slope and intercept of curvature. Would prefer that it just depended on curvature. A better avenue would be to solve our ode with a less stringent assumption than ri =const. It may be that solving the ode numerically may be the best.

BTW you mentioned having a good job. Where are you working?

David, I overlooked the PDF you sent. Did you just solve the problem I asked? Incredible! Will have to take a longer look. But I am impressed big time.

From:  Me

Think so, going to look at Ri curvature for other formulations. Been reading papers, found a few interesting ones. Also derived a quadratic expansion for SBL in z/L.  Not sure if can get analytic expressions in terms of Ri yet. Should be able to repeat curvature analysis for at least four formulations. Also thinking about other planets and Titan, might offer insights about Earth’s polar regions.

Curvature just depended on parameters of Businger-Dyer profiles (1-bz/L)^-a for heat and momentum.  If L constant then scales by L^2 for z.

Guess I was just looking at near neutral case. See now that there are z/L and F(z/L) terms in general that will change.
