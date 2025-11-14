
In our paper we try to find a correction to the stability function f(Ri). Because in the discretized  Ri=g0zVzVz) = (g0)V2∆z  becomes bigger with ∆z  which means fs(Ri)=e-γRiRic becomes smaller (too little mixing). We then say we want to find a correction function
dd∆z(fsRifcRi,∆z)=0
We then solve the ODE  dd∆z-∆zre-g0V2∆zRicfcRi,∆z-∆zr =0
By assuming V2  is constant which is equivalent to Ri is a constant. My first question is can you solve this with a less restrictive form – e.g. V2   = a quadratic i.e.  . V2=c+ m∆z2.  If so would the final form have a form that contained c+ m∆z2.
Because we assumed that V2  is a constant we said  to account that this is not generally true we said we will define a tunable parameter D.   Since we know that no correction is needed if V2   is a linear function,  V2=c+ m∆z,    we said D must depend on curvature. So we said
D = D0 + M* |(d2 (Ri)/dz2|  and found through trials that
D = 0.3 +300 * |(d2 (Ri)/dz2|  with D < 0.7
So this is why we need curvature Ri or  2∂z2Ri.  We are assuming since there is no need for correction if V2=c+ m∆z  then the correction must depend on 2∂z2Ri.      But our form D = 0.3 +300 * |(d2 (Ri)/dz2|  with D < 0.7 seems to underestimate the correction near the surface where Ri is concave down but d2 (Ri)/dz2   is difficult to quantify in coarse grids. I think it over estimates the correction aloft because of our D0 = 0.3.  I believe D should only depend on 2∂z2Ri  i.e. D =  M* |(d2 (Ri)/dz2|  . Thus if I can get better estimates of curvature near the surface using similarity we may can make D0=0.
That is the reason that now I need to have a good formula for 2∂z2Ri.  Here is what I have from my analytics
2∂z2Ri=kg0 *u*2(0.94L)(1(1+0.47(z/L) )3)
You said in your email that  “I got essentially the same thing, have more powers of L in the denominator than you,”  Do you mean you have 0.94L2
If you can please send me your expression in the same notation as above.  I need to have everything correct before I start coding. Also, check my  u* and theta * below.









Ri=(g0)(∂θ/∂z) / (∂V/∂z)2.                                        (1)
Use M-O non-dimensional shear and temperature gradient
                       m=(kzu*)∂V∂z                                                                (2)
                       h=(kz*)∂θ∂z                                                               (3)
Or

 ∂V∂z = mu*kz
∂θ∂z = h*kz

Use (2) and (3) in (1)

Ri=(g0)((h*kz/ (mu*kz)2.
Assume Pr = 1
Then
m= h
Ri=(g0)kz*u*2(m)-1
Use Businger forms
m=1+0.47
h=1+0.47 ζ
  = z/L
Ri=(g0)(1+0.47 ζ)-1
Then
Ri=(g0)(kz*u*2(1+0.47 ζ)-1

Next take derivatives to find curvature in Ri  - curvature = 2∂z2Ri

∂Ri∂z= (g0)∂z(kz*u*2)(1+0.47 (z/L))-1
∂Ri∂z= g0 *u*2∂z(kz(1+0.47 (z/L))-1
or
∂Ri∂z= g0 *u*2((k)(1+0.47 (z/L))-1)+kz(-1 (0.47L) (1+0.47(z/L) )-2

Let g0 *u*2= j
ChatGPT simplifies by taking common denominator
∂Ri∂z=  jk(1+0.47z/L)-0.47z/L)(1+0.47 (z/L) )2​
∂Ri∂z=  =jk(1+0.47(z/L) )2

Take a second derivative to get curvature in Ri
2∂z2Ri=∂z(jk(1+0.47(z/L) )2)
Using chain rule
2∂z2Ri=jk(-2)0.47L(1(1+0.47(z/L) )3)
So
2∂z2Ri=kg0 *u*2(0.94L)(1(1+0.47(z/L) )3)

Use
    L=u*2/k(g/0) *

And from England and McNider 1995


u*2=k2fmRi1/2  |V1|2/ (lnz1z0)2,                               51a
u**=k2fhRi1/2V11-G /(lnz1z0)2              51b
For fm=fh
*=k  fmRi1/21/2V11-G/ (lnz1z0)

I plan to calculate this curvature in my visual basic code in Excel that simulates the GABLS 1 cases







Chatgpt
https://chatgpt.com/s/t_6912a565042081919b176d56cf899449
First Derivative
https://chatgpt.com/s/t_6912a565042081919b176d56cf899449
Chatgpt
Second derivative
https://chatgpt.com/s/t_6912aaf6fb248191b4d29c385dde2f64

