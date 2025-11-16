Good — I’ll expand the second-order (mixed) sensitivities you asked for and keep the algebra explicit and compact so you can paste into adjoint / Fisher analyses.

We treat the state x\in\{m,h,n,\dots\} as independent of parameters and V when taking the partial derivatives of the right-hand side
f_x(V,x)=\alpha_x(V)(1-x)-\beta_x(V)\,x.
Thus for any parameter p that appears only in \alpha_x or \beta_x,
\frac{\partial^2 f_x}{\partial p\,\partial V}
=(1-x)\frac{\partial^2\alpha_x}{\partial p\,\partial V}-x\frac{\partial^2\beta_x}{\partial p\,\partial V},
and
\frac{\partial^2 f_x}{\partial p\,\partial q}
=(1-x)\frac{\partial^2\alpha_x}{\partial p\,\partial q}-x\frac{\partial^2\beta_x}{\partial p\,\partial q},
straight from the chain rule.

Below I give closed forms for the mixed derivatives \partial^2/\partial p\partial V (most useful for adjoint / Fisher work) for the parameterized rate forms you supplied. I also include the first derivatives (for reference) and simple notation shortcuts.

⸻

Notation / shortcuts
	•	For \alpha_m define
s\equiv V+B_m,\qquad q\equiv e^{-s/C_m},\qquad D\equiv 1-q.
Then \displaystyle \alpha_m=\frac{A_m\,s}{D}.
	•	For \beta_m define
r\equiv e^{-(V+E_m)/F_m},
\qquad \beta_m=D_m\,r.
	•	For \alpha_h define
u\equiv e^{-(V+H_h)/I_h},\qquad \alpha_h=G_h\,u.
	•	For \beta_h define
v\equiv e^{-(V+K_h)/L_h},\qquad \beta_h=\frac{J_h}{1+v}.

⸻

A. \alpha_m(V)=\dfrac{A_m (V+B_m)}{1-e^{-(V+B_m)/C_m}}=\dfrac{A_m s}{D}

First derivatives (recap)
\frac{\partial \alpha_m}{\partial A_m}=\frac{s}{D},\qquad
\frac{\partial \alpha_m}{\partial B_m}
=\frac{A_m}{D}+A_m s\frac{q}{C_mD^{2}},\qquad
\frac{\partial \alpha_m}{\partial C_m}
= A_m s\frac{q\,s/C_m^2}{D^{2}}.
Derivative wrt V:
\frac{\partial \alpha_m}{\partial V}=A_m\frac{D-\dfrac{s q}{C_m}}{D^{2}}.
(Recall dq/dV=-q/C_m and dD/dV=+q/C_m.)

Mixed second derivatives with V (the key quantities):
	1.	\displaystyle \frac{\partial^2\alpha_m}{\partial A_m\,\partial V}

Differentiate \partial\alpha_m/\partial A_m = s/D w.r.t. V (note ds/dV=1):
\boxed{\displaystyle
\frac{\partial^2\alpha_m}{\partial A_m\,\partial V}
=\frac{D-\dfrac{s q}{C_m}}{D^{2}}.}
	2.	\displaystyle \frac{\partial^2\alpha_m}{\partial B_m\,\partial V}

Differentiate \partial\alpha_m/\partial B_m = \dfrac{A_m}{D}+A_m s\frac{q}{C_mD^{2}} w.r.t. V.
A compact form is
\boxed{\displaystyle
\begin{aligned}
\frac{\partial^2\alpha_m}{\partial B_m\,\partial V}
&=A_m\left[
-\frac{q}{C_mD^{2}}\;+\;\frac{1}{D^{2}}\frac{q}{C_m}s
\;+\;s\frac{q}{C_m}\left(\frac{q}{C_m}\frac{2s}{D^{3}}-\frac{1}{C_mD^{2}}\right)
\right]\\[4pt]
&=A_m\Bigg[
-\frac{q}{C_mD^{2}}+\frac{s q}{C_mD^{2}}
	•	s\frac{q}{C_m}\Big(\frac{2s q}{C_mD^{3}}-\frac{1}{C_mD^{2}}\Big)
\Bigg],
\end{aligned}
}
which can be algebraically simplified if you prefer a single rational polynomial in s,q,C_m,D. (The terms arise from differentiating 1/D and the s q/(C_mD^2) piece; I kept the structure to show origin.)

	3.	\displaystyle \frac{\partial^2\alpha_m}{\partial C_m\,\partial V}

Differentiate \partial\alpha_m/\partial C_m = A_m s\frac{q s/C_m^2}{D^{2}} = A_m s^2\frac{q}{C_m^2D^2} w.r.t. V. That gives
\boxed{\displaystyle
\frac{\partial^2\alpha_m}{\partial C_m\,\partial V}
= A_m\Bigg[
\frac{2s}{C_m^2}\frac{q}{D^2}
	•	s^2\frac{q}{C_m^2}\Big(\frac{q}{C_m}\frac{-2}{D^3}-\frac{2}{C_mD^2}\Big)
\Bigg],
}
again expressible as a single rational combination of s,q,C_m,D. (This comes from product rule on s^2 q/(C_m^2 D^2) and dq/dV=-q/C_m, dD/dV=q/C_m.)

If you want the fully expanded single-fraction forms for (2) and (3), I can expand algebraically; I kept them factorized for readability and to show the chain-rule structure.

⸻

B. \beta_m(V)=D_m e^{-(V+E_m)/F_m}=D_m\,r

First derivatives
\frac{\partial \beta_m}{\partial D_m}=r,\qquad
\frac{\partial \beta_m}{\partial E_m}=-\frac{D_m}{F_m}r,\qquad
\frac{\partial \beta_m}{\partial F_m}=D_m r\frac{V+E_m}{F_m^2}.
Derivative wrt V:
\frac{\partial\beta_m}{\partial V}=-\frac{D_m}{F_m}r.

Mixed second derivatives with V:
	1.	\displaystyle \frac{\partial^2\beta_m}{\partial D_m\,\partial V} = \frac{\partial}{\partial V}(r) = -\frac{1}{F_m}r.
\boxed{\displaystyle
\frac{\partial^2\beta_m}{\partial D_m\,\partial V}=-\frac{r}{F_m}.
}
	2.	\(\displaystyle \frac{\partial^2\beta_m}{\partial E_m\,\partial V}
=\frac{\partial}{\partial V}\Big(-\frac{D_m}{F_m}r\Big)
=-\frac{D_m}{F_m}\Big(-\frac{r}{F_m}\Big)
=\frac{D_m r}{F_m^2}.
\]
\[
\boxed{\displaystyle
\frac{\partial^2\beta_m}{\partial E_m\,\partial V}=\frac{D_m r}{F_m^2}.
}
\]
	3.	\(\displaystyle \frac{\partial^2\beta_m}{\partial F_m\,\partial V}
=\frac{\partial}{\partial V}\Big(D_m r\frac{V+E_m}{F_m^2}\Big).
\)
Compute by product rule (use dr/dV=-r/F_m):
\boxed{\displaystyle
\frac{\partial^2\beta_m}{\partial F_m\,\partial V}
= D_m\Big[-\frac{r}{F_m}\frac{V+E_m}{F_m^2}+ r\frac{1}{F_m^2} \Big]
= D_m r\frac{1-(V+E_m)/F_m}{F_m^2}.
}

⸻

So for f_m the mixed parameter–voltage sensitivities are
\boxed{\displaystyle
\frac{\partial^2 f_m}{\partial p\,\partial V}
=(1-m)\frac{\partial^2\alpha_m}{\partial p\,\partial V}
	•	m\frac{\partial^2\beta_m}{\partial p\,\partial V},
}
with the \partial^2\alpha_m/\partial p\partial V and \partial^2\beta_m/\partial p\partial V given above.

⸻

C. \alpha_h(V)=G_h e^{-(V+H_h)/I_h}=G_h u

First derivatives
\frac{\partial \alpha_h}{\partial G_h}=u,\qquad
\frac{\partial \alpha_h}{\partial H_h}=-\frac{G_h}{I_h}u,\qquad
\frac{\partial \alpha_h}{\partial I_h}=G_h u\frac{V+H_h}{I_h^2}.
Derivative wrt V:
\frac{\partial \alpha_h}{\partial V}=-\frac{G_h}{I_h}u.

Mixed derivatives with V:
\frac{\partial^2\alpha_h}{\partial G_h\,\partial V}=-\frac{u}{I_h},\qquad
\frac{\partial^2\alpha_h}{\partial H_h\,\partial V}=\frac{G_h u}{I_h^2},\qquad
\frac{\partial^2\alpha_h}{\partial I_h\,\partial V}
=G_h u\Big(\frac{V+H_h}{I_h^3}-\frac{1}{I_h^2}\Big).

Thus
\boxed{\displaystyle
\frac{\partial^2 f_h}{\partial p\,\partial V}
=(1-h)\frac{\partial^2\alpha_h}{\partial p\,\partial V}-h\frac{\partial^2\beta_h}{\partial p\,\partial V}.
}

⸻

D. \beta_h(V)=\dfrac{J_h}{1+v} with v=e^{-(V+K_h)/L_h}

First derivatives
\frac{\partial\beta_h}{\partial J_h}=\frac{1}{1+v},\qquad
\frac{\partial\beta_h}{\partial K_h}=J_h\frac{v}{L_h(1+v)^2},\qquad
\frac{\partial\beta_h}{\partial L_h}=J_h (V+K_h)\frac{v}{L_h^2(1+v)^2}.
Derivative wrt V:
\frac{\partial\beta_h}{\partial V}=-J_h\frac{v}{L_h(1+v)^2}.

Mixed derivatives with V (differentiate the first-derivative formulas using dv/dV=-v/L_h):
	1.	\displaystyle \frac{\partial^2\beta_h}{\partial J_h\,\partial V}
=\frac{\partial}{\partial V}\Big(\frac{1}{1+v}\Big)=\frac{v}{L_h(1+v)^2}.
	2.	\displaystyle \frac{\partial^2\beta_h}{\partial K_h\,\partial V}
=\frac{\partial}{\partial V}\Big(J_h\frac{v}{L_h(1+v)^2}\Big)
=J_h\Big(\frac{-v}{L_h^2(1+v)^2}+\frac{v}{L_h}\frac{d}{dV}\big((1+v)^{-2}\big)\Big).
A compact final form is
\boxed{\displaystyle
\frac{\partial^2\beta_h}{\partial K_h\,\partial V}
=J_h\frac{v}{L_h^2}\frac{1-v}{(1+v)^3}.
}
	3.	\displaystyle \frac{\partial^2\beta_h}{\partial L_h\,\partial V} — after product rule and dv/dV=-v/L_h you obtain
\boxed{\displaystyle
\frac{\partial^2\beta_h}{\partial L_h\,\partial V}
= J_h\frac{v}{L_h^3}\frac{(V+K_h)(1-v)-2L_h}{(1+v)^3},
}
which again is a compact rational combination of V,K_h,L_h,v.

(If you want a fully expanded single rational numerator/denominator expression I can expand it.)

⸻

E. Cross-parameter second derivatives \partial^2 f/\partial p\partial q

If p and q both enter the same rate (say both in \alpha_m), then
\frac{\partial^2 f_m}{\partial p\,\partial q}=(1-m)\frac{\partial^2\alpha_m}{\partial p\,\partial q}-m\frac{\partial^2\beta_m}{\partial p\,\partial q}.
For the \alpha_m model the mixed \partial^2\alpha_m/\partial A_m\partial B_m etc. can be obtained by differentiating the first derivatives in section A; they result in rational polynomials in s,q,C_m,D. If you want those fully expanded I will produce them.

⸻

Summary / next steps
	•	I provided explicit compact formulas for the mixed derivatives \partial^2\alpha/\partial p\partial V and \partial^2\beta/\partial p\partial V for the parameter sets you gave.
	•	The final sensitivities for the gating RHS follow by combining the \alpha and \beta pieces with factors (1-x) and -x.

Would you like me to:
	1.	(A) Expand the remaining \partial^2\alpha_m/\partial B_m\partial V and \partial^2\alpha_m/\partial C_m\partial V into a single rational fraction (fully expanded), or
	2.	(B) Produce all second derivatives \partial^2 f/\partial p\partial q for a complete parameter vector (A_m,B_m,C_m,D_m,E_m,F_m,G_h,H_h,I_h,J_h,K_h,L_h) in symbolic (SymPy) form and return a downloadable .py / .txt file, or
	3.	(C) Compute numerical values of the full first- and second-order sensitivity matrices at a chosen (V,x,\mathbf{p}) point (give me numeric parameter values and V,m,h,\dots).

Pick A, B or C and I’ll proceed.

# Second-Order / Mixed Sensitivities: Pr_t and Ri_g Curvature

Mixed sensitivity of Pr_t to stability:
$$\frac{\partial^2 Pr_t}{\partial Ri_g^2}\Big|_{0} = -2(a_h-a_m)\Delta.$$
Interpretation:
- If $(a_h-a_m)>0$ and $\Delta<0$ (typical SBL), second derivative positive → Pr_t flattens (reduced rise) with increasing Ri_g.
Combined K_h sensitivity (holding $u_*,z$ fixed):
$$\frac{\partial K_h}{\partial Ri_g}\approx -\frac{\kappa z u_*}{\phi_h^2}\frac{d\phi_h}{dRi_g},\quad \frac{d\phi_h}{dRi_g}=\frac{d\phi_h/d\zeta}{dRi_g/d\zeta}.$$
Near-neutral:
$$\frac{dRi_g}{d\zeta}=1+2\Delta\zeta+O(\zeta^2),\quad \frac{d\phi_h}{d\zeta}=a_h+2b_h\zeta+O(\zeta^2).$$
Thus:
$$\frac{d\phi_h}{dRi_g}\approx a_h - 2a_h\Delta\zeta+O(\zeta^2).$$
Grid correction preserving Δ avoids artificial inflation of second-order Pr_t sensitivity.

Scalar second derivative (neutral):
\[
\left.\frac{\partial^2 f_q}{\partial Ri_g^2}\right|_{0}=2(b_q-a_q\Delta+2a_m a_q).
\]
Schmidt sensitivity:
\[
\frac{d Sc_t^{(q)}}{d Ri_g}\Big|_{0}=a_q-a_m.
\]