This file contains tiny reference snippets (VB and Fortran77) to compute the power‑law (SOL‑PL) and exponential (SOL‑EXP) forms of f_c described in the overview. Use these as drop‑in examples for Dick (VB) and Arastoo (Fortran77).

Visual Basic (VB6 / VBA) — exponential & power‑law
```vb
' filepath: /Users/davidengland/Documents/GitHub/ABL/implementations/fc_examples.md
' VB example: compute fc (exponential and power-law variants)
Function fc_exp(Dz As Double, zeta As Double, B As Double, _
                alpha As Double, Dz_ref As Double, zeta_ref As Double, _
                p As Double, q As Double, fc_min As Double) As Double
    Dim exponent As Double, ratio As Double
    If B <= 1.05 Then
        fc_exp = 1#
        Exit Function
    End If
    ratio = (Dz / Dz_ref) ^ p
    exponent = -alpha * (B - 1#) * ratio * ( (zeta / zeta_ref) ^ q )
    fc_exp = Exp(exponent)
    If fc_exp < fc_min Then fc_exp = fc_min
End Function

Function fc_power(Dz As Double, zeta As Double, B As Double, _
                  alpha As Double, Dz_ref As Double, zeta_ref As Double, _
                  q As Double, fc_min As Double) As Double
    ' SOL-PL exact power-law form
    If B <= 1.05 Then
        fc_power = 1#
        Exit Function
    End If
    fc_power = (Dz / Dz_ref) ^ ( -alpha * (B - 1#) * ( (zeta / zeta_ref) ^ q ) )
    If fc_power < fc_min Then fc_power = fc_min
End Function
```

Fortran 77 — exponential & power‑law
```fortran
! filepath: /Users/davidengland/Documents/GitHub/ABL/implementations/fc_examples.md
C Fortran77 example: compute fc (exponential and power-law variants)
      DOUBLE PRECISION FUNCTION FC_EXP(DZ, ZETA, B, ALPHA, DZ_REF, ZETA_REF, P, Q, FC_MIN)
      DOUBLE PRECISION DZ, ZETA, B, ALPHA, DZ_REF, ZETA_REF, P, Q, FC_MIN
      DOUBLE PRECISION RATIO, EXPNT
      IF (B .LE. 1.05D0) THEN
         FC_EXP = 1.0D0
         RETURN
      END IF
      RATIO = (DZ / DZ_REF) ** P
      EXPNT = -ALPHA * (B - 1.0D0) * RATIO * ( (ZETA / ZETA_REF) ** Q )
      FC_EXP = DEXP(EXPNT)
      IF (FC_EXP .LT. FC_MIN) FC_EXP = FC_MIN
      RETURN
      END

      DOUBLE PRECISION FUNCTION FC_POWER(DZ, ZETA, B, ALPHA, DZ_REF, ZETA_REF, Q, FC_MIN)
      DOUBLE PRECISION DZ, ZETA, B, ALPHA, DZ_REF, ZETA_REF, Q, FC_MIN
      IF (B .LE. 1.05D0) THEN
         FC_POWER = 1.0D0
         RETURN
      END IF
      FC_POWER = (DZ / DZ_REF) ** ( -ALPHA * (B - 1.0D0) * ( (ZETA / ZETA_REF) ** Q ) )
      IF (FC_POWER .LT. FC_MIN) FC_POWER = FC_MIN
      RETURN
      END
```

Notes
- The VB code is VBA/VB6 style; it can be pasted into Excel VBA modules or legacy VB code.
- The Fortran77 example uses standard intrinsic DEXP and double precision; adapt naming conventions to your build.
- Both snippets use a simple B threshold (1.05) and apply a floor fc_min; adjust constants to your tuning.
