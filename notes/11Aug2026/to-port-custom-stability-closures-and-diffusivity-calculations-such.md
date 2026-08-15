To port custom stability closures and diffusivity calculations (such as the MOSTGeneral formulation and topological $d_{\text{fold}}$ manifold correction) into existing atmospheric models like WRF (Weather Research and Forecasting model), you must interface with WRF's Fortran physics hierarchy and registry system.  
  
In WRF, vertical mixing is split across two primary physics sub-drivers in phys/:  
  
1. **Surface Layer (SL) Schemes** (phys/module_sf_*.F): Compute diagnostic exchange coefficients ($C_m, C_h$), friction velocity ($u_*$), and surface fluxes near $z_1$.   
2. **Planetary Boundary Layer (PBL) Schemes** (phys/module_bl_*.F): Compute column eddy diffusivities ($K_m, K_h$) and vertical turbulent transport across layers $k = 1 \dots N_z$.   
# Step-by-Step Implementation Strategy in WRF  
## 1. Register New Variables and Namelist Switches (Registry/Registry.EM_COMMON)  
To expose new parameters ($\beta_m, \beta_h$) or 3D fields ($d_{\text{fold}}$) to WRF's input files (namelist.input) and state arrays:  
  
Plaintext  
  
# Add a custom PBL option flag to namelist &physics  
rship   custom_beta_m   namelist,physics   1   5.0   h   "custom_beta_m"   "Linear MOST beta_m parameter" ""  
rship   custom_beta_h   namelist,physics   1   5.0   h   "custom_beta_h"   "Linear MOST beta_h parameter" ""  
  
# Add 3D auxiliary topological field d_fold to state memory (if dynamic or read from input)  
state   real   d_fold   ikj   dyn_em   1   -   rh   "D_FOLD"   "Topological fold distance"   "m"  
## 2. Fortran Implementation of the Closure (phys/module_bl_custom.F)  
Create a custom PBL module that encapsulates the numerically stable quadratic root solver and column diffusivity calculation:  
  
Fortran  
  
MODULE module_bl_custom  
   USE module_model_constants, ONLY: g, karman  
  
CONTAINS  
  
   SUBROUTINE custom_pbl( &  
        u3d, v3d, th3d, z3d, dz8w, &  
        d_fold, km, kh, pr_t, &  
        beta_m, beta_h, delta_0, &  
        ids, ide, jds, jde, kds, kde, &  
        ims, ime, jms, jme, kms, kme, &  
        its, ite, jts, jte, kts, kte &  
   )  
      IMPLICIT NONE  
  
      ! Array dimensions & bounds  
      INTEGER, INTENT(IN) :: ids, ide, jds, jde, kds, kde, &  
                             ims, ime, jms, jme, kms, kme, &  
                             its, ite, jts, jte, kts, kte  
  
      ! State fields  
      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(IN)  :: u3d, v3d, th3d, z3d, dz8w, d_fold  
      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(OUT) :: km, kh, pr_t  
      REAL, INTENT(IN) :: beta_m, beta_h, delta_0  
  
      ! Local variables  
      INTEGER :: i, j, k  
      REAL :: dz, du, dv, dth, shear_sq, shear, l_mix, N2, Ri_g, Ri_c  
      REAL :: bm, bh, disc, denom, zeta, phi_m, phi_h, fm, fh, C_fold, base_K  
      REAL, PARAMETER :: l0 = 30.0, shear2_min = 1.0e-12  
  
      Ri_c = beta_h / (beta_m**2)  
  
      DO j = jts, jte  
         DO i = its, ite  
            DO k = kts+1, kte-1  
               ! 1. Grid gradients  
               dz  = 0.5 * (dz8w(i,k+1,j) + dz8w(i,k,j))  
               du  = (u3d(i,k+1,j) - u3d(i,k-1,j)) / (2.0 * dz)  
               dv  = (v3d(i,k+1,j) - v3d(i,k-1,j)) / (2.0 * dz)  
               dth = (th3d(i,k+1,j) - th3d(i,k-1,j)) / (2.0 * dz)  
  
               shear_sq = MAX(du**2 + dv**2, shear2_min)  
               shear    = SQRT(shear_sq)  
  
               ! 2. Blackadar mixing length  
               l_mix = (karman * z3d(i,k,j)) / (1.0 + (karman * z3d(i,k,j)) / l0)  
  
               ! 3. Gradient Richardson number  
               N2   = (g / th3d(i,1,j)) * dth  
               Ri_g = MAX(0.0, N2 / shear_sq)  
  
               ! 4. Evaluate MOSTGeneral via Citardauq formulation  
               IF (Ri_g >= Ri_c) THEN  
                  fm   = 0.0  
                  fh   = 0.0  
                  pr_t(i,k,j) = 1.0  
               ELSE  
                  disc  = MAX(0.0, 1.0 + 4.0 * (beta_h - beta_m) * Ri_g)  
                  denom = (1.0 - 2.0 * beta_m * Ri_g) + SQRT(disc)  
                    
                  IF (denom <= 0.0) THEN  
                     fm = 0.0  
                     fh = 0.0  
                     pr_t(i,k,j) = 1.0  
                  ELSE  
                     zeta  = (2.0 * Ri_g) / denom  
                     phi_m = 1.0 + beta_m * zeta  
                     phi_h = 1.0 + beta_h * zeta  
                       
                     fm = 1.0 / (phi_m**2)  
                     fh = 1.0 / (phi_m * phi_h)  
                     pr_t(i,k,j) = phi_h / phi_m  
                  END IF  
               END IF  
  
               ! 5. Diffusivities & Geometric Suppression Factor  
               base_K = (l_mix**2) * shear  
               C_fold = TANH(d_fold(i,k,j) / delta_0)  
  
               km(i,k,j) = base_K * fm * C_fold  
               kh(i,k,j) = base_K * fh * C_fold  
            END DO  
  
            ! Boundary conditions at bottom/top boundaries  
            km(i,1,j)   = km(i,2,j)  
            kh(i,1,j)   = kh(i,2,j)  
            pr_t(i,1,j) = pr_t(i,2,j)  
  
            km(i,kte,j)   = km(i,kte-1,j)  
            kh(i,kte,j)   = kh(i,kte-1,j)  
            pr_t(i,kte,j) = pr_t(i,kte-1,j)  
         END DO  
      END DO  
   END SUBROUTINE custom_pbl  
  
END MODULE module_bl_custom  
## 3. Hook Module into WRF Drivers & Build System  
1. **Add to phys/Makefile**: Include module_bl_custom.o in the physics compilation target list so make compiles it before the driver.   
2. **Hook into phys/module_pbl_driver.F**:   
    * Add USE module_bl_custom, ONLY: custom_pbl at top.   
    * Add a CASE statement inside pbl_driver mapped to your chosen bl_pbl_physics option integer:  Fortran  CASE (CUSTOMPBLSCHEME) ! e.g., bl_pbl_physics = 99  
    *    CALL custom_pbl( &  
    *         u3d=u_phy, v3d=v_phy, th3d=th_phy, z3d=z, dz8w=dz8w, &  
    *         d_fold=d_fold, km=rublbl, kh=rthblbl, pr_t=pr_t, &  
    *         beta_m=config_flags%custom_beta_m, beta_h=config_flags%custom_beta_h, &  
    *         delta_0=1.0e-3, &  
    *         ids=ids, ide=ide, jds=jds, jde=jde, kds=kds, kde=kde, &  
    *         ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme, &  
    *         its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte &  
    *    )  
    *      
3. **Recompile WRF**: Run ./compile em_real from the root directory.   
## Alternative Modeling Frameworks  
If modifying WRF's legacy Fortran codebase presents too high of an integration overhead for rapid iteration:  
  
* **MPAS-Atmosphere**: WRF's modern Fortran successor (uses dynamic unstructured Voronoi meshes). Physics routines share similar structure to WRF but feature cleaner module interfaces.   
* **ClimaAtmos.jl / Oceananigans.jl**: Julia-native geophysical fluid dynamics engines. If you already have your closures written in Julia, these frameworks allow direct zero-cost function dispatch of your ClosureState routines without Fortran translation or driver wrapper overhead.   
Are you looking to modify an existing WRF scheme directly (e.g., tweaking MYNN or YSU in-place), or are you building a standalone custom physics module?  
