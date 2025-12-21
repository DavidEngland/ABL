England and McNider’s 1995 stability equations were corrected primarily because they contained two distinct mathematical errors for the stable boundary layer: a **"tipping error"** in the stability parameter equation and a **"wrongly derived relation"** for the stability function for heat.,,

The specific reasons for these corrections are detailed below:

### **1. The Stability Parameter Tipping Error**
The original 1995 publication included a quadratic equation for the non-dimensional stability parameter ($\zeta$), but the published solution was **missing a minus sign** before the square root.,
*   **Reason for Correction:** Without the negative root, the equation would produce a **nonphysical solution** where $\zeta$ does not equal zero when the Richardson number ($Ri$) is zero.
*   **Context:** It is believed this was a simple typographical (tipping) error, as the authors elsewhere in the paper discussed the necessity of taking the negative root to maintain physical consistency.,

### **2. The Heat Stability Function Derivation Error**
The second correction addressed an error in the derivation of the stability function for heat ($f_h$) for stable stratification when the Prandtl number is not equal to one.,
*   **Reason for Correction:** The original derivation was found to produce **nonphysical negative values** for the heat stability function within certain Richardson number limits.,,
*   **Consequence:** In reality, the value for $f_h$ must not fall below zero; therefore, a corrected analytical version was required to ensure the model produces realistic results.,

### **3. Generalization for Modern Application**
Beyond these two errors, modern researchers corrected and generalized the equations to account for cases where the **Prandtl number does not equal one**.,,
*   **Empirical Realities:** The 1995 equations were derived under the assumption that the empirical constants for momentum ($a_m$) and heat ($a_h$) were equal.,,
*   **Updated Physics:** Subsequent re-evaluations of similarity theory (such as by Högström in 1996) showed that these constants **do not have the same value** in the real atmosphere., Correcting these equations allowed for the development of **generalized formulae** that remain stable and accurate even when the transport of heat and momentum are not assumed to be identical.,

**Analogy**
Correcting these equations is like **fixing a GPS formula that forgot to subtract the elevation**. Without the correction, the GPS might tell you that you are at the top of the mountain while you are still standing at the trailhead ($Ri = 0$). Fixing the "tipping error" and "derivation error" ensures the math actually lands the user at the **correct physical location** rather than giving impossible, "below-sea-level" (negative) coordinates.,