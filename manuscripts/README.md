# Manuscript Drafts — Coordination Guide

## Active Manuscripts

## Program-Level Planning (Cross-Paper)

- Paper-series roadmap (not a manuscript draft): `program/PaperSeries_Gegenbauer_Roadmap_v01.md`

### 1. Grid-Curvature Correction for SBL (Primary)
**Working Title:** "Grid-Dependent Corrections to Stable Boundary Layer Turbulence Closures: A Curvature-Aware Approach"

**Target Journal:** Boundary-Layer Meteorology (first choice) or JAMC  
**Lead Author:** England  
**Co-authors:** McNider, Biazar  
**Status:** First draft in progress  
**File:** `Grid_Curvature_SBL_v01.md`

**Key Sections Status:**
- [x] Abstract (draft)
- [x] Introduction
- [x] Theory (curvature derivation)
- [ ] Methods (needs McNider input on GABLS1 setup)
- [ ] Results (needs Biazar validation cases)
- [ ] Discussion
- [ ] Conclusions

**Action Items:**
- McNider: Review dynamic Ri_c* formulation (Section 5)
- Biazar: Provide Dallas tower comparison (Figure 6)
- England: Complete sensitivity analysis (Section 4.3)

---

### 2. Arctic Amplification & Boundary Layer Physics (Secondary)
**Working Title:** "The Role of Near-Surface Stability Misrepresentation in Arctic Climate Model Biases"

**Target Journal:** Journal of Climate or GRL  
**Lead Author:** McNider (TBD)  
**Co-authors:** England, Biazar  
**Status:** Outline only  
**File:** `Arctic_SBL_Climate_v01.md`

**Focus:** Link grid-curvature bias to Arctic amplification errors in CMIP6 models

---

### 3. WRF Integration Paper (Implementation-Focused)
**Working Title:** "Grid-Aware Ri-Curvature and Dynamic Critical Richardson Number Integration in WRF Stable Boundary Layer Schemes"

**Target Journal:** Monthly Weather Review (primary) or JAMC  
**Lead Author:** England (proposed)  
**Co-authors:** McNider, Biazar  
**Status:** Outline drafted  
**File:** `WRF_Ri_Curvature_Integration_Outline_v01.md`

**Focus:** WRF-ready implementation, SCM+3D validation, and operationally safe defaults

---

## File Naming Convention

- Use semantic versioning: `Title_v01.md`, `Title_v02.md`, etc.
- Date major revisions in filename: `Title_v01_2025-01-15.md`
- Keep previous versions in `archive/` folder
- Final submission version: `Title_FINAL_SUBMITTED_YYYY-MM-DD.md`

---

## Collaboration Workflow

1. **Draft Development (Markdown)**
   - Use Markdown for easy version control and commenting
   - Inline comments: `<!-- McNider: Please review this claim -->`
   - Action items: `**[ACTION-McNider]** Verify SHEBA case setup`

2. **Figure Preparation**
   - Store figures in `manuscripts/figures/`
   - Use descriptive names: `fig01_bias_ratio_vs_grid.png`
   - Keep raw data/scripts in `manuscripts/data/` and `manuscripts/scripts/`

3. **LaTeX Conversion (Final Stage)**
   - Convert Markdown → LaTeX when ready for submission
   - Use journal-specific templates in `manuscripts/templates/`

4. **Review Cycles**
   - McNider reviews: tag with `<!-- RTM-COMMENT: ... -->`
   - Biazar reviews: tag with `<!-- APB-COMMENT: ... -->`
   - Resolve comments before incrementing version number

---

## Shared Resources

### Literature Database
- `manuscripts/references.bib` — Shared BibTeX database
- Key papers already included:
  - England & McNider (1995)
  - Businger et al. (1971)
  - GABLS intercomparison papers
  - Arctic amplification references

### Data Access
- GABLS1 LES data: `data/GABLS1/`
- ARM NSA tower: `data/ARM_NSA/`
- SHEBA: `data/SHEBA/`
- Dallas tower (Biazar): Request access via APB

### Code Repositories
- Fortran correction module: `implementations/McNider_1DBLM_fc_module.f90`
- Python analysis tools: `notebooks/`
- Validation scripts: `manuscripts/scripts/`

---

## Journal-Specific Guidelines

### Boundary-Layer Meteorology
- Format: LaTeX (Springer template)
- Length: ~8000 words (flexible)
- Figures: Max 10–12
- Supplementary material: Encouraged
- Submission: Editorial Manager

### JAMC
- Format: LaTeX (AMS template)
- Length: ~6000 words
- Figures: Max 8–10
- Supplementary material: Limited
- Submission: BAMS system

### Journal of Climate
- Format: LaTeX (AMS template)
- Length: ~7000 words (broader scope)
- Figures: Max 12–15
- Supplementary material: Encouraged
- Focus: Climate model implications

---

## Communication Protocol

### Weekly Check-ins
- **Frequency:** Every Monday 10 AM CST
- **Platform:** Zoom or email thread
- **Agenda template:**
  1. Progress since last meeting
  2. Blockers/questions
  3. Action items for next week
  4. Timeline adjustments

### Email Subject Lines
- Use prefix: `[ABL-MS1]` for primary manuscript
- Use prefix: `[ABL-MS2]` for Arctic paper
- Example: `[ABL-MS1] Draft Section 3 ready for review`

### Urgent Issues
- Direct phone/text for time-sensitive matters
- Tag in document: `**[URGENT-RTM]**` or `**[URGENT-APB]**`

---

## Timeline (Tentative)

### Primary Manuscript (Grid-Curvature SBL)
- **Jan 31, 2025:** First complete draft to co-authors
- **Feb 15, 2025:** Incorporate feedback, finalize figures
- **Feb 28, 2025:** Submit to BLM
- **Apr 2025:** Target acceptance (2-month review typical)

### Secondary Manuscript (Arctic Climate)
- **Mar 15, 2025:** Complete outline and literature review
- **Apr 30, 2025:** First draft
- **Jun 2025:** Submit to J. Climate

---

## Contact Information

**David E. England**  
Email: david.england@uah.edu  
Phone: [redacted]  
GitHub: @DavidEngland

**Richard T. McNider**  
Email: mcnider@uah.edu  
Office: NSSTC 3051

**Arastoo P. Biazar**  
Email: biazar@uah.edu  
Office: NSSTC 3053

---

**Last Updated:** January 2025
