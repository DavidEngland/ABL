"""Placeholder: MOST profile families and Newton inversion (to be extracted/refactored).

Will move and refactor code/profiles.py into this module with:
  - make_profile(tag, pars) for 8+ MOST families
  - zeta_from_ri_gradient() for Newton inversion
  - profile_catalog() for canonical IDs (BD71, BH91, CB05, GF97)
"""


def make_profile(tag: str, pars: dict):
    """Placeholder: Return (phi_m, phi_h) callables for MOST family."""
    raise NotImplementedError("Profile module pending extraction from code/profiles.py (Milestone 0)")


def zeta_from_ri_gradient(ri_target: float, profile_tag: str, pars: dict) -> float:
    """Placeholder: Invert Ri_g → ζ using Newton iteration."""
    raise NotImplementedError("Profile inversion pending extraction (Milestone 0)")


def profile_catalog() -> dict:
    """Placeholder: Return canonical profile registry.

    Returns
    -------
    dict
        Keys: 'BD71', 'BH91', 'CB05', 'GF97'
        Values: (tag, pars) tuples
    """
    raise NotImplementedError("Profile catalog pending implementation (Milestone 0)")
