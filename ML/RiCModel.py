class RiCModel:
    def __init__(self, coeff, meta):
        self.coeff = coeff
        self.meta = meta  # ModelVersion, DatasetVersion, Thresholds

    def infer(self, state):
        raw = self._raw(state)
        ric = clamp(raw)
        flags = self._constraints_check(state, raw, ric)
        self._log(state, raw, ric, flags)
        if flags["invalid_input"] or flags["hard_violation"]:
            return self.meta["fallback_Ri_c"]
        return ric

    # Implement _raw, _constraints_check, and _log with guards and telemetry
