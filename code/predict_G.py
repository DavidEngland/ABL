def predict_G(session, features: list[float]) -> Optional[float]:
    """Run ONNX model; return scalar G or None if session is None or error.

    Parameters
    ----------
    session : onnxruntime.InferenceSession or None
        Loaded ONNX session. If None, returns None.
    features : list[float]
        Single feature vector (flat).

    Returns
    -------
    float or None
        Predicted G value (scalar). Returns None if model not available or fails.
    """
    if session is None:
        return None

    try:
        # onnxruntime expects float32 2D input [N=1, F]
        import numpy as _np
        X = _np.asarray([features], dtype=_np.float32)

        inp_name = session.get_inputs()[0].name
        out_name = session.get_outputs()[0].name

        out = session.run([out_name], {inp_name: X})
        y = out[0]

        # Expect shape (1,1) or (1,) depending on model
        val = float(y.reshape(-1)[0])
        return val

    except Exception:
        return None