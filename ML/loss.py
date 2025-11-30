def loss(pred_Ric, true_Ric, pred_events, true_events, raw_unclamped, monotone_penalty_w):
    mae = np.mean(np.abs(pred_Ric - true_Ric))
    precision, recall = compute_precision_recall(pred_events, true_events)
    f1 = 2 * precision * recall / max(precision + recall, 1e-6)
    event_term = 1.0 - f1

    out_of_bounds = np.mean(np.maximum(0, RI_MIN - raw_unclamped) + np.maximum(0, raw_unclamped - RI_MAX))
    curvature = np.mean(np.abs(np.gradient(np.gradient(raw_unclamped))))
    monotone_pen = monotone_penalty_w * compute_monotone_penalty()

    return mae + event_term + 0.1 * out_of_bounds + 0.01 * curvature + monotone_pen
