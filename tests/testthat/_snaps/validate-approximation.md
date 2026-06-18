# print.fb_approximation_validation snapshot

    Code
      print(v)
    Output
      <fb_approximation_validation>
        scheme:    low_rank_smooth
        threshold: Frobenius capture >= 0.99
        verdict:   FAIL
          XX s(x): rank 4/9  capture 0.95  (bias bound 0.05)
        fallback:  Re-fit with a higher rank, or drop the approximation to fit the exact smooth. If the basis dimension k is large only because of a high default k in s(x, k = ...), reducing k directly is the exact alternative to truncation.

