Inflection modeling notes:

- Piecewise linear: breakpoint x0 estimated via non-linear least squares; AIC with k=4 parameters.
- Logistic: L/(1+exp(-k(x-x0))) + b; inflection at x0; AIC with k=4.
- We require at least 3 distinct age points per tissue; otherwise report NA.
