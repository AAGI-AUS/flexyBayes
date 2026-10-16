# Security policy

## Supported versions

Security fixes are made on the 0.18.x line (0.18.1 and later 0.18.x
releases) and on the latest CRAN release once there is one. Older versions
are not patched.

| Version | Supported |
|---------|-----------|
| 0.18.x  | yes |
| 0.10.x to 0.17.x | best effort |
| 0.3.x to 0.9.x | best effort (never released publicly) |
| < 0.3   | no |

## Reporting a vulnerability

Please **do not open a public GitHub issue** for a security problem.
Email the maintainer privately instead:

> **max.moldovan@adelaide.edu.au**

Use a subject line such as `[flexyBayes SECURITY] <one-line description>`.

The maintainer will:

1. Acknowledge the report within 7 days.
2. Reply in substance within 30 days, with a fix plan or an explanation of
   why the issue will not be fixed.
3. Coordinate the fix and any public disclosure, normally with a 90-day
   embargo from the first acknowledgement.

If you have had no response within these times, open a public GitHub issue
at <https://github.com/AAGI-AUS/flexyBayes/issues>, without sensitive
details, asking whether the report arrived.

## Scope

flexyBayes is written in R and generates INLA and Stan model code. Areas of
likely security interest are:

- input that can reach `eval()` or `parse()` during formula parsing or code
  generation
- file paths passed to `R CMD INSTALL` or to vignette build scripts
- user-supplied numeric data in code generation.

Problems in INLA, brms or Stan themselves should be reported to those
projects.
