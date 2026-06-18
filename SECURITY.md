# Security policy

## Supported versions

Security fixes are issued against the current development release
(`0.8.3`, on the 0.8.x line) and the latest CRAN release (when present).
Older released versions are not back-patched.

| Version | Supported |
|----|----|
| 0.8.x | yes (current development line) |
| 0.3.x – 0.7.x | best-effort (local release tarballs only; no public release) |
| \< 0.3 | no |

## Reporting a vulnerability

Please **do not file a public GitHub issue** for security problems.
Instead, email the maintainer privately:

> **<max.moldovan@adelaide.edu.au>**

Use a clear subject line such as
`[flexyBayes SECURITY] <one-line description>`.

We commit to:

1.  Acknowledge receipt within 7 days.
2.  Provide a substantive reply within 30 days, including either a patch
    / fix plan or an explanation of why the issue is not actionable.
3.  Coordinate a fix and (where appropriate) public disclosure on a
    timeline that protects users — typically a 90-day embargo from first
    acknowledgement.

If you do not receive a response within the windows above, please follow
up with a public GitHub issue (without sensitive details) so we can
confirm the report was received.

## Scope

The package is pure R glue plus generated greta / INLA code. Likely
areas of security interest:

- inputs that can reach [`eval()`](https://rdrr.io/r/base/eval.html) /
  [`parse()`](https://rdrr.io/r/base/parse.html) paths during formula
  parsing or code generation;
- file paths supplied to `R CMD INSTALL` or to vignette compilation
  scripts;
- handling of user-supplied numeric data in code-generation paths.

Reports outside the package — e.g. issues in upstream `greta`, `INLA`,
or `brms` — should be directed to those projects’ security channels.
