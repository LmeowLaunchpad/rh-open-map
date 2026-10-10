# Contributing

1. **Pick a problem** from open-problems/ and open an issue titled with its number (e.g. "P3: ...") saying what you will try.
2. **Check the graveyard first.** If your route is listed there, say why your version escapes the stated reason.
3. **Label every claim**: [proved] (complete argument), [proved, computer-assisted] (certified interval arithmetic, with the script and
   its output), [checked] (numerics), [heuristic], [conjecture]. A statement equivalent to RH is a restatement, not a step.
4. **Kill tests.** Any proposed mechanism must fail in the fake worlds (Davenport-Heilbronn, zeta(s-a)zeta(s+a), Epstein, Gonek's
   zeta_X) where RH fails; if it also "works" there, it carries no RH content.
5. **Reproducibility.** Certificates: python-flint (Arb) with all error terms bounded; give the exact command line. Numerics: state the
   basis size and precision and show convergence (the window problems need basis N >= 2 x log x).
6. **Negative results are welcome.** A clean refutation of a route, with the reason, goes into the graveyard and is credited.
7. **Verification.** Headline claims get an independent re-run before they are marked verified.
