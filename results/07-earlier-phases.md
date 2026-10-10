# 07. Earlier phases: the common-signal class and its wall

Starting point: github.com/openai/math (719 AI-generated manuscripts), in particular family 003 (a quasi-RH zero-free strip
Re s > 7/8 for a class of L-functions, assembled in Lean through one proposition).

## The common-signal class and its law [proved/checked]
- Mechanism: a probe with a direct bound plus a bootstrap against sup beta* (the supremum of zero real parts); 7/8 = 11/16 + 3/16.
- Law: sigma_0 = 1 - h/N + b/(2N) = 11/12 - l/4 at N = 6, with caps 7/8 (exact if a reflected-modulus cap M' >= 1/2 is structural),
  13/15, 5/6; at N = 4: 13/16, 4/5, 3/4. The Hasse identity (mu times quadratic Gauss sum = cubic Gauss sum = Patterson theta
  coefficient) forces sextic density; Gauss sums carry Moebius fractions mu^0, mu^{1/3}, mu^{1/2}, mu^{2/3} at orders 2, 3, 4, 6.
- Four requirements any escape needs: (K) a square-root-dense common kernel of the row characters, (M) a Moebius source decoupled from
  the pole variable, (F) Moebius hidden in bounded automorphic phases, (S) self-referential rows.
- Theorems: kernel lemma (power-residue families: kernel = N-th powers x S-units); trace-function theorem ((M) and (S) are incompatible
  for l-adic trace functions on A^1 at every order and rank; Lefschetz-sign principle); principal-count law beta* <= e - g_P;
  integrality/divisibility engines closed at Mertens scale; Landau-Deligne positivity class void over Q.
- Jacobi-sum Moebius source: sum_{x mod c} (x(1-x)/c) = mu(c) chi_{-4}(c), failing at (S) (rows are Kloosterman-Selberg zeta functions
  with tempered poles on the principal row's line).

## By-products
- Knauf (1998) Ramanujan conjecture refuted with integer certificates (d = 29 graphs; d = 433 lumped chain); Knauf (2013)
  reconciliation hope contradicted to n = 390625. [checked, certificates in the original run]
- Finding 004 (published in LmeowLaunchpad/Legendary): every real abelian number field of class number one is Euclidean, conditional
  on two openai/math manuscripts (families 029 and 011).

## Why this phase ended
Every common-signal variant hits the same law; the named missing objects (a Moebius-twisted quadratic theta series with a direct
bound, a Chebotarev-aware Archimedean height, a positive singularity other than s = 1, inverse-L rows at order 2, a depth-multiplying
spectral object) were not found. The attack then moved to Weil positivity on windows (results 01-06).
