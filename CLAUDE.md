# CLAUDE.md — bkstautau

Running notes for Claude Code sessions in this repo. Keep this file current:
whenever you change the generator, card, notebook or produce a new sample,
add a dated entry to the **Work log** and, if a result changes our
understanding, to **Conclusions**. The README.md is the user-facing manual
(build/run instructions, CSV column table); do not duplicate it here, link
to it.

## What this study is

Generator-level (Pythia 8.318 + FastJet 3.5.1, no detector simulation) study
of `pp -> b bbar` at 13.6 TeV where one b-hadron, the anti-B0, is forced to
decay `anti-B0 -> anti-K*0 (-> K- pi+) tau+ tau-`, with one tau -> mu nu nu
and the other tau -> hadrons. The opposite-side b-hadron decays generically.
We look at:

* kinematics (pT, eta) of the signal-side visible decay products: K, pi,
  mu, visible hadronic tau;
* the anti-kT R=0.4 gen-jets of the signal B and of the recoiling b-hadron
  (with and without neutrinos), gen-MET;
* all of this **within the CMS acceptance**: the analysis-level requirement
  is |eta| < 2.4 on the four visible signal objects (tracker / muon
  coverage). No acceptance cut is applied at generation; the cut lives in
  the notebook so the full sample stays available.

Physics context: B -> K* tau tau is a lepton-flavour-universality test
channel; the sample is meant to tell us what the signal looks like at
CMS (soft K/pi/mu from a ~90 GeV B inside a b-jet, large neutrino energy
loss, MET aligned with the signal jet) and what the recoiling b-jet looks
like for tagging/triggering.

## Layout (details in README.md)

```
gen/gen_bkstautau.cc   generator main: decay forcing, tau-mode steering, ghost-tagged jets, CSV writer
gen/bkstautau.cmnd     Pythia card: hard-scatter bbbar, pTHat > 100 (decay table edits, jet params, seed)
gen/bkstautau_incl.cmnd  card for the inclusive sample: HardQCD:all, pTHat > 10
gen/Makefile           builds against ../install (rpath set; nothing to source)
gen/bkstautau.csv      hard-bbbar sample, 15105 rows (one per kept event), from 40k generated events; 93 columns
gen/gen.log            Pythia log of that run
gen/bkstautau_incl.csv inclusive sample, 16361 rows from 20 seeds x 40k events (seeds 1-20); same columns
gen/incl_logs/         Pythia logs of the 20 inclusive runs
gen/bkstautau_incl100.csv  inclusive sample at pTHat > 100 (HardQCD:all), 10896 rows from 20 seeds (101-120) x 10k
gen/incl100_logs/      its Pythia logs (card = bkstautau_incl.cmnd with pTHatMin = 100, numberOfEvents = 10000)
ana/plots.ipynb        all plots; kernel "Python 3.12 (bkstautau)" = ./env
install/               Pythia8 + FastJet install prefix
src/                   Pythia/FastJet sources + build logs
env/                   micromamba Python 3.12 env (numpy, pandas, matplotlib, jupyter)
venv/                  stale Python 3.9 venv, not used; ignore
```

## How to work here

* Build/run: `cd gen && make && ./gen_bkstautau bkstautau.cmnd bkstautau.csv`.
  ~11 s per 1000 generated events on one core; a 40k run takes ~7-8 min.
  Use a different `Random:seed` per parallel run and concatenate CSVs
  (recipe in README). The machine has 128 cores; 20 parallel 40k runs of
  the inclusive card (~9 ms/event) finish in ~7 min wall time.
* Python: always `./env/bin/python` (or the registered kernel). Do not use
  `venv/`.
* Execute the notebook headless with
  `cd ana && ../env/bin/jupyter nbconvert --to notebook --execute --inplace plots.ipynb`.
* Quick checks from a shell: `env/bin/python -c "import pandas as pd; ..."`
  on `gen/bkstautau.csv`. pT/eta/phi are not stored; derive them from px,py,pz.
* The generator is a single C++ file; keep it that way unless it grows a lot.
  New observables go in as new CSV columns (append at the end of the header
  and the row, and add them to the README column table).
* `Jets:R`, `Jets:pTMin` are user settings registered in the C++, read from
  the card.
* The repo is https://github.com/kap01/bkstautau (remote `origin`, SSH).
  Commit and push after each session's changes.
* Do not commit/regenerate `bkstautau.csv` casually: 8 MB, and the notebook
  numbers below refer to this exact sample (seed 12345, 40k events). The run
  is reproducible: re-running with the same seed after adding output columns
  gave bit-identical values in all pre-existing columns (checked 2026-10-05),
  so new columns can be added without changing the sample.

## Generator design decisions (why things are the way they are)

* **Decay forcing via the decay table, not a user hook.** `511:onMode = 2`
  keeps generic decays for the particle (B0) only; `511:addChannel = 3 ...`
  opens the signal channel for the antiparticle (anti-B0) only. The flavour
  that matters is the one at decay time, i.e. after B0 oscillation, so
  every anti-B0 at decay is signal. Events where *both* b-hadrons are
  anti-B0 at decay (two signal decays, ~10 % of signal events) are
  dropped so each row has exactly one signal B and one generic
  opposite-side b-hadron.
* **Tau modes toggled per event in C++** (`setTauModes`) using Pythia's
  particle-only/antiparticle-only `onMode` codes 2/3, redrawn at random
  each event so the muon charge is 50/50. Electron modes are off.
* **B -> K* tau tau is flat phase space** (`meMode 0`): Pythia has no
  matrix element for it. q2 and angular distributions are *not* physical;
  reweight in q2 or couple EvtGen if that matters.
* **Jets:** anti-kT R=0.4, pT > 3 GeV, on all final-state particles,
  clustered twice (with/without neutrinos). b-hadrons are ghost-associated
  (momentum x 1e-18 added to the input list). No eta cut on jets.
* **"Other" b-hadron** = highest-pT weakly decaying b-hadron that is not
  the signal B. `n_oth_bhad` > 1 signals extra g -> bb splittings.
* **Production mechanism** (`mech`, `oth_mech`, `bMechanism()` in the C++)
  from the status of the earliest b-quark ancestor of the b-hadron, after
  climbing through B* -> B gamma / oscillation copies to the hadron made by
  string fragmentation: 1 pair creation (status 23, no incoming b),
  2 flavour excitation (21, 41-44), 3 gluon splitting (51), 4 MPI (31, 33),
  5 other (beam-remnant b, 61/63: the PDF companion of a flavour-excited
  b that stayed in the remnant, FE-like). Known limitation: ISR companions
  of an MPI incoming b are labelled FE (7 % of the inclusive sample has
  mech=FE with oth_mech=MPI or vice versa).
* **Inclusive sample** (`bkstautau_incl.cmnd`): `HardQCD:all` includes
  the massive gg/qq -> bbbar processes, so pair creation is in. pTHat > 10
  GeV chosen (user request) to be as inclusive as practical; sigma(gen)
  = 8.3 mb, 2.0 % of generated events give a signal decay.
* **Vertices** are written in mm: PV = signal-B production vertex (always
  the origin: `Beams:allowVertexSpread` is off), B decay vertex, and the
  decay vertices + flight distances of both taus. Impact parameters are
  computed in the notebook (straight lines, no B field, CMS `dxy` sign
  convention), not in the C++.
* **pTHatMin = 100 GeV** so the b-jets are at CMS-trigger-like pT. The
  cross section is then ~5.4 nb vs ~0.37 mb inclusive; the sample is *not*
  representative of the inclusive bb spectrum and notebook axis ranges
  assume this hard sample.

## Conclusions so far (sample: seed 12345, 40k generated, 15105 kept)

Unless stated otherwise the numbers refer to the hard-bbbar sample; the
inclusive sample has its own block at the end of this section.

Numbers are from `ana/plots.ipynb` executed on the current CSV, cross-checked
directly from the CSV on 2026-10-05.

**Sample composition**
* 40000 generated; 27563 have a B0/anti-B0; 16715 have >=1 signal decay;
  1610 have two (dropped); 15105 written. sigma(gen) = 5.40e-6 mb.
* 14.2 % of signal B's were produced as B0 and oscillated. Muon charge is
  50.6 % mu+ / 49.4 % mu-. Hadronic tau prongs: 1-prong 74.6 %, 3-prong
  25.1 %, 5-prong 0.3 %.
* Opposite-side b-hadron: B+ 39 %, B0 32 %, Bs 9 %, anti-Lambda_b 3 %
  (right-sign: the signal anti-B0 = (b dbar) carries the b, so the recoil
  carries the bbar). Wrong-sign B- 11 %, anti-B0 2 %, anti-Bs 2 %,
  Lambda_b 1 % (~16 % total) come from oscillated signal B's (14.2 %) plus
  g -> bb. 91 % of events have exactly one other b-hadron; 9 % have 3+.

**CMS acceptance**
* Each visible object (K, pi, mu, tau_h) is within |eta| < 2.4 in ~87.7 %
  of events; all four together in **85.9 %** (12978 / 15105). The objects
  are highly correlated in eta because they come from one ~90 GeV B, so
  the combined acceptance is barely below the single-object one.
* Median pT in acceptance: K 11 GeV, pi 6 GeV, mu 9.5 GeV, tau_h visible
  20 GeV. The muon is above 3 GeV in 81 % and above 5 GeV in 70 % of
  in-acceptance events. The pion is the softest object.
* Signal B median pT 89 GeV, other b-hadron 91 GeV (set by pTHatMin).

**Jets**
* Signal B is inside a jet (pT > 3 GeV) in 98.7 % of events, the other
  b-hadron in 99.9 %; both in the same jet in 0.09 %. Median
  dR(b-hadron, jet axis) ~0.014. ~18.6 jets/event above 3 GeV.
* Signal jet: median pT 114 GeV with neutrinos, 81 GeV without, vs B pT
  89 GeV. Median response jet/B = 1.23 with nu, 0.88 without; the
  no-neutrino jet keeps a median **72 %** of the with-neutrino jet pT.
  The signal b-jet loses roughly a quarter of its energy to the four tau
  neutrinos.
* Other b-jet: median pT 116 GeV with / 110 GeV without neutrinos,
  response 1.23 / 1.18, no-nu/with-nu median 1.00 (neutrinos only in the
  ~20 % semileptonic decays). Response > 1 for both jets because the jet
  picks up the rest of the b fragmentation + underlying event.
* 88 % of signal jets and 89 % of other-b jets are within |eta| < 2.4.

**Transverse balance of the two b-jets** (14780 events with two distinct
b-jets; numbers identical within 1 % when both jets are required in |eta| < 2.4)
* With neutrinos the pair balances: pT(sig jet)/pT(other jet) 16/50/84 % =
  0.72 / 0.99 / 1.35, asymmetry A = (s-o)/(s+o) mean -0.02, width 0.25.
* Without neutrinos the signal jet is systematically softer: ratio
  0.48 / 0.74 / 1.11, A mean -0.15 (a 15 % shift, same width 0.26). This is
  the tau-neutrino loss; the recoil jet is almost unaffected.
* Mixed case, signal jet *with* neutrinos vs recoil jet *without*: ratio
  0.75 / 1.05 / 1.51, A mean +0.02, width 0.26. The median overshoots 1 and
  the upper tail grows (1.35 -> 1.51) because the recoil jet is now the one
  missing energy: 37 % of recoil jets contain > 1 GeV of neutrino pT
  (semileptonic b/c decays, median 15 GeV when present).
* "No-neutrino signal jet + signal-chain gen-MET" is numerically the same
  as the mixed case: the four tau neutrinos land inside the R=0.4 signal
  jet cone (neutrino content of the jet equals MET_sig exactly in 80 % of
  events, within 5 GeV in 97.5 %). So at gen level, correcting the signal
  jet with MET_sig is equivalent to clustering with neutrinos.
* |dphi| between the two jets is the same with/without neutrinos (median
  2.93 rad, 67 % above 2.8): neutrinos change the magnitude, not the
  direction. The remaining imbalance and acoplanarity come from ISR/FSR
  and the third jet, as the with-neutrino width (0.25) shows.
* Vector sum |pT(sig)+pT(oth)| / HT median 0.18 with nu, 0.25 without.
* Plots in `ana/plots.ipynb`, section "Transverse balance of the two b-jets".

**MET**
* Gen-MET median 28 GeV, mean 31 GeV; the signal-chain neutrinos alone give
  median 31 GeV and dominate (|MET_sig|/|MET| median 1.00). MET > 50 GeV in
  15.6 % of events, > 100 GeV in 0.9 %. The MET points along the signal
  jet (and away from the other b-jet).

**B0 proper time (open point resolved 2026-10-05)**
* Mean true proper time 1.417 ps vs Pythia tau(B0) = 1.530 ps (0.4587 mm)
  is *not* a bug. The sample is selected on flavour at decay, so unmixed
  B's follow e^{-t/tau}(1+cos dm t)/2 (mean 0.71 tau = 1.09 ps; observed
  1.08) and mixed ones e^{-t/tau}(1-cos dm t)/2 (mean 2.26 tau = 3.44 ps;
  observed 3.47). With the unbiased mixed fraction chi_d = 18.8 % the two
  would average to exactly tau, but the sample has only 14.2 % mixed:
  the veto on events with two signal decays preferentially drops mixed
  signal B's (their recoil is a wrong-sign B that is anti-B0 at decay
  ~80 % of the time vs ~20 % for the unmixed case). 14.2 % mixed gives an
  expected mean of 1.423 ps, matching the 1.417 observed. Plot and check
  in the notebook, "Proper time split by oscillation".

**Impact parameters w.r.t. the PV** (straight-line, true momenta, mm)
* Tau flight distance: median 1.86 mm (muonic) / 1.89 mm (hadronic), mean
  4.6 mm, tails to tens of mm (tau0 = 87 um, beta*gamma ~ 20).
* Muon |dxy| w.r.t. PV: 16/50/84 % = 0.038 / 0.193 / 0.71 mm, mean 0.42
  mm; |dz| median 0.38 mm (0.31 in acceptance); 3D IP median 0.35 mm.
  Numbers are the same in and out of acceptance.
* The muon is about twice as displaced as the K (|dxy| median 0.088 mm)
  and pi (0.122 mm) from the B vertex. The tau flight alone (muon |dxy|
  w.r.t. the B decay vertex) contributes median 0.027 mm, mean 0.07 mm:
  most of the muon displacement is the B flight, the tau adds a tail and
  smears the lifetime sign (96.5 % of muons have positive lifetime-signed
  dxy w.r.t. the B direction, vs 100 % for K and pi).
* Fraction of in-acceptance muons with pT > 3 GeV above |dxy| thresholds:
  20 um 90 %, 50 um 78 %, 100 um 64 %, 200 um 44 %, 500 um 18 %, 1 mm 6 %
  (kaon with pT > 3: 79 / 63 / 46 / 28 / 9 / 2 %). So a ~50 um |dxy|
  requirement on the muon keeps ~80 % of signal at gen level; a displaced
  muon trigger needing > 0.5-1 mm keeps < 20 %.
* Plots in `ana/plots.ipynb`, section "Impact parameters w.r.t. the
  primary vertex" (overlays of |dxy|, lifetime-signed dxy, |dz|, 3D IP for
  mu/K/pi; muon decomposition B-flight vs tau-flight; tau flight distances;
  |dxy| vs muon pT; cumulative threshold curves).

**Inclusive QCD sample vs hard bbbar (production mechanism)**
* Hard-bbbar sample is 95.2 % pair creation (PC); the rest (FE 1.2 %, GS
  2.9 %, MPI 0.6 %) are extra b pairs from the shower/MPI in which the
  signal B is soft (median pT 6-9 GeV) and the recoil is the 100 GeV
  hard b.
* Inclusive sample (800k generated, 29864 with a B0, 17810 with a signal
  decay, 1449 with two, 16361 kept): PC 5.6 %, flavour excitation 43.5 %,
  gluon splitting 33.1 %, MPI 16.4 %, beam remnant 1.4 %. So the hard-bbbar
  sample represents ~6 % of inclusive signal events at this pTHat. The
  sample is soft: signal B median pT 6 GeV (PC 12, FE 6.3, GS 5.4, MPI
  5.7), 45 % of events have all four signal objects in |eta| < 2.4 (86 %
  in the hard sample); only 120 events have B pT > 30 GeV.
* Recoil b-hadron: median pT 6.2 GeV; within |eta| < 2.4 in 53 % (PC 66 %,
  MPI 45 %); in a jet above 3 GeV in only 66 % (PC 96 %, FE 71 %, GS 62 %,
  MPI 51 %) vs 99.9 % in the hard sample. Both b-hadrons in distinct jets
  in 40 % of events (hard: 98.5 %).
* Topology: median |dphi|(B, other b) 1.82 rad overall (PC 2.74, FE 1.86,
  GS 1.28, MPI 2.31) vs 2.92 in the hard sample; the back-to-back peak is
  specific to PC. Median dR 2.44 (GS 1.78). Only 2 % of pairs are within
  dR < 0.4 / in the same jet (GS 5 %): at pTHat 10 the gluon-splitting
  pairs are not collinear because the opening angle ~ m_bb/pT is large;
  expect this to change at CMS-trigger pT.
* Transverse balance (both jets, inclusive): ratio sig/other 16/50/84 % =
  0.42 / 0.93 / 2.13 with neutrinos (hard: 0.72 / 0.99 / 1.35), 0.40 /
  0.86 / 1.76 without; A mean -0.03 width 0.36 with nu (hard -0.02 / 0.25),
  -0.08 / 0.33 without. The pair is still balanced on average but the
  spread is 45 % wider and the upper tail (signal jet harder than the
  recoil jet) is much longer, from FE/GS/MPI where the second b is soft.
  Per-mechanism shapes in the notebook. Requiring B pT > 30 GeV in the
  inclusive sample selects a strongly asymmetric configuration (ratio
  median 3.1, A = +0.44): a hard signal B in a pTHat ~ 10 event means the
  recoil b is soft, so the balance seen in the hard sample is a property
  of the PC topology, not of signal events with a hard B in general.
* Plots in `ana/plots.ipynb`, section "Inclusive QCD sample".

**Inclusive sample at pTHat > 100 GeV** (`bkstautau_incl100.csv`, 200k
generated, 10896 signal, sigma(gen) = 1.43e-3 mb, 5.4 % signal/event; not
yet in the notebook, numbers from a shell check)
* Mechanism fractions of signal events: gluon splitting 61 %, flavour
  excitation 29 %, MPI 7 %, **pair creation 2.4 %**, remnant 0.5 %.
  Cross-check: sigma x P(signal) of the hard-bbbar sample (5.40 ub x
  0.378) over the inclusive one (1.43 mb x 0.054) = 2.6 %, consistent.
  So the hard-bbbar sample is ~2.5 % of signal events at this pTHat.
* Selecting a hard signal B raises the PC share but it stays small: B pT
  > 30 GeV: PC 7 %, FE 34 %, GS 59 %; > 50: 9 / 40 / 51; > 80: 13 / 50 /
  38. Signal jet pT > 100: PC 8 %, FE 33 %, GS 58 %.
* Signal B median pT: PC 88 GeV, FE 17, GS 17, MPI 6. In GS the two
  b-hadrons are within dR < 0.4 / in the same jet in 23-24 % of events
  (vs 2 % at pTHat 10): at this scale gluon splitting does give the
  collinear pair.
* **Mechanism vs b-hadron pT** (both b-hadrons per event counted, each
  sample used only well above its pTHat cut where the cut does not sculpt
  the mix): pair creation is ~15-17 % of b-hadrons for every pT bin from
  10 to 250 GeV, flavour excitation ~50-57 %, gluon splitting ~25-38 %,
  MPI 10 % at 10-15 GeV and negligible above 20 GeV. The 2.4 % PC share of
  all signal events at pTHat > 100 is low only because that cut admits
  soft GS/FE b's while PC b's sit at ~pTHat; at fixed b pT the PC share is
  flat. The hard-bbbar cross section as a fraction of all 2->2 scatters is
  also flat (0.32 % at pTHat > 10 with sigma = 26.4 ub, 0.38 % at pTHat >
  100).
* CMS B-parking-like tag (muon pT > 9, |eta| < 1.5, tag B ~ 30 GeV): the
  other b-hadron is in |eta| < 2.4 with pT > 5 GeV in ~80-85 % of events,
  median pT 15-25 GeV, within dR < 0.4 of the tag in ~5 % at 30 GeV (18
  events only in the pTHat > 10 sample) and 17 % at ~85 GeV. Tag-and-probe
  works because the probe exists and is in acceptance, not because it
  balances the tag.
* **Selecting pair creation in data with topology cuts** (pTHat > 100
  sample, two distinct no-nu b-jets with pT > 30 in |eta| < 2.4: 9 % of
  signal events, PC purity 20 %): |dphi| > 2.8 -> purity 52 %, PC
  efficiency 51 % (of all PC events); adding |A| < 0.3 -> 59 % / 40 %;
  |dphi| > 3.0 & |A| < 0.2 -> 60 % / 16 %. The residual is FE ~20-25 %
  and GS ~20 %, so back-to-back b-jet pairs are at best ~60 % pair
  creation. Same purities at pTHat > 10 (jets > 10 GeV) but with only
  2 % of signal events passing the two-b-jet requirement.

**Neutrino collinearity in the tau decays** (shell check 2026-10-06, not
in the notebook)
* The nu pair from tau -> mu nu nu follows the muon to ~1.5 x m_tau/p_tau:
  angle(mu, nu pair) median 0.05 rad (dR 0.10) at tau pT ~ 32 GeV (hard
  sample), 0.12 rad (dR 0.46) at tau pT ~ 6 GeV (inclusive pTHat > 100),
  0.25 rad (dR 1.2) at tau pT ~ 2 GeV (pTHat > 10). Same for tau_h vis vs
  its neutrino. The neutrino momentum perpendicular to the muon direction
  is boost-independent, median 1.7 GeV (16/84 %: 0.75 / 3.6 GeV).
* Direction is not the problem, magnitude is: the nu pair carries a
  median 67 % of the tau momentum (16/84 %: 38 / 90 %); tau_h visible
  keeps 68 % (38 / 90 %).
* The two taus are nearly parallel: angle(tau+, tau-) median 0.02 rad at
  B ~ 90 GeV (dR(mu, tau_h vis) 0.08), 0.05 rad at pTHat > 100. The
  di-tau collinear approximation (two fractions from two MET components)
  is therefore ill-conditioned here; MET only measures the sum along the
  common direction. Constraints have to come from vertexing: B line of
  flight (PV -> K* vertex), B mass, tau masses, and the tau_h decay vertex
  for 3-prong decays (25 %); with a 3-prong tau_h the system is
  overconstrained by one, with a 1-prong it is one short.
* Signal-jet pT > 50 GeV selection (no-nu jet, |eta| < 2.4, four objects
  in acceptance): keeps ~0.1 % of all signal events (10 of 16361 in the
  pTHat > 10 sample, so +-30 %), 26.5 % of the pTHat > 100 sample, 76 % of
  the hard sample. In the selection the muonic tau has pT ~ 20-24 GeV,
  angle(mu, nu pair) median 0.08 rad (dR 0.13; 16/84 %: 0.05 / 0.34),
  tau opening angle 0.03 rad, nu fraction unchanged (67 %). Mechanism mix
  of the selected events: GS 62 %, FE 31 %, PC 7 % (pTHat > 100 sample).
* **Expected yield in 1 ab-1 with no-nu signal jet pT > 100, |eta| < 2.4**
  (pTHat > 100 sample): anti-B0-at-decay per generated event = (11952 +
  1056) / 200000 = 0.065 (dropped two-signal events counted twice, since
  with a real BR both B's contribute); fraction of signal events passing
  the jet cut 10.1 % (9.7 % with the four objects in acceptance);
  sigma(cand, BR = 1) = 1.43 ub x 0.065 x 0.101 = 9.4 nb -> 9.4e9 anti-B0
  with such a jet per ab-1. Times BR(B0 -> K*0 tau tau)_SM ~ 1e-7, BR(K*0
  -> K+ pi-) = 2/3, 2 x BR(tau -> mu nu nu) x BR(tau -> had) = 0.226:
  **~140 events / ab-1** before trigger and reconstruction (~135 with the
  four objects in |eta| < 2.4; ~380 for jet pT > 50; ~270 if the cut were
  on the with-neutrino jet). Uncertainties: LO Pythia cross section
  (factor ~2), SM BR (~20 %), pTHat > 100 slightly undercounts jets > 100
  from softer scatters.
* Same chain started from the pTHat > 10 sample (user request): 8.3 mb x
  0.0241 anti-B0/event = 200 ub = 2.0e14 anti-B0 per ab-1; the 4-object
  acceptance is 45 % overall but ~100 % for events with a 100 GeV jet, so
  it must not be factorised. Fraction with no-nu jet > 100, |eta| < 2.4,
  in acceptance: 2 events / 16361 = 1.2e-4 (Poisson 68 %: 0.4-2.8e-4);
  power-law extrapolation of the jet tail above 30 GeV (index 3.2):
  7.1e-5; from the pTHat > 100 sample 9.4 nb / 200 ub = 4.7e-5. All
  consistent; the pTHat > 100 value (1099 events) is the one to use ->
  ~140 events / ab-1 at SM BR. A pTHat > 50 bridge sample would remove the
  stitching assumption.

**q2**
* Flat-phase-space envelope: min 12.6 GeV2 = (2 m_tau)2, mean 15.7, max
  20.6 GeV2. The upper edge exceeds (m_B - m_K*)2 = 19.2 GeV2 because the
  K* mass is Breit-Wigner distributed. Not a physical q2 spectrum.

**Open points / things to check**
* Impact parameters are straight-line from true vertices and momenta: no
  B-field curvature, no track smearing. For CMS-like numbers apply a dxy
  resolution (~20-50 um for these pT) and compute significances.
* The two-signal veto biases the mixed fraction (14.2 % vs chi_d 18.8 %)
  and hence the proper-time distribution; if a flavour-unbiased sample is
  needed, keep both-signal events and pick one, or weight.
* The flat-phase-space decay model should be replaced or reweighted before
  any q2- or angle-dependent conclusion.
* The inclusive sample is at pTHat > 10 GeV, far below any CMS b-jet /
  tau trigger threshold; for a CMS-relevant inclusive picture generate
  `HardQCD:all` with a higher pTHat (or pTHat bins) and select on the
  signal B / signal jet pT, where gluon splitting is expected to give
  collinear b pairs inside one jet. The mechanism column makes this easy
  to split.
* No pT thresholds are applied in the "acceptance" definition yet; a
  realistic CMS selection would add e.g. muon pT > 3-5 GeV, track pT > ~0.5-1
  GeV, and a tau_h / jet pT threshold.

## Work log

* **2026-10-03 / 10-04** (before Claude sessions): Pythia 8.318 and FastJet
  3.5.1 built into `install/`; generator written; 40k-event run with seed
  12345 produced `gen/bkstautau.csv` (15105 rows) and `gen/gen.log`;
  `ana/plots.ipynb` written and executed; README.md written.
* **2026-10-05**: Created this CLAUDE.md. Read the generator, card, notebook
  and run log; recomputed all headline numbers directly from the CSV and
  confirmed they match the notebook outputs. Recorded design decisions,
  conclusions and open points above. No code changes.
* **2026-10-05**: Studied the transverse balance of signal vs recoil b-jet
  with and without neutrinos (user request). Added a notebook section
  (two cells after "Neutrino energy fraction") and executed the notebook.
  Conclusions recorded above. No generator changes.
* **2026-10-05**: Added the mixed case (signal jet with nu, recoil jet
  without) to the balance plots (user request); notebook re-executed.
  Found it is identical to "no-nu signal jet + MET_sig" because the tau
  neutrinos are always inside the signal jet cone; corrected the earlier
  (wrong) explanation of the widened ratio in Conclusions.
* **2026-10-05**: Generator: added 14 vertex columns at the end of the CSV
  (`pv_*`, `Bdec_*`, `taumu_fd`, `taumu_dec_*`, `tauh_fd`, `tauh_dec_*`;
  helpers `writeV3`, `decayVertex`). Rebuilt and re-ran 40k events with
  seed 12345: all 76 old columns bit-identical to the previous CSV, same
  15105 rows and summary counts. README column table updated. Notebook:
  new section "Impact parameters w.r.t. the primary vertex" (user request:
  muon IP w.r.t. the B's PV; K and pi added for comparison) with dxy,
  lifetime-signed dxy, dz, 3D IP, B-flight vs tau-flight decomposition and
  threshold-efficiency curves; new subsection "Proper time split by
  oscillation" that resolves the open point on the B0 lifetime (selection
  on decay flavour + two-signal veto, not a vertex bug). Notebook executed
  headless; conclusions above updated.
* **2026-10-05**: Production mechanism and inclusive sample (user
  question: LHCb-style gluon splitting vs hard-scatter bbbar). Generator:
  `bMechanism()` + columns `mech`, `mech_stat`, `oth_mech` (93 columns
  now); summary line in the log. New card `bkstautau_incl.cmnd`
  (`HardQCD:all`, pTHat > 10 GeV, user choice). Ran 20 seeds x 40k in
  parallel -> `gen/bkstautau_incl.csv` (16361 rows), logs in
  `gen/incl_logs/`; re-ran the hard sample with seed 12345 (old columns
  bit-identical again). Notebook: new section "Inclusive QCD sample: b
  production mechanism and b-jet topology" (composition table per
  mechanism, topology plots, balance per mechanism vs hard sample),
  inserted before the MET section; executed headless. README: card, parallel
  recipe, new columns. Conclusions and open points updated.
* **2026-10-05**: User asked for the pair-creation fraction at pTHat ~ 100.
  Generated 20 x 10k inclusive events with `bkstautau_incl.cmnd` and
  pTHatMin = 100 (seeds 101-120) -> `gen/bkstautau_incl100.csv`, logs in
  `gen/incl100_logs/`. Mechanism fractions recorded in Conclusions. No
  notebook changes yet.
* **2026-10-06**: Put the project on GitHub (user request):
  https://github.com/kap01/bkstautau, branch `main`, pushed over SSH
  (no `gh` on this machine; the repo was created by the user on the web).
  `.gitignore` excludes `install/`, `src/`, `env/`, `venv/`, `.vscode/`
  and the compiled binary; the samples, logs and executed notebook are
  tracked. README gained a "Setting up from a clone" section. Commit
  with `git add -A && git commit && git push` after changes; keep the
  CSVs in the repo only while they stay at the current ~30 MB total.
