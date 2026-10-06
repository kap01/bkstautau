# bkstautau — B0 -> K*0 tau+ tau- with Pythia 8

Self-contained Pythia 8.318 install plus a generator for hard-scatter
`pp -> b bbar` events in which at least one B0 (or anti-B0) decays to
`K*0 tau+ tau-`, with `K*0 -> K+ pi-` and one muonic + one hadronic tau.
Only the anti-B0 is forced; B0 decays generically; exactly one signal
decay per event. Anti-kT gen-jets (FastJet) of both b-hadrons are written,
with and without neutrinos.

```
~/bkstautau/
  src/pythia8318/      Pythia source (configure/build logs in src/)
  src/fastjet-3.5.1/   FastJet source
  install/             Pythia + FastJet installed here (include/, lib/, share/, bin/)
  gen/                 the generator
    gen_bkstautau.cc   main program
    bkstautau.cmnd     Pythia settings card (hard-scatter b bbar, pTHat > 100 GeV)
    bkstautau_incl.cmnd  card for the inclusive QCD sample (HardQCD:all, pTHat > 10 GeV)
    Makefile
    bkstautau.csv      output of the last hard-bbbar run (40k events, seed 12345)
    gen.log            Pythia log of that run
    bkstautau_incl.csv inclusive sample, 20 seeds x 40k events concatenated
    incl_logs/         Pythia logs of the 20 inclusive runs
    bkstautau_incl100.csv  same card with pTHatMin = 100, 20 seeds x 10k events
    incl100_logs/      Pythia logs of those runs
  ana/
    plots.ipynb        generator-level plots from the CSV
  env/                 Python 3.12 env (micromamba, conda-forge): numpy, pandas, matplotlib, jupyter
```

## Setting up from a clone

The repository holds the generator, cards, samples, logs, notebook and
notes; the Pythia/FastJet builds and the Python environment are local and
ignored by git. To recreate them (paths relative to the clone):

```bash
mkdir -p src install
# Pythia 8.318
cd src && curl -O https://pythia.org/download/pythia83/pythia8318.tgz && tar xzf pythia8318.tgz
cd pythia8318 && ./configure --prefix=$PWD/../../install && make -j8 install && cd ../..
# FastJet 3.5.1
cd src && curl -O https://fastjet.fr/repo/fastjet-3.5.1.tar.gz && tar xzf fastjet-3.5.1.tar.gz
cd fastjet-3.5.1 && ./configure --prefix=$PWD/../../install && make -j8 install && cd ../..
# Python environment for the notebook
micromamba create -y -p ./env -c conda-forge python=3.12 numpy pandas matplotlib jupyter ipykernel
./env/bin/python -m ipykernel install --user --name bkstautau --display-name "Python 3.12 (bkstautau)"
```

`CLAUDE.md` is the running log of what has been done and concluded, kept
for Claude Code sessions and for humans alike.

## Build and run

```bash
cd ~/bkstautau/gen
make                        # builds ./gen_bkstautau against ../install
./gen_bkstautau bkstautau.cmnd bkstautau.csv
```

Edit `Main:numberOfEvents`, `Random:seed`, `Beams:eCM` etc. in the card.
Nothing needs to be sourced: the binary carries an rpath to
`../install/lib` and Pythia finds its XML data via the compiled-in prefix.
The Makefile picks up FastJet flags from `../install/bin/fastjet-config`.
~1000 events / 11 s on one core; parallelise by running several seeds.

The inclusive sample is made from `bkstautau_incl.cmnd` (`HardQCD:all`,
`PhaseSpace:pTHatMin = 10.`, ~9 ms/event, ~2 % of events give a signal
decay) by running one copy per seed and concatenating:

```bash
cd ~/bkstautau/gen && mkdir -p incl_logs
for s in $(seq 1 20); do
  sed "s/^Random:seed = .*/Random:seed = $s/" bkstautau_incl.cmnd > /tmp/seed$s.cmnd
  ./gen_bkstautau /tmp/seed$s.cmnd /tmp/seed$s.csv > incl_logs/seed$s.log 2>&1 &
done; wait
head -1 /tmp/seed1.csv > bkstautau_incl.csv; tail -q -n +2 /tmp/seed*.csv >> bkstautau_incl.csv
```

## Plots

The notebook kernel is `Python 3.12 (bkstautau)` (registered in
`~/.local/share/jupyter/kernels/bkstautau`, points at `~/bkstautau/env`), so
VS Code / JupyterLab list it directly.  Nothing needs sourcing; to get
`python`/`jupyter` on PATH in a shell: `micromamba activate ~/bkstautau/env`
(or just call `~/bkstautau/env/bin/python`).

```bash
cd ~/bkstautau/ana
../env/bin/jupyter notebook plots.ipynb         # interactive
../env/bin/jupyter nbconvert --to notebook --execute --inplace plots.ipynb  # batch
```

## What the generator does

* Process: `HardQCD:hardbbbar` (gg -> bb, qq -> bb) at 13.6 TeV with
  `PhaseSpace:pTHatMin = 100.` (sigma ~ 5.4 nb; set to 0 for inclusive bb,
  sigma ~ 0.37 mb, much softer spectra -- the notebook axis ranges are
  tuned for the 100 GeV sample).
* `511:onMode = 2` keeps the full generic B0 decay table open for B0 only;
  `511:addChannel = 3 1.0 0 313 15 -15` adds `anti-B0 -> anti-K*0 tau+ tau-`
  open for anti-B0 only (flat phase space, `meMode 0`). The flavour that
  counts is the one at decay time, i.e. after B0 oscillation: every anti-B0
  at decay is signal, every B0 at decay is generic. Events where both
  b-hadrons are anti-B0 at decay (two signal decays) are dropped, so each
  written event has exactly one signal B and one generically decaying
  opposite-side b-hadron.
* `313:onMode = off` / `313:onIfMatch = 321 211` restricts K*0 -> K+ pi-
  (comment them out to also get K0 pi0 and K0 gamma).
* Tau decays: `setTauModes()` in the C++ sets the tau decay channels so that
  one charge decays only to `mu nu nu` and the other only hadronically
  (electron modes off), using Pythia's particle-only / antiparticle-only
  `onMode` codes 2 / 3. Which charge is muonic is redrawn at random for
  every event, so the sample is charge-symmetric.
* Events without a signal decay are skipped; one CSV row per kept event.
* Production mechanism (`mech`): the b quark of the signal B is found among
  the string partons of the first-fragmented b-hadron ancestor (climbing
  through B* -> B gamma and oscillation copies), its mother chain is
  followed while the mother is the same b quark, and the status of the
  earliest b decides: outgoing from the hardest 2->2 with no incoming b =
  pair creation; incoming to the hardest 2->2 or produced by ISR = flavour
  excitation; produced by FSR = gluon splitting; from a secondary 2->2 =
  MPI. ISR companions of an MPI incoming b are counted as flavour
  excitation (status codes do not distinguish the ISR of the two systems).
* Vertices: the primary vertex (signal-B production vertex), the B decay
  vertex and both tau decay vertices are written in mm, so that impact
  parameters of the muon, K and pi w.r.t. the PV can be computed offline
  (done in the notebook, straight-line tracks, no B field).
* Gen-jets: anti-kT (`Jets:R = 0.4`, `Jets:pTMin = 3.`, CMS ak4GenJets-like)
  clustered with FastJet from all final-state particles, twice: including
  neutrinos (`*jet_*`) and excluding them (`*jetnonu_*`). The signal B and
  the "other" b-hadron (highest-pT weakly decaying b-hadron that is not the
  signal; extras come from g -> bb) are ghost-associated: each is added to
  the clustering with momentum scaled by 1e-18, and the jet that contains
  the ghost is its jet. No eta cut is applied to jets.

## CSV columns

| column | meaning |
|---|---|
| `event` | event number |
| `Bid` | B flavour at production (511 = B0, -511 = anti-B0) |
| `B_osc` | 1 if the B oscillated before decaying (decay flavour is always anti-B0, so `Bid = 511` implies `B_osc = 1`) |
| `B_fd` | true flight distance, production -> decay vertex [mm] |
| `B_*`, `Kst_*` | 4-vectors (px, py, pz, e) [GeV] |
| `K_id`, `K_*`, `pi_id`, `pi_*` | K*0 daughters |
| `taup_*`, `taum_*` | true tau+ / tau- 4-vectors |
| `mu_id`, `mu_*` | the muon (13 = mu- from tau-, -13 = mu+ from tau+) |
| `tauh_id`, `tauh_nprong`, `tauh_*` | hadronic tau: id, charged-prong count, visible 4-vector (all final-state descendants except neutrinos) |
| `q2` | m^2(tau+ tau-) [GeV^2] |
| `oth_id`, `oth_fd`, `oth_*` | other b-hadron: PDG id, flight distance [mm], 4-vector (0/nan if none) |
| `n_oth_bhad` | number of weakly decaying b-hadrons other than the signal B (1 normally; 3, 5 with g -> bb) |
| `njets`, `njets_nonu` | number of jets with pT > `Jets:pTMin`, with / without neutrinos |
| `jets_shared` | 1 if signal B and other b-hadron sit in the same (with-neutrino) jet |
| `sigjet_*`, `sigjetnonu_*` | jet containing the signal B, with / without neutrinos: px, py, pz, e, nconst (nan/0 if the B is in no jet above threshold) |
| `othjet_*`, `othjetnonu_*` | same for the other b-hadron |
| `met_px`, `met_py` | gen-MET: vector sum of all final-state neutrinos (CMS genMetTrue-like) [GeV] |
| `metsig_px`, `metsig_py` | the part of gen-MET from neutrinos in the signal B decay chain (tau decays) |
| `pv_x`, `pv_y`, `pv_z` | primary vertex = production vertex of the signal B [mm] (the origin unless `Beams:allowVertexSpread` is on) |
| `Bdec_x`, `Bdec_y`, `Bdec_z` | signal B decay vertex (= K*0 production vertex) [mm] |
| `taumu_fd`, `taumu_dec_*` | muonic tau: flight distance [mm] and decay vertex (= muon production vertex) [mm] |
| `tauh_fd`, `tauh_dec_*` | hadronic tau: flight distance [mm] and decay vertex (= production vertex of the prongs) [mm] |
| `mech`, `mech_stat` | production mechanism of the signal B's b quark: 1 pair creation, 2 flavour excitation, 3 gluon splitting, 4 MPI, 5 other; and the Pythia status of its earliest b-quark ancestor (see `bMechanism` in the C++) |
| `oth_mech` | same for the other b-hadron (0 if none) |

## Notes / possible next steps

* The B decay is flat phase space. For a physical `q2` / angular distribution
  (e.g. the `BTOSLLBALL` model you use for B+ -> K e e) couple EvtGen via
  `include/Pythia8Plugins/EvtGen.h`, or reweight the CSV in `q2`.
* The warning `TauDecays::decay: unknown correlated tau production` is
  expected: Pythia has no matrix element for B -> K* tau tau, so the tau
  pair is decayed as if from an unpolarised photon (spin correlations
  approximate).
* For HepMC output build HepMC3 and reconfigure Pythia with
  `--with-hepmc3=<path>`; `Pythia8Plugins/HepMC3.h` then gives a writer.
* `PhaseSpace:pTHatMin = 0` is also fine for the massive bb final state if
  you want an inclusive sample.
