// gen_bkstautau.cc
//
// pp -> b bbar (hard QCD) with anti-B0 -> anti-K*0 tau+ tau- forced via the
// decay table (B0 keeps its generic decays), K*0 -> K pi, and the tau pair
// forced to one muonic + one hadronic decay (which tau takes which role is
// chosen at random per event).  Only events with exactly one such signal
// decay are kept, one CSV row per event.
//
// For every kept event, anti-kT gen-jets are clustered with FastJet from all
// final-state particles, once including neutrinos and once without, and the
// jets containing the signal B and the other (opposite-side) b-hadron are
// identified by ghost association.  Gen-MET (all neutrinos) is written too.
// Production/decay vertices [mm] of the signal B and the two taus are written
// so that impact parameters of the muon, K and pi w.r.t. the primary vertex
// can be computed offline.
//
// Usage:  ./gen_bkstautau [card.cmnd] [output.csv]

#include "Pythia8/Pythia.h"
#include "fastjet/ClusterSequence.hh"
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace Pythia8;

namespace {

// First descendant of `i` whose |id| is in `ids`, searching recursively.
int findDesc(const Event& ev, int i, const std::vector<int>& ids) {
  for (int d : ev[i].daughterList()) {
    for (int id : ids) if (ev[d].idAbs() == id) return d;
    int r = findDesc(ev, d, ids);
    if (r >= 0) return r;
  }
  return -1;
}

// K*0 -> K+ pi- is forced in the card, but accept neutral modes too.
const std::vector<int> kaonIds = {321, 311};
const std::vector<int> pionIds = {211, 111};

bool isNeutrino(const Particle& p) {
  int id = p.idAbs();
  return id == 12 || id == 14 || id == 16;
}

// Number of charged final-state descendants of `i` (tau prong count).
int nCharged(const Event& ev, int i) {
  int n = 0;
  for (int d : ev[i].daughterList()) {
    if (ev[d].isFinal()) { if (ev[d].isCharged()) ++n; }
    else n += nCharged(ev, d);
  }
  return n;
}

// Sum of final-state, non-neutrino descendants of `i` (visible tau "jet").
Vec4 visibleP4(const Event& ev, int i) {
  Vec4 sum;
  for (int d : ev[i].daughterList()) {
    if (ev[d].isFinal()) { if (!isNeutrino(ev[d])) sum += ev[d].p(); }
    else sum += visibleP4(ev, d);
  }
  return sum;
}

// Sum of final-state neutrino descendants of `i`.
Vec4 neutrinoP4(const Event& ev, int i) {
  Vec4 sum;
  for (int d : ev[i].daughterList()) {
    if (ev[d].isFinal()) { if (isNeutrino(ev[d])) sum += ev[d].p(); }
    else sum += neutrinoP4(ev, d);
  }
  return sum;
}

// Does tau `i` decay directly to a muon?
bool isMuonicTau(const Event& ev, int i) {
  for (int d : ev[i].daughterList()) if (ev[d].idAbs() == 13) return true;
  return false;
}

// Is particle `i` a B0/anti-B0 that decays directly to K*0 tau+ tau- ?
bool isSignalB(const Event& ev, int i) {
  if (ev[i].idAbs() != 511) return false;
  std::vector<int> dl = ev[i].daughterList();
  if (dl.size() != 3) return false;
  bool kst = false, taup = false, taum = false;
  for (int d : dl) {
    if (ev[d].idAbs() == 313) kst  = true;
    if (ev[d].id()    == -15) taup = true;   // tau+
    if (ev[d].id()    ==  15) taum = true;   // tau-
  }
  return kst && taup && taum;
}

// Open-beauty hadron (B mesons incl. excited states, b baryons); excludes
// b quarks, b diquarks and bb-bar quarkonia.
bool isBHadron(int id) {
  id = std::abs(id) % 10000;              // strip excitation digits
  if (id < 100) return false;
  if (id < 1000) return (id / 100) % 10 == 5 && (id / 10) % 10 != 5;
  return (id / 1000) % 10 == 5 && (id / 10) % 10 != 0;
}

// Weakly decaying b-hadron: a b-hadron none of whose daughters is one
// (this also skips intermediate record copies of the same hadron).
bool isWeakBHadron(const Event& ev, int i) {
  if (!isBHadron(ev[i].id())) return false;
  for (int d : ev[i].daughterList()) if (isBHadron(ev[d].id())) return false;
  return true;
}

// Production mechanism of the b quark inside b-hadron `iHad`, from the
// ancestry of the b quark with the hadron's flavour.  Walk up the mother
// chain as long as the mother is a b quark of the same id; the status of
// the earliest such b classifies it (Pythia status codes):
//   1 pair creation      earliest b is outgoing from the hardest 2->2 (23)
//                        with no incoming b
//   2 flavour excitation earliest b is incoming to the hardest 2->2 (21),
//                        an ISR incoming b (41, 42) or the ISR companion
//                        from a backward g -> b bbar branching (43, 44)
//   3 gluon splitting    earliest b made in the final-state shower (51-59)
//   4 MPI                earliest b in a secondary 2->2 (31-39)
//   5 other              beam remnant etc.
// `stat` receives the raw status of the earliest b for cross-checks.
int bMechanism(const Event& ev, int iHad, int& stat) {
  int id = ev[iHad].id();
  bool baryon = (std::abs(id) % 10000) / 1000 != 0;
  int bId = (baryon ? 5 : -5) * (id > 0 ? 1 : -1);   // b quark in the hadron
  // Climb through b-hadron ancestors (B* -> B gamma, oscillation copies,
  // record copies) to the hadron made directly by string fragmentation.
  int had = iHad;
  for (;;) {
    int up = -1;
    for (int m : ev[had].motherList()) if (m > 0 && isBHadron(ev[m].id())) { up = m; break; }
    if (up < 0) break;
    had = up;
  }
  // The b quark among that hadron's mothers (string partons).
  int cur = -1;
  for (int m : ev[had].motherList()) if (m > 0 && ev[m].id() == bId) { cur = m; break; }
  if (cur < 0)   // fall back to any b
    for (int m : ev[had].motherList()) if (m > 0 && ev[m].idAbs() == 5) { cur = m; break; }
  if (cur < 0) { stat = 0; return 5; }
  for (;;) {
    int next = -1;
    for (int m : ev[cur].motherList()) if (m > 0 && ev[m].id() == ev[cur].id()) { next = m; break; }
    if (next < 0) break;
    cur = next;
  }
  stat = ev[cur].statusAbs();
  if (stat == 23) {
    for (int m : ev[cur].motherList()) if (m > 0 && ev[m].idAbs() == 5) return 2;
    return 1;
  }
  if (stat == 21 || (stat >= 41 && stat <= 49)) return 2;
  if (stat >= 51 && stat <= 59) return 3;
  if (stat >= 31 && stat <= 39) return 4;
  return 5;
}

// Flight distance [mm] of particle `i`: production vertex to the production
// vertex of its first daughter.
double flightDistance(const Event& ev, int i) {
  int d = ev[i].daughter1();
  if (d <= 0) return NAN;
  return (ev[d].vProd() - ev[i].vProd()).pAbs();
}

// Steer the tau decay table so that one charge decays only muonically and
// the other only hadronically.  DecayChannel::onMode codes:
//   0 = off, 1 = on, 2 = on for tau- only, 3 = on for tau+ only.
// Electron modes are switched off entirely.  Pythia re-reads the channel
// list at every decay, so this can be called between events.
void setTauModes(Pythia& pythia, bool tauMinusIsMuonic) {
  ParticleDataEntryPtr tau = pythia.particleData.particleDataEntryPtr(15);
  int muMode  = tauMinusIsMuonic ? 2 : 3;
  int hadMode = tauMinusIsMuonic ? 3 : 2;
  for (int i = 0; i < tau->sizeChannels(); ++i) {
    DecayChannel& ch = tau->channel(i);
    bool hasMu = false, hasE = false;
    for (int j = 0; j < ch.multiplicity(); ++j) {
      int id = std::abs(ch.product(j));
      if (id == 13) hasMu = true;
      if (id == 11) hasE  = true;
    }
    ch.onMode(hasE ? 0 : (hasMu ? muMode : hadMode));
  }
}

// --- jets -----------------------------------------------------------------

// Ghost tags stored in PseudoJet::user_index (real particles carry their
// event-record index, which is >= 0).
const int ghostSig = -1, ghostOth = -2;

fastjet::PseudoJet makeGhost(const Particle& p, int tag) {
  const double eps = 1e-18;
  fastjet::PseudoJet g(p.px() * eps, p.py() * eps, p.pz() * eps, p.e() * eps);
  g.set_user_index(tag);
  return g;
}

// Index of the jet whose constituents include the ghost `tag`, or -1.
int ghostJet(const std::vector<fastjet::PseudoJet>& jets, int tag) {
  for (size_t j = 0; j < jets.size(); ++j)
    for (const auto& c : jets[j].constituents())
      if (c.user_index() == tag) return j;
  return -1;
}

// Number of real (non-ghost) constituents.
int nConst(const fastjet::PseudoJet& jet) {
  int n = 0;
  for (const auto& c : jet.constituents()) if (c.user_index() >= 0) ++n;
  return n;
}

// --- output helpers -------------------------------------------------------

void writeP4(std::ofstream& o, const Vec4& p) {
  o << "," << p.px() << "," << p.py() << "," << p.pz() << "," << p.e();
}
void writeP4(std::ofstream& o, const Particle& p) { writeP4(o, p.p()); }
void writeNaN4(std::ofstream& o) { o << ",nan,nan,nan,nan"; }

// Vertex columns x,y,z [mm].
void writeV3(std::ofstream& o, const Vec4& v) {
  o << "," << v.px() << "," << v.py() << "," << v.pz();
}

// Decay vertex of particle `i` = production vertex of its first daughter.
Vec4 decayVertex(const Event& ev, int i) {
  int d = ev[i].daughter1();
  if (d <= 0) return Vec4(NAN, NAN, NAN, NAN);
  return ev[d].vProd();
}

// Jet columns: px,py,pz,e,nconst (nan/0 if no jet).
void writeJet(std::ofstream& o, const std::vector<fastjet::PseudoJet>& jets,
  int j) {
  if (j < 0) { writeNaN4(o); o << ",0"; return; }
  const auto& jet = jets[j];
  o << "," << jet.px() << "," << jet.py() << "," << jet.pz() << "," << jet.e()
    << "," << nConst(jet);
}

} // namespace

int main(int argc, char* argv[]) {

  std::string card = argc > 1 ? argv[1] : "bkstautau.cmnd";
  std::string outf = argc > 2 ? argv[2] : "bkstautau.csv";

  Pythia pythia;
  // User settings for the jet clustering, readable from the card.
  pythia.settings.addParm("Jets:R",     0.4, true, true, 0.1, 1.5);
  pythia.settings.addParm("Jets:pTMin", 3.0, true, false, 0.0, 0.0);
  pythia.readFile(card);
  int    nEvents = pythia.mode("Main:numberOfEvents");
  int    nAbort  = pythia.mode("Main:timesAllowErrors");
  double jetR    = pythia.parm("Jets:R");
  double jetPtMin = pythia.parm("Jets:pTMin");
  if (!pythia.init()) return 1;
  Event& ev = pythia.event;

  fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, jetR);
  std::cout << "\n" << jetDef.description()
            << ", pT > " << jetPtMin << " GeV\n\n";

  std::ofstream out(outf);
  out << "event,Bid,B_osc,B_fd"
      << ",B_px,B_py,B_pz,B_e"
      << ",Kst_px,Kst_py,Kst_pz,Kst_e"
      << ",K_id,K_px,K_py,K_pz,K_e"
      << ",pi_id,pi_px,pi_py,pi_pz,pi_e"
      << ",taup_px,taup_py,taup_pz,taup_e"
      << ",taum_px,taum_py,taum_pz,taum_e"
      << ",mu_id,mu_px,mu_py,mu_pz,mu_e"
      << ",tauh_id,tauh_nprong,tauh_px,tauh_py,tauh_pz,tauh_e"
      << ",q2"
      << ",oth_id,oth_fd,oth_px,oth_py,oth_pz,oth_e,n_oth_bhad"
      << ",njets,njets_nonu,jets_shared"
      << ",sigjet_px,sigjet_py,sigjet_pz,sigjet_e,sigjet_nconst"
      << ",sigjetnonu_px,sigjetnonu_py,sigjetnonu_pz,sigjetnonu_e,sigjetnonu_nconst"
      << ",othjet_px,othjet_py,othjet_pz,othjet_e,othjet_nconst"
      << ",othjetnonu_px,othjetnonu_py,othjetnonu_pz,othjetnonu_e,othjetnonu_nconst"
      << ",met_px,met_py,metsig_px,metsig_py"
      << ",pv_x,pv_y,pv_z"
      << ",Bdec_x,Bdec_y,Bdec_z"
      << ",taumu_fd,taumu_dec_x,taumu_dec_y,taumu_dec_z"
      << ",tauh_fd,tauh_dec_x,tauh_dec_y,tauh_dec_z"
      << ",mech,mech_stat,oth_mech"
      << "\n";

  long nMech[6] = {0, 0, 0, 0, 0, 0};
  long nGen = 0, nWithB0 = 0, nSig = 0, nMulti = 0, nCand = 0, nBad = 0,
       nWrongFlav = 0, nNoOther = 0, nNoSigJet = 0, nNoOthJet = 0, iAbort = 0;
  for (int iEv = 0; iEv < nEvents; ++iEv) {

    // Randomly assign the muonic role to tau- or tau+ for this event.
    setTauModes(pythia, pythia.rndm.flat() < 0.5);

    if (!pythia.next()) {
      if (++iAbort < nAbort) continue;
      std::cerr << "Too many event-generation errors, stopping.\n";
      break;
    }
    ++nGen;

    bool hasB0 = false;
    std::vector<int> sig;
    for (int i = 0; i < ev.size(); ++i) {
      if (ev[i].idAbs() == 511) hasB0 = true;
      if (isSignalB(ev, i)) sig.push_back(i);
    }
    if (hasB0) ++nWithB0;
    if (sig.empty()) continue;
    ++nSig;
    // Exactly one signal B per event: drop the rare events where both
    // b-hadrons are anti-B0 at decay time (via oscillation).
    if (sig.size() > 1) { ++nMulti; continue; }
    {
      int i = sig[0];
      ++nCand;
      // Sanity check: the decaying B must be an anti-B0 (-> anti-K*0).
      if (ev[findDesc(ev, i, {313})].id() != -313) ++nWrongFlav;

      int iKst = findDesc(ev, i, {313});
      int iTp = -1, iTm = -1;
      for (int d : ev[i].daughterList()) {
        if (ev[d].id() == -15) iTp = d;
        if (ev[d].id() ==  15) iTm = d;
      }
      int iK  = findDesc(ev, iKst, kaonIds);
      int iPi = findDesc(ev, iKst, pionIds);

      // Which tau is muonic?  (Both or neither would mean the forcing
      // failed; count it and skip.)
      bool muP = isMuonicTau(ev, iTp), muM = isMuonicTau(ev, iTm);
      if (muP == muM) { ++nBad; continue; }
      int iTmu = muP ? iTp : iTm;
      int iTh  = muP ? iTm : iTp;
      int iMu  = findDesc(ev, iTmu, {13});

      // Oscillation flag (daughters of an oscillated B carry status 92/94)
      // and flight distance [mm] from production to decay vertex.
      int st = ev[iKst].statusAbs();
      int osc = (st == 92 || st == 94) ? 1 : 0;
      double fd = (ev[iKst].vProd() - ev[i].vProd()).pAbs();

      Vec4 q = ev[iTp].p() + ev[iTm].p();

      // The other b-hadron: highest-pT weakly decaying b-hadron that is not
      // the signal B (extra ones come from gluon splitting g -> bb).
      int iOth = -1, nOth = 0;
      for (int k = 0; k < ev.size(); ++k) {
        if (k == i || !isWeakBHadron(ev, k)) continue;
        ++nOth;
        if (iOth < 0 || ev[k].pT() > ev[iOth].pT()) iOth = k;
      }
      if (iOth < 0) ++nNoOther;

      // Gen-jets: all final-state particles (+ ghosts), with and without
      // neutrinos.
      // Gen-MET: vector sum of all final-state neutrinos (like CMS
      // genMetTrue), and the part from the signal B decay chain alone.
      std::vector<fastjet::PseudoJet> inAll, inNoNu;
      Vec4 met;
      for (int k = 0; k < ev.size(); ++k) {
        if (!ev[k].isFinal()) continue;
        fastjet::PseudoJet pj(ev[k].px(), ev[k].py(), ev[k].pz(), ev[k].e());
        pj.set_user_index(k);
        inAll.push_back(pj);
        if (isNeutrino(ev[k])) met += ev[k].p(); else inNoNu.push_back(pj);
      }
      Vec4 metSig = neutrinoP4(ev, i);
      inAll.push_back(makeGhost(ev[i], ghostSig));
      inNoNu.push_back(makeGhost(ev[i], ghostSig));
      if (iOth >= 0) {
        inAll.push_back(makeGhost(ev[iOth], ghostOth));
        inNoNu.push_back(makeGhost(ev[iOth], ghostOth));
      }
      fastjet::ClusterSequence csAll(inAll, jetDef), csNoNu(inNoNu, jetDef);
      auto jetsAll  = fastjet::sorted_by_pt(csAll.inclusive_jets(jetPtMin));
      auto jetsNoNu = fastjet::sorted_by_pt(csNoNu.inclusive_jets(jetPtMin));

      int jSig   = ghostJet(jetsAll,  ghostSig);
      int jSigNN = ghostJet(jetsNoNu, ghostSig);
      int jOth   = iOth >= 0 ? ghostJet(jetsAll,  ghostOth) : -1;
      int jOthNN = iOth >= 0 ? ghostJet(jetsNoNu, ghostOth) : -1;
      if (jSig < 0) ++nNoSigJet;
      if (iOth >= 0 && jOth < 0) ++nNoOthJet;
      int shared = (jSig >= 0 && jSig == jOth) ? 1 : 0;

      out << iEv << "," << ev[i].id() << "," << osc << "," << fd;
      writeP4(out, ev[i]);
      writeP4(out, ev[iKst]);
      out << "," << (iK  >= 0 ? ev[iK].id()  : 0);
      if (iK  >= 0) writeP4(out, ev[iK]);  else writeNaN4(out);
      out << "," << (iPi >= 0 ? ev[iPi].id() : 0);
      if (iPi >= 0) writeP4(out, ev[iPi]); else writeNaN4(out);
      writeP4(out, ev[iTp]);
      writeP4(out, ev[iTm]);
      out << "," << ev[iMu].id();
      writeP4(out, ev[iMu]);
      out << "," << ev[iTh].id() << "," << nCharged(ev, iTh);
      writeP4(out, visibleP4(ev, iTh));
      out << "," << q.m2Calc();
      // other b-hadron
      if (iOth >= 0) {
        out << "," << ev[iOth].id() << "," << flightDistance(ev, iOth);
        writeP4(out, ev[iOth]);
      } else { out << ",0,nan"; writeNaN4(out); }
      out << "," << nOth;
      // jets
      out << "," << jetsAll.size() << "," << jetsNoNu.size() << "," << shared;
      writeJet(out, jetsAll,  jSig);
      writeJet(out, jetsNoNu, jSigNN);
      writeJet(out, jetsAll,  jOth);
      writeJet(out, jetsNoNu, jOthNN);
      out << "," << met.px() << "," << met.py()
          << "," << metSig.px() << "," << metSig.py();
      // Vertices [mm]: primary vertex = production vertex of the signal B
      // (the origin unless Beams:allowVertexSpread is on), B decay vertex,
      // and the decay vertices of the muonic and hadronic tau (= production
      // vertices of the muon and of the tau_h prongs), plus tau flight
      // distances.
      writeV3(out, ev[i].vProd());
      writeV3(out, ev[iKst].vProd());
      out << "," << flightDistance(ev, iTmu);
      writeV3(out, decayVertex(ev, iTmu));
      out << "," << flightDistance(ev, iTh);
      writeV3(out, decayVertex(ev, iTh));
      // Production mechanism of the signal B and of the other b-hadron.
      int mechStat = 0, othStat = 0;
      int mech = bMechanism(ev, i, mechStat);
      int othMech = iOth >= 0 ? bMechanism(ev, iOth, othStat) : 0;
      out << "," << mech << "," << mechStat << "," << othMech;
      ++nMech[mech];
      out << "\n";
    }
  }

  pythia.stat();
  std::cout << "\n==== bkstautau summary ====\n"
            << "generated events          : " << nGen    << "\n"
            << "events with a B0/B0bar    : " << nWithB0 << "\n"
            << "events with >=1 signal B  : " << nSig    << "\n"
            << "  of which 2 signal B (dropped): " << nMulti << "\n"
            << "events written (1 signal) : " << nCand - nBad << "\n"
            << "candidates w/o mu+had taus: " << nBad    << "\n"
            << "signal B not anti-B0 (BUG): " << nWrongFlav << "\n"
            << "events w/o other b-hadron : " << nNoOther << "\n"
            << "signal B not in a jet     : " << nNoSigJet << "\n"
            << "other b not in a jet      : " << nNoOthJet << "\n"
            << "signal B mechanism (PC/FE/GS/MPI/other): " << nMech[1] << " / "
            << nMech[2] << " / " << nMech[3] << " / " << nMech[4] << " / "
            << nMech[5] << "\n"
            << "sigma(gen) [mb]           : " << pythia.info.sigmaGen() << "\n"
            << "output                    : " << outf    << "\n";
  return 0;
}
