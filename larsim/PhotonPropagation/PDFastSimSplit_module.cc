// LArSoft libraries
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardataobj/Simulation/OpDetBacktrackerRecord.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/Simulation/SimPhotons.h"
#include "larsim/PhotonPropagation/PhotonVisibilityService.h"
#include "larsim/PhotonPropagation/PropagationTimeModel.h"
#include "larsim/PhotonPropagation/ScintTimeTools/ScintTime.h"
#include "larsim/Simulation/LArG4Parameters.h"

#include "nurandom/RandomUtils/NuRandomService.h"

// Art libraries
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Utilities/make_tool.h"
#include "canvas/Utilities/Exception.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "CLHEP/Random/RandPoissonQ.h"
#include "CLHEP/Units/SystemOfUnits.h"

#include <iostream>
#include <memory>
#include <utility>
#include <vector>

namespace phot {
  class PDFastSimSplit : public art::EDProducer {
  public:
    explicit PDFastSimSplit(fhicl::ParameterSet const&);
    void produce(art::Event&) override;
    void AddOpDetBTR(std::vector<sim::OpDetBacktrackerRecord>& opbtr,
                     std::vector<int>& ChannelMap,
                     const sim::OpDetBacktrackerRecord& btr) const;

  private:
    art::ServiceHandle<PhotonVisibilityService const> fPVS;
    const bool fDoFastComponent;
    const bool fDoSlowComponent;
    const bool fIncludePropTime;
    const size_t fNOpChannels;
    const art::InputTag fSimTag;
    CLHEP::HepRandomEngine& fPhotonEngine;
    std::unique_ptr<CLHEP::RandPoissonQ> fRandPoissPhot;
    CLHEP::HepRandomEngine& fScintTimeEngine;
    std::unique_ptr<ScintTime> fScintTime; // Tool to retrive timing of scintillation
    std::unique_ptr<PropagationTimeModel> fPropTimeModel;
  };

  //......................................................................
  PDFastSimSplit::PDFastSimSplit(fhicl::ParameterSet const& pset)
    : art::EDProducer{pset}
    , fPVS(art::ServiceHandle<PhotonVisibilityService const>())
    , fDoFastComponent(pset.get<bool>("DoFastComponent", true))
    , fDoSlowComponent(pset.get<bool>("DoSlowComponent", true))
    , fIncludePropTime(pset.get<bool>("IncludePropTime", false))
    , fNOpChannels(fPVS->NOpChannels())
    , fSimTag{pset.get<art::InputTag>("SimulationLabel")}
    , fPhotonEngine(art::ServiceHandle<rndm::NuRandomService> {}
                      ->createEngine(*this, "HepJamesRandom", "photon", pset, "SeedPhoton"))
    , fRandPoissPhot(std::make_unique<CLHEP::RandPoissonQ>(fPhotonEngine))
    , fScintTimeEngine(art::ServiceHandle<rndm::NuRandomService>()->createEngine(*this,
                                                                                 "HepJamesRandom",
                                                                                 "scinttime",
                                                                                 pset,
                                                                                 "SeedScintTime"))
    , fScintTime{art::make_tool<phot::ScintTime>(pset.get<fhicl::ParameterSet>("ScintTimeTool"))}
  {
    mf::LogInfo("PDFastSimSplit") << "Initializing PDFastSimSplit." << std::endl;
    fhicl::ParameterSet VUVTimingParams;
    fhicl::ParameterSet VISTimingParams;
    // Validate configuration options
    if (fIncludePropTime &&
        !pset.get_if_present<fhicl::ParameterSet>("VUVTiming", VUVTimingParams)) {
      throw art::Exception(art::errors::Configuration)
        << "Propagation time simulation requested, but VUVTiming not specified."
        << "\n";
    }
    if (fIncludePropTime && false/*fStoreReflected*/ &&
        !pset.get_if_present<fhicl::ParameterSet>("VISTiming", VISTimingParams)) {
      throw art::Exception(art::errors::Configuration)
        << "Reflected light propagation time simulation requested, but VISTiming not specified."
        << "\n";
    }
    // Initialise the Scintillation Time
    fScintTime->initRand(fScintTimeEngine);

    // propagation time model
    if (fIncludePropTime) {
      fPropTimeModel = std::make_unique<PropagationTimeModel>(
        VUVTimingParams, VISTimingParams, fScintTimeEngine, false/*fStoreReflected*/);
    }

    mf::LogInfo("PDFastSimSplit") << "Use Lite Photon." << std::endl;
    produces<std::vector<sim::SimPhotonsLite>>();
    produces<std::vector<sim::OpDetBacktrackerRecord>>();

    mf::LogInfo("PDFastSimSplit") << "PDFastSimSplit Initialization finish" << std::endl;
  }

  //......................................................................
  void PDFastSimSplit::produce(art::Event& event)
  {
    std::vector<int> PDChannelToSOCMapDirect(fNOpChannels, -1);  // Where each OpChan is.

    // SimPhotonsLite
    auto phlit = std::make_unique<std::vector<sim::SimPhotonsLite>>();
    auto opbtr = std::make_unique<std::vector<sim::OpDetBacktrackerRecord>>();
    auto& dir_phlitcol(*phlit);
    dir_phlitcol.resize(fNOpChannels);
    for (unsigned int i = 0; i < fNOpChannels; i++) {
      dir_phlitcol[i].OpChannel = i;
    }

    art::Handle<std::vector<sim::SimEnergyDeposit>> edepHandle;
    if (!event.getByLabel(fSimTag, edepHandle)) {
      mf::LogWarning("PDFastSimSplit") << "PDFastSimSplit Module Cannot getByLabel: " << fSimTag;
      return;
    }

    auto const& edeps = edepHandle;
    std::cout << edepHandle->size() << std::endl;

    size_t n_split = 1;
    for (auto s : fPVS->GetSplitVals()) n_split *= s;
    std::cout << "Split vals: " << n_split << std::endl;

    std::vector<std::vector<const sim::SimEnergyDeposit*>> split_edeps(n_split);
    for (auto const & edepi : *edeps) {
      auto const & prt = edepi.MidPoint();
      auto coords = fPVS->GetSplitRegion(prt);
      size_t index = fPVS->GetTreeIndex(coords[0], coords[1], coords[2]);
      //std::array<size_t, 3> index_to_region = fPVS->TreeIndexToRegion(index);
      split_edeps[index].push_back(&edepi);
      //std::cout << prt.X() << " " << prt.Y() << " " << prt.Z() << " " <<
      //             coords[0] << " " << coords[1] << " " << coords[2] << " " <<
      //             index_to_region[0] << " " << index_to_region[1] << " " << index_to_region[2] << " " <<
      //             std::endl;
    }

    for (auto & edep_v : split_edeps) {
      std::cout << edep_v.size() << std::endl;
    }


    for (size_t i = 0; i < split_edeps.size(); ++i) {
      auto & edep_v = split_edeps[i];

      std::array<size_t, 3> region = fPVS->TreeIndexToRegion(i);
      fPVS->LoadSplitLibrary(region[0], region[1], region[2]);

      for (auto const * edepi : edep_v) {

        int nphot_fast = edepi->NumFPhotons();
        int nphot_slow = edepi->NumSPhotons();
        if (!((nphot_fast > 0 && fDoFastComponent) ||
              (nphot_slow > 0 && fDoSlowComponent))) continue;

        auto const& prt = edepi->MidPoint();

        MappedCounts_t Visibilities = fPVS->GetAllVisibilities(prt);

        if (!Visibilities) {
          mf::LogWarning("PDFastSimSplit")
            << "There is no entry in the PhotonLibrary for this position in space. Position: "
            << edepi->MidPoint() << "\n Move to next point";

          mf::LogWarning("PDFastSimSplit") <<
            "Index: " << i << " Region: " << "(" << region[0] << ", " <<
            region[1] << ", " << region[2] << ")" << std::endl;

          continue;
        }

        int trackID = edepi->TrackID();
        int nphot = edepi->NumPhotons();
        double edeposit = edepi->Energy() / nphot;
        double pos[3] = {edepi->MidPointX(), edepi->MidPointY(), edepi->MidPointZ()};

        geo::Point_t const ScintPoint = {pos[0], pos[1], pos[2]};

        // loop through direct photons then reflected photons cases
        for (unsigned int channel = 0; channel < fNOpChannels; ++channel) {

          double visibleFraction = Visibilities[channel];
          if (visibleFraction < 1e-9) continue; // voxel is not visible at this optical channel

          // number of detected photons
          int ndetected_fast =
            (nphot_fast > 0) ? fRandPoissPhot->fire(nphot_fast * visibleFraction) : 0;
          int ndetected_slow =
            (nphot_slow > 0) ? fRandPoissPhot->fire(nphot_slow * visibleFraction) : 0;
          if (!((ndetected_fast > 0 && fDoFastComponent) ||
                (ndetected_slow > 0 && fDoSlowComponent)))
            continue;

          // calculate propagation times if included, does not matter whether fast or slow photon
          std::vector<double> transport_time;
          if (fIncludePropTime) {
            transport_time.resize(ndetected_fast + ndetected_slow);
            fPropTimeModel->propagationTime(transport_time, ScintPoint, channel, 0);
          }

          // SimPhotonsLite case
          sim::OpDetBacktrackerRecord tmpbtr(channel);
          if (ndetected_fast > 0 && fDoFastComponent) {
            for (int i = 0; i < ndetected_fast; ++i) {
              // calculate the time at which each photon is seen
              double dtime = edepi->StartT() + fScintTime->fastScintTime();
              if (fIncludePropTime) dtime += transport_time[i];
              int time = static_cast<int>(std::round(dtime));
              ++dir_phlitcol[channel].DetectedPhotons[time];
              tmpbtr.AddScintillationPhotons(trackID, time, 1, pos, edeposit);
            }
          }
          if ((ndetected_slow > 0) && fDoSlowComponent) {
            for (int i = 0; i < ndetected_slow; ++i) {
              // calculate the time at which each photon is seen
              double dtime = edepi->StartT() + fScintTime->slowScintTime();
              if (fIncludePropTime) dtime += transport_time[ndetected_fast + i];
              int time = static_cast<int>(std::round(dtime));
              ++dir_phlitcol[channel].DetectedPhotons[time];
              tmpbtr.AddScintillationPhotons(trackID, time, 1, pos, edeposit);
            }
          }
          AddOpDetBTR(*opbtr, PDChannelToSOCMapDirect, tmpbtr);
        }
      }
    }

    event.put(move(phlit));
    event.put(move(opbtr));

    return;
  }

  //......................................................................
  void PDFastSimSplit::AddOpDetBTR(std::vector<sim::OpDetBacktrackerRecord>& opbtr,
                                 std::vector<int>& ChannelMap,
                                 const sim::OpDetBacktrackerRecord& btr) const
  {
    int iChan = btr.OpDetNum();
    if (ChannelMap[iChan] < 0) {
      ChannelMap[iChan] = opbtr.size();
      opbtr.emplace_back(std::move(btr));
    }
    else {
      size_t idtest = ChannelMap[iChan];
      auto const& timePDclockSDPsMap = btr.timePDclockSDPsMap();
      for (auto const& timePDclockSDP : timePDclockSDPsMap) {
        for (auto const& sdp : timePDclockSDP.second) {
          double xyz[3] = {sdp.x, sdp.y, sdp.z};
          opbtr.at(idtest).AddScintillationPhotons(
            sdp.trackID, timePDclockSDP.first, sdp.numPhotons, xyz, sdp.energy);
        }
      }
    }
  }

} // namespace

DEFINE_ART_MODULE(phot::PDFastSimSplit)
