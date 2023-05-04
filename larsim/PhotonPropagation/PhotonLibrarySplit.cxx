#include "canvas/Utilities/Exception.h"

#include "cetlib_except/exception.h"
#include "larsim/PhotonPropagation/PhotonLibrarySplit.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "TTree.h"
#include "TString.h"
#include "TFile.h"

#include <cassert>

namespace phot {
  std::string const PhotonLibrarySplit::OpChannelBranchName = "OpChannel";

  //------------------------------------------------------------
  PhotonLibrarySplit::PhotonLibrarySplit(sim::PhotonVoxelDef vox_def,
                                         size_t n_split_x,
                                         size_t n_split_y,
                                         size_t n_split_z,
                                         size_t n_op_chans)
    : fNOpChannels(n_op_chans),
      fVoxelDef(vox_def),
      fNSplitX(n_split_x),
      fNSplitY(n_split_y),
      fNSplitZ(n_split_z) {
    //Each sub-library will be made of
    //(nx/nx_split)*(ny/ny_split)*(nz/nz_split) voxels
    //
    //This will be used to resize the vectors later
    auto steps = fVoxelDef.GetSteps();
    fNx = (steps[0]/fNSplitX);
    fNy = (steps[1]/fNSplitY);
    fNz = (steps[2]/fNSplitZ);
    fNVoxels = (fNx*fNy*fNz);
  }
  //------------------------------------------------------------
 
  size_t PhotonLibrarySplit::GetTreeIndex(size_t i, size_t j, size_t k) const {
    return i + fNSplitX*j + fNSplitY*fNSplitX*k;
  }

  std::array<size_t, 3> PhotonLibrarySplit::TreeIndexToRegion(size_t i) const {
    size_t ix = i % fNSplitX;
    size_t iy = ((i - ix)/fNSplitX) % fNSplitY;
    size_t iz = (i - ix - iy*fNSplitX)/(fNSplitX*fNSplitY);
    return {ix, iy, iz};
  }

  void PhotonLibrarySplit::OpenFile(std::string LibraryFile) {
    mf::LogInfo("PhotonLibrarySplit") << "Opening photon library input file: "
                                 << LibraryFile.c_str() << std::endl;

    fInputFile = TFile::Open(LibraryFile.c_str());
    fInputTrees.resize(fNSplitX*fNSplitY*fNSplitZ);

    TString library_name;
    try {
      for (size_t i = 0; i < fNSplitX; ++i) {
        for (size_t j = 0; j < fNSplitY; ++j) {
          for (size_t k = 0; k < fNSplitZ; ++k) {
            library_name.Form("PhotonLibraryData_%zu_%zu_%zu", i, j, k);
            size_t index = GetTreeIndex(i, j, k);
            fInputTrees[index] = (TTree*)fInputFile->Get(library_name);
            //fInputTrees[index]->SetDirectory(0);
            mf::LogInfo("PhotonLibrarySplit") << "Added Tree " <<
                fInputTrees[index]->GetName() << " " << i << " " << j << " " <<
                k << " " << index << std::endl;
          }
        }
      }
    }
    catch (...) {
      throw cet::exception("PhotonLibrarySplit")
        << "Error in ttree load, reading photon library: " << LibraryFile << "\n";
    }
  }

  void PhotonLibrarySplit::LoadLibraryFromFile(size_t ix, size_t iy, size_t iz) {

    mf::LogInfo("PhotonLibrarySplit") <<
        "Loading Library from Tree" << std::endl;

    fLookupTable.clear();

    //TODO-- check that ix, iy, iz in bounds

    /*TFile* f = nullptr;
    TTree* tt = nullptr;
    TDirectory* pSrcDir = nullptr;*/

    /*
    TString library_name;
    library_name.Form("PhotonLibrarySplitData_%i_%i_%i", ix, iy, iz);
    try {
      tt = (TTree*)f->Get(library_name);

      if (tt) { pSrcDir = f; }
      else { // Library not in the top directory
        TKey* key = f->FindKeyAny(library_name);
        if (key) {
          tt = (TTree*)key->ReadObj();
          pSrcDir = key->GetMotherDir();
        }
        else {
          mf::LogError("PhotonLibrarySplit") << "PhotonLibrarySplitData not found in file" << LibraryFile;
        }
      }
    }
    catch (...) {
      throw cet::exception("PhotonLibrarySplit")
        << "Error in ttree load, reading photon library: " << LibraryFile << "\n";
    }*/

    int Voxel;
    int OpChannel;
    float Visibility;

    TTree * tt = fInputTrees[GetTreeIndex(ix, iy, iz)];
    std::cout << "Tree: " << tt << std::endl;
    tt->SetBranchAddress("Voxel", &Voxel);
    tt->SetBranchAddress("OpChannel", &OpChannel);
    tt->SetBranchAddress("Visibility", &Visibility);

    // with STL vectors, where `resize()` directly controls the allocation of
    // memory, reserving the space is redundant; not so with `util::LazyVector`,
    // where `resize()` never increases the memory; `data_init()` allocates
    // all the storage we need at once, effectively suppressing the laziness
    // of the vector (by design, that was only relevant in `CreateEmptyLibrary()`)
    fLookupTable.resize(LibrarySize());
    fLookupTable.data_init(LibrarySize());

    SetCurrentRegion(ix, iy, iz);

    std::cout << fLookupTable.size() << std::endl;

    size_t NEntries = tt->GetEntries();
    for (size_t i = 0; i < NEntries; ++i) {
      tt->GetEntry(i);
      //std::cout << Voxel << " " << OpChannel << " " << Visibility << std::endl;
      // Set the visibility at this optical channel
      auto coords = fVoxelDef.GetVoxelCoords(Voxel);

      if ((size_t(coords[0]) < ix*fNx || size_t(coords[0]) >= (ix+1)*fNx) ||
          (size_t(coords[1]) < iy*fNy || size_t(coords[1]) >= (iy+1)*fNy) ||
          (size_t(coords[2]) < iz*fNz || size_t(coords[2]) >= (iz+1)*fNz)) {
        std::cout << "Error " << coords[0] << " " << coords[1] << " " << coords[2] << std::endl;
        std::cout << "\t" << ix << " " << iy << " " << iz << std::endl;
        break;
      }

      uncheckedAccess(Voxel, OpChannel) = Visibility;
    } // for entries

  }

  void PhotonLibrarySplit::CloseFile() {
    try {
      fInputFile->Close();
    }
    catch (...) {
      mf::LogError("PhotonLibrarySplit") << "Error in closing file : " <<
          fInputFile->GetName();
    }
  }
  float PhotonLibrarySplit::GetCount(size_t Voxel, size_t OpChannel) const {
    //if (((Voxel - GetBaseVoxel(Voxel)) >= fNVoxels) || (OpChannel >= fNOpChannels))
    if ((TransformVoxel(Voxel) >= fNVoxels) || (OpChannel >= fNOpChannels))
      return 0;
    else
      return uncheckedAccess(Voxel, OpChannel);
  }

  float const* PhotonLibrarySplit::GetCounts(size_t Voxel) const {
    if (TransformVoxel(Voxel) >= fNVoxels)
      return nullptr;
    else
      return fLookupTable.data_address(uncheckedIndex(Voxel, 0));
  }
}

