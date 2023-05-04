////# PhotonLibrarySplit.h header file
////#
////# Jake Calcutt, MIT, 2023
#ifndef PHOTONLIBRARYSPLIT_H
#define PHOTONLIBRARYSPLIT_H

#include "larsim/PhotonPropagation/IPhotonLibrary.h"

#include "larsim/Simulation/PhotonVoxels.h"

#include "TF1.h"
#include "TFile.h"
#include "TTree.h"
class TTree;

#include "lardataobj/Utilities/LazyVector.h"

#include <limits> // std::numeric_limits
#include <optional>

namespace art {
  class TFileDirectory;
}

namespace phot {

  class PhotonLibrarySplit : public IPhotonLibrary {
  public:
    PhotonLibrarySplit(sim::PhotonVoxelDef vox_def,
                       size_t n_split_x,
                       size_t n_split_y,
                       size_t n_split_z,
                       size_t n_op_chans);

    void OpenFile(std::string LibraryFile);
    void CloseFile();

    virtual float GetCount(size_t Voxel, size_t OpChannel) const override;
    //void SetCount(size_t Voxel, size_t OpChannel, float Count);

    /// Don't implement reflected light
    virtual bool hasReflected() const override { return false; }
    virtual const float* GetReflCounts(size_t Voxel) const override { return 0; }
    virtual float GetReflCount(size_t Voxel, size_t OpChannel) const override { return 0; }

    /// Don't implement reflected light timing
    virtual bool hasReflectedT0() const override { return false; }
    virtual const float* GetReflT0s(size_t Voxel) const override { return 0; }
    virtual float GetReflT0(size_t Voxel, size_t OpChannel) const override { return 0; }


    /// Returns a pointer to NOpChannels() visibility values, one per channel
    virtual float const* GetCounts(size_t Voxel) const override;
    void LoadLibraryFromFile(size_t ix, size_t iy, size_t iz);

    // --- BEGIN --- Metadata: voxel information -------------------------------
    /// @name Metadata: voxel information
    /// @{

    /// Returns whether voxel metadata is available.
    //bool hasVoxelDef() const { return fVoxelDef.has_value(); }

    /// Returns the current voxel metadata (undefined behaviour if none).
    /// @see `hasVoxelDef()`
    /*sim::PhotonVoxelDef const& GetVoxelDef() const
    {
      assert(fVoxelDef);
      return *fVoxelDef;
    }*/

    ///
    void SetCurrentRegion(size_t ix, size_t iy, size_t iz) {
      fCurrentIx = ix;
      fCurrentIy = iy;
      fCurrentIz = iz;
    }

    size_t GetTreeIndex(size_t i, size_t j, size_t k) const;
    std::array<size_t, 3> TreeIndexToRegion(size_t i) const;

    std::array<size_t, 3> GetSplitRegion(size_t i, size_t j, size_t k) const {
      return {i/fNx, j/fNy, k/fNz};
    }

    /// @}
    // --- END --- Metadata: voxel information ---------------------------------

  private:
    virtual int NOpChannels() const override { return fNOpChannels; }
    virtual int NVoxels() const override { return fNVoxels; }

    virtual bool isVoxelValid(size_t Voxel) const override { return isVoxelValidImpl(Voxel); }

    /*bool fHasReflected = false; ///< Whether the current library deals with reflected light counts.
    bool fHasReflectedT0 =
      false; ///< Whether the current library deals with reflected light timing.

    size_t fHasTiming =
      0; ///< Whether the current library deals with time propagation distribution.*/

    // fLookupTable[unchecked_index(Voxel, OpChannel)] = Count
    // for each voxel, all NChannels() channels are stored in sequence
    util::LazyVector<float> fLookupTable;
    /*util::LazyVector<float> fReflLookupTable;
    util::LazyVector<float> fReflTLookupTable;
    util::LazyVector<std::vector<float>> fTimingParLookupTable;
    util::LazyVector<TF1> fTimingParTF1LookupTable;
    std::string fTimingParFormula;
    size_t fTimingParNParameters;*/

    size_t fNOpChannels;
    size_t fNVoxels;

    /// Voxel definition loaded from library metadata.
    sim::PhotonVoxelDef fVoxelDef;
    size_t fNSplitX, fNSplitY, fNSplitZ;
    size_t fNx, fNy, fNz;
    size_t fCurrentIx = 0, fCurrentIy = 0, fCurrentIz = 0;

    TFile * fInputFile;
    std::vector<TTree*> fInputTrees;

    bool isVoxelValidImpl(size_t Voxel) const { return Voxel < fNVoxels; }

    size_t TransformVoxel(size_t Voxel) const {
      auto coords = fVoxelDef.GetVoxelCoords(Voxel);
      size_t i = coords[0] - fCurrentIx*fNx;
      size_t j = coords[1] - fCurrentIy*fNy;
      size_t k = coords[2] - fCurrentIz*fNz;

      //std::cout << coords[0] << " " << coords[1] << " " << coords[2] << std::endl;
      //std::cout << i << " " << j << " " << k << std::endl;
      return i + j*fNx + k*fNy*fNz;
    }

    /*size_t GetBaseVoxel(size_t Voxel) const {
      auto coords = fVoxelDef.GetVoxelCoords(Voxel);
      size_t base_voxel = ((fNx*(fNSplitX - 1)*fCurrentIx) +
                           (fNx*fNSplitX*coords[1]) +
                           (fNx*fNSplitX*fNy*fNSplitY*coords[2]));

                            



      return base_voxel;
    }*/

    /// Returns the index of visibility of specified voxel and cell
    size_t uncheckedIndex(size_t Voxel, size_t OpChannel) const
    {
      //std::cout << "\t" << GetBaseVoxel(Voxel) << std::endl;
      return TransformVoxel(Voxel)*fNOpChannels + OpChannel;
    }

    /// Unchecked access to a visibility datum
    float uncheckedAccess(size_t Voxel, size_t OpChannel) const
    {
      return fLookupTable[uncheckedIndex(Voxel, OpChannel)];
    }

    /// Unchecked access to a visibility datum
    float& uncheckedAccess(size_t Voxel, size_t OpChannel)
    {
      size_t index = uncheckedIndex(Voxel, OpChannel);
      //std::cout << "\t" << fLookupTable.size() << " " <<  index << std::endl;
      return fLookupTable[index];
    }

    /// Name of the optical channel number in the input tree
    static std::string const OpChannelBranchName;

    /// Returns the number of optical channels in the specified tree
    static size_t ExtractNOpChannels(TTree* tree);

  };

}

#endif
