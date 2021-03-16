//----------------------------------*-C++-*----------------------------------//
// Copyright 2020 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file SeltzerBergerInteractor.hh
//---------------------------------------------------------------------------//
#pragma once

#include "base/Macros.hh"
#include "base/Types.hh"
#include "physics/base/Interaction.hh"
#include "physics/base/ParticleTrackView.hh"
#include "physics/base/Secondary.hh"
#include "base/StackAllocator.hh"
#include "physics/base/Units.hh"
#include "physics/material/ElementView.hh"
#include "SeltzerBerger.hh"

namespace celeritas
{
namespace detail
{
//---------------------------------------------------------------------------//
/*!
 * Seltzer-Berger model for electron/positron bremsstrahlung.
 *
 * Provides energy loss of electrons and positrons due to radiative photons in
 * the field of a nucleus (details in the Geant4 Physics Reference 10.6).
 * The Seltzer-Berger model implements cross sections based on interpolation
 * of pre-calculated tables.
 *
 * \note This interactor performs the sampling routines from the Geant4
 * G4SeltzerBergerModel, documented in section 10.2.1 of the Geant4 Physics
 * Reference (release 10.6), covering electron and positron energies from
 * 1 keV to 10 GeV.
 */
class SeltzerBergerInteractor
{
  public:
    //!@{
    //! Type aliases
    using SeltzerBergerPointers
        = SeltzerBergerData<Ownership::const_reference, MemSpace::native>;
    using Energy   = units::MevEnergy;
    using EnergySq = Quantity<UnitProduct<units::Mev, units::Mev>>;
    //!@}

  public:
    //! Construct sampler from shared and state data
    inline CELER_FUNCTION
    SeltzerBergerInteractor(const SeltzerBergerPointers& shared,
                            const ParticleTrackView&     particle,
                            const Real3&                 inc_direction,
                            StackAllocator<Secondary>&   allocate,
                            const ElementView&           element);

    // Sample an interaction with the given RNG
    template<class Engine>
    inline CELER_FUNCTION Interaction operator()(Engine& rng);

  private:
    // Sample brems-gamma energy. The analogous G4 sampling routine uses a
    // Ranlux engine for random number generation; here we use the canonical
    // generation.
    template<class Engine>
    inline CELER_FUNCTION real_type
    sample_energy_transfer(real_type kinetic_energy_min,
                           real_type kinetic_energy_max,
                           Engine&   rng);

    // Device data for the model
    const SeltzerBergerPointers& shared_;

    // Incident electron (or positron) energy
    const real_type inc_energy_;

    // Incident direction
    const Real3& inc_direction_;

    // Allocate space for a secondary particle
    StackAllocator<Secondary>& allocate_;

    // Element properties for calculating screening functions and variables
    const ElementView& element_;

    const real_type migdal = 4.0 * constants::pi * constants::re_electron
                             * ipow<2>(constants::lambda_compton);

    // SBParamsData sb_tables;
    // class SBSampler (const SBParamsData& tables, const ElementView & el,
    //                  real_type inc_energy, particleID?)
};

//---------------------------------------------------------------------------//
} // namespace detail
} // namespace celeritas

#include "SeltzerBergerInteractor.i.hh"
