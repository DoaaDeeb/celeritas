//----------------------------------*-C++-*----------------------------------//
// Copyright 2020 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file BremAngularSampler.hh
//---------------------------------------------------------------------------//
#pragma once

#include "physics/base/Units.hh"

namespace celeritas
{
namespace detail
{
//---------------------------------------------------------------------------//
/*!
 * Angular distribution sampler for bremsstrahlung processes.
 *
 * \note This sampler performs similarly the analogous routines from the Geant4
 * G4ModifiedTsai class; suggested by L.Urban (Geant3 manual (1993) Phys211)
 * and derived from Tsai distribution (Rev Mod Phys 49,421(1977)).
 */
class BremAngularSampler
{
    //!@{
    //! Type aliases
    using Energy = units::MevEnergy;
    //!@}
  public:
    // Construct sampler. Note that `inc_energy` is the kinetic energy of the
    // incident particle.
    inline CELER_FUNCTION BremAngularSampler(real_type inc_mass,
                                             real_type inc_energy,
                                             Real3     inc_direction);

    // Sample direction with the given RNG
    template<class Engine>
    inline CELER_FUNCTION Real3 operator()(Engine& rng);

  private:
    const real_type inc_mass_;
    const real_type inc_energy_;
    const Real3     inc_direction_;

    BremAngularSampler& operator=(const BremAngularSampler& right) = delete;
    BremAngularSampler(const BremAngularSampler&)                  = delete;
};

//---------------------------------------------------------------------------//
} // namespace detail
} // namespace celeritas

#include "BremAngularSampler.i.hh"
