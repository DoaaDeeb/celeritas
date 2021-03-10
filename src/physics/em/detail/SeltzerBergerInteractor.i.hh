//----------------------------------*-C++-*----------------------------------//
// Copyright 2020 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file SeltzerBergerInteractor.i.hh
//---------------------------------------------------------------------------//

#include "base/ArrayUtils.hh"
#include "base/Constants.hh"
#include "random/distributions/BernoulliDistribution.hh"
#include "random/distributions/GenerateCanonical.hh"
#include "random/distributions/UniformRealDistribution.hh"

namespace celeritas
{
namespace detail
{
//---------------------------------------------------------------------------//
/*!
 * Construct with shared and state data.
 *
 * The incident particle must be above the energy threshold: this should be
 * handled in code *before* the interactor is constructed.
 */
SeltzerBergerInteractor::SeltzerBergerInteractor(
    const SeltzerBergerPointers& shared,
    const ParticleTrackView&     particle,
    const Real3&                 inc_direction,
    SecondaryAllocatorView&      allocate,
    const ElementView&           element)
    : shared_(shared)
    , inc_energy_(particle.energy().value())
    , inc_direction_(inc_direction)
    , allocate_(allocate)
    , element_(element)
{
    CELER_EXPECT(particle.particle_id() == shared_.electron_id);

    static_assert(sizeof(real_type) == sizeof(double),
                  "Embedded constants are hardcoded to double precision.");
}

//---------------------------------------------------------------------------//
/*!
 * e-/e+ bremsstrahlung using the Seltzer-Berger model.
 *
 * See section 6.5 of the Geant physics reference 10.6.
 */
template<class Engine>
CELER_FUNCTION Interaction SeltzerBergerInteractor::operator()(Engine& rng)
{
    // Allocate space for the pair-produced electrons
    Secondary* secondaries = this->allocate_(2);
    if (secondaries == nullptr)
    {
        // Failed to allocate space for a secondary
        return Interaction::from_failure();
    }

    // Construct interaction for change to primary (incident) particle (gamma)
    Interaction result;
    result.action = Action::spawned;

    return result;
}

//---------------------------------------------------------------------------//
} // namespace detail
} // namespace celeritas
