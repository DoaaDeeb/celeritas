//----------------------------------*-C++-*----------------------------------//
// Copyright 2020 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file SeltzerBergerInteractor.i.hh
//---------------------------------------------------------------------------//
#include "BremAngularSampler.hh"

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
    StackAllocator<Secondary>&   allocate,
    const ElementView&           element)
    : shared_(shared)
    , inc_energy_(particle.energy().value())
    , inc_direction_(inc_direction)
    , allocate_(allocate)
    , element_(element)
{
    CELER_EXPECT(particle.particle_id() == shared_.ids.electron
                 || particle.particle_id() == shared_.ids.positron);
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
    // Allocate space for the brem-gamma
    Secondary* brems_gamma = this->allocate_(1);
    if (brems_gamma == nullptr)
    {
        // Failed to allocate space for a secondary
        return Interaction::from_failure();
    }

    // Data to be loaded from SBParamsData for SBSampler, based on material
    real_type cut_energy, max_energy;
    // Determine thresholds based on SB data.
    real_type kinetic_energy_min = celeritas::min(cut_energy, inc_energy_);
    real_type kinetic_energy_max = celeritas::min(max_energy, inc_energy_);
    CELER_ASSERT(cut_energy > kinetic_energy_max);

    /// TODO: How to get material density factor?
    // real_type density_factor       = mat.electron_density() * migdal;
    Energy primary_total_energy = Energy{
        inc_energy_ + shared_.electron_mass * ipow<2>(constants::c_light)};

    // Brem-gamma energy sampled either by rejection or sampling tables.
    // In G4, a flag is used to select sampling the method.
    // SBEnergyDistribution(shared_,
    //                      inc_energy_,
    //                      element_,
    //                      density_correction, // EnergySq
    //                      min_gamma_energy);  // Energy
    real_type gamma_energy = this->sample_energy_transfer(
        kinetic_energy_min, kinetic_energy_max, rng);

    // This should never happen under normal conditions but protect anyway
    CELER_ASSERT(gamma_energy > 0.0);

    BremAngularSampler sample_gamma_direction{
        shared_.electron_mass, inc_energy_, inc_direction_};
    Real3 gamma_direction = sample_gamma_direction(rng);

    // Construct interaction
    Interaction result;
    result.action      = Action::spawned;
    result.secondaries = {brems_gamma, 1};

    // Outgoing secondary is a photon
    brems_gamma[0].particle_id = shared_.ids.gamma;
    /// TODO: sample energy
    brems_gamma[0].energy    = units::MevEnergy{0.0};
    brems_gamma[0].direction = gamma_direction;

    // If brems-gamma is highly energetic, G4 stops tracking the
    // incoming/primary particle and creates new secondary e-/e+

    // Compute post-interaction kinematics of e-/e+

    return result;
}

template<class Engine>
CELER_FUNCTION real_type SeltzerBergerInteractor::sample_energy_transfer(
    real_type kinetic_energy_min, real_type kinetic_energy_max, Engine& rng)
{
    /// TODO:
    // density_factor == MigdalConstant * material_electron_density
    // density_correction == density_factor * inc_total_energy^2
    return 1.0;
}

//---------------------------------------------------------------------------//
} // namespace detail
} // namespace celeritas
