//----------------------------------*-C++-*----------------------------------//
// Copyright 2020 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file BremAngularSampler.hh
//---------------------------------------------------------------------------//
#include "BremAngularSampler.hh"

#include <cmath>
#include "base/ArrayUtils.hh"
#include "base/Algorithms.hh"
#include "random/distributions/GenerateCanonical.hh"

namespace celeritas
{
namespace detail
{
BremAngularSampler::BremAngularSampler(real_type inc_mass,
                                       real_type inc_energy,
                                       Real3     inc_direction)
    : inc_mass_{inc_mass}
    , inc_energy_{inc_energy}
    , inc_direction_{inc_direction}
{
    CELER_EXPECT(inc_mass_ > 0.0);
    CELER_EXPECT(inc_energy_ > 0.0);
}

template<class Engine>
Real3 BremAngularSampler::operator()(Engine& rng)
{
    real_type umax = 2.0 * (1.0 * inc_energy_ / inc_mass_);
    real_type u;

    do
    {
        u = -std::log(generate_canonical(rng) * generate_canonical(rng));
        u /= 0.25 > generate_canonical(rng) ? 0.625 : 1.875;
    } while (u > umax);

    real_type cost = 1.0 - 2.0 * ipow<2>(u) / ipow<2>(umax);
    real_type sint = std::sqrt(1 - ipow<2>(cost));
    real_type phi  = 2.0 * constants::pi * generate_canonical(rng);

    return rotate({sint * std::cos(phi), sint * std::sin(phi), cost},
                  inc_direction_);
}

//---------------------------------------------------------------------------//
} // namespace detail
} // namespace celeritas
