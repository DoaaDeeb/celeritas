//----------------------------------*-C++-*----------------------------------//
// Copyright 2021 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file CutoffView.hh
//---------------------------------------------------------------------------//
#pragma once

#include "CutoffInterface.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Access invariant material- and particle-dependent cutoff values.
 *
 * \c CutoffParamsData is defined in \c CutoffInterface and constructed by
 * \c CutoffParams .
 *
 * \code
 * CutoffParams cutoffs(input);
 * CutoffView cutoff_view(cutoffs.host_pointers(), particle_id, material_id);
 * cutoff_view.energy();
 * cutoff_view.range();
 * \endcode
 */
class CutoffView
{
  public:
    //!@{
    //! Type aliases
    using CutoffPointers
        = CutoffParamsData<Ownership::const_reference, MemSpace::native>;
    using Energy = units::MevEnergy;
    //!@}

  public:
    // Construct for the given particle and material ids
    inline CELER_FUNCTION CutoffView(const CutoffPointers& params,
                                     ParticleId            particle,
                                     MaterialId            material);

    //! Return energy cutoff value
    CELER_FORCEINLINE_FUNCTION Energy energy() const { return cutoff_.energy; }
    //! Return range cutoff value
    CELER_FORCEINLINE_FUNCTION real_type range() const
    {
        return cutoff_.range;
    }

  private:
    ParticleCutoff cutoff_;
};

//---------------------------------------------------------------------------//
} // namespace celeritas

#include "CutoffView.i.hh"