//----------------------------------*-C++-*----------------------------------//
// Copyright 2021 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file SeltzerBergerInterface.hh
//---------------------------------------------------------------------------//
#pragma once

#include "base/Collection.hh"
#include "base/Types.hh"
#include "physics/grid/TwodGridInterface.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Persistent shared bremsstrahlung DCS data.
 *
 * \sa SeltzerBergerParams (owns the pointed-to data)
 */
template<Ownership W, MemSpace M>
struct SeltzerBergerParamsData
{
    template<class T>
    using Data = Collection<T, W, M>;

    // Backend storage
    Data<real_type>    reals;
    Data<TwodGridData> value_grids; //!< One grid per element

    //// MEMBER FUNCTIONS ////

    //! Whether the data is assigned
    explicit inline CELER_FUNCTION operator bool() const
    {
        return !value_grids.empty();
    }

    //! Assign from another set of data
    template<Ownership W2, MemSpace M2>
    SeltzerBergerParamsData&
    operator=(const SeltzerBergerParamsData<W2, M2>& other)
    {
        CELER_EXPECT(other);
        reals       = other.reals;
        value_grids = other.value_grids;
        return *this;
    }
};

//---------------------------------------------------------------------------//
} // namespace celeritas
