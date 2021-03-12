//----------------------------------*-C++-*----------------------------------//
// Copyright 2021 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file SeltzerBergerParams.cc
//---------------------------------------------------------------------------//
#include "SeltzerBergerParams.hh"

#include "base/CollectionBuilder.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Construct from a vector of cross section data.
 */
SeltzerBergerParams::SeltzerBergerParams(const Input& inp)
{
    CELER_EXPECT(!inp.empty());

    // Build data on host
    HostValue host_data;
    auto      reals       = make_builder(&host_data.reals);
    auto      value_grids = make_builder(&host_data.value_grids);
    value_grids.reserve(inp.size());

    for (const auto& el_grid : inp)
    {
        CELER_ASSERT(!el_grid.value.empty()
                     && el_grid.value.size()
                            == el_grid.x.size() * el_grid.y.size());
        TwodGridData grid;

        // TODO: we could probably use a single x and y grid for all elements.
        // Only Z = 100 has different energy grids.
        // Incident charged particle log energy grid
        grid.x = reals.insert_back(el_grid.x.begin(), el_grid.x.end());

        // Photon reduced energy grid
        grid.y = reals.insert_back(el_grid.y.begin(), el_grid.y.end());

        // 2D scaled DCS grid
        grid.values
            = reals.insert_back(el_grid.value.begin(), el_grid.value.end());

        CELER_ASSERT(grid);
        value_grids.push_back(grid);
    }

    // Move to mirrored data, copying to device
    data_ = CollectionMirror<SeltzerBergerParamsData>{std::move(host_data)};
    CELER_ENSURE(this->data_);
}

//---------------------------------------------------------------------------//
} // namespace celeritas
