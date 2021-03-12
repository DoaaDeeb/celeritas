//----------------------------------*-C++-*----------------------------------//
// Copyright 2020 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file SeltzerBergerModel.cc
//---------------------------------------------------------------------------//
#include "SeltzerBergerModel.hh"

#include "base/Assert.hh"
#include "base/CollectionBuilder.hh"
#include "base/Range.hh"
#include "physics/base/ParticleParams.hh"
#include "physics/base/PDGNumber.hh"
#include "physics/material/MaterialParams.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Construct from model ID and other necessary data.
 */
SeltzerBergerModel::SeltzerBergerModel(ModelId               id,
                                       const ParticleParams& particles,
                                       const MaterialParams& materials,
                                       ReadData              load_sb_table)
{
    CELER_EXPECT(id);
    CELER_EXPECT(load_sb_table);

    detail::SeltzerBergerData<Ownership::value, MemSpace::host> host_data;

    // Save IDs
    host_data.ids.model    = id;
    host_data.ids.electron = particles.find(pdg::electron());
    host_data.ids.positron = particles.find(pdg::positron());
    host_data.ids.gamma    = particles.find(pdg::gamma());
    CELER_VALIDATE(host_data.ids,
                   "Electron, positron and gamma particles must be enabled to "
                   "use the Seltzer-Berger Model.");

    // Save particle properties
    host_data.electron_mass
        = particles.get(host_data.ids.electron).mass().value();

    // Load differential cross sections
    make_builder(&host_data.differential_xs.grids)
        .reserve(materials.num_elements());
    for (auto el_id : range(ElementId{materials.num_elements()}))
    {
        AtomicNumber z_number = materials.get(el_id).atomic_number();
        this->append_table(load_sb_table(z_number), &host_data.differential_xs);
    }
    CELER_ASSERT(host_data.differential_xs.grids.size()
                 == materials.num_elements());

    // Move to mirrored data, copying to device
    data_ = CollectionMirror<detail::SeltzerBergerData>{std::move(host_data)};

    CELER_ENSURE(this->data_);
}

//---------------------------------------------------------------------------//
/*!
 * Particle types and energy ranges that this model applies to.
 */
auto SeltzerBergerModel::applicability() const -> SetApplicability
{
    Applicability electron_applic;
    electron_applic.particle = this->host_pointers().ids.electron;
    electron_applic.lower    = units::MevEnergy{1};
    electron_applic.upper    = units::MevEnergy{1e5};

    Applicability positron_applic;
    positron_applic.particle = this->host_pointers().ids.positron;
    positron_applic.lower    = units::MevEnergy{1};
    positron_applic.upper    = units::MevEnergy{1e5};

    return {electron_applic, positron_applic};
}

//---------------------------------------------------------------------------//
/*!
 * Apply the interaction kernel.
 */
void SeltzerBergerModel::interact(
    CELER_MAYBE_UNUSED const ModelInteractPointers& pointers) const
{
#if CELERITAS_USE_CUDA
    detail::seltzer_berger_interact(this->device_pointers(), pointers);
#else
    CELER_ASSERT_UNREACHABLE();
#endif
}

//---------------------------------------------------------------------------//
/*!
 * Get the model ID for this model.
 */
ModelId SeltzerBergerModel::model_id() const
{
    return this->host_pointers().ids.model;
}

//---------------------------------------------------------------------------//
/*!
 * Construct differential cross section tables for a single element.
 */
void SeltzerBergerModel::append_table(const ImportSBTable& imported,
                                      HostXsTables*        tables) const
{
    auto reals = make_builder(&tables->reals);

    CELER_ASSERT(!imported.value.empty()
                 && imported.value.size()
                        == imported.x.size() * imported.y.size());
    TwodGridData grid;

    // TODO: we could probably use a single x and y grid for all elements.
    // Only Z = 100 has different energy grids.
    // Incident charged particle log energy grid
    grid.x = reals.insert_back(imported.x.begin(), imported.x.end());

    // Photon reduced energy grid
    grid.y = reals.insert_back(imported.y.begin(), imported.y.end());

    // 2D scaled DCS grid
    grid.values
        = reals.insert_back(imported.value.begin(), imported.value.end());

    CELER_ASSERT(grid);
    make_builder(&tables->grids).push_back(grid);
}

//---------------------------------------------------------------------------//
} // namespace celeritas
