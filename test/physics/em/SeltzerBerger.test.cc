//----------------------------------*-C++-*----------------------------------//
// Copyright 2020 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file SeltzerBerger.test.cc
//---------------------------------------------------------------------------//
#include "physics/em/detail/SeltzerBergerInteractor.hh"

#include "physics/material/ElementView.hh"
#include "physics/material/Types.hh"
#include "physics/base/Units.hh"
#include "physics/grid/TwodGridCalculator.hh"
#include "physics/em/SeltzerBergerModel.hh"
#include "physics/em/GammaConversionProcess.hh"
#include "physics/material/MaterialTrackView.hh"

#include "celeritas_test.hh"
#include "gtest/Main.hh"
#include "base/ArrayUtils.hh"
#include "base/Range.hh"
#include "physics/base/Units.hh"
#include "physics/em/SeltzerBergerParams.hh"
#include "../InteractorHostTestBase.hh"
#include "../InteractionIO.hh"

using celeritas::ElementId;
using celeritas::GammaConversionProcess;
using celeritas::SeltzerBergerParams;
using celeritas::SeltzerBergerReader;
using celeritas::detail::SeltzerBergerInteractor;
using celeritas::units::AmuMass;
namespace constants = celeritas::constants;
namespace pdg       = celeritas::pdg;

//---------------------------------------------------------------------------//
// TEST HARNESS
//---------------------------------------------------------------------------//

class SeltzerBergerInteractorTest
    : public celeritas_test::InteractorHostTestBase
{
    using Base = celeritas_test::InteractorHostTestBase;

  protected:
    void set_sb_params(SeltzerBergerParams::Input inp)
    {
        CELER_EXPECT(!inp.empty());
        sb_params_ = std::make_shared<SeltzerBergerParams>(std::move(inp));
    }

    void SetUp() override
    {
        using celeritas::MatterState;
        using celeritas::ParticleDef;
        using namespace celeritas::units;
        using namespace celeritas::constants;
        constexpr auto zero   = celeritas::zero_quantity();
        constexpr auto stable = ParticleDef::stable_decay_constant();

        // Set up shared particle data
        Base::set_particle_params(
            {{"electron",
              pdg::electron(),
              MevMass{0.5109989461},
              ElementaryCharge{-1},
              stable},
             {"positron",
              pdg::positron(),
              MevMass{0.5109989461},
              ElementaryCharge{1},
              stable},
             {"gamma", pdg::gamma(), zero, zero, stable}});

        // Set up shared material data
        MaterialParams::Input mat_inp;
        mat_inp.elements  = {{29, AmuMass{63.546}, "Cu"}};
        mat_inp.materials = {
            {1.0 * na_avogadro,
             293.0,
             MatterState::solid,
             {{ElementId{0}, 1.0}},
             "Cu"},
        };
        this->set_material_params(mat_inp);
        this->set_material("Cu");

        // Set up Seltzer-Berger cross section data
        std::string         data_path = this->test_data_path("physics/em", "");
        SeltzerBergerReader read_element_data(data_path);
        std::vector<SeltzerBergerReader::result_type> sb_inp;
        sb_inp.push_back(read_element_data(29));
        this->set_sb_params(sb_inp);
    }

    void sanity_check(const Interaction& interaction) const
    {
        ASSERT_TRUE(interaction);
    }

  protected:
    celeritas::detail::SeltzerBergerPointers pointers_;
    std::shared_ptr<SeltzerBergerParams>     sb_params_;
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(SeltzerBergerInteractorTest, basic) {}

TEST_F(SeltzerBergerInteractorTest, stress_test) {}

TEST_F(SeltzerBergerInteractorTest, model) {}
