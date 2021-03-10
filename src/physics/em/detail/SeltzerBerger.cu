//---------------------------------*-CUDA-*----------------------------------//
// Copyright 2020 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file SeltzerBerger.cu
//---------------------------------------------------------------------------//
#include "SeltzerBerger.hh"

#include "base/Assert.hh"
#include "base/KernelParamCalculator.cuda.hh"
#include "random/cuda/RngEngine.hh"
#include "physics/base/ModelInterface.hh"
#include "physics/base/ParticleTrackView.hh"
#include "physics/base/PhysicsTrackView.hh"
#include "physics/base/SecondaryAllocatorView.hh"
#include "physics/material/MaterialTrackView.hh"
#include "SeltzerBergerInteractor.hh"

namespace celeritas
{
namespace detail
{
namespace
{
//---------------------------------------------------------------------------//
// KERNELS
//---------------------------------------------------------------------------//
/*!
 * Interact using the Seltzer-Berger model on applicable tracks.
 */
__global__ void seltzer_berger_interact_kernel(const SeltzerBergerPointers sb,
                                               const ModelInteractPointers ptrs)
{
    auto tid = celeritas::KernelParamCalculator::thread_id();
    if (tid.get() >= ptrs.states.size())
        return;

    SecondaryAllocatorView allocate_secondaries(ptrs.secondaries);
    ParticleTrackView particle(ptrs.params.particle, ptrs.states.particle, tid);

    // Setup for ElementView access
    MaterialTrackView material(ptrs.params.material, ptrs.states.material, tid);
    // Cache the associated MaterialView as function calls to MaterialTrackView
    // are expensive
    MaterialView material_view = material.material_view();

    PhysicsTrackView physics(ptrs.params.physics,
                             ptrs.states.physics,
                             particle.particle_id(),
                             material.material_id(),
                             tid);

    // This interaction only applies if the Seltzer-Berger model was selected
    if (physics.model_id() != sb.model_id)
        return;

    // Assume only a single element in the material, for now
    CELER_ASSERT(material_view.num_elements() == 1);
    SeltzerBergerInteractor interact(
        sb,
        particle,
        ptrs.states.direction[tid.get()],
        allocate_secondaries,
        material_view.element_view(celeritas::ElementComponentId{0}));

    RngEngine rng(ptrs.states.rng, tid);
    ptrs.result[tid.get()] = interact(rng);
    CELER_ENSURE(ptrs.result[tid.get()]);
}

} // namespace

//---------------------------------------------------------------------------//
// LAUNCHERS
//---------------------------------------------------------------------------//
/*!
 * Launch the Seltzer-Berger interaction.
 */
void seltzer_berger_interact(const SeltzerBergerPointers& sb,
                             const ModelInteractPointers& model)
{
    CELER_EXPECT(sb);
    CELER_EXPECT(model);

    static const KernelParamCalculator calc_kernel_params(
        seltzer_berger_interact_kernel, "seltzer_berger_interact");
    auto params = calc_kernel_params(model.states.size());
    seltzer_berger_interact_kernel<<<params.grid_size, params.block_size>>>(
        sb, model);
    CELER_CUDA_CHECK_ERROR();
}

//---------------------------------------------------------------------------//
} // namespace detail
} // namespace celeritas
