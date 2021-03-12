//----------------------------------*-C++-*----------------------------------//
// Copyright 2021 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file SeltzerBergerParams.hh
//---------------------------------------------------------------------------//
#pragma once

#include <vector>
#include "base/CollectionMirror.hh"
#include "io/SeltzerBergerReader.hh"
#include "physics/base/Units.hh"
#include "SeltzerBergerInterface.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Data management for the Seltzer-Berger bremsstrahlung differential cross
 * sections.
 *
 * The total scaled bremsstrahlung differential cross section for an element Z
 * is defined as
 * \f[
 *   \chi(Z,E,\kappa) = \frac{\beta^2}{Z^2} k \frac{d \sigma}{dk},
 * \f]
 * where \f$ \kappa = k / E \f$ is the ratio of the emitted photon energy to
 * the incident charged particle energy, \f$ \beta \f$ is the ratio of the
 * charged particle velocity to the speed of light, and \f$ \frac{d \sigma}{dk}
 * \f$ is the bremsstrahlung differential cross section.
 *
 * Seltzer and Berger have tabulated the scaled DCS (in mb) for elements Z = 1
 * - 100 and for incident charged particle energies from 1 keV to 10 GeV
 * (reported in MeV) in Seltzer S.M. and M.J. Berger (1986), "Bremsstrahlung
 * energy spectra from electrons with kinetic energy 1 keV–10 GeV incident on
 * screened nuclei and orbital electrons of neutral atoms with Z = 1–100", At.
 * Data Nucl. Data Tables 35, 345–418.
 */
class SeltzerBergerParams
{
  public:
    //!@{
    //! Type aliases
    using EnergyUnits = units::Mev;
    using XsUnits     = units::Millibarn;
    using HostRef
        = SeltzerBergerParamsData<Ownership::const_reference, MemSpace::host>;
    using DeviceRef
        = SeltzerBergerParamsData<Ownership::const_reference, MemSpace::device>;
    using Input = std::vector<SeltzerBergerReader::result_type>;
    //!@}

  public:
    // Construct from a vector of cross section data
    explicit SeltzerBergerParams(const Input& inp);

    //! Access SB data on the host
    const HostRef& host_pointers() const { return data_.host(); }

    //! Access SB data on the device
    const DeviceRef& device_pointers() const { return data_.device(); }

  private:
    using HostValue = SeltzerBergerParamsData<Ownership::value, MemSpace::host>;

    // Host/device storage and reference
    CollectionMirror<SeltzerBergerParamsData> data_;
};

//---------------------------------------------------------------------------//
} // namespace celeritas
