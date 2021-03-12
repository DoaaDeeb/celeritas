//----------------------------------*-C++-*----------------------------------//
// Copyright 2021 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file SeltzerBergerReader.hh
//---------------------------------------------------------------------------//
#pragma once

#include <string>
#include <vector>
#include "base/Types.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Read Seltzer-Berger data from Geant4's $G4LEDATA files.
 * Use \c operator() to retrieve data for different atomic numbers.
 *
 * \code
    SeltzerBergerReader sb_reader();
    auto sb_data_vector = sb_reader(1); // Hydrogen
   \endcode
 */
class SeltzerBergerReader
{
  public:
    struct result_type
    {
        std::vector<real_type> x; //!< Log energy of incident particle [MeV]
        std::vector<real_type> y; //!< Ratio of gamma energy to incident energy
        std::vector<real_type> value; //!< Scaled DCS [mb]
    };

  public:
    //! Construct using $G4LEDATA
    SeltzerBergerReader();
    //! Construct from a user defined path
    explicit SeltzerBergerReader(std::string folder_path);

    //! Read data from ascii and return result_type data
    result_type operator()(unsigned int atomic_number) const;

  private:
    std::string path_to_file_;
};

//---------------------------------------------------------------------------//
} // namespace celeritas
