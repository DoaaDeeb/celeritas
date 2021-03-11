//----------------------------------*-C++-*----------------------------------//
// Copyright 2021 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file SeltzerBergerReader.cc
//---------------------------------------------------------------------------//
#include "io/SeltzerBergerReader.hh"
#include "base/Range.hh"
#include "base/Assert.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Construct using environmental variable $G4LEDATA.
 */
SeltzerBergerReader::SeltzerBergerReader()
{
    CELER_VALIDATE(std::getenv("G4LEDATA"),
                   "Environment variable G4LEDATA is not defined.");
    path_to_file_ = std::getenv("G4LEDATA");
}

//---------------------------------------------------------------------------//
/*!
 * Construct using a user defined path to the folder containing the data.
 * The path should point to the files that are usually stored in
 *
 * [Geant4-install]/share/Geant4-10.7.0/data/G4EMLOW7.12/brem_SB/.
 */
SeltzerBergerReader::SeltzerBergerReader(std::string folder_path)
    : path_to_file_(folder_path)
{
    CELER_ENSURE(path_to_file_.size());
}

//---------------------------------------------------------------------------//
/*!
 * Fetch data for a given atomic number.
 */
SeltzerBergerReader::result_type
SeltzerBergerReader::operator()(unsigned int atomic_number)
{
    // Standard data files encompass Z = [1, 99]
    CELER_EXPECT(atomic_number > 0 && atomic_number < 100);

    // Open file for given atomic number
    std::string file = path_to_file_ + "/brem_SB/br"
                       + std::to_string(atomic_number);
    std::ifstream input_stream(file.c_str());
    CELER_VALIDATE(input_stream, "Could not open file " << file << ".");

    result_type data_vector = this->retrieve(input_stream);
    input_stream.close();

    CELER_ENSURE(data_vector.x.size() && data_vector.y.size());
    return data_vector;
}

//---------------------------------------------------------------------------//
// PRIVATE
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
/*!
 * Parse selected input data and store it into a result_type structure.
 * This method does *NOT* closes the input stream.
 */
SeltzerBergerReader::result_type
SeltzerBergerReader::retrieve(std::ifstream& input_stream)
{
    result_type vector;

    // Fetch binning information
    unsigned int g4_physics_vector_type; // Not used
    unsigned int x_size;
    unsigned int y_size;

    input_stream >> g4_physics_vector_type >> x_size >> y_size;

    CELER_VALIDATE(input_stream, "Could not read file.");
    CELER_VALIDATE(x_size >= 2 && y_size >= 2, "Number of bins is too small.");

    // Resize result_type
    vector.x.resize(x_size);
    vector.y.resize(y_size);
    vector.value.resize(x_size * y_size);

    // Fetch content
    for (auto i : range(x_size))
    {
        CELER_ASSERT(input_stream);
        input_stream >> vector.x[i];
    }

    for (auto i : range(y_size))
    {
        CELER_ASSERT(input_stream);
        input_stream >> vector.y[i];
    }

    for (auto i : range(x_size))
    {
        for (auto j : range(y_size))
        {
            CELER_ASSERT(input_stream);
            input_stream >> vector.value[i * y_size + j];
        }
    }

    CELER_ENSURE(vector.x.size() == x_size && vector.y.size() == y_size
                 && vector.value.size() == x_size * y_size);
    return vector;
}

//---------------------------------------------------------------------------//
} // namespace celeritas
