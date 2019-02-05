/*
 * MIT License
 *
 * Copyright (c) 2018-2019 Benjamin Köhler
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include "ImporterScientific.h"

#include <algorithm>
#include <cstdint>
#include <filesystem>
#include <iostream>

#include <bk/StringUtils>

//====================================================================================================
//===== CONSTRUCTORS & DESTRUCTOR
//====================================================================================================
ImporterScientific::ImporterScientific() = default;
//ImporterScientific::ImporterScientific(const ImporterScientific&) = default;
ImporterScientific::ImporterScientific(ImporterScientific&&) = default;
ImporterScientific::~ImporterScientific() = default;

//====================================================================================================
//===== GETTER
//====================================================================================================
std::string ImporterScientific::result() const
{ return _res.str(); }

//====================================================================================================
//===== SETTER
//====================================================================================================
//[[maybe_unused]] ImporterScientific& ImporterScientific::operator=(const ImporterScientific&) = default;
[[maybe_unused]] ImporterScientific& ImporterScientific::operator=(ImporterScientific&&) = default;

void ImporterScientific::set_dir(std::string_view dir)
{ _dir = dir; }

//====================================================================================================
//===== FUNCTIONS
//====================================================================================================
void ImporterScientific::_read_nd_scalar_image_in_sparse_matrix_style(std::ifstream& file)
{
    /*
     *       [1] x [uint32] : numDims
     * [numDims] x [uint32] : size per dimension
     * [numDims] x [double] : scale per dimension
     *      [16] x [double] : world matrix (4x4 from dicom)
     *      [16] x [double] : inverse world matrix (4x4 from dicom)
     *      [25] x [double] : world matrix (5x5 including time in 4th row/col)
     *      [25] x [double] : inverse world matrix (5x5 including time in 4th row/col)
     *       [1] x [uint32] : numNonZeroValues
     *       for numNonZeroValues:
     *               [numDims] x [uint32] : xyzt grid id
     *                     [1] x [double] : value
     */

    /*
     * num dimensions
     */
    std::uint32_t numDims = 0;
    file.read(reinterpret_cast<char*>(&numDims), sizeof(std::uint32_t));
    _res << "\t\t- num. dimensions: " << numDims << std::endl;

    /*
     * grid size
     */
    std::vector<std::uint32_t> gridsize(numDims);
    file.read(reinterpret_cast<char*>(gridsize.data()), numDims * sizeof(std::uint32_t));

    _res << "\t\t- grid size: ";
    for (unsigned int i = 0; i < numDims; ++i)
    {
        _res << gridsize[i];
        if (i < numDims - 1)
        { _res << " x "; }
    }
    _res << std::endl;

    /*
     * voxel scale
     */
    std::vector<double> voxelscale(numDims);
    file.read(reinterpret_cast<char*>(voxelscale.data()), numDims * sizeof(double));

    _res << "\t\t- voxel scale: ";
    for (unsigned int i = 0; i < numDims; ++i)
    {
        _res << voxelscale[i];
        if (i < numDims - 1)
        { _res << " x "; }
    }
    _res << std::endl;

    /*
     * - world matrix (4x4 from dicom)
     * - inverse world matrix (4x4 from dicom)
     *
     * - world matrix (5x5 including time in 4th row/col)
     * - inverse world matrix (5x5 including time in 4th row/col)
     */
    std::vector<double> dbuffer(16);
    file.read(reinterpret_cast<char*>(dbuffer.data()), dbuffer.size() * sizeof(double));

    unsigned int cnt = 0;
    _res << "\t\t- world matrix:" << std::endl;
    for (unsigned int rowid = 0; rowid < 4; ++rowid)
    {
        _res << "\t\t\t";
        for (unsigned int colid = 0; colid < 4; ++colid)
        { _res << dbuffer[cnt++] << " "; }

        _res << std::endl;
    }

    //dbuffer.resize(16);
    file.read(reinterpret_cast<char*>(dbuffer.data()), dbuffer.size() * sizeof(double));

    cnt = 0;
    _res << "\t\t- inverse world matrix:" << std::endl;
    for (unsigned int rowid = 0; rowid < 4; ++rowid)
    {
        _res << "\t\t\t";
        for (unsigned int colid = 0; colid < 4; ++colid)
        { _res << dbuffer[cnt++] << " "; }

        _res << std::endl;
    }

    dbuffer.resize(25);
    file.read(reinterpret_cast<char*>(dbuffer.data()), dbuffer.size() * sizeof(double));

    cnt = 0;
    _res << "\t\t- world matrix with time:" << std::endl;
    for (unsigned int rowid = 0; rowid < 5; ++rowid)
    {
        _res << "\t\t\t";
        for (unsigned int colid = 0; colid < 5; ++colid)
        { _res << dbuffer[cnt++] << " "; }

        _res << std::endl;
    }

    //dbuffer.resize(25);
    file.read(reinterpret_cast<char*>(dbuffer.data()), dbuffer.size() * sizeof(double));

    cnt = 0;
    _res << "\t\t- inverse world matrix with time:" << std::endl;
    for (unsigned int rowid = 0; rowid < 5; ++rowid)
    {
        _res << "\t\t\t";
        for (unsigned int colid = 0; colid < 5; ++colid)
        { _res << dbuffer[cnt++] << " "; }

        _res << std::endl;
    }

    //------------------------------------------------------------------------------------------------------
    // non-zero values
    //------------------------------------------------------------------------------------------------------
    std::uint32_t numNonZeroValues = 0;
    file.read(reinterpret_cast<char*>(&numNonZeroValues), sizeof(std::uint32_t));
    _res << "\t\t- num. non-zero values: " << numNonZeroValues << std::endl;

    std::vector<std::uint32_t> gridpos(numDims);

    for (unsigned int i = 0; i < numNonZeroValues; ++i)
    {
        //------------------------------------------------------------------------------------------------------
        // grid pos
        //------------------------------------------------------------------------------------------------------
        file.read(reinterpret_cast<char*>(gridpos.data()), numDims * sizeof(std::uint32_t));

        //------------------------------------------------------------------------------------------------------
        // value
        //------------------------------------------------------------------------------------------------------
        double val = 0;
        file.read(reinterpret_cast<char*>(&val), sizeof(double));

        if (i < NUM_DEMO)
        {
            _res << "\t\t\t- " << i << ": [";
            for (unsigned int k = 0; k < numDims; ++k)
            {
                _res << gridpos[k];
                if (k < numDims - 1)
                { _res << ", "; }
            }
            _res << "] = " << val << std::endl;
        }
    }
    _res << "\t\t\t- ..." << std::endl;
}

bool ImporterScientific::read_mesh(std::string_view filepath)
{
    /*
     *                           [1] x [uint32] : numPoints
     *               [numPoints * 3] x [double] : list of points
     *               [numPoints * 3] x [double] : list of normals per point
     *                           [1] x [uint32] : numTriangles
     *            [numTriangles * 3] x [uint32] : list of triangles (indices of points per triangle)
     *            [numTriangles * 3] x [double] : list of normals per triangle
     *                           [1] x [uint32] : numTimes
     *        [numPoints * numTimes] x [double] : wall shear stress per point over time
     *        [numPoints * numTimes] x [double] : wall shear stress per point over time : axial
     *        [numPoints * numTimes] x [double] : wall shear stress per point over time : circumferential
     *    [numPoints * numTimes * 3] x [double] : wall shear stress VECTOR per point over time
     *    [numPoints * numTimes * 3] x [double] : wall shear stress VECTOR per point over time : axial
     *    [numPoints * numTimes * 3] x [double] : wall shear stress VECTOR per point over time : circumferential
     *                   [numPoints] x [double] : mean wall shear stress
     *                   [numPoints] x [double] : mean wall shear stress : axial
     *                   [numPoints] x [double] : mean wall shear stress : circumferential
     *                   [numPoints] x [double] : oscillatory shear index
     *                   [numPoints] x [double] : oscillatory shear index : axial
     *                   [numPoints] x [double] : oscillatory shear index : circumferential
     *               [numPoints * 3] x [double] : mean wss vector
     *               [numPoints * 3] x [double] : mean wss vector : axial
     *               [numPoints * 3] x [double] : mean wss vector : circumferential
     */

    if (!std::filesystem::exists(filepath.data()))
    {
        _res << "\t- vessel has no mesh (path \"" << filepath.data() << "\")" << std::endl;
        return false;
    }

    _res << "\t- reading mesh (path \"" << filepath.data() << "\")" << std::endl;

    std::ifstream file(filepath.data(), std::ios_base::in | std::ios_base::binary);

    if (!file.good())
    {
        _res << "\t\tFAILED! Could not open file!" << std::endl;
        return false;
    }

    //------------------------------------------------------------------------------------------------------
    // num points
    //------------------------------------------------------------------------------------------------------
    std::uint32_t numPoints = 0;
    file.read(reinterpret_cast<char*>(&numPoints), sizeof(std::uint32_t));
    _res << "\t\t- num. points: " << numPoints << std::endl;

    //------------------------------------------------------------------------------------------------------
    // list of points
    //------------------------------------------------------------------------------------------------------
    std::vector<double> dbuffer(3 * numPoints);
    file.read(reinterpret_cast<char*>(dbuffer.data()), dbuffer.size() * sizeof(double));

    for (unsigned int pointid = 0; pointid < NUM_DEMO; ++pointid)
    {
        const unsigned int off = pointid * 3;
        _res << "\t\t\t- point" << pointid << ": [" << dbuffer[off] << ", " << dbuffer[off + 1] << ", " << dbuffer[off + 2] << "]" << std::endl;
    }
    _res << "\t\t\t- ..." << std::endl;

    //------------------------------------------------------------------------------------------------------
    // list of point normals
    //------------------------------------------------------------------------------------------------------
    //dbuffer.resize(3 * numPoints);
    file.read(reinterpret_cast<char*>(dbuffer.data()), dbuffer.size() * sizeof(double));

    for (unsigned int pointid = 0; pointid < NUM_DEMO; ++pointid)
    {
        const unsigned int off = pointid * 3;
        _res << "\t\t\t- normal" << pointid << ": [" << dbuffer[off] << ", " << dbuffer[off + 1] << ", " << dbuffer[off + 2] << "]" << std::endl;
    }
    _res << "\t\t\t- ..." << std::endl;

    //------------------------------------------------------------------------------------------------------
    // num triangles
    //------------------------------------------------------------------------------------------------------
    std::uint32_t numTriangles = 0;
    file.read(reinterpret_cast<char*>(&numTriangles), sizeof(std::uint32_t));

    //------------------------------------------------------------------------------------------------------
    // list of triangles
    //------------------------------------------------------------------------------------------------------
    std::vector<std::uint32_t> ui32buffer(3 * numTriangles);
    file.read(reinterpret_cast<char*>(ui32buffer.data()), ui32buffer.size() * sizeof(std::uint32_t));
    _res << "\t\t- num. triangles: " << numTriangles << std::endl;

    for (unsigned int cellid = 0; cellid < NUM_DEMO; ++cellid)
    {
        const unsigned int off = cellid * 3;
        _res << "\t\t\t- triangle" << cellid << ": [" << ui32buffer[off] << ", " << ui32buffer[off + 1] << ", " << ui32buffer[off + 2] << "]" << std::endl;
    }
    _res << "\t\t\t- ..." << std::endl;

    //------------------------------------------------------------------------------------------------------
    // list of triangle normals
    //------------------------------------------------------------------------------------------------------
    dbuffer.resize(3 * numTriangles);
    file.read(reinterpret_cast<char*>(dbuffer.data()), dbuffer.size() * sizeof(double));

    for (unsigned int cellid = 0; cellid < NUM_DEMO; ++cellid)
    {
        const unsigned int off = cellid * 3;
        _res << "\t\t\t- normal" << cellid << ": [" << dbuffer[off] << ", " << dbuffer[off + 1] << ", " << dbuffer[off + 2] << "]" << std::endl;
    }
    _res << "\t\t\t- ..." << std::endl;

    //------------------------------------------------------------------------------------------------------
    // num temporal positions
    //------------------------------------------------------------------------------------------------------
    std::uint32_t numTimes = 0;
    file.read(reinterpret_cast<char*>(&numTimes), sizeof(std::uint32_t));
    _res << "\t\t- num. temporal positions: " << numTimes << std::endl;

    //------------------------------------------------------------------------------------------------------
    //  wall shear stress per point over time
    //------------------------------------------------------------------------------------------------------
    dbuffer.resize(numPoints * numTimes);
    file.read(reinterpret_cast<char*>(dbuffer.data()), dbuffer.size() * sizeof(double));
    _res << "\t\t- WSS per point per time:" << std::endl;

    for (unsigned int pointid = 0; pointid < NUM_DEMO; ++pointid)
    {
        _res << "\t\t\t- point" << pointid << ": ";

        for (unsigned int timeid = 0; timeid < NUM_DEMO; ++timeid)
        {
            const unsigned int off = pointid * numTimes + timeid;
            _res << dbuffer[off] << ", ";
        }

        _res << "..." << std::endl;
    }

    //------------------------------------------------------------------------------------------------------
    // wall shear stress per point over time : AXIAL
    //------------------------------------------------------------------------------------------------------
    //dbuffer.resize(numPoints * numTimes);
    file.read(reinterpret_cast<char*>(dbuffer.data()), dbuffer.size() * sizeof(double));
    _res << "\t\t- Axial WSS per point per time:" << std::endl;

    for (unsigned int pointid = 0; pointid < NUM_DEMO; ++pointid)
    {
        _res << "\t\t\t- point" << pointid << ": ";

        for (unsigned int timeid = 0; timeid < NUM_DEMO; ++timeid)
        {
            const unsigned int off = pointid * numTimes + timeid;
            _res << dbuffer[off] << ", ";
        }

        _res << "..." << std::endl;
    }
    _res << "\t\t\t- ..." << std::endl;

    //------------------------------------------------------------------------------------------------------
    // wall shear stress per point over time : CIRCUMFERENTIAL
    //------------------------------------------------------------------------------------------------------
    //dbuffer.resize(numPoints * numTimes);
    file.read(reinterpret_cast<char*>(dbuffer.data()), dbuffer.size() * sizeof(double));
    _res << "\t\t- Circumferential WSS per point per time:" << std::endl;

    for (unsigned int pointid = 0; pointid < NUM_DEMO; ++pointid)
    {
        _res << "\t\t\t- point" << pointid << ": ";

        for (unsigned int timeid = 0; timeid < NUM_DEMO; ++timeid)
        {
            const unsigned int off = pointid * numTimes + timeid;
            _res << dbuffer[off] << ", ";
        }

        _res << "..." << std::endl;
    }
    _res << "\t\t\t- ..." << std::endl;

    //------------------------------------------------------------------------------------------------------
    // wall shear stress vector per point over time
    //------------------------------------------------------------------------------------------------------
    dbuffer.resize(numPoints * numTimes * 3);
    file.read(reinterpret_cast<char*>(dbuffer.data()), dbuffer.size() * sizeof(double));
    _res << "\t\t- WSS vector per point per time:" << std::endl;

    for (unsigned int pointid = 0; pointid < NUM_DEMO; ++pointid)
    {
        _res << "\t\t\t- point" << pointid << ": ";

        for (unsigned int timeid = 0; timeid < NUM_DEMO; ++timeid)
        {
            const unsigned int off = pointid * numTimes + timeid * 3;
            _res << "[" << dbuffer[off] << ", " << dbuffer[off + 1] << ", " << dbuffer[off + 2] << "], ";
        }

        _res << "..." << std::endl;
    }
    _res << "\t\t\t- ..." << std::endl;

    //------------------------------------------------------------------------------------------------------
    // wall shear stress vector per point over time : AXIAL
    //------------------------------------------------------------------------------------------------------
    //dbuffer.resize(numPoints * numTimes * 3);
    file.read(reinterpret_cast<char*>(dbuffer.data()), dbuffer.size() * sizeof(double));
    _res << "\t\t- Axial WSS vector per point per time:" << std::endl;

    for (unsigned int pointid = 0; pointid < NUM_DEMO; ++pointid)
    {
        _res << "\t\t\t- point" << pointid << ": ";

        for (unsigned int timeid = 0; timeid < NUM_DEMO; ++timeid)
        {
            const unsigned int off = pointid * numTimes + timeid * 3;
            _res << "[" << dbuffer[off] << ", " << dbuffer[off + 1] << ", " << dbuffer[off + 2] << "], ";
        }

        _res << "..." << std::endl;
    }
    _res << "\t\t\t- ..." << std::endl;

    //------------------------------------------------------------------------------------------------------
    // wall shear stress vector per point over time : CIRCUMFERENTIAL
    //------------------------------------------------------------------------------------------------------
    //dbuffer.resize(numPoints * numTimes * 3);
    file.read(reinterpret_cast<char*>(dbuffer.data()), dbuffer.size() * sizeof(double));
    _res << "\t\t- Circumferential WSS vector per point per time:" << std::endl;

    for (unsigned int pointid = 0; pointid < NUM_DEMO; ++pointid)
    {
        _res << "\t\t\t- point" << pointid << ": ";

        for (unsigned int timeid = 0; timeid < NUM_DEMO; ++timeid)
        {
            const unsigned int off = pointid * numTimes + timeid * 3;
            _res << "[" << dbuffer[off] << ", " << dbuffer[off + 1] << ", " << dbuffer[off + 2] << "], ";
        }

        _res << "..." << std::endl;
    }
    _res << "\t\t\t- ..." << std::endl;

    //------------------------------------------------------------------------------------------------------
    // mean wall shear stress per point
    //------------------------------------------------------------------------------------------------------
    dbuffer.resize(numPoints);
    file.read(reinterpret_cast<char*>(dbuffer.data()), dbuffer.size() * sizeof(double));
    _res << "\t\t- Mean WSS per point:" << std::endl;

    for (unsigned int pointid = 0; pointid < NUM_DEMO; ++pointid)
    { _res << "\t\t\t- point" << pointid << ": " << dbuffer[pointid] << std::endl; }
    _res << "\t\t\t- ..." << std::endl;

    //------------------------------------------------------------------------------------------------------
    // mean wall shear stress per point : AXIAL
    //------------------------------------------------------------------------------------------------------
    //dbuffer.resize(numPoints);
    file.read(reinterpret_cast<char*>(dbuffer.data()), dbuffer.size() * sizeof(double));
    _res << "\t\t- Mean axial WSS per point:" << std::endl;

    for (unsigned int pointid = 0; pointid < NUM_DEMO; ++pointid)
    { _res << "\t\t\t- point" << pointid << ": " << dbuffer[pointid] << std::endl; }
    _res << "\t\t\t- ..." << std::endl;

    //------------------------------------------------------------------------------------------------------
    // mean wall shear stress per point : CIRCUMFERENTIAL
    //------------------------------------------------------------------------------------------------------
    //dbuffer.resize(numPoints);
    file.read(reinterpret_cast<char*>(dbuffer.data()), dbuffer.size() * sizeof(double));
    _res << "\t\t- Mean circumferential WSS per point:" << std::endl;

    for (unsigned int pointid = 0; pointid < NUM_DEMO; ++pointid)
    { _res << "\t\t\t- point" << pointid << ": " << dbuffer[pointid] << std::endl; }
    _res << "\t\t\t- ..." << std::endl;

    //------------------------------------------------------------------------------------------------------
    // oscillatory shear index (OSI) per point
    //------------------------------------------------------------------------------------------------------
    //dbuffer.resize(numPoints);
    file.read(reinterpret_cast<char*>(dbuffer.data()), dbuffer.size() * sizeof(double));
    _res << "\t\t- OSI per point:" << std::endl;

    for (unsigned int pointid = 0; pointid < NUM_DEMO; ++pointid)
    { _res << "\t\t\t- point" << pointid << ": " << dbuffer[pointid] << std::endl; }
    _res << "\t\t\t- ..." << std::endl;

    //------------------------------------------------------------------------------------------------------
    // oscillatory shear index (OSI) per point : AXIAL
    //------------------------------------------------------------------------------------------------------
    //dbuffer.resize(numPoints);
    file.read(reinterpret_cast<char*>(dbuffer.data()), dbuffer.size() * sizeof(double));
    _res << "\t\t- Axial OSI per point:" << std::endl;

    for (unsigned int pointid = 0; pointid < NUM_DEMO; ++pointid)
    { _res << "\t\t\t- point" << pointid << ": " << dbuffer[pointid] << std::endl; }
    _res << "\t\t\t- ..." << std::endl;

    //------------------------------------------------------------------------------------------------------
    // oscillatory shear index (OSI) per point : CIRCUMFERENTIAL
    //------------------------------------------------------------------------------------------------------
    //dbuffer.resize(numPoints);
    file.read(reinterpret_cast<char*>(dbuffer.data()), dbuffer.size() * sizeof(double));
    _res << "\t\t- Circumferential OSI per point:" << std::endl;

    for (unsigned int pointid = 0; pointid < NUM_DEMO; ++pointid)
    { _res << "\t\t\t- point" << pointid << ": " << dbuffer[pointid] << std::endl; }
    _res << "\t\t\t- ..." << std::endl;

    //------------------------------------------------------------------------------------------------------
    // mean wall shear stress vector per point
    //------------------------------------------------------------------------------------------------------
    dbuffer.resize(numPoints * 3);
    file.read(reinterpret_cast<char*>(dbuffer.data()), dbuffer.size() * sizeof(double));
    _res << "\t\t- Mean WSS vector per point:" << std::endl;

    for (unsigned int pointid = 0; pointid < NUM_DEMO; ++pointid)
    {
        const unsigned int off = pointid * 3;
        _res << "\t\t\t- point" << pointid << ": [" << dbuffer[off] << ", " << dbuffer[off + 1] << ", " << dbuffer[off + 2] << "]" << std::endl;
    }
    _res << "\t\t\t- ..." << std::endl;

    //------------------------------------------------------------------------------------------------------
    //  mean wall shear stress vector per point : AXIAL
    //------------------------------------------------------------------------------------------------------
    //dbuffer.resize(numPoints * 3);
    file.read(reinterpret_cast<char*>(dbuffer.data()), dbuffer.size() * sizeof(double));
    _res << "\t\t- Mean axial WSS vector per point:" << std::endl;

    for (unsigned int pointid = 0; pointid < NUM_DEMO; ++pointid)
    {
        const unsigned int off = pointid * 3;
        _res << "\t\t\t- point" << pointid << ": [" << dbuffer[off] << ", " << dbuffer[off + 1] << ", " << dbuffer[off + 2] << "]" << std::endl;
    }
    _res << "\t\t\t- ..." << std::endl;

    //------------------------------------------------------------------------------------------------------
    // mean wall shear stress vector per point : CIRCUMFERENTIAL
    //------------------------------------------------------------------------------------------------------
    //dbuffer.resize(numPoints * 3);
    file.read(reinterpret_cast<char*>(dbuffer.data()), dbuffer.size() * sizeof(double));
    _res << "\t\t- Mean circumferential WSS vector per point:" << std::endl;

    for (unsigned int pointid = 0; pointid < NUM_DEMO; ++pointid)
    {
        const unsigned int off = pointid * 3;
        _res << "\t\t\t- point" << pointid << ": [" << dbuffer[off] << ", " << dbuffer[off + 1] << ", " << dbuffer[off + 2] << "]" << std::endl;
    }
    _res << "\t\t\t- ..." << std::endl;

    file.close();

    return true;
}

bool ImporterScientific::read_centerlines(std::string_view filepath)
{
    /*
     * [1] x [uint32] : numCenterlines
     *
     * for numCenterlines:
     *                          [1] x [uint32] : numPoints
     *              [numPoints * 3] x [double] : list of points
     *                  [numPoints] x [double] : vessel radius estimation per point
     *          [numPoints * 3 * 3] x [double] : local coordinate system (x,y,z vector) per point
     */

    if (!std::filesystem::exists(filepath.data()))
    {
        _res << "\t- vessel has no centerlines (path \"" << filepath.data() << "\")" << std::endl;
        return false;
    }

    _res << "\t- reading centerlines (path \"" << filepath.data() << "\")" << std::endl;

    std::ifstream file(filepath.data(), std::ios_base::in | std::ios_base::binary);

    if (!file.good())
    {
        _res << "\t\tFAILED! Could not open file!" << std::endl;
        return false;
    }

    std::vector<double> dbuffer;

    //------------------------------------------------------------------------------------------------------
    // num centerlines
    //------------------------------------------------------------------------------------------------------
    std::uint32_t numCenterlines = 0;
    file.read(reinterpret_cast<char*>(&numCenterlines), sizeof(std::uint32_t));
    _res << "\t- num. centerlines: " << numCenterlines << std::endl;

    for (unsigned int clid = 0; clid < numCenterlines; ++clid)
    {
        //------------------------------------------------------------------------------------------------------
        // num points
        //------------------------------------------------------------------------------------------------------
        std::uint32_t numPoints = 0;
        file.read(reinterpret_cast<char*>(&numPoints), sizeof(std::uint32_t));
        if (clid < NUM_DEMO)
        { _res << "\t\t- num. points of centerline " << clid << ": " << numPoints << std::endl; }

        //------------------------------------------------------------------------------------------------------
        // list of points
        //------------------------------------------------------------------------------------------------------
        dbuffer.resize(numPoints * 3);
        file.read(reinterpret_cast<char*>(dbuffer.data()), dbuffer.size() * sizeof(double));

        if (clid < NUM_DEMO)
        {
            for (unsigned int pointid = 0; pointid < std::min(NUM_DEMO, numPoints); ++pointid)
            {
                const unsigned int off = pointid * 3;
                _res << "\t\t\t- point" << pointid << ": [" << dbuffer[off] << ", " << dbuffer[off + 1] << ", " << dbuffer[off + 2] << "]" << std::endl;
            }
            _res << "\t\t\t- ..." << std::endl;
        }

        //------------------------------------------------------------------------------------------------------
        // vessel radius estimation per point
        //------------------------------------------------------------------------------------------------------
        dbuffer.resize(numPoints);
        file.read(reinterpret_cast<char*>(dbuffer.data()), dbuffer.size() * sizeof(double));

        if (clid < NUM_DEMO)
        {
            for (unsigned int pointid = 0; pointid < std::min(NUM_DEMO, numPoints); ++pointid)
            { _res << "\t\t\t- point" << pointid << " vessel radius [mm]: " << dbuffer[pointid] << std::endl; }
            _res << "\t\t\t- ..." << std::endl;
        }

        //------------------------------------------------------------------------------------------------------
        // local coordinate system (xyz vector) per point
        // - x/y are vectors in vessel's cross-section
        // - z is parallel to the centerline tangent
        //------------------------------------------------------------------------------------------------------
        dbuffer.resize(numPoints * 9);
        file.read(reinterpret_cast<char*>(dbuffer.data()), dbuffer.size() * sizeof(double));

        if (clid < NUM_DEMO)
        {
            for (unsigned int pointid = 0; pointid < std::min(NUM_DEMO, numPoints); ++pointid)
            {
                const unsigned int off = pointid * 9;
                _res << "\t\t\t- LCS at point" << pointid << ": ";
                _res << "X=[" << dbuffer[off + 0] << ", " << dbuffer[off + 1] << ", " << dbuffer[off + 2] << "], ";
                _res << "Y=[" << dbuffer[off + 3] << ", " << dbuffer[off + 4] << ", " << dbuffer[off + 5] << "], ";
                _res << "Z=[" << dbuffer[off + 6] << ", " << dbuffer[off + 7] << ", " << dbuffer[off + 8] << "]" << std::endl;
            }
            _res << "\t\t\t- ..." << std::endl;
        }
    } // for clid: num centerlines

    file.close();

    return true;
}

bool ImporterScientific::read_landmark_measuring_planes(std::string_view filepath)
{
    /*
     * [1] x [uint32] : num measuring planes
     * [1] x [uint32] : num measuring planes of land marks
     * for num measuring planes
     *      [measuring plane]
     * for num measuring planes of land marks
     *      [1] x [uint32] : land mark semantic
     *      [measuring plane]
     *
     * measuring plane:
     * [3] x [uint32] : grid size xyt
     * [3] x [double] : scale xyt
     * [3] x [double] : plane center xyz (world coordinates)
     * [3] x [double] : local coordinate system x axis (world coordinates)
     * [3] x [double] : local coordinate system y axis (world coordinates)
     * [3] x [double] : local coordinate system z axis (world coordinates) == normal
     * [1] x [double] : vessel diameter in mm
     * [size x * size y * size t * 3] x [double] : flow vector (rotated for usage in world space)
     * [size x * size y] x [uint8] : cross-section segmentation
     * [size x * size y * size t] x [double] : axial velocity
     * [size x * size y * size t] x [double] : circumferential velocity
     * [1] x [double] : min_flow_rate_per_time
     * [1] x [double] : max_flow_rate_per_time
     * [1] x [double] : mean_flow_rate_per_time
     * [1] x [double] : median_flow_rate_per_time
     * [1] x [double] : forward_flow_volume
     * [1] x [double] : backward_flow_volume
     * [1] x [double] : net_flow_volume
     * [1] x [double] : percentaged_back_flow_volume
     * [1] x [double] : cardiac_output
     * [1] x [double] : max_velocity
     * [1] x [double] : min_velocity
     * [1] x [double] : mean_velocity
     * [1] x [double] : median_velocity
     * [1] x [double] : min_velocity_axial
     * [1] x [double] : max_velocity_axial
     * [1] x [double] : mean_velocity_axial
     * [1] x [double] : median_velocity_axial
     * [1] x [double] : min_velocity_circumferential
     * [1] x [double] : max_velocity_circumferential
     * [1] x [double] : mean_velocity_circumferential
     * [1] x [double] : median_velocity_circumferential
     * [1] x [double] : area_mm2
     * [size t] x [double] : flow_rate_per_time
     * [size t] x [double] : areal_mean_velocity_per_time
     * [size t] x [double] : areal_mean_velocity_axial_per_time
     * [size t] x [double] : areal_mean_velocity_circumferential_per_time
     * [size t] x [double] : flow_jet_angle_per_time
     * [size t] x [double] : flow_jet_displacement_per_time
     * [size t] x [double] : flow_jet_high_velocity_area_percent_per_time
     * [1] x [double] : mp.max_flow_jet_angle_per_time
     * [1] x [double] : min_flow_jet_angle_per_time
     * [1] x [double] : mean_flow_jet_angle_per_time
     * [1] x [double] : median_flow_jet_angle_per_time
     * [1] x [double] : flow_jet_angle_at_fastest_time
     * [1] x [double] : mean_flow_jet_angle_velocity_weighted
     * [1] x [double] : min_flow_jet_displacement_per_time
     * [1] x [double] : max_flow_jet_displacement_per_time
     * [1] x [double] : mean_flow_jet_displacement_per_time
     * [1] x [double] : median_flow_jet_displacement_per_time
     * [1] x [double] : flow_jet_displacement_at_fastest_time
     * [1] x [double] : mean_flow_jet_displacement_velocity_weighted
     * [1] x [double] : min_flow_jet_high_velocity_area_percent_per_time
     * [1] x [double] : max_flow_jet_high_velocity_area_percent_per_time
     * [1] x [double] : mean_flow_jet_high_velocity_area_percent_per_time
     * [1] x [double] : median_flow_jet_high_velocity_area_percent_per_time
     * [1] x [double] : flow_jet_high_velocity_at_fastest_time
     * [1] x [double] : mean_flow_jet_high_velocity_velocity_weighted
     * [size t * 3] x [double] : flow_jet_position_per_time
     * [3] x [uint32] : numSamples for uncertainty calculation
     * [numSamples] x [double] : samples_net_flow_volume
     * [numSamples] x [double] : samples_forward_flow_volume
     * [numSamples] x [double] : samples_backward_flow_volume
     * [numSamples] x [double] : samples_percentaged_backward_flow_volume
     * [numSamples] x [double] : samples_cardiac_output
     */

    std::uint8_t ui8temp = 0;
    double dtemp = 0;

    const auto readMeasuringPlane = [&](std::ifstream& file)
    {
        //------------------------------------------------------------------------------------------------------
        // vessel id
        //------------------------------------------------------------------------------------------------------
        ui8temp = 0;
        file.read(reinterpret_cast<char*>(&ui8temp), sizeof(std::uint8_t));
        _res << "\t\t\t- vessel id: " << static_cast<int>(ui8temp) << std::endl;

        //------------------------------------------------------------------------------------------------------
        // grid size
        //   - x/y [0,1] in the plane + time steps [2]
        //------------------------------------------------------------------------------------------------------
        std::vector<std::uint32_t> gridsize(3);
        file.read(reinterpret_cast<char*>(gridsize.data()), gridsize.size() * sizeof(std::uint32_t));
        _res << "\t\t\t- grid size: [" << gridsize[0] << ", " << gridsize[1] << ", " << gridsize[2] << "]" << std::endl;

        //------------------------------------------------------------------------------------------------------
        // voxel scale
        //    - x/y [0,1] in the plane + temporal resolution [2]
        //------------------------------------------------------------------------------------------------------
        std::vector<double> dbuffer(3);
        file.read(reinterpret_cast<char*>(dbuffer.data()), dbuffer.size() * sizeof(double));
        _res << "\t\t\t- voxel scale: " << dbuffer[0] << " x " << dbuffer[1] << " [mm] / " << dbuffer[2] << " [ms]" << std::endl;

        //------------------------------------------------------------------------------------------------------
        // center
        //------------------------------------------------------------------------------------------------------
        // dbuffer.resize(3);
        file.read(reinterpret_cast<char*>(dbuffer.data()), dbuffer.size() * sizeof(double));
        _res << "\t\t\t- center: [" << dbuffer[0] << ", " << dbuffer[1] << ", " << dbuffer[2] << "]" << std::endl;

        //------------------------------------------------------------------------------------------------------
        // local coordinate system
        //   - orthonormal
        //   - nx/ny in the plane; nz is normal
        //------------------------------------------------------------------------------------------------------
        // dbuffer.resize(3);
        file.read(reinterpret_cast<char*>(dbuffer.data()), dbuffer.size() * sizeof(double));
        _res << "\t\t\t- LCS X: [" << dbuffer[0] << ", " << dbuffer[1] << ", " << dbuffer[2] << "]" << std::endl;

        // dbuffer.resize(3);
        file.read(reinterpret_cast<char*>(dbuffer.data()), dbuffer.size() * sizeof(double));
        _res << "\t\t\t- LCS Y: [" << dbuffer[0] << ", " << dbuffer[1] << ", " << dbuffer[2] << "]" << std::endl;

        // dbuffer.resize(3);
        file.read(reinterpret_cast<char*>(dbuffer.data()), dbuffer.size() * sizeof(double));
        _res << "\t\t\t- LCS Z: [" << dbuffer[0] << ", " << dbuffer[1] << ", " << dbuffer[2] << "]" << std::endl;

        //------------------------------------------------------------------------------------------------------
        // vessel diameter in mm
        //------------------------------------------------------------------------------------------------------
        file.read(reinterpret_cast<char*>(&dtemp), sizeof(double));
        _res << "\t\t\t- vessel diameter: " << dtemp << std::endl;

        //------------------------------------------------------------------------------------------------------
        // velocity vector per grid point
        //    - already rotated for use in world space and venc-scaled
        //------------------------------------------------------------------------------------------------------
        dbuffer.resize(gridsize[0] * gridsize[1] * gridsize[2] * 3);
        file.read(reinterpret_cast<char*>(dbuffer.data()), dbuffer.size() * sizeof(double));

        unsigned int cnt = 0;
        for (unsigned int x = 0; x < gridsize[0] && cnt < NUM_DEMO; ++x)
        {
            for (unsigned int y = 0; y < gridsize[1] && cnt < NUM_DEMO; ++y)
            {
                for (unsigned int t = 0; t < gridsize[2] && cnt < NUM_DEMO; ++t, ++cnt)
                {
                    const unsigned int off = x * gridsize[1] * gridsize[2] * 3 + y * gridsize[2] * 3 + t * 3;
                    _res << "\t\t\t- flow vector " << cnt << ": [" << dbuffer[off] << ", " << dbuffer[off+1] << ", " << dbuffer[off+2] << "]" << std::endl;
                } // for t
            } // for y
        } // for x

        //------------------------------------------------------------------------------------------------------
        // segmentation
        //    - static seg.; not time-dependent
        //------------------------------------------------------------------------------------------------------
        std::vector<std::uint8_t> ui8buffer(gridsize[0] * gridsize[1]);
        file.read(reinterpret_cast<char*>(ui8buffer.data()), ui8buffer.size() * sizeof(std::uint8_t));

        cnt = 0;
        for (unsigned int x = 0; x < gridsize[0] && cnt < NUM_DEMO; ++x)
        {
            for (unsigned int y = 0; y < gridsize[1] && cnt < NUM_DEMO; ++y)
            {
                const unsigned int off = x * gridsize[1] * +y;
                _res << "\t\t\t- seg value " << cnt << ": " << ui8buffer[off] << std::endl;
            } // for y
        } // for x

        //------------------------------------------------------------------------------------------------------
        // axial velocity per grid point
        //------------------------------------------------------------------------------------------------------
        dbuffer.resize(gridsize[0] * gridsize[1] * gridsize[2]);
        file.read(reinterpret_cast<char*>(dbuffer.data()), dbuffer.size() * sizeof(double));

        cnt = 0;
        for (unsigned int x = 0; x < gridsize[0] && cnt < NUM_DEMO; ++x)
        {
            for (unsigned int y = 0; y < gridsize[1] && cnt < NUM_DEMO; ++y)
            {
                for (unsigned int t = 0; t < gridsize[2] && cnt < NUM_DEMO; ++t, ++cnt)
                {
                    const unsigned int off = x * gridsize[1] * gridsize[2] + y * gridsize[2] + t;
                    _res << "\t\t\t- axial velocity " << cnt << ": " << dbuffer[off] << std::endl;
                } // for t
            } // for y
        } // for x

        //------------------------------------------------------------------------------------------------------
        // circumferential velocity per grid point
        //------------------------------------------------------------------------------------------------------
        dbuffer.resize(gridsize[0] * gridsize[1] * gridsize[2]);
        file.read(reinterpret_cast<char*>(dbuffer.data()), dbuffer.size() * sizeof(double));

        cnt = 0;
        for (unsigned int x = 0; x < gridsize[0] && cnt < NUM_DEMO; ++x)
        {
            for (unsigned int y = 0; y < gridsize[1] && cnt < NUM_DEMO; ++y)
            {
                for (unsigned int t = 0; t < gridsize[2] && cnt < NUM_DEMO; ++t, ++cnt)
                {
                    const unsigned int off = x * gridsize[1] * gridsize[2] + y * gridsize[2] + t;
                    _res << "\t\t\t- circumferential velocity " << cnt << ": " << dbuffer[off] << std::endl;
                } // for t
            } // for y
        } // for x

        //------------------------------------------------------------------------------------------------------
        // stats
        //------------------------------------------------------------------------------------------------------
        file.read(reinterpret_cast<char*>(&dtemp), sizeof(double));
        _res << "\t\t\t- min flow rate per time" << std::endl;

        file.read(reinterpret_cast<char*>(&dtemp), sizeof(double));
        _res << "\t\t\t- max flow rate per time" << std::endl;

        file.read(reinterpret_cast<char*>(&dtemp), sizeof(double));
        _res << "\t\t\t- mean flow rate per time" << std::endl;

        file.read(reinterpret_cast<char*>(&dtemp), sizeof(double));
        _res << "\t\t\t- median flow rate per time" << std::endl;

        file.read(reinterpret_cast<char*>(&dtemp), sizeof(double));
        _res << "\t\t\t- forward flow volume" << std::endl;

        file.read(reinterpret_cast<char*>(&dtemp), sizeof(double));
        _res << "\t\t\t- backward flow volume" << std::endl;

        file.read(reinterpret_cast<char*>(&dtemp), sizeof(double));
        _res << "\t\t\t- net flow volume" << std::endl;

        file.read(reinterpret_cast<char*>(&dtemp), sizeof(double));
        _res << "\t\t\t- percentaged back flow volume" << std::endl;

        file.read(reinterpret_cast<char*>(&dtemp), sizeof(double));
        _res << "\t\t\t- cardiac output" << std::endl;

        file.read(reinterpret_cast<char*>(&dtemp), sizeof(double));
        _res << "\t\t\t- max velocity" << std::endl;

        file.read(reinterpret_cast<char*>(&dtemp), sizeof(double));
        _res << "\t\t\t- min velocity" << std::endl;

        file.read(reinterpret_cast<char*>(&dtemp), sizeof(double));
        _res << "\t\t\t- mean velocity" << std::endl;

        file.read(reinterpret_cast<char*>(&dtemp), sizeof(double));
        _res << "\t\t\t- median velocity" << std::endl;

        file.read(reinterpret_cast<char*>(&dtemp), sizeof(double));
        _res << "\t\t\t- min velocity axial" << std::endl;

        file.read(reinterpret_cast<char*>(&dtemp), sizeof(double));
        _res << "\t\t\t- max velocity axial" << std::endl;

        file.read(reinterpret_cast<char*>(&dtemp), sizeof(double));
        _res << "\t\t\t- mean velocity axial" << std::endl;

        file.read(reinterpret_cast<char*>(&dtemp), sizeof(double));
        _res << "\t\t\t- median velocity axial" << std::endl;

        file.read(reinterpret_cast<char*>(&dtemp), sizeof(double));
        _res << "\t\t\t- min velocity circumferential" << std::endl;

        file.read(reinterpret_cast<char*>(&dtemp), sizeof(double));
        _res << "\t\t\t- max velocity circumferential" << std::endl;

        file.read(reinterpret_cast<char*>(&dtemp), sizeof(double));
        _res << "\t\t\t- mean velocity circumferential" << std::endl;

        file.read(reinterpret_cast<char*>(&dtemp), sizeof(double));
        _res << "\t\t\t- median velocity circumferential" << std::endl;

        file.read(reinterpret_cast<char*>(&dtemp), sizeof(double));
        _res << "\t\t\t- area mm2" << std::endl;

        dbuffer.resize(gridsize[2]);

        file.read(reinterpret_cast<char*>(dbuffer.data()), dbuffer.size() * sizeof(double));
        _res << "\t\t\t- flow rate per time" << std::endl;

        file.read(reinterpret_cast<char*>(dbuffer.data()), dbuffer.size() * sizeof(double));
        _res << "\t\t\t- areal mean velocity per time" << std::endl;

        file.read(reinterpret_cast<char*>(dbuffer.data()), dbuffer.size() * sizeof(double));
        _res << "\t\t\t- areal mean velocity axial per time" << std::endl;

        file.read(reinterpret_cast<char*>(dbuffer.data()), dbuffer.size() * sizeof(double));
        _res << "\t\t\t- areal mean velocity circumferential per time" << std::endl;

        /*
         * flow jet
         */
        file.read(reinterpret_cast<char*>(dbuffer.data()), dbuffer.size() * sizeof(double));
        _res << "\t\t\t- flow jet angle per time" << std::endl;

        file.read(reinterpret_cast<char*>(dbuffer.data()), dbuffer.size() * sizeof(double));
        _res << "\t\t\t- flow jet displacement per time" << std::endl;

        file.read(reinterpret_cast<char*>(dbuffer.data()), dbuffer.size() * sizeof(double));
        _res << "\t\t\t- flow jet high velocity area percent per time" << std::endl;

        file.read(reinterpret_cast<char*>(&dtemp), sizeof(double));
        _res << "\t\t\t- max flow jet angle per time" << std::endl;

        file.read(reinterpret_cast<char*>(&dtemp), sizeof(double));
        _res << "\t\t\t- min flow jet angle per time" << std::endl;

        file.read(reinterpret_cast<char*>(&dtemp), sizeof(double));
        _res << "\t\t\t- mean flow jet angle per time" << std::endl;

        file.read(reinterpret_cast<char*>(&dtemp), sizeof(double));
        _res << "\t\t\t- median flow jet angle per time" << std::endl;

        file.read(reinterpret_cast<char*>(&dtemp), sizeof(double));
        _res << "\t\t\t- flow jet angle at fastest time" << std::endl;

        file.read(reinterpret_cast<char*>(&dtemp), sizeof(double));
        _res << "\t\t\t- mean flow jet angle velocity weighted" << std::endl;

        file.read(reinterpret_cast<char*>(&dtemp), sizeof(double));
        _res << "\t\t\t- min flow jet displacement per time" << std::endl;

        file.read(reinterpret_cast<char*>(&dtemp), sizeof(double));
        _res << "\t\t\t- max flow jet displacement per time" << std::endl;

        file.read(reinterpret_cast<char*>(&dtemp), sizeof(double));
        _res << "\t\t\t- mean flow jet displacement per time" << std::endl;

        file.read(reinterpret_cast<char*>(&dtemp), sizeof(double));
        _res << "\t\t\t- median flow jet displacement per time" << std::endl;

        file.read(reinterpret_cast<char*>(&dtemp), sizeof(double));
        _res << "\t\t\t- flow jet displacement at fastest time" << std::endl;

        file.read(reinterpret_cast<char*>(&dtemp), sizeof(double));
        _res << "\t\t\t- mean flow jet displacement velocity weighted" << std::endl;

        file.read(reinterpret_cast<char*>(&dtemp), sizeof(double));
        _res << "\t\t\t- min flow jet high velocity area percent per time" << std::endl;

        file.read(reinterpret_cast<char*>(&dtemp), sizeof(double));
        _res << "\t\t\t- max flow jet high velocity area percent per time" << std::endl;

        file.read(reinterpret_cast<char*>(&dtemp), sizeof(double));
        _res << "\t\t\t- mean flow jet high velocity area percent per time" << std::endl;

        file.read(reinterpret_cast<char*>(&dtemp), sizeof(double));
        _res << "\t\t\t- median flow jet high velocity area percent per time" << std::endl;

        file.read(reinterpret_cast<char*>(&dtemp), sizeof(double));
        _res << "\t\t\t- flow jet high velocity at fastest time" << std::endl;

        file.read(reinterpret_cast<char*>(&dtemp), sizeof(double));
        _res << "\t\t\t- mean flow jet high velocity velocity weighted" << std::endl;

        dbuffer.resize(gridsize[2] * 3);
        file.read(reinterpret_cast<char*>(dbuffer.data()), dbuffer.size() * sizeof(double));

        for (unsigned int t = 0; t < gridsize[2]; ++t)
        {
            if (t < NUM_DEMO)
            {
                const unsigned int off = t * 3;
                _res << "\t\t\t- flow jet position per time " << t << ": [" << dbuffer[off] << ", " << dbuffer[off + 1] << ", " << dbuffer[off + 2] << "]" << std::endl;
            }
        }
        _res << "\t\t\t- ..." << std::endl;

        //------------------------------------------------------------------------------------------------------
        // uncertainty
        //------------------------------------------------------------------------------------------------------
        /*
         * numSamples
         */
        std::uint32_t numSamples = 0;
        file.read(reinterpret_cast<char*>(&numSamples), sizeof(std::uint32_t));

        dbuffer.resize(numSamples);

        /*
         * samples: net flow volume
         */
        file.read(reinterpret_cast<char*>(dbuffer.data()), dbuffer.size() * sizeof(double));

        _res << "\t\t\t- samples net flow volume: ";
        for (unsigned int i = 0; i < NUM_DEMO; ++i)
        { _res << dbuffer[i] << ", "; }
        _res << "..." << std::endl;

        /*
         * samples: forward flow volume
         */
        file.read(reinterpret_cast<char*>(dbuffer.data()), dbuffer.size() * sizeof(double));

        _res << "\t\t\t- samples forward flow volume: ";
        for (unsigned int i = 0; i < NUM_DEMO; ++i)
        { _res << dbuffer[i] << ", "; }
        _res << "..." << std::endl;

        /*
         * samples: backward flow volume
         */
        file.read(reinterpret_cast<char*>(dbuffer.data()), dbuffer.size() * sizeof(double));

        _res << "\t\t\t- samples backward flow volume: ";
        for (unsigned int i = 0; i < NUM_DEMO; ++i)
        { _res << dbuffer[i] << ", "; }
        _res << "..." << std::endl;

        /*
         * samples: percentaged back flow volume
         */
        file.read(reinterpret_cast<char*>(dbuffer.data()), dbuffer.size() * sizeof(double));

        _res << "\t\t\t- samples percentaged backward flow volume: ";
        for (unsigned int i = 0; i < NUM_DEMO; ++i)
        { _res << dbuffer[i] << ", "; }
        _res << "..." << std::endl;

        /*
         * samples: cardiac output
         */
        file.read(reinterpret_cast<char*>(dbuffer.data()), dbuffer.size() * sizeof(double));

        _res << "\t\t\t- samples cardiac output: ";
        for (unsigned int i = 0; i < NUM_DEMO; ++i)
        { _res << dbuffer[i] << ", "; }
        _res << "..." << std::endl;
    }; // readMeasuringPlane()


    if (!std::filesystem::exists(filepath.data()))
    {
        _res << "\t- vessel has no land marks of measuring planes (path \"" << filepath.data() << "\")" << std::endl;
        _res.flush();
        return false;
    }

    _res << "\t- reading land marks of measuring planes (path \"" << filepath.data() << "\")" << std::endl;

    std::ifstream file(filepath.data(), std::ios_base::in | std::ios_base::binary);

    if (!file.good())
    {
        _res << "\t\tFAILED! Could not open file!" << std::endl;
        return false;
    }

    // num measuring planes
    std::uint32_t numMeasuringPlanes = 0;
    file.read(reinterpret_cast<char*>(&numMeasuringPlanes), sizeof(std::uint32_t));
    std::cout << "\t\t- num. measuring planes: " << numMeasuringPlanes << std::endl;

    // num measuring planes of land marks
    std::uint32_t numMeasuringPlanesOfLandMarks = 0;
    file.read(reinterpret_cast<char*>(&numMeasuringPlanesOfLandMarks), sizeof(std::uint32_t));
    std::cout << "\t\t- num. measuring planes of landmarks: " << numMeasuringPlanesOfLandMarks << std::endl;

    // measuring planes
    for (unsigned int i = 0; i < numMeasuringPlanes; ++i)
    {
        std::cout << "\t\t- measuring plane " <<i<<": "<< std::endl;
        readMeasuringPlane(file);
    }

    // measuring planes of land marks
    for (unsigned int i = 0; i < numMeasuringPlanesOfLandMarks; ++i)
    {
        std::cout << "\t\t- measuring plane " <<i<<" of landmarks: "<< std::endl;

        std::uint32_t semantic = 0;
        file.read(reinterpret_cast<char*>(&semantic), sizeof(std::uint32_t));

        std::cout << "\t\t\t- semantic: " << semantic << " (";
        switch (semantic)
        {
            case 1: std::cout << "LandMarkSemantic_Aorta_AboveAorticValve"; break;
            case 2: std::cout << "LandMarkSemantic_Aorta_MidAscendingAorta"; break;
            case 3: std::cout << "LandMarkSemantic_Aorta_BeforeBrachiocephalicArtery"; break;
            case 4: std::cout << "LandMarkSemantic_Aorta_BetweenLeftCommonCarotid_and_LeftSubclavianArtery"; break;
            case 5: std::cout << "LandMarkSemantic_Aorta_DistalToLeftSubclavianArtery"; break;
            case 6: std::cout << "LandMarkSemantic_Aorta_MidDescendingAorta"; break;
            case 7: std::cout << "LandMarkSemantic_PulmonaryArtery_AbovePulmonaryValve"; break;
            case 8: std::cout << "LandMarkSemantic_PulmonaryArtery_BeforeJunction"; break;
            case 9: std::cout << "LandMarkSemantic_PulmonaryArtery_LeftPulmonaryArtery_Begin"; break;
            case 10: std::cout << "LandMarkSemantic_PulmonaryArtery_RightPulmonaryArtery_Begin"; break;
            default: std::cout << "None";
        }
        std::cout << ")" << std::endl;

        readMeasuringPlane(file);
    }

    file.close();

    return true;
}

bool ImporterScientific::read_pathlines(std::string_view filepath)
{
    /*
     *             [1] x [uint32] : numPathlines
     *             
     *  for numPathlines
     *                  [1] x [uint32] : numPoints
     *      [numPoints * 4] x [double] : list of points (xyz+time)
     *          [numPoints] x [double] : attribute per point : relative pressure
     *          [numPoints] x [double] : attribute per point : cos(angle) between pathline and centerline tangent
     *          [numPoints] x [double] : attribute per point : rotation direction
     *          [numPoints] x [double] : attribute per point : velocity
     *          [numPoints] x [double] : attribute per point : axial velocity
     *                  [1] x [double] : attribute : length
     */

    if (!std::filesystem::exists(filepath.data()))
    {
        _res << "\t- vessel has no pathlines (path \"" << filepath.data() << "\")" << std::endl;
        _res.flush();
        return false;
    }

    _res << "\t- reading pathlines (path \"" << filepath.data() << "\")" << std::endl;

    std::ifstream file(filepath.data(), std::ios_base::in | std::ios_base::binary);

    if (!file.good())
    {
        _res << "\t\tFAILED! Could not open file!" << std::endl;
        return false;
    }

    std::vector<double> dbuffer;

    //------------------------------------------------------------------------------------------------------
    // num pathlines
    //------------------------------------------------------------------------------------------------------
    std::uint32_t numPathines = 0;
    file.read(reinterpret_cast<char*>(&numPathines), sizeof(std::uint32_t));
    _res << "\t- num. pathlines: " << numPathines << std::endl;

    for (unsigned int plid = 0; plid < numPathines; ++plid)
    {
        //------------------------------------------------------------------------------------------------------
        // num points
        //------------------------------------------------------------------------------------------------------
        std::uint32_t numPoints = 0;
        file.read(reinterpret_cast<char*>(&numPoints), sizeof(std::uint32_t));
        if (plid < NUM_DEMO)
        { _res << "\t\t- num. points of pathline" << plid << ": " << numPoints << std::endl; }

        //------------------------------------------------------------------------------------------------------
        // list of points (xyz + time)
        //------------------------------------------------------------------------------------------------------
        dbuffer.resize(numPoints * 4); // x y z time
        file.read(reinterpret_cast<char*>(dbuffer.data()), dbuffer.size() * sizeof(double));

        if (plid < NUM_DEMO)
        {
            for (unsigned int pointid = 0; pointid < std::min(NUM_DEMO, numPoints); ++pointid)
            {
                const unsigned int off = pointid * 4;
                _res << "\t\t\t- point" << pointid << ": [" << /*x=*/dbuffer[off] << ", " << /*y=*/dbuffer[off + 1] << ", " << /*z=*/dbuffer[off + 2] << ", " << /*t=*/dbuffer[off + 3] << "]"
                     << std::endl;
            }
            _res << "\t\t\t- ..." << std::endl;
        }

        //------------------------------------------------------------------------------------------------------
        // attribute: relative pressure per point
        //------------------------------------------------------------------------------------------------------
        dbuffer.resize(numPoints);
        file.read(reinterpret_cast<char*>(dbuffer.data()), dbuffer.size() * sizeof(double));

        if (plid < NUM_DEMO)
        {
            for (unsigned int pointid = 0; pointid < std::min(NUM_DEMO, numPoints); ++pointid)
            { _res << "\t\t\t- relative pressure [mmHg] of point" << pointid << ": " << dbuffer[pointid] << std::endl; }
            _res << "\t\t\t- ..." << std::endl;
        }

        //------------------------------------------------------------------------------------------------------
        // attribute: cos(angle) between pathline tangent and centerline per point
        //------------------------------------------------------------------------------------------------------
        //dbuffer.resize(numPoints);
        file.read(reinterpret_cast<char*>(dbuffer.data()), dbuffer.size() * sizeof(double));

        if (plid < NUM_DEMO)
        {
            for (unsigned int pointid = 0; pointid < std::min(NUM_DEMO, numPoints); ++pointid)
            { _res << "\t\t\t- cos(angle) pathline/centerline tangent of point" << pointid << ": " << dbuffer[pointid] << std::endl; }
            _res << "\t\t\t- ..." << std::endl;
        }

        //------------------------------------------------------------------------------------------------------
        // attribute: rotation direction per point
        //------------------------------------------------------------------------------------------------------
        //dbuffer.resize(numPoints);
        file.read(reinterpret_cast<char*>(dbuffer.data()), dbuffer.size() * sizeof(double));

        if (plid < NUM_DEMO)
        {
            for (unsigned int pointid = 0; pointid < std::min(NUM_DEMO, numPoints); ++pointid)
            { _res << "\t\t\t- rotation direction of point" << pointid << ": " << dbuffer[pointid] << std::endl; }
            _res << "\t\t\t- ..." << std::endl;
        }

        //------------------------------------------------------------------------------------------------------
        // attribute: velocity per point
        //------------------------------------------------------------------------------------------------------
        //dbuffer.resize(numPoints);
        file.read(reinterpret_cast<char*>(dbuffer.data()), dbuffer.size() * sizeof(double));

        if (plid < NUM_DEMO)
        {
            for (unsigned int pointid = 0; pointid < std::min(NUM_DEMO, numPoints); ++pointid)
            { _res << "\t\t\t- velocity [m/s] at point" << pointid << ": " << dbuffer[pointid] << std::endl; }
            _res << "\t\t\t- ..." << std::endl;
        }

        //------------------------------------------------------------------------------------------------------
        // attribute: axial velocity per point
        //------------------------------------------------------------------------------------------------------
        //dbuffer.resize(numPoints);
        file.read(reinterpret_cast<char*>(dbuffer.data()), dbuffer.size() * sizeof(double));

        if (plid < NUM_DEMO)
        {
            for (unsigned int pointid = 0; pointid < std::min(NUM_DEMO, numPoints); ++pointid)
            { _res << "\t\t\t- axial velocity [m/s] at point" << pointid << ": " << dbuffer[pointid] << std::endl; }
            _res << "\t\t\t- ..." << std::endl;
        }

        //------------------------------------------------------------------------------------------------------
        // pathline length (spatial; temporal component is ignored)
        //------------------------------------------------------------------------------------------------------
        double length = 0;
        file.read(reinterpret_cast<char*>(&length), sizeof(double));

        if (plid < NUM_DEMO)
        { _res << "\t\t\t- spatial length [mm]: " << length << std::endl; }
    } // for plid: num pathlines

    file.close();

    return true;
}

bool ImporterScientific::read_flowfield(std::string_view filepath)
{
    /*
     *                                 [4] x [uint32] : size x y z t
     *                                 [4] x [double] : scale x y z t
     *                                [16] x [double] : world matrix (4x4 from dicom)
     *                                [16] x [double] : inverse world matrix (4x4 from dicom)
     *                                [25] x [double] : world matrix (5x5 including time in 4th row/col)
     *                                [25] x [double] : inverse world matrix (5x5 including time in 4th row/col)
     *                                 [9] x [double] : rotational part of world matrix (3x3; used below to transform the velocity vectors to world space)
     *                                 [9] x [double] : inverse rotational part of world matrix (3x3)
     * [sizeX * sizeY * sizeZ * sizeT * 3] x [double] : flow vectors (rotated in world coordinates)
     */

    if (!std::filesystem::exists(filepath.data()))
    {
        _res << "\t- no flow field (path \"" << filepath.data() << "\")" << std::endl;
        return false;
    }

    _res << "\t- reading flow field (path \"" << filepath.data() << "\")" << std::endl;

    std::ifstream file(filepath.data(), std::ios_base::in | std::ios_base::binary);

    if (!file.good())
    {
        _res << "\t\tFAILED! Could not open file!" << std::endl;
        return false;
    }

    //------------------------------------------------------------------------------------------------------
    // grid size x y z t
    //------------------------------------------------------------------------------------------------------
    std::vector<std::uint32_t> gridsize(4);
    file.read(reinterpret_cast<char*>(gridsize.data()), gridsize.size() * sizeof(std::uint32_t));
    _res << "\t\t- grid size: " << gridsize[0] << " x " << gridsize[1] << " x " << gridsize[2] << " x " << gridsize[3] << std::endl;

    //------------------------------------------------------------------------------------------------------
    // voxel scale x y z t
    //------------------------------------------------------------------------------------------------------
    std::vector<double> dbuffer(4);
    file.read(reinterpret_cast<char*>(dbuffer.data()), dbuffer.size() * sizeof(double));
    _res << "\t\t- voxel scale: " << dbuffer[0] << " x " << dbuffer[1] << " x " << dbuffer[2] << " mm / " << dbuffer[3] << " ms" << std::endl;

    //------------------------------------------------------------------------------------------------------
    // world matrix (4x4 from dicom)
    //------------------------------------------------------------------------------------------------------
    dbuffer.resize(16); // 4 x 4
    file.read(reinterpret_cast<char*>(dbuffer.data()), dbuffer.size() * sizeof(double));

    unsigned int cnt = 0;
    _res << "\t\t- world matrix:" << std::endl;
    for (unsigned int rowid = 0; rowid < 4; ++rowid)
    {
        _res << "\t\t\t";
        for (unsigned int colid = 0; colid < 4; ++colid)
        { _res << dbuffer[cnt++] << " "; }

        _res << std::endl;
    }

    //------------------------------------------------------------------------------------------------------
    // inverse world matrix
    //------------------------------------------------------------------------------------------------------
    //dbuffer.resize(16); // 4 x 4
    file.read(reinterpret_cast<char*>(dbuffer.data()), dbuffer.size() * sizeof(double));

    cnt = 0;
    _res << "\t\t- inverse world matrix:" << std::endl;
    for (unsigned int rowid = 0; rowid < 4; ++rowid)
    {
        _res << "\t\t\t";
        for (unsigned int colid = 0; colid < 4; ++colid)
        { _res << dbuffer[cnt++] << " "; }

        _res << std::endl;
    }

    //------------------------------------------------------------------------------------------------------
    // world matrix with time (5x5 including time in 4th row/col)
    //------------------------------------------------------------------------------------------------------
    dbuffer.resize(25); // 5 x 5
    file.read(reinterpret_cast<char*>(dbuffer.data()), dbuffer.size() * sizeof(double));

    cnt = 0;
    _res << "\t\t- world matrix with time:" << std::endl;
    for (unsigned int rowid = 0; rowid < 5; ++rowid)
    {
        _res << "\t\t\t";
        for (unsigned int colid = 0; colid < 5; ++colid)
        { _res << dbuffer[cnt++] << " "; }

        _res << std::endl;
    }

    //------------------------------------------------------------------------------------------------------
    // inverse world matrix with time
    //------------------------------------------------------------------------------------------------------
    //dbuffer.resize(25); // 5 x 5
    file.read(reinterpret_cast<char*>(dbuffer.data()), dbuffer.size() * sizeof(double));

    cnt = 0;
    _res << "\t\t- inverse world matrix with time:" << std::endl;
    for (unsigned int rowid = 0; rowid < 5; ++rowid)
    {
        _res << "\t\t\t";
        for (unsigned int colid = 0; colid < 5; ++colid)
        { _res << dbuffer[cnt++] << " "; }

        _res << std::endl;
    }

    //------------------------------------------------------------------------------------------------------
    // rotational part of world matrix (3x3; used to transform velocity vectors to world space)
    //------------------------------------------------------------------------------------------------------
    dbuffer.resize(9); // 3 x 3
    file.read(reinterpret_cast<char*>(dbuffer.data()), dbuffer.size() * sizeof(double));

    cnt = 0;
    _res << "\t\t- rotational part of world matrix:" << std::endl;
    for (unsigned int rowid = 0; rowid < 3; ++rowid)
    {
        _res << "\t\t\t";
        for (unsigned int colid = 0; colid < 3; ++colid)
        { _res << dbuffer[cnt++] << " "; }

        _res << std::endl;
    }

    //------------------------------------------------------------------------------------------------------
    // inverse rotational part of world matrix
    //------------------------------------------------------------------------------------------------------
    //dbuffer.resize(9); // 3 x 3
    file.read(reinterpret_cast<char*>(dbuffer.data()), dbuffer.size() * sizeof(double));

    cnt = 0;
    _res << "\t\t- inverse rotational part of world matrix:" << std::endl;
    for (unsigned int rowid = 0; rowid < 3; ++rowid)
    {
        _res << "\t\t\t";
        for (unsigned int colid = 0; colid < 3; ++colid)
        { _res << dbuffer[cnt++] << " "; }

        _res << std::endl;
    }

    //------------------------------------------------------------------------------------------------------
    // flow vectors
    // - already rotated in world coordinates
    // - already venc-scaled 
    //------------------------------------------------------------------------------------------------------
    dbuffer.resize(gridsize[0] * gridsize[1] * gridsize[2] * gridsize[3] * 3);
    file.read(reinterpret_cast<char*>(dbuffer.data()), dbuffer.size() * sizeof(double));

    cnt = 0;
    unsigned int cntDemo = 0;
    for (unsigned int x = 0; x < gridsize[0] && cntDemo < NUM_DEMO; ++x)
    {
        for (unsigned int y = 0; y < gridsize[1] && cntDemo < NUM_DEMO; ++y)
        {
            for (unsigned int z = 0; z < gridsize[2] && cntDemo < NUM_DEMO; ++z)
            {
                for (unsigned int t = 0; t < gridsize[3] && cntDemo < NUM_DEMO; ++t, ++cntDemo)
                {
                    const unsigned int off = x * gridsize[1] * gridsize[2] * gridsize[3] * 3 + y * gridsize[2] * gridsize[3] * 3 + z * gridsize[3] * 3 + t * 3;
                    _res << "\t\t- flow vector" << (cnt++) << " [m/s]: " << "[" << dbuffer[off] << ", " << dbuffer[off + 1] << ", " << dbuffer[off + 2] << "]" << std::endl;
                }
            }
        }
    }
    _res << "\t\t- ..." << std::endl;

    file.close();

    return true;
}

bool ImporterScientific::read_pressure_map(std::string_view filepath)
{
    if (!std::filesystem::exists(filepath.data()))
    {
        _res << "\t- no pressure map (path \"" << filepath.data() << "\")" << std::endl;
        return false;
    }

    _res << "\t- reading pressure map (path \"" << filepath.data() << "\")" << std::endl;

    std::ifstream file(filepath.data(), std::ios_base::in | std::ios_base::binary);

    if (!file.good())
    {
        _res << "\t\tFAILED! Could not open file!" << std::endl;
        return false;
    }

    _read_nd_scalar_image_in_sparse_matrix_style(file);

    file.close();

    return true;
}

bool ImporterScientific::read_rotation_direction_map(std::string_view filepath)
{

    if (!std::filesystem::exists(filepath.data()))
    {
        _res << "\t- no rotation direction map (path \"" << filepath.data() << "\")" << std::endl;
        return false;
    }

    _res << "\t- reading rotation direction map (path \"" << filepath.data() << "\")" << std::endl;

    std::ifstream file(filepath.data(), std::ios_base::in | std::ios_base::binary);

    if (!file.good())
    {
        _res << "\t\tFAILED! Could not open file!" << std::endl;
        return false;
    }

    _read_nd_scalar_image_in_sparse_matrix_style(file);

    file.close();

    return true;
}

bool ImporterScientific::read_axial_velocity_map(std::string_view filepath)
{

    if (!std::filesystem::exists(filepath.data()))
    {
        _res << "\t- no axial velocity map (path \"" << filepath.data() << "\")" << std::endl;
        return false;
    }

    _res << "\t- reading axial velocity map (path \"" << filepath.data() << "\")" << std::endl;

    std::ifstream file(filepath.data(), std::ios_base::in | std::ios_base::binary);

    if (!file.good())
    {
        _res << "\t\tFAILED! Could not open file!" << std::endl;
        return false;
    }

    _read_nd_scalar_image_in_sparse_matrix_style(file);

    file.close();

    return true;
}

bool ImporterScientific::read_cos_angle_to_centerline_map(std::string_view filepath)
{

    if (!std::filesystem::exists(filepath.data()))
    {
        _res << "\t- no cos(angle) to centerline (path \"" << filepath.data() << "\")" << std::endl;
        return false;
    }

    _res << "\t- reading cos(angle) to centerline (path \"" << filepath.data() << "\")" << std::endl;

    std::ifstream file(filepath.data(), std::ios_base::in | std::ios_base::binary);

    if (!file.good())
    {
        _res << "\t\tFAILED! Could not open file!" << std::endl;
        return false;
    }

    _read_nd_scalar_image_in_sparse_matrix_style(file);

    file.close();

    return true;
}

bool ImporterScientific::read_turbulent_kinetic_energy_map(std::string_view filepath)
{

    if (!std::filesystem::exists(filepath.data()))
    {
        _res << "\t- no turbulent kinetic energy map (path \"" << filepath.data() << "\")" << std::endl;
        return false;
    }

    _res << "\t- reading turbulent kinetic energy map (path \"" << filepath.data() << "\")" << std::endl;

    std::ifstream file(filepath.data(), std::ios_base::in | std::ios_base::binary);

    if (!file.good())
    {
        _res << "\t\tFAILED! Could not open file!" << std::endl;
        return false;
    }

    _read_nd_scalar_image_in_sparse_matrix_style(file);

    file.close();

    return true;
}

bool ImporterScientific::read_flow_jet(std::string_view filepath)
{
    /*
     * [1] x [uint32] : numFlowjets
     *
     * for numFlowjets:
     *         [1] x [uint32] : numPoints
     *         [1] x [uint32] : numTimes
     *
     *         for numPoints:
     *                 for numTimes:
     *                         [3] x [double] : position of the flow jet tube, which represents the peak velocity position in the cross-section
     *                         [1] x [double] : peak velocity in the cross-section
     *                         [3] x [double] : center of the magenta high-velocity area net
     *                         [3] x [double] : direction 0 of the high-velocity area
     *                         [1] x [double] : radius 0 of the high-velocity area
     *                         [3] x [double] : direction 1 of the high-velocity area
     *                         [1] x [double] : radius 1 of the high-velocity area
     *
     *                 [3] x [double] : vessel center (centerline position) for this cross-section
     *                 [1] x [double] : vessel radius
     *                 [3] x [double] : x direction of centerline's local coordinate system
     *                 [3] x [double] : y direction of centerline's local coordinate system
     */

    if (!std::filesystem::exists(filepath.data()))
    {
        _res << "\t- no flow jet (path \"" << filepath.data() << "\")" << std::endl;
        return false;
    }

    _res << "\t- reading flow jet (path \"" << filepath.data() << "\")" << std::endl;

    std::ifstream file(filepath.data(), std::ios_base::in | std::ios_base::binary);

    if (!file.good())
    {
        _res << "\t\tFAILED! Could not open file!" << std::endl;
        return false;
    }

    //------------------------------------------------------------------------------------------------------
    // num flowjets
    //------------------------------------------------------------------------------------------------------
    std::uint32_t numFlowjets = 0;
    file.read(reinterpret_cast<char*>(&numFlowjets), sizeof(std::uint32_t));
    _res << "\t\t- num. flow jets: " << numFlowjets << std::endl;

    for (unsigned int fjid = 0; fjid < numFlowjets; ++fjid)
    {
        _res << "\t\t- flow jet " << fjid << ":" << std::endl;

        //------------------------------------------------------------------------------------------------------
        // num points
        //------------------------------------------------------------------------------------------------------
        std::uint32_t numPoints = 0;
        file.read(reinterpret_cast<char*>(&numPoints), sizeof(std::uint32_t));
        _res << "\t\t\t- num points: " << numPoints << std::endl;

        //------------------------------------------------------------------------------------------------------
        // num times
        //------------------------------------------------------------------------------------------------------
        std::uint32_t numTimes = 0;
        file.read(reinterpret_cast<char*>(&numTimes), sizeof(std::uint32_t));
        _res << "\t\t\t- num times: " << numTimes << std::endl;

        //------------------------------------------------------------------------------------------------------
        // list of points per time
        //------------------------------------------------------------------------------------------------------
        std::vector<double> dbuffer;
        unsigned int cnt = 0;

        for (unsigned int pointid = 0; pointid < numPoints; ++pointid)
        {
            for (unsigned int timeid = 0; timeid < numTimes; ++timeid)
            {
                if (cnt < NUM_DEMO && timeid < NUM_DEMO)
                { _res << "\t\t\t- point " << pointid << " time " << timeid << std::endl; }

                //------------------------------------------------------------------------------------------------------
                // position of the flow jet tube, which represents the peak velocity position in the cross-section
                //------------------------------------------------------------------------------------------------------
                dbuffer.resize(3);
                file.read(reinterpret_cast<char*>(dbuffer.data()), dbuffer.size() * sizeof(double));
                if (cnt < NUM_DEMO && timeid < NUM_DEMO)
                { _res << "\t\t\t\t- peak velocity position: [" << dbuffer[0] << ", " << dbuffer[1] << ", " << dbuffer[2] << "]" << std::endl; }

                //------------------------------------------------------------------------------------------------------
                // peak velocity in the cross-section
                //------------------------------------------------------------------------------------------------------
                dbuffer.resize(1);
                file.read(reinterpret_cast<char*>(dbuffer.data()), dbuffer.size() * sizeof(double));
                if (cnt < NUM_DEMO && timeid < NUM_DEMO)
                { _res << "\t\t\t\t- peak velocity [m/s]: " << dbuffer[0] << std::endl; }

                //------------------------------------------------------------------------------------------------------
                // center of the magenta high-velocity area net
                //------------------------------------------------------------------------------------------------------
                dbuffer.resize(3);
                file.read(reinterpret_cast<char*>(dbuffer.data()), dbuffer.size() * sizeof(double));
                if (cnt < NUM_DEMO && timeid < NUM_DEMO)
                { _res << "\t\t\t\t- area center: [" << dbuffer[0] << ", " << dbuffer[1] << ", " << dbuffer[2] << "]" << std::endl; }

                //------------------------------------------------------------------------------------------------------
                // direction 0 + radius of the high-velocity area
                //------------------------------------------------------------------------------------------------------
                dbuffer.resize(3);
                file.read(reinterpret_cast<char*>(dbuffer.data()), dbuffer.size() * sizeof(double));
                if (cnt < NUM_DEMO && timeid < NUM_DEMO)
                { _res << "\t\t\t\t- area dir0: [" << dbuffer[0] << ", " << dbuffer[1] << ", " << dbuffer[2] << "]" << std::endl; }

                dbuffer.resize(1);
                file.read(reinterpret_cast<char*>(dbuffer.data()), dbuffer.size() * sizeof(double));
                if (cnt < NUM_DEMO && timeid < NUM_DEMO)
                { _res << "\t\t\t\t- area radius0 [mm]: " << dbuffer[0] << std::endl; }

                //------------------------------------------------------------------------------------------------------
                // direction 1 + radius of the high-velocity area
                //------------------------------------------------------------------------------------------------------
                dbuffer.resize(3);
                file.read(reinterpret_cast<char*>(dbuffer.data()), dbuffer.size() * sizeof(double));
                if (cnt < NUM_DEMO && timeid < NUM_DEMO)
                { _res << "\t\t\t\t- area dir1: [" << dbuffer[0] << ", " << dbuffer[1] << ", " << dbuffer[2] << "]" << std::endl; }

                dbuffer.resize(1);
                file.read(reinterpret_cast<char*>(dbuffer.data()), dbuffer.size() * sizeof(double));
                if (cnt < NUM_DEMO && timeid < NUM_DEMO)
                { _res << "\t\t\t\t- area radius1 [mm]: " << dbuffer[0] << std::endl; }
            } // for times

            if (cnt < NUM_DEMO)
            { _res << "\t\t\t- ..." << std::endl; }

            //------------------------------------------------------------------------------------------------------
            // vessel center (centerline position) + radius for this cross-section
            //------------------------------------------------------------------------------------------------------
            dbuffer.resize(3);
            file.read(reinterpret_cast<char*>(dbuffer.data()), dbuffer.size() * sizeof(double));
            if (cnt < NUM_DEMO)
            { _res << "\t\t\t- vessel center: [" << dbuffer[0] << ", " << dbuffer[1] << ", " << dbuffer[2] << "]" << std::endl; }

            dbuffer.resize(1);
            file.read(reinterpret_cast<char*>(dbuffer.data()), dbuffer.size() * sizeof(double));
            if (cnt < NUM_DEMO)
            { _res << "\t\t\t- vessel radius [mm]: " << dbuffer[0] << std::endl; }

            //------------------------------------------------------------------------------------------------------
            // x direction of centerline's local coordinate system
            //------------------------------------------------------------------------------------------------------
            dbuffer.resize(3);
            file.read(reinterpret_cast<char*>(dbuffer.data()), dbuffer.size() * sizeof(double));
            if (cnt < NUM_DEMO)
            { _res << "\t\t\t- x direction of local coordinate system: [" << dbuffer[0] << ", " << dbuffer[1] << ", " << dbuffer[2] << "]" << std::endl; }

            //------------------------------------------------------------------------------------------------------
            // y direction of centerline's local coordinate system
            //------------------------------------------------------------------------------------------------------
            dbuffer.resize(3);
            file.read(reinterpret_cast<char*>(dbuffer.data()), dbuffer.size() * sizeof(double));
            if (cnt < NUM_DEMO)
            { _res << "\t\t\t- y direction of local coordinate system: [" << dbuffer[0] << ", " << dbuffer[1] << ", " << dbuffer[2] << "]" << std::endl; }

            ++cnt;
        } // for points
    } // for fjid: num numFlowjets

    file.close();

    return true;
}

bool ImporterScientific::read_ivsd(std::string_view filepath)
{
    if (!std::filesystem::exists(filepath.data()))
    {
        _res << "\t- no ivsd (path \"" << filepath.data() << "\")" << std::endl;
        return false;
    }

    _res << "\t- reading ivsd (path \"" << filepath.data() << "\")" << std::endl;

    std::ifstream file(filepath.data(), std::ios_base::in | std::ios_base::binary);

    if (!file.good())
    {
        _res << "\t\tFAILED! Could not open file!" << std::endl;
        return false;
    }

    _read_nd_scalar_image_in_sparse_matrix_style(file);

    file.close();

    return true;
}

bool ImporterScientific::read_magnitude_tmip(std::string_view filepath)
{
    if (!std::filesystem::exists(filepath.data()))
    {
        _res << "\t- no mag tmip (path \"" << filepath.data() << "\")" << std::endl;
        return false;
    }

    _res << "\t- reading mag tmip (path \"" << filepath.data() << "\")" << std::endl;

    std::ifstream file(filepath.data(), std::ios_base::in | std::ios_base::binary);

    if (!file.good())
    {
        _res << "\t\tFAILED! Could not open file!" << std::endl;
        return false;
    }

    _read_nd_scalar_image_in_sparse_matrix_style(file);

    file.close();

    return true;
}

bool ImporterScientific::read_anatomical_images()
{
    _res << "\t- reading anatomical images (path \"" << _dir << "\")" << std::endl;

    //------------------------------------------------------------------------------------------------------
    // find 3d anatomical images
    //------------------------------------------------------------------------------------------------------
    std::vector<std::string> anatomicalImage3DNames;

    for (auto& dirIt : std::filesystem::directory_iterator(_dir))
    {
        if (!dirIt.is_regular_file())
        { continue; }

        if (bk::string_utils::contains(dirIt.path().filename().string(), "3d_anatomical_image", false))
        { anatomicalImage3DNames.emplace_back(dirIt.path().filename().string()); }
    }

    _res << "\t\t- found " << anatomicalImage3DNames.size() << " 3D anatomical images: ";
    for (unsigned int i = 0; i < anatomicalImage3DNames.size(); ++i)
    {
        _res << anatomicalImage3DNames[i];
        if (static_cast<int>(i) < static_cast<int>(anatomicalImage3DNames.size()) - 1)
        { _res << ", "; }
    }
    _res << std::endl;

    //------------------------------------------------------------------------------------------------------
    // read 3d anatomical images
    //------------------------------------------------------------------------------------------------------
    for (const std::string& imgName: anatomicalImage3DNames)
    {
        const std::string filepath = _dir + "/" + imgName;
        std::ifstream file(filepath, std::ios_base::in | std::ios_base::binary);

        _read_nd_scalar_image_in_sparse_matrix_style(file);

        file.close();
    }

    //------------------------------------------------------------------------------------------------------
    // find 3d+t anatomical images
    //------------------------------------------------------------------------------------------------------
    std::vector<std::string> anatomicalImage3DTNames;

    for (auto& dirIt : std::filesystem::directory_iterator(_dir))
    {
        if (!dirIt.is_regular_file())
        { continue; }

        if (bk::string_utils::contains(dirIt.path().filename().string(), "3dt_anatomical_image", false))
        { anatomicalImage3DTNames.emplace_back(dirIt.path().filename().string()); }
    }

    _res << "\t\t- found " << anatomicalImage3DTNames.size() << " 3D+T anatomical images: ";
    for (unsigned int i = 0; i < anatomicalImage3DTNames.size(); ++i)
    {
        _res << anatomicalImage3DTNames[i];
        if (static_cast<int>(i) < static_cast<int>(anatomicalImage3DTNames.size()) - 1)
        { _res << ", "; }
    }
    _res << std::endl;

    //------------------------------------------------------------------------------------------------------
    // read 3d+t anatomical images
    //------------------------------------------------------------------------------------------------------
    for (const std::string& imgName: anatomicalImage3DTNames)
    {
        const std::string filepath = _dir + "/" + imgName;
        std::ifstream file(filepath, std::ios_base::in | std::ios_base::binary);

        _read_nd_scalar_image_in_sparse_matrix_style(file);

        file.close();
    }

    return true;
}

bool ImporterScientific::read_flow2dt_images()
{
    _res << "\t- searching 2D+T flow images in \"" << _dir << "\"" << std::endl;

    //------------------------------------------------------------------------------------------------------
    // find 2d+t flow images
    //------------------------------------------------------------------------------------------------------
    std::vector<std::string> flowImage2DTNames;

    for (auto& dirIt : std::filesystem::directory_iterator(_dir))
    {
        if (!dirIt.is_regular_file())
        { continue; }

        if (bk::string_utils::contains(dirIt.path().filename().string(), "flowfield_2dt", false))
        { flowImage2DTNames.emplace_back(dirIt.path().filename().string()); }
    }

    std::sort(flowImage2DTNames.begin(), flowImage2DTNames.end());

    _res << "\t\t- found " << flowImage2DTNames.size() << " 2D+T flow images: ";
    for (unsigned int i = 0; i < flowImage2DTNames.size(); ++i)
    {
        _res << flowImage2DTNames[i];
        if (static_cast<int>(i) < static_cast<int>(flowImage2DTNames.size()) - 1)
        { _res << ", "; }
    }
    _res << std::endl;

    //------------------------------------------------------------------------------------------------------
    // read 2d+t flow images
    //------------------------------------------------------------------------------------------------------
    for (const std::string& imgName: flowImage2DTNames)
    {
        /*
         *                     [3] x [uint32] : size x y t
         *                     [3] x [double] : scale x y t
         *                    [16] x [double] : world matrix (4x4 from dicom)
         *                    [16] x [double] : inverse world matrix (4x4 from dicom)
         *                    [25] x [double] : world matrix (5x5 including time in 4th row/col)
         *                    [25] x [double] : inverse world matrix (5x5 including time in 4th row/col)
         * [sizeX * sizeY * sizeT] x [double] : flow velocities
         */

        const std::string filepath = _dir + "/" + imgName;
        std::ifstream file(filepath, std::ios_base::in | std::ios_base::binary);
        _res << "\t\t- image " << imgName << ":" << std::endl;

        //------------------------------------------------------------------------------------------------------
        // grid size
        //------------------------------------------------------------------------------------------------------
        std::vector<std::uint32_t> gridsize(3);
        file.read(reinterpret_cast<char*>(gridsize.data()), gridsize.size() * sizeof(std::uint32_t));
        _res << "\t\t\t- grid size (xyt): " << gridsize[0] << " x " << gridsize[1] << " x " << gridsize[2] << std::endl;

        //------------------------------------------------------------------------------------------------------
        // voxel scale
        //------------------------------------------------------------------------------------------------------
        std::vector<double> voxelscale(3);
        file.read(reinterpret_cast<char*>(voxelscale.data()), voxelscale.size() * sizeof(double));
        _res << "\t\t\t- voxel scale (xyt): " << voxelscale[0] << " x " << voxelscale[1] << " [mm] / " << voxelscale[2] << " [ms]" << std::endl;

        //------------------------------------------------------------------------------------------------------
        // world matrix
        //------------------------------------------------------------------------------------------------------
        std::vector<double> dbuffer(16);
        file.read(reinterpret_cast<char*>(dbuffer.data()), dbuffer.size() * sizeof(double));
        unsigned int cnt = 0;
        _res << "\t\t\t- world matrix:" << std::endl;
        for (unsigned int rowid = 0; rowid < 4; ++rowid)
        {
            _res << "\t\t\t\t";
            for (unsigned int colid = 0; colid < 4; ++colid)
            { _res << dbuffer[cnt++] << " "; }

            _res << std::endl;
        }

        //------------------------------------------------------------------------------------------------------
        // inverse world matrix
        //------------------------------------------------------------------------------------------------------
        //dbuffer.resize(16);
        file.read(reinterpret_cast<char*>(dbuffer.data()), dbuffer.size() * sizeof(double));
        cnt = 0;
        _res << "\t\t\t- inverse world matrix:" << std::endl;
        for (unsigned int rowid = 0; rowid < 4; ++rowid)
        {
            _res << "\t\t\t\t";
            for (unsigned int colid = 0; colid < 4; ++colid)
            { _res << dbuffer[cnt++] << " "; }

            _res << std::endl;
        }

        //------------------------------------------------------------------------------------------------------
        // world matrix with time
        //------------------------------------------------------------------------------------------------------
        dbuffer.resize(25);
        file.read(reinterpret_cast<char*>(dbuffer.data()), dbuffer.size() * sizeof(double));
        cnt = 0;
        _res << "\t\t\t- world matrix with time:" << std::endl;
        for (unsigned int rowid = 0; rowid < 5; ++rowid)
        {
            _res << "\t\t\t\t";
            for (unsigned int colid = 0; colid < 5; ++colid)
            { _res << dbuffer[cnt++] << " "; }

            _res << std::endl;
        }

        //------------------------------------------------------------------------------------------------------
        // inverse world matrix with time
        //------------------------------------------------------------------------------------------------------
        //dbuffer.resize(25);
        file.read(reinterpret_cast<char*>(dbuffer.data()), dbuffer.size() * sizeof(double));
        cnt = 0;
        _res << "\t\t\t- inverse world matrix with time:" << std::endl;
        for (unsigned int rowid = 0; rowid < 5; ++rowid)
        {
            _res << "\t\t\t\t";
            for (unsigned int colid = 0; colid < 5; ++colid)
            { _res << dbuffer[cnt++] << " "; }

            _res << std::endl;
        }

        //------------------------------------------------------------------------------------------------------
        // flow velocities
        //------------------------------------------------------------------------------------------------------
        dbuffer.resize(gridsize[0] * gridsize[1] * gridsize[2]);
        file.read(reinterpret_cast<char*>(dbuffer.data()), dbuffer.size() * sizeof(double));

        cnt = 0;
        for (unsigned int x = 0; x < gridsize[0] && cnt < NUM_DEMO; ++x)
        {
            for (unsigned int y = 0; y < gridsize[1] && cnt < NUM_DEMO; ++y)
            {
                for (unsigned int t = 0; t < gridsize[2] && cnt < NUM_DEMO; ++t, ++cnt)
                {
                    const unsigned int off = x * gridsize[1] * gridsize[2] + y * gridsize[2] + t;
                    _res << "\t\t\t- velocity at grid pos [" << x << ", " << y << ", " << t << "] = " << dbuffer[off] << std::endl;
                }
            }
        }

        file.close();
    }

    return true;
}

bool ImporterScientific::read_flow_statistics(std::string_view filepath)
{
    if (!std::filesystem::exists(filepath.data()))
    {
        _res << "\t- no flow statistics (path \"" << filepath.data() << "\")" << std::endl;
        return false;
    }

    _res << "\t- reading flow statistics (path \"" << filepath.data() << "\")" << std::endl;

    std::ifstream file(filepath.data(), std::ios_base::in | std::ios_base::binary);

    if (!file.good())
    {
        _res << "\t\tFAILED! Could not open file!" << std::endl;
        return false;
    }

    /*
     * numTimes
     */
    std::uint32_t numTimes = 0;
    file.read(reinterpret_cast<char*>(&numTimes), sizeof(std::uint32_t));

    /*
     * helper functions
     */
    const auto read_and_print_value = [&](std::string_view name)
    {
        double dtemp = 0;
        file.read(reinterpret_cast<char*>(&dtemp), sizeof(double));
        _res << "\t\t- " << name << ": " << dtemp << std::endl;
    };

    const auto read_and_print_vector = [&](std::string_view name)
    {
        std::vector<double> dbuffer(numTimes);
        file.read(reinterpret_cast<char*>(dbuffer.data()), numTimes * sizeof(double));
        _res << "\t\t- " << name << ": ";

        for (unsigned int i = 0; i < NUM_DEMO; ++i)
        { _res << dbuffer[i] << ", "; }
        _res << "..." << std::endl;
    };

    read_and_print_value("vortex pressure threshold");
    read_and_print_value("volume total in ml");
    read_and_print_value("section volume in ml");
    read_and_print_value("section volume in percent");

    //------------------------------------------------------------------------------------------------------
    // diameter
    //------------------------------------------------------------------------------------------------------
    read_and_print_value("min diameter in mm");
    read_and_print_value("max diameter in mm");
    read_and_print_value("mean diameter in mm");
    read_and_print_value("median diameter in mm");

    //------------------------------------------------------------------------------------------------------
    // cross-sectional area
    //------------------------------------------------------------------------------------------------------
    read_and_print_value("min cross sectional area in mm2");
    read_and_print_value("max cross sectional area in mm2");
    read_and_print_value("mean cross sectional area in mm2");
    read_and_print_value("median cross sectional area in mm2");

    //------------------------------------------------------------------------------------------------------
    // vortex volume
    //------------------------------------------------------------------------------------------------------
    read_and_print_vector("vortex volume in ml per time");
    read_and_print_vector("vortex volume in percent per time");

    read_and_print_value("max vortex volume in ml");
    read_and_print_value("max vortex volume in percent");
    read_and_print_value("max vortex volume time in ms");
    read_and_print_value("mean vortex volume in ml");
    read_and_print_value("mean vortex volume in percent");
    read_and_print_value("median vortex volume in ml");
    read_and_print_value("median vortex volume in percent");
    read_and_print_value("systolic max vortex volume in ml");
    read_and_print_value("systolic max vortex volume in percent");
    read_and_print_value("systolic max vortex volume time in ms");
    read_and_print_value("systolic mean vortex volume in ml");
    read_and_print_value("systolic mean vortex volume in percent");
    read_and_print_value("systolic median vortex volume in ml");
    read_and_print_value("systolic median vortex volume in percent");
    read_and_print_value("diastolic max vortex volume in ml");
    read_and_print_value("diastolic max vortex volume in percent");
    read_and_print_value("diastolic max vortex volume time in ms");
    read_and_print_value("diastolic mean vortex volume in ml");
    read_and_print_value("diastolic mean vortex volume in percent");
    read_and_print_value("diastolic median vortex volume in ml");
    read_and_print_value("diastolic median vortex volume in percent");

    //------------------------------------------------------------------------------------------------------
    // vortex coverage
    //------------------------------------------------------------------------------------------------------
    read_and_print_value("vortex coverage in ml");
    read_and_print_value("vortex coverage in percent");
    read_and_print_value("systolic vortex coverage in ml");
    read_and_print_value("systolic vortex coverage in percent");
    read_and_print_value("diastolic vortex coverage in ml");
    read_and_print_value("diastolic vortex coverage in percent");

    //------------------------------------------------------------------------------------------------------
    // velocity
    //------------------------------------------------------------------------------------------------------
    read_and_print_vector("max velocity per time");
    read_and_print_vector("max axial velocity per time");
    read_and_print_vector("max circumferential velocity per time");
    read_and_print_vector("mean velocity per time");
    read_and_print_vector("mean axial velocity per time");
    read_and_print_vector("mean circumferential velocity per time");
    read_and_print_vector("median velocity per time");
    read_and_print_vector("median axial velocity per time");
    read_and_print_vector("median circumferential velocity per time");

    read_and_print_value("max mean velocity");
    read_and_print_value("max mean velocity time in ms");
    read_and_print_value("max mean axial velocity");
    read_and_print_value("max mean axial velocity time in ms");
    read_and_print_value("max mean circumferential velocity");
    read_and_print_value("max mean circumferential velocity time in ms");
    read_and_print_value("mean mean velocity");
    read_and_print_value("mean mean axial velocity");
    read_and_print_value("mean mean circumferential velocity");
    read_and_print_value("median mean velocity");
    read_and_print_value("median mean axial velocity");
    read_and_print_value("median mean circumferential velocity");
    read_and_print_value("max overall velocity");
    read_and_print_value("max overall velocity time in ms");
    read_and_print_value("max overall velocity q99");
    read_and_print_value("max overall velocity q99 time in ms");
    read_and_print_value("max overall axial velocity");
    read_and_print_value("max overall axial velocity time in ms");
    read_and_print_value("max overall axial velocity q99");
    read_and_print_value("max overall axial velocity q99 time in ms");
    read_and_print_value("max overall circumferential velocity");
    read_and_print_value("max overall circumferential velocity time in ms");
    read_and_print_value("max overall circumferential velocity q99");
    read_and_print_value("max overall circumferential velocity q99 time in ms");
    read_and_print_value("systolic max mean velocity");
    read_and_print_value("systolic max mean velocity time in ms");
    read_and_print_value("systolic max mean axial velocity");
    read_and_print_value("systolic max mean axial velocity time in ms");
    read_and_print_value("systolic max mean circumferential velocity");
    read_and_print_value("systolic max mean circumferential velocity time in ms");
    read_and_print_value("systolic mean mean velocity");
    read_and_print_value("systolic mean mean axial velocity");
    read_and_print_value("systolic mean mean circumferential velocity");
    read_and_print_value("systolic median mean velocity");
    read_and_print_value("systolic median mean axial velocity");
    read_and_print_value("systolic median mean circumferential velocity");
    read_and_print_value("systolic max overall velocity");
    read_and_print_value("systolic max overall velocity time in ms");
    read_and_print_value("systolic max overall velocity q99");
    read_and_print_value("systolic max overall velocity q99 time in ms");
    read_and_print_value("systolic max overall axial velocity");
    read_and_print_value("systolic max overall axial velocity time in ms");
    read_and_print_value("systolic max overall axial velocity q99");
    read_and_print_value("systolic max overall axial velocity q99 time in ms");
    read_and_print_value("systolic max overall circumferential velocity");
    read_and_print_value("systolic max overall circumferential velocity time in ms");
    read_and_print_value("systolic max overall circumferential velocity q99");
    read_and_print_value("systolic max overall circumferential velocity q99 time in ms");
    read_and_print_value("diastolic max mean velocity");
    read_and_print_value("diastolic max mean velocity time in ms");
    read_and_print_value("diastolic max mean axial velocity");
    read_and_print_value("diastolic max mean axial velocity time in ms");
    read_and_print_value("diastolic max mean circumferential velocity");
    read_and_print_value("diastolic max mean circumferential velocity time in ms");
    read_and_print_value("diastolic mean mean velocity");
    read_and_print_value("diastolic mean mean axial velocity");
    read_and_print_value("diastolic mean mean circumferential velocity");
    read_and_print_value("diastolic median mean velocity");
    read_and_print_value("diastolic median mean axial velocity");
    read_and_print_value("diastolic median mean circumferential velocity");
    read_and_print_value("diastolic max overall velocity");
    read_and_print_value("diastolic max overall velocity time in ms");
    read_and_print_value("diastolic max overall velocity q99");
    read_and_print_value("diastolic max overall velocity q99 time in ms");
    read_and_print_value("diastolic max overall axial velocity");
    read_and_print_value("diastolic max overall axial velocity time in ms");
    read_and_print_value("diastolic max overall axial velocity q99");
    read_and_print_value("diastolic max overall axial velocity q99 time in ms");
    read_and_print_value("diastolic max overall circumferential velocity");
    read_and_print_value("diastolic max overall circumferential velocity time in ms");
    read_and_print_value("diastolic max overall circumferential velocity q99");
    read_and_print_value("diastolic max overall circumferential velocity q99 time in ms");

    //------------------------------------------------------------------------------------------------------
    // rotation
    //------------------------------------------------------------------------------------------------------
    read_and_print_vector("left rotation volume in ml per time");
    read_and_print_vector("left rotation volume in percent per time");

    read_and_print_value("max left rotation volume in ml");
    read_and_print_value("max left rotation volume in percent");
    read_and_print_value("max left rotation volume time in ms");
    read_and_print_value("mean left rotation volume in ml");
    read_and_print_value("mean left rotation volume in percent");
    read_and_print_value("median left rotation volume in ml");
    read_and_print_value("median left rotation volume in percent");
    read_and_print_value("systolic max left rotation volume in ml");
    read_and_print_value("systolic max left rotation volume in percent");
    read_and_print_value("systolic max left rotation volume time in ms");
    read_and_print_value("systolic mean left rotation volume in ml");
    read_and_print_value("systolic mean left rotation volume in percent");
    read_and_print_value("systolic median left rotation volume in ml");
    read_and_print_value("systolic median left rotation volume in percent");
    read_and_print_value("diastolic max left rotation volume in ml");
    read_and_print_value("diastolic max left rotation volume in percent");
    read_and_print_value("diastolic max left rotation volume time in ms");
    read_and_print_value("diastolic mean left rotation volume in ml");
    read_and_print_value("diastolic mean left rotation volume in percent");
    read_and_print_value("diastolic median left rotation volume in ml");
    read_and_print_value("diastolic median left rotation volume in percent");

    read_and_print_vector("right rotation volume in ml per time");
    read_and_print_vector("right rotation volume in percent per time");

    read_and_print_value("max right rotation volume in ml");
    read_and_print_value("max right rotation volume in percent");
    read_and_print_value("max right rotation volume time in ms");
    read_and_print_value("mean right rotation volume in ml");
    read_and_print_value("mean right rotation volume in percent");
    read_and_print_value("median right rotation volume in ml");
    read_and_print_value("median right rotation volume in percent");
    read_and_print_value("systolic max right rotation volume in ml");
    read_and_print_value("systolic max right rotation volume in percent");
    read_and_print_value("systolic max right rotation volume time in ms");
    read_and_print_value("systolic mean right rotation volume in ml");
    read_and_print_value("systolic mean right rotation volume in percent");
    read_and_print_value("systolic median right rotation volume in ml");
    read_and_print_value("systolic median right rotation volume in percent");
    read_and_print_value("diastolic max right rotation volume in ml");
    read_and_print_value("diastolic max right rotation volume in percent");
    read_and_print_value("diastolic max right rotation volume time in ms");
    read_and_print_value("diastolic mean right rotation volume in ml");
    read_and_print_value("diastolic mean right rotation volume in percent");
    read_and_print_value("diastolic median right rotation volume in ml");
    read_and_print_value("diastolic median right rotation volume in percent");

    //------------------------------------------------------------------------------------------------------
    // pressure
    //------------------------------------------------------------------------------------------------------
    /*
     * mean pressure per time
     */
    read_and_print_vector("mean pressure per time");

    /*
     * min max mean median
     */
    read_and_print_value("min mean pressure");
    read_and_print_value("min mean pressure time in ms");
    read_and_print_value("max mean pressure");
    read_and_print_value("max mean pressure time in ms");
    read_and_print_value("mean mean pressure");
    read_and_print_value("median mean pressure");

    /*
     * min max mean median : systolic
     */
    read_and_print_value("systolic min mean pressure");
    read_and_print_value("systolic min mean pressure time in ms");
    read_and_print_value("systolic max mean pressure");
    read_and_print_value("systolic max mean pressure time in ms");
    read_and_print_value("systolic mean mean pressure");
    read_and_print_value("systolic median mean pressure");

    /*
     * min max mean median : diastolic
     */
    read_and_print_value("diastolic min mean pressure");
    read_and_print_value("diastolic min mean pressure time in ms");
    read_and_print_value("diastolic max mean pressure");
    read_and_print_value("diastolic max mean pressure time in ms");
    read_and_print_value("diastolic mean mean pressure");
    read_and_print_value("diastolic median mean pressure");

    /*
     * mean pressure per time in vortex region
     */
    read_and_print_vector("mean pressure in vortex region per time");

    /*
     * min max mean median in vortex region
     */
    read_and_print_value("min mean pressure in vortex region");
    read_and_print_value("min mean pressure in vortex region time in ms");
    read_and_print_value("max mean pressure in vortex region");
    read_and_print_value("max mean pressure in vortex region time in ms");
    read_and_print_value("mean mean pressure in vortex region");
    read_and_print_value("median mean pressure in vortex region");

    /*
     * min max mean median in vortex region : systolic
     */
    read_and_print_value("systolic min mean pressure in vortex region");
    read_and_print_value("systolic min mean pressure in vortex region time in ms");
    read_and_print_value("systolic max mean pressure in vortex region");
    read_and_print_value("systolic max mean pressure in vortex region time in ms");
    read_and_print_value("systolic mean mean pressure in vortex region");
    read_and_print_value("systolic median mean pressure in vortex region");

    /*
     * min max mean median in vortex region : diastolic
     */
    read_and_print_value("diastolic min mean pressure in vortex region");
    read_and_print_value("diastolic min mean pressure in vortex region time in ms");
    read_and_print_value("diastolic max mean pressure in vortex region");
    read_and_print_value("diastolic max mean pressure in vortex region time in ms");
    read_and_print_value("diastolic mean mean pressure in vortex region");
    read_and_print_value("diastolic median mean pressure in vortex region");

    //------------------------------------------------------------------------------------------------------
    // flow displacement
    //------------------------------------------------------------------------------------------------------
    /*
     * min max mean median displacement (velocity-weighted)
     */
    read_and_print_value("max flow jet displacement velocity weighted");
    read_and_print_value("min flow jet displacement velocity weighted");
    read_and_print_value("mean flow jet displacement velocity weighted");
    read_and_print_value("median flow jet displacement velocity weighted");

    /*
     * min max mean median angle to centerline (velocity-weighted)
     */
    read_and_print_value("max flow jet angle velocity weighted");
    read_and_print_value("min flow jet angle velocity weighted");
    read_and_print_value("mean flow jet angle velocity weighted");
    read_and_print_value("median flow jet angle velocity weighted");

    /*
     * min max mean median high-velocity area percent of cross-section (velocity-weighted)
     */
    read_and_print_value("max flow jet high velocity area percent velocity weighted");
    read_and_print_value("min flow jet high velocity area percent velocity weighted");
    read_and_print_value("mean flow jet high velocity area percent velocity weighted");
    read_and_print_value("median flow jet high velocity area percent velocity weighted");

    file.close();

    return true;
}

bool ImporterScientific::read_segmentation(std::string_view filepath)
{
    if (!std::filesystem::exists(filepath.data()))
    {
        _res << "\t- no segmentation (path \"" << filepath.data() << "\")" << std::endl;
        return false;
    }

    _res << "\t- reading segmentation (path \"" << filepath.data() << "\")" << std::endl;

    std::ifstream file(filepath.data(), std::ios_base::in | std::ios_base::binary);

    if (!file.good())
    {
        _res << "\t\tFAILED! Could not open file!" << std::endl;
        return false;
    }

    _read_nd_scalar_image_in_sparse_matrix_style(file);

    file.close();

    return true;
}

bool ImporterScientific::read_segmentation_info(std::string_view filepath)
{
    if (!std::filesystem::exists(filepath.data()))
    {
        _res << "\t- no segmentation info (path \"" << filepath.data() << "\")" << std::endl;
        return false;
    }

    _res << "\t- reading segmentation info (path \"" << filepath.data() << "\")" << std::endl;

    std::ifstream file(filepath.data(), std::ios_base::in);

    if (!file.good())
    {
        _res << "\t\tFAILED! Could not open file!" << std::endl;
        return false;
    }

    std::string line;
    while (!file.eof())
    {
        std::getline(file, line);

        if (line.empty())
        { continue; }

        _res << "\t\t-> \"" << line << "\"" << std::endl;
    }

    // result is one of the following:
    //
    // "The segmentation was performed on the magnitude images' TMIP."
    // "The segmentation was performed on the LPC."
    // "The segmentation was performed on 3D anatomical image <id>."
    // "The segmentation was performed on 3D+T anatomical image <id>."
    // "The segmentation was performed on the signal intensity image's TMIP."
    // "The segmentation was performed on the IVSD."

    file.close();

    return true;
}

bool ImporterScientific::read_segmentation_graphcut_inside_outside_ids(std::string_view filepath)
{
    if (!std::filesystem::exists(filepath.data()))
    {
        _res << "\t- no segmentation graph cut inside/outside ids (path \"" << filepath.data() << "\")" << std::endl;
        return false;
    }

    _res << "\t- reading segmentation graph cut inside/outside ids (path \"" << filepath.data() << "\")" << std::endl;

    std::ifstream file(filepath.data(), std::ios_base::in | std::ios_base::binary);

    if (!file.good())
    {
        _res << "\t\tFAILED! Could not open file!" << std::endl;
        return false;
    }

    //------------------------------------------------------------------------------------------------------
    // numInsideIds
    //------------------------------------------------------------------------------------------------------
    std::uint32_t numInsideIds = 0;
    file.read(reinterpret_cast<char*>(&numInsideIds), sizeof(std::uint32_t));
    _res << "\t\t- num. inside ids: " << numInsideIds << std::endl;

    //------------------------------------------------------------------------------------------------------
    // numOutsideIds
    //------------------------------------------------------------------------------------------------------
    std::uint32_t numOutsideIds = 0;
    file.read(reinterpret_cast<char*>(&numOutsideIds), sizeof(std::uint32_t));
    _res << "\t\t- num. outside ids: " << numOutsideIds << std::endl;

    //------------------------------------------------------------------------------------------------------
    // list of inside id triplets (xyz)
    //------------------------------------------------------------------------------------------------------
    std::vector<std::uint32_t> ui32buffer(numInsideIds * 3);
    file.read(reinterpret_cast<char*>(ui32buffer.data()), ui32buffer.size() * sizeof(std::uint32_t));

    for (unsigned int i = 0; i < std::min(NUM_DEMO, numInsideIds); ++i)
    {
        const unsigned int off = i * 3;
        _res << "\t\t- inside grid pos " << i << ": [" << ui32buffer[off] << ", " << ui32buffer[off + 1] << ", " << ui32buffer[off + 2] << "]" << std::endl;
    }
    _res << "\t\t- ..." << std::endl;

    //------------------------------------------------------------------------------------------------------
    // list of outside id triplets (xyz)
    //------------------------------------------------------------------------------------------------------
    ui32buffer.clear();
    ui32buffer.resize(numOutsideIds * 3, 0);
    file.read(reinterpret_cast<char*>(ui32buffer.data()), ui32buffer.size() * sizeof(std::uint32_t));

    for (unsigned int i = 0; i < std::min(NUM_DEMO, numOutsideIds); ++i)
    {
        const unsigned int off = i * 3;
        _res << "\t\t- outside grid pos " << i << ": [" << ui32buffer[off] << ", " << ui32buffer[off+1] << ", " << ui32buffer[off+2] << "]" << std::endl;
    }
    if (numOutsideIds != 0)
    { _res << "\t\t- ..." << std::endl; }
    else
    { _res << "\t\t- no outside ids specified" << std::endl; }

    file.close();

    return true;
}

bool ImporterScientific::read_segmentation_in_flowfield_size(std::string_view filepath)
{
    if (!std::filesystem::exists(filepath.data()))
    {
        _res << "\t- no segmentation in flow field size (path \"" << filepath.data() << "\")" << std::endl;
        return false;
    }

    _res << "\t- reading segmentation in flow field size (path \"" << filepath.data() << "\")" << std::endl;

    std::ifstream file(filepath.data(), std::ios_base::in | std::ios_base::binary);

    if (!file.good())
    {
        _res << "\t\tFAILED! Could not open file!" << std::endl;
        return false;
    }

    _read_nd_scalar_image_in_sparse_matrix_style(file);

    file.close();

    return true;
}

bool ImporterScientific::read_vessel_section_segmentation_in_flowfield_size(std::string_view filepath)
{
    if (!std::filesystem::exists(filepath.data()))
    {
        _res << "\t- no vessel section segmentation in flow field size (path \"" << filepath.data() << "\")" << std::endl;
        return false;
    }

    _res << "\t- reading vessel section segmentation in flow field size (path \"" << filepath.data() << "\")" << std::endl;

    std::ifstream file(filepath.data(), std::ios_base::in | std::ios_base::binary);

    if (!file.good())
    {
        _res << "\t\tFAILED! Could not open file!" << std::endl;
        return false;
    }

    //------------------------------------------------------------------------------------------------------
    // numSections
    //------------------------------------------------------------------------------------------------------
    std::uint32_t numSections = 0;
    file.read(reinterpret_cast<char*>(&numSections), sizeof(std::uint32_t));
    _res << "\t\tnum. sections: " << numSections << std::endl;

    //------------------------------------------------------------------------------------------------------
    // section segmentations sparse matrix style
    //------------------------------------------------------------------------------------------------------
    for (unsigned int sectionid = 0; sectionid < numSections; ++sectionid)
    {
        _res << "\t\t- section " << sectionid << " of " << numSections << ":" << std::endl;
        _read_nd_scalar_image_in_sparse_matrix_style(file);
    }

    file.close();

    return true;
}

bool ImporterScientific::read_vessel_section_segmentation_semantics(std::string_view filepath)
{
    if (!std::filesystem::exists(filepath.data()))
    {
        _res << "\t- no vessel section segmentation semantics (path \"" << filepath.data() << "\")" << std::endl;
        return false;
    }

    _res << "\t- reading vessel section segmentation semantics (path \"" << filepath.data() << "\")" << std::endl;

    std::ifstream file(filepath.data(), std::ios_base::in);

    if (!file.good())
    {
        _res << "\t\tFAILED! Could not open file!" << std::endl;
        return false;
    }

    std::string line;
    while (!file.eof())
    {
        std::getline(file, line);

        if (line.empty())
        { continue; }

        _res << "\t\t" << line << std::endl;
    }

    file.close();

    return true;
}

bool ImporterScientific::read_centerline_start_end_ids_on_mesh(std::string_view filepath)
{
    /*
     * [1] x [uint32] : seedId
     * [1] x [uint32] : numTargetIds
     * [numTargetIds] x [uint32] : targetIds
     */

    if (!std::filesystem::exists(filepath.data()))
    {
        _res << "\t- no centerline start/end ids (path \"" << filepath.data() << "\")" << std::endl;
        return false;
    }

    _res << "\t- reading centerline start/end ids (path \"" << filepath.data() << "\")" << std::endl;

    std::ifstream file(filepath.data(), std::ios_base::in | std::ios_base::binary);

    if (!file.good())
    {
        _res << "\t\tFAILED! Could not open file!" << std::endl;
        return false;
    }

    //------------------------------------------------------------------------------------------------------
    // seedId
    //------------------------------------------------------------------------------------------------------
    std::uint32_t seedId = 0;
    file.read(reinterpret_cast<char*>(&seedId), sizeof(std::uint32_t));
    _res << "\t\t- seed id: " << seedId << std::endl;

    //------------------------------------------------------------------------------------------------------
    // numTargetIds
    //------------------------------------------------------------------------------------------------------
    std::uint32_t numTargetIds = 0;
    file.read(reinterpret_cast<char*>(&numTargetIds), sizeof(std::uint32_t));

    //------------------------------------------------------------------------------------------------------
    // targetIds
    //------------------------------------------------------------------------------------------------------
    std::vector<std::uint32_t> targetIds(numTargetIds);
    file.read(reinterpret_cast<char*>(targetIds.data()), numTargetIds * sizeof(std::uint32_t));

    _res << "\t\t- target ids: ";
    for (std::uint32_t tid : targetIds)
    { _res << tid << " "; }
    _res << std::endl;

    return true;
}

bool ImporterScientific::read_static_tissue_mask(std::string_view filepath)
{
    if (!std::filesystem::exists(filepath.data()))
    {
        _res << "\t- no static tissue mask (path \"" << filepath.data() << "\")" << std::endl;
        return false;
    }

    _res << "\t- reading static tissue mask (path \"" << filepath.data() << "\")" << std::endl;

    std::ifstream file(filepath.data(), std::ios_base::in | std::ios_base::binary);

    if (!file.good())
    {
        _res << "\t\tFAILED! Could not open file!" << std::endl;
        return false;
    }

    _read_nd_scalar_image_in_sparse_matrix_style(file);

    return true;
}

bool ImporterScientific::read_static_tissue_ivsd_thresholds(std::string_view filepath)
{
    if (!std::filesystem::exists(filepath.data()))
    {
        _res << "\t- no static tissue ivsd thresholds (path \"" << filepath.data() << "\")" << std::endl;
        return false;
    }

    _res << "\t- reading static tissue ivsd thresholds (path \"" << filepath.data() << "\")" << std::endl;

    std::ifstream file(filepath.data(), std::ios_base::in | std::ios_base::binary);

    if (!file.good())
    {
        _res << "\t\tFAILED! Could not open file!" << std::endl;
        return false;
    }

    double lowerThreshold = 0;
    file.read(reinterpret_cast<char*>(&lowerThreshold), sizeof(double));

    double upperThreshold = 0;
    file.read(reinterpret_cast<char*>(&upperThreshold), sizeof(double));

    _res << "\t\t- lower threshold: " << lowerThreshold << std::endl;
    _res << "\t\t- upper threshold: " << upperThreshold << std::endl;

    return true;
}

bool ImporterScientific::read_dataset_filter_tags(std::string_view filepath)
{
    if (!std::filesystem::exists(filepath.data()))
    {
        _res << "\t- no filter tags (path \"" << filepath.data() << "\")" << std::endl;
        return false;
    }

    _res << "\t- reading filter tags (path \"" << filepath.data() << "\")" << std::endl;

    std::ifstream file(filepath.data(), std::ios_base::in);

    if (!file.good())
    {
        _res << "\t\tFAILED! Could not open file!" << std::endl;
        return false;
    }

    std::string csv;
    std::getline(file, csv);

    std::vector<std::string> tags = bk::string_utils::split(csv, ";");
    tags.erase(std::remove_if(tags.begin(), tags.end(), [](const std::string& s)
    { return s.empty(); }));

    _res << "\t\t- " << tags.size() << " filter tags: ";

    for (const std::string& s: tags)
    { _res << s << " "; }
    _res << std::endl;

    return true;
}

bool ImporterScientific::read_phase_wrapped_voxels(std::string_view filepath)
{
    if (!std::filesystem::exists(filepath.data()))
    {
        _res << "\t- no phase wraps (path \"" << filepath.data() << "\")" << std::endl;
        return false;
    }

    _res << "\t- reading phase wraps (path \"" << filepath.data() << "\")" << std::endl;

    std::ifstream file(filepath.data(), std::ios_base::in | std::ios_base::binary);

    if (!file.good())
    {
        _res << "\t\tFAILED! Could not open file!" << std::endl;
        return false;
    }

    for (unsigned int dimid = 0; dimid < 3; ++dimid)
    {
        //------------------------------------------------------------------------------------------------------
        // numWrappedVoxels
        //------------------------------------------------------------------------------------------------------
        std::uint32_t numWrappedVoxels = 0;
        file.read(reinterpret_cast<char*>(&numWrappedVoxels), sizeof(std::uint32_t));
        _res << "\t\tnum. wrapped voxels of 3D+T flow image " << dimid << ": " << numWrappedVoxels << std::endl;

        for (unsigned int i = 0; i < numWrappedVoxels; ++i)
        {
            //------------------------------------------------------------------------------------------------------
            // x y z t gridPos
            //------------------------------------------------------------------------------------------------------
            std::vector<std::uint32_t> gridpos(4);
            file.read(reinterpret_cast<char*>(gridpos.data()), 4 * sizeof(std::uint32_t));
            if (i < NUM_DEMO)
            { _res << "\t\t\t- " << i << ": grid pos [" << gridpos[0] << ", " << gridpos[1] << ", " << gridpos[2] << ", " << gridpos[3] << "]"; }

            /*
             * wrapFactor
             * - x was corrected via:   x += factor * 2 * venc
             */
            std::int8_t wrapFactor = 0;
            file.read(reinterpret_cast<char*>(&wrapFactor), sizeof(std::int8_t));
            if (i < NUM_DEMO)
            { _res << " is wrapped " << static_cast<int>(wrapFactor) << "x" << std::endl; }
        } // for i
        _res << "\t\t\t- ..." << std::endl;
    } // for dimid

    return true;
}

bool ImporterScientific::read_velocity_offset_correction_3dt(std::string_view filepath)
{
    if (!std::filesystem::exists(filepath.data()))
    {
        _res << "\t- no 3D+T flow images' velocity offset correction (path \"" << filepath.data() << "\")" << std::endl;
        return false;
    }

    _res << "\t- reading flow images' velocity offset correction (path \"" << filepath.data() << "\")" << std::endl;

    std::ifstream file(filepath.data(), std::ios_base::in | std::ios_base::binary);

    if (!file.good())
    {
        _res << "\t\tFAILED! Could not open file!" << std::endl;
        return false;
    }

    //------------------------------------------------------------------------------------------------------
    // end diastolic time point
    //------------------------------------------------------------------------------------------------------
    std::uint32_t end_diastolic_time_id = 0;
    file.read(reinterpret_cast< char*>(&end_diastolic_time_id), sizeof(std::uint32_t));
    _res << "\t\t- end diastolic time point id: " << end_diastolic_time_id << std::endl;

    //------------------------------------------------------------------------------------------------------
    // ivsd threshold
    //------------------------------------------------------------------------------------------------------
    double ivsd_static_tissue_threshold = 0;
    file.read(reinterpret_cast< char*>(&ivsd_static_tissue_threshold), sizeof(double));
    _res << "\t\t- ivsd static tissue threshold: " << ivsd_static_tissue_threshold << std::endl;

    std::vector<double> planeCoeffs(3);

    for (unsigned int v = 0; v < 3; ++v)
    {
        //------------------------------------------------------------------------------------------------------
        // number of slices
        //------------------------------------------------------------------------------------------------------
        std::uint32_t numSlices = 0;
        file.read(reinterpret_cast< char*>(&numSlices), sizeof(std::uint32_t));
        _res << "\t\t\t- num. slices in flow image " << v << ": " << numSlices << std::endl;

        std::vector<double> buf(numSlices * 3);
        file.read(reinterpret_cast<char*>(buf.data()), buf.size() * sizeof(double));

        for (unsigned int z = 0; z < NUM_DEMO/*numSlices*/; ++z)
        {
            const unsigned int off = 3 * z;
            _res << "\t\t\t- plane coeffs of slice " << z << ": " << planeCoeffs[off] << ", " << planeCoeffs[off + 1] << ", " << planeCoeffs[off + 2] << std::endl;
        } // for z : numSlices
        _res << "\t\t\t- ..." << std::endl;
    } // for v

    file.close();

    return true;
}

bool ImporterScientific::read_dicom_tags(std::string_view filepath)
{
    if (!std::filesystem::exists(filepath.data()))
    {
        _res << "\t- no dicom tags (path \"" << filepath.data() << "\")" << std::endl;
        return false;
    }

    _res << "\t- reading dicom tags (path \"" << filepath.data() << "\")" << std::endl;

    std::ifstream file(filepath.data(), std::ios_base::in | std::ios_base::binary);

    if (!file.good())
    {
        _res << "\t\tFAILED! Could not open file!" << std::endl;
        return false;
    }

    const auto read_string = [&]() -> std::string
    {
        std::uint16_t len = 0;
        file.read(reinterpret_cast<char*>(&len), sizeof(std::uint16_t));

        std::string s;
        s.resize(len);
        file.read(s.data(), s.length());

        return s;
    };

    std::uint16_t numImages = 0;
    file.read(reinterpret_cast<char*>(&numImages), sizeof(std::uint16_t));

    _res << "\t\t- DICOM tags of " << numImages << " images: " << std::endl;

    for (unsigned int i = 0; i < numImages; ++i)
    {
        std::uint16_t dcmImgId = 0;
        file.read(reinterpret_cast<char*>(&dcmImgId), sizeof(std::uint16_t));
        _res << "\t\t- DICOM image ID: " << dcmImgId << std::endl;

        std::uint16_t nDimensions = 0;
        file.read(reinterpret_cast<char*>(&nDimensions), sizeof(std::uint16_t));
        _res << "\t\t\t- nDimensions: " << nDimensions << std::endl;

        std::uint16_t Rows = 0;
        file.read(reinterpret_cast<char*>(&Rows), sizeof(std::uint16_t));
        _res << "\t\t\t- Rows: " << Rows << std::endl;

        std::uint16_t Columns = 0;
        file.read(reinterpret_cast<char*>(&Columns), sizeof(std::uint16_t));
        _res << "\t\t\t- Columns: " << Columns << std::endl;

        std::uint16_t Slices = 0;
        file.read(reinterpret_cast<char*>(&Slices), sizeof(std::uint16_t));
        _res << "\t\t\t- Slices: " << Slices << std::endl;

        std::uint16_t TemporalPositions = 0;
        file.read(reinterpret_cast<char*>(&TemporalPositions), sizeof(std::uint16_t));
        _res << "\t\t\t- TemporalPositions: " << TemporalPositions << std::endl;

        std::uint16_t NumberOfFrames = 0;
        file.read(reinterpret_cast<char*>(&NumberOfFrames), sizeof(std::uint32_t));
        _res << "\t\t\t- NumberOfFrames: " << NumberOfFrames << std::endl;

        double RowSpacing = 0;
        file.read(reinterpret_cast<char*>(&RowSpacing), sizeof(double));
        _res << "\t\t\t- RowSpacing: " << RowSpacing << std::endl;

        double ColSpacing = 0;
        file.read(reinterpret_cast<char*>(&ColSpacing), sizeof(double));
        _res << "\t\t\t- ColSpacing: " << ColSpacing << std::endl;

        double SliceSpacing = 0;
        file.read(reinterpret_cast<char*>(&SliceSpacing), sizeof(double));
        _res << "\t\t\t- SliceSpacing: " << SliceSpacing << std::endl;

        double TemporalResolution = 0;
        file.read(reinterpret_cast<char*>(&TemporalResolution), sizeof(double));
        _res << "\t\t\t- TemporalResolution: " << TemporalResolution << std::endl;

        const std::string PatientName = read_string();
        _res << "\t\t\t- PatientName: " << PatientName << std::endl;

        const std::string PatientID = read_string();
        _res << "\t\t\t- PatientID: " << PatientID << std::endl;

        const std::string PatientSex = read_string();
        _res << "\t\t\t- PatientSex: " << PatientSex << std::endl;

        std::uint8_t PatientAge = 0;
        file.read(reinterpret_cast<char*>(&PatientAge), sizeof(std::uint8_t));
        _res << "\t\t\t- PatientAge: " << static_cast<int>(PatientAge) << std::endl;

        double PatientWeight = 0;
        file.read(reinterpret_cast<char*>(&PatientWeight), sizeof(double));
        _res << "\t\t\t- PatientWeight: " << PatientWeight << std::endl;

        const std::string PatientBirthDate = read_string();
        _res << "\t\t\t- PatientBirthDate: " << PatientBirthDate << std::endl;

        const std::string SequenceName = read_string();
        _res << "\t\t\t- SequenceName: " << SequenceName << std::endl;

        const std::string SequenceName_Private = read_string();
        _res << "\t\t\t- SequenceName_Private: " << SequenceName_Private << std::endl;

        const std::string PatientPosition = read_string();
        _res << "\t\t\t- PatientPosition: " << PatientPosition << std::endl;

        const std::string StudyDescription = read_string();
        _res << "\t\t\t- StudyDescription: " << StudyDescription << std::endl;

        const std::string SeriesDescription = read_string();
        _res << "\t\t\t- SeriesDescription: " << SeriesDescription << std::endl;

        const std::string SeriesInstanceUID = read_string();
        _res << "\t\t\t- SeriesInstanceUID: " << SeriesInstanceUID << std::endl;

        const std::string StudyInstanceUID = read_string();
        _res << "\t\t\t- StudyInstanceUID: " << StudyInstanceUID << std::endl;

        const std::string ProtocolName = read_string();
        _res << "\t\t\t- ProtocolName: " << ProtocolName << std::endl;

        const std::string Modality = read_string();
        _res << "\t\t\t- Modality: " << Modality << std::endl;

        std::uint8_t SamplesPerPixel = 0;
        file.read(reinterpret_cast<char*>(&SamplesPerPixel), sizeof(std::uint8_t));
        _res << "\t\t\t- SamplesPerPixel: " << static_cast<int>(SamplesPerPixel) << std::endl;

        std::uint32_t LargestImagePixelValue = 0;
        file.read(reinterpret_cast<char*>(&LargestImagePixelValue), sizeof(std::uint32_t));
        _res << "\t\t\t- LargestImagePixelValue: " << LargestImagePixelValue << std::endl;

        std::uint8_t BitsAllocated = 0;
        file.read(reinterpret_cast<char*>(&BitsAllocated), sizeof(std::uint8_t));
        _res << "\t\t\t- BitsAllocated: " << static_cast<int>(BitsAllocated) << std::endl;

        std::uint8_t BitsStored = 0;
        file.read(reinterpret_cast<char*>(&BitsStored), sizeof(std::uint8_t));
        _res << "\t\t\t- BitsStored: " << static_cast<int>(BitsStored) << std::endl;

        std::uint8_t HighBit = 0;
        file.read(reinterpret_cast<char*>(&HighBit), sizeof(std::uint8_t));
        _res << "\t\t\t- HighBit: " << static_cast<int>(HighBit) << std::endl;

        const std::string AcquisitionDate = read_string();
        _res << "\t\t\t- AcquisitionDate: " << AcquisitionDate << std::endl;

        const std::string InstitutionName = read_string();
        _res << "\t\t\t- InstitutionName: " << InstitutionName << std::endl;

        std::vector<double> ImageOrientationPatientX(3);
        file.read(reinterpret_cast<char*>(ImageOrientationPatientX.data()), 3 * sizeof(double));
        _res << "\t\t\t- ImageOrientationPatientX: [" << ImageOrientationPatientX[0] << ", " << ImageOrientationPatientX[1] << ", " << ImageOrientationPatientX[2] << "]" << std::endl;

        std::vector<double> ImageOrientationPatientY(3);
        file.read(reinterpret_cast<char*>(ImageOrientationPatientY.data()), 3 * sizeof(double));
        _res << "\t\t\t- ImageOrientationPatientY: [" << ImageOrientationPatientY[0] << ", " << ImageOrientationPatientY[1] << ", " << ImageOrientationPatientY[2] << "]" << std::endl;

        std::vector<double> worldMatrix(16);
        file.read(reinterpret_cast<char*>(worldMatrix.data()), 16 * sizeof(double));
        unsigned int cnt = 0;
        _res << "\t\t\t- world matrix:" << std::endl;
        for (unsigned int rowid = 0; rowid < 4; ++rowid)
        {
            _res << "\t\t\t\t";
            for (unsigned int colid = 0; colid < 4; ++colid)
            { _res << worldMatrix[cnt++] << " "; }

            _res << std::endl;
        }
    } // for image ids

    file.close();

    return true;
}

bool ImporterScientific::read_cardiac_cycle_definition(std::string_view filepath)
{
    /*
     * [1] x [uint32] : numTimes
     * [1] x [uint32] : idSystoleBegin
     * [1] x [double] : msSystoleBegin
     * [1] x [uint32] : idSystoleEnd
     * [1] x [double] : msSystoleEnd
     * [1] x [uint32] : numVessels
     * for numVessels
     *      [numTimes] x [double] : axial velocity per time in vessel
     */

    if (!std::filesystem::exists(filepath.data()))
    {
        _res << "\t- no cardiac cycle definition (path \"" << filepath.data() << "\")" << std::endl;
        return false;
    }

    _res << "\t- reading cardiac cycle definition (path \"" << filepath.data() << "\")" << std::endl;

    std::ifstream file(filepath.data(), std::ios_base::in | std::ios_base::binary);

    if (!file.good())
    {
        _res << "\t\tFAILED! Could not open file!" << std::endl;
        return false;
    }

    //------------------------------------------------------------------------------------------------------
    // numTimes
    //------------------------------------------------------------------------------------------------------
    std::uint32_t numTimes = 0;
    file.read(reinterpret_cast< char*>(&numTimes), sizeof(std::uint32_t));
    _res << "\t\t- num. times: " << numTimes << std::endl;

    //------------------------------------------------------------------------------------------------------
    // idSystoleBegin
    //------------------------------------------------------------------------------------------------------
    std::uint32_t idSystoleBegin = 0;
    file.read(reinterpret_cast<char*>(&idSystoleBegin), sizeof(std::uint32_t));

    //------------------------------------------------------------------------------------------------------
    // msSystoleBegin
    //------------------------------------------------------------------------------------------------------
    double msSystoleBegin = 0;
    file.read(reinterpret_cast<char*>(&msSystoleBegin), sizeof(double));

    _res << "\t\t- systole begin (= diastole end): " << msSystoleBegin << " [ms] (time point id " << idSystoleBegin << ")" << std::endl;

    //------------------------------------------------------------------------------------------------------
    // idSystoleEnd
    //------------------------------------------------------------------------------------------------------
    std::uint32_t idSystoleEnd = 0;
    file.read(reinterpret_cast<char*>(&idSystoleEnd), sizeof(std::uint32_t));

    //------------------------------------------------------------------------------------------------------
    // msSystoleEnd
    //------------------------------------------------------------------------------------------------------
    double msSystoleEnd = 0;
    file.read(reinterpret_cast<char*>(&msSystoleEnd), sizeof(double));

    _res << "\t\t- systole end (= diastole begin): " << msSystoleEnd << " [ms] (time point id " << idSystoleEnd << ")" << std::endl;

    //------------------------------------------------------------------------------------------------------
    // numVessels
    //------------------------------------------------------------------------------------------------------
    std::uint32_t numVessels = 0;
    file.read(reinterpret_cast<char*>(&numVessels), sizeof(std::uint32_t));

    _res << "\t\t- num. vessels: " << numVessels << std::endl;

    //------------------------------------------------------------------------------------------------------
    // mean axial velocity per vessel
    //------------------------------------------------------------------------------------------------------
    std::vector<double> dbuffer(numVessels * numTimes);
    file.read(reinterpret_cast<char*>(dbuffer.data()), dbuffer.size() * sizeof(double));

    for (unsigned int vid = 0; vid < numVessels; ++vid)
    {
        const unsigned int off = vid * numTimes;
        _res << "\t\t- mean axial velocity [m/s] per time in vessel " << vid << ": ";

        for (unsigned int t = 0; t < NUM_DEMO; ++t)
        { _res << dbuffer[off + t] << ", "; }
        _res << "..." << std::endl;
    } // for vid : num vessels

    file.close();

    return true;
}

bool ImporterScientific::read_venc(std::string_view filepath)
{
    /*
     * [1] x [uint16] : dicom image id of 3D+T flow image X (LR)
     * [1] x [double] : venc 3D+T flow image X (LR)
     * [1] x [uint16] : dicom image id of 3D+T flow image Y (AP)
     * [1] x [double] : venc 3D+T flow image Y (AP)
     * [1] x [uint16] : dicom image id of 3D+T flow image Z (FH)
     * [1] x [double] : venc 3D+T flow image Z (FH)
     * [1] x [uint8] : num2DTFlowImages
     * for num2DTFlowImages
     *      [1] x [uint16] : dicom image id
     *      [1] x [double] : venc
     */

    if (!std::filesystem::exists(filepath.data()))
    {
        _res << "\t- no venc (path \"" << filepath.data() << "\")" << std::endl;
        return false;
    }

    _res << "\t- reading venc (path \"" << filepath.data() << "\")" << std::endl;

    std::ifstream file(filepath.data(), std::ios_base::in | std::ios_base::binary);

    if (!file.good())
    {
        _res << "\t\tFAILED! Could not open file!" << std::endl;
        return false;
    }

    //------------------------------------------------------------------------------------------------------
    // vencs of 3D+T flow images
    //------------------------------------------------------------------------------------------------------
    _res << "\t\t- VENCs of 3D+T flow images:" << std::endl;

    std::vector<std::uint16_t> dcmImgIdsFlow3DT(3);
    std::vector<double> vencsFlow3DT(3);

    for (unsigned int i = 0; i < 3; ++i)
    {
        file.read(reinterpret_cast< char*>(&dcmImgIdsFlow3DT[i]), sizeof(std::uint16_t));
        file.read(reinterpret_cast< char*>(&vencsFlow3DT[i]), sizeof(double));
    }

    _res << "\t\t\t- X (LR) image (ID " << dcmImgIdsFlow3DT[0] << "): " << vencsFlow3DT[0] << " [m/s]" << std::endl;
    _res << "\t\t\t- Y (AP) image (ID " << dcmImgIdsFlow3DT[1] << "): " << vencsFlow3DT[1] << " [m/s]" << std::endl;
    _res << "\t\t\t- Z (FH) image (ID " << dcmImgIdsFlow3DT[2] << "): " << vencsFlow3DT[2] << " [m/s]" << std::endl;

    //------------------------------------------------------------------------------------------------------
    // num2DTFlowImages
    //------------------------------------------------------------------------------------------------------
    std::uint8_t num2DTFlowImages = 0;
    file.read(reinterpret_cast<char*>(&num2DTFlowImages), sizeof(std::uint8_t));
    _res << "\t\t- num. 2D+T flow images: " << static_cast<int>(num2DTFlowImages) << std::endl;

    //------------------------------------------------------------------------------------------------------
    // vencs of 2D+T flow images
    //------------------------------------------------------------------------------------------------------
    std::vector<std::uint16_t> dcmImgIdsFlow2DT(num2DTFlowImages);
    std::vector<double> vencsFlow2DT(num2DTFlowImages);

    for (unsigned int i = 0; i < num2DTFlowImages; ++i)
    {
        file.read(reinterpret_cast< char*>(&dcmImgIdsFlow2DT[i]), sizeof(std::uint16_t));
        file.read(reinterpret_cast< char*>(&vencsFlow2DT[i]), sizeof(double));
    }

    for (unsigned int i = 0; i < num2DTFlowImages; ++i)
    { _res << "\t\t\t- Image " << i << " (DICOM image ID" << dcmImgIdsFlow2DT[i] << ") : VENC " << vencsFlow2DT[i] << " [m/s]" << std::endl; }

    file.close();

    return true;
}

std::string ImporterScientific::read_all()
{
    /*
     * clean up
     */
    _res.str("");
    _res << std::fixed;
    _res.precision(2);

    _vessel_names.clear();

    /*
     * clear dir path
     */
    _dir = bk::string_utils::replace(_dir, "\\", "/");

    if (bk::string_utils::ends_with(_dir, "/"))
    { _dir = bk::string_utils::chop_back(_dir, 1); }

    _res << "Reading directory \"" << _dir << "\"" << std::endl;

    /*
     * iterate directory -> find available vessels
     */
    for (auto& dirIt : std::filesystem::directory_iterator(_dir))
    {
        if (dirIt.is_directory())
        { _vessel_names.emplace_back(dirIt.path().filename().string()); }
    }

    _res << "\t- found " << _vessel_names.size() << " vessel(s): ";
    for (std::string_view vname: _vessel_names)
    { _res << "\"" << vname << "\" "; }
    _res << std::endl;

    /*
     * read dataset
     */
    read_dataset_filter_tags(_dir + "/dataset_tags.txt");
    read_dicom_tags(_dir + "/dicom_tags_3dt_flow");
    read_venc(_dir + "/venc");
    read_cardiac_cycle_definition(_dir + "/cardiac_cycle");
    read_static_tissue_mask(_dir + "/static_tissue_mask_in_flowfield_size");
    read_static_tissue_ivsd_thresholds(_dir + "/static_tissue_ivsd_thresholds");
    read_phase_wrapped_voxels(_dir + "/phase_wraps_3dt");
    read_flowfield(_dir + "/flowfield");
    read_velocity_offset_correction_3dt(_dir + "/velocity_offset_correction_3dt.voc");
    read_flow2dt_images();
    read_magnitude_tmip(_dir + "/magnitude3dt_tmip");
    read_anatomical_images();
    read_pressure_map(_dir + "/pressuremap");
    read_rotation_direction_map(_dir + "/rotationdirection");
    read_axial_velocity_map(_dir + "/axialvelocity");
    read_cos_angle_to_centerline_map(_dir + "/cosangletocenterline");
    read_turbulent_kinetic_energy_map(_dir + "/tke");
    read_ivsd(_dir + "/ivsd");
    //read_dicom_tags(_dir + "/dicom_tags_3dt_anatomical");
    //read_dicom_tags(_dir + "/dicom_tags_3dt_magnitude");
    //read_dicom_tags(_dir + "/dicom_tags_3dt_signal_intensity");
    //read_dicom_tags(_dir + "/dicom_tags_3d_anatomical");
    //read_dicom_tags(_dir + "/dicom_tags_2dt_flow");
    //read_dicom_tags(_dir + "/dicom_tags_2dt_anatomical");
    //read_dicom_tags(_dir + "/dicom_tags_2d_anatomical");
    read_flow_statistics(_dir + "/flow_stats");

    /*
     * read vessels
     */
    for (std::string_view vname: _vessel_names)
    {
        const std::string vesselPath = _dir + "/" + vname.data() + "/";

        _res << "-----------------------------------------------------------------------------------------------------------------------" << std::endl;
        _res << "-----------------------------------------------------------------------------------------------------------------------" << std::endl;
        _res << "Reading vessel \"" << vname << "\" (path \"" << vesselPath << "\")" << std::endl;

        read_mesh(vesselPath + "mesh");
        read_centerline_start_end_ids_on_mesh(vesselPath + "centerline_seed_target_ids_on_mesh");
        read_centerlines(vesselPath + "centerlines");
        read_flow_jet(vesselPath + "flowjets");
        read_pathlines(vesselPath + "pathlines");
        read_landmark_measuring_planes(vesselPath + "measuring_planes");
        read_segmentation(vesselPath + "segmentation");
        read_segmentation_info(vesselPath + "segmentation_info.txt");
        read_segmentation_graphcut_inside_outside_ids(vesselPath + "graphcut_segmentation_inside_outside_ids");
        read_segmentation_in_flowfield_size(vesselPath + "segmentation_in_flowfield_size");
        read_vessel_section_segmentation_in_flowfield_size(vesselPath + "vessel_section_segmentation_in_flowfield_size");
        read_vessel_section_segmentation_semantics(vesselPath + "vessel_section_info.txt");
    } // for vesselNames

    return _res.str();
}