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

#ifndef BLOODLINE_IMPORTERSCIENTIFIC_H
#define BLOODLINE_IMPORTERSCIENTIFIC_H

#include <fstream>
#include <memory>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

class ImporterScientific
{
    //====================================================================================================
    //===== DEFINITIONS
    //====================================================================================================
    static constexpr unsigned int NUM_DEMO = 3;

    //====================================================================================================
    //===== MEMBERS
    //====================================================================================================
    std::string _dir;
    std::vector<std::string> _vessel_names;
    std::stringstream _res;

    //====================================================================================================
    //===== CONSTRUCTORS & DESTRUCTOR
    //====================================================================================================
  public:
    ImporterScientific();
    ImporterScientific(const ImporterScientific&) = delete;
    ImporterScientific(ImporterScientific&&);

    ~ImporterScientific();

    //====================================================================================================
    //===== GETTER
    //====================================================================================================
    [[nodiscard]] std::string result() const;

    //====================================================================================================
    //===== SETTER
    //====================================================================================================
    [[maybe_unused]] ImporterScientific& operator=(const ImporterScientific&) = delete;
    [[maybe_unused]] ImporterScientific& operator=(ImporterScientific&&);

    void set_dir(std::string_view dir);

    //====================================================================================================
    //===== FUNCTIONS
    //====================================================================================================
    void _read_nd_scalar_image_in_sparse_matrix_style(std::ifstream& file);

    [[maybe_unused]] bool read_mesh(std::string_view filepath);
    [[maybe_unused]] bool read_centerlines(std::string_view filepath);
    [[maybe_unused]] bool read_landmark_measuring_planes(std::string_view filepath);
    [[maybe_unused]] bool read_pathlines(std::string_view filepath);
    [[maybe_unused]] bool read_flowfield(std::string_view filepath);
    [[maybe_unused]] bool read_pressure_map(std::string_view filepath);
    [[maybe_unused]] bool read_rotation_direction_map(std::string_view filepath);
    [[maybe_unused]] bool read_axial_velocity_map(std::string_view filepath);
    [[maybe_unused]] bool read_cos_angle_to_centerline_map(std::string_view filepath);
    [[maybe_unused]] bool read_turbulent_kinetic_energy_map(std::string_view filepath);
    [[maybe_unused]] bool read_flow_jet(std::string_view filepath);
    [[maybe_unused]] bool read_ivsd(std::string_view filepath);
    [[maybe_unused]] bool read_magnitude_tmip(std::string_view filepath);
    [[maybe_unused]] bool read_anatomical_images();
    [[maybe_unused]] bool read_flow2dt_images();
    [[maybe_unused]] bool read_flow_statistics(std::string_view filepath);
    [[maybe_unused]] bool read_segmentation(std::string_view filepath);
    [[maybe_unused]] bool read_segmentation_info(std::string_view filepath);
    [[maybe_unused]] bool read_segmentation_graphcut_inside_outside_ids(std::string_view filepath);
    [[maybe_unused]] bool read_segmentation_in_flowfield_size(std::string_view filepath);
    [[maybe_unused]] bool read_vessel_section_segmentation_in_flowfield_size(std::string_view filepath);
    [[maybe_unused]] bool read_vessel_section_segmentation_semantics(std::string_view filepath);
    [[maybe_unused]] bool read_centerline_start_end_ids_on_mesh(std::string_view filepath);
    [[maybe_unused]] bool read_static_tissue_mask(std::string_view filepath);
    [[maybe_unused]] bool read_static_tissue_ivsd_thresholds(std::string_view filepath);
    [[maybe_unused]] bool read_dataset_filter_tags(std::string_view filepath);
    [[maybe_unused]] bool read_phase_wrapped_voxels(std::string_view filepath);
    [[maybe_unused]] bool read_velocity_offset_correction_3dt(std::string_view filepath);
    [[maybe_unused]] bool read_dicom_tags(std::string_view filepath);
    [[maybe_unused]] bool read_cardiac_cycle_definition(std::string_view filepath);
    [[maybe_unused]] bool read_venc(std::string_view filepath);

    [[maybe_unused]] std::string read_all();

}; // class ImporterScientific

#endif //BLOODLINE_IMPORTERSCIENTIFIC_H
