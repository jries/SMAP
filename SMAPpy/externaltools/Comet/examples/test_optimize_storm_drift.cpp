#include <iostream>
#include <chrono>
#include <random>
#include <string>
#include "../gpu_optimize_storm_drift.h"



void test_3d_example()
{
    int n_molecules = 20;
    int n_loc_per_molecule = 40;
    int n_timepoints = 100;

    float loc_precision = 0.05f;

    float xd_min = 0.0f;
    float xd_size = 10.0f;
    float yd_min = 0.0f;
    float yd_size = 10.0f;
    float zd_min = 0.0f;
    float zd_size = 10.0f;

    float xd_max = xd_min + xd_size;
    float yd_max = yd_min + yd_size;
    float zd_max = zd_min + zd_size;

    // random localizations
    std::chrono::high_resolution_clock::time_point time_now
        = std::chrono::high_resolution_clock::now();

    unsigned int time_now_in_ms
        = static_cast<unsigned int>
        (std::chrono::duration_cast<std::chrono::milliseconds>(time_now.time_since_epoch()).count());

    std::mt19937 rng;
    rng.seed(time_now_in_ms);

    std::uniform_real_distribution< float > uniform_dist(0.0f, 1.0f);
    std::normal_distribution< float > normal_dist(0.0f, 1.0f);

    std::vector<float> mol_x_coords(n_molecules);
    std::vector<float> mol_y_coords(n_molecules);
    std::vector<float> mol_z_coords(n_molecules);

    for (int i = 0; i < n_molecules; i++)
    {
        mol_x_coords[i] = uniform_dist(rng) * xd_size + xd_min;
        mol_y_coords[i] = uniform_dist(rng) * yd_size + yd_min;
        mol_z_coords[i] = uniform_dist(rng) * zd_size + zd_min;
    }

    int n_locs = n_molecules * n_loc_per_molecule;

    std::vector<int> loc_times(n_locs);
    std::vector<float> loc_errors_x(n_locs);
    std::vector<float> loc_errors_y(n_locs);
    std::vector<float> loc_errors_z(n_locs);

    for (int i = 0; i < n_locs; i++) 
    { 
        loc_times[i] = (int)floorf(uniform_dist(rng) * n_timepoints);
        loc_errors_x[i] = normal_dist(rng) * loc_precision;
        loc_errors_y[i] = normal_dist(rng) * loc_precision;
        loc_errors_z[i] = normal_dist(rng) * loc_precision;
    }

    float drift_scale_x = 2 * loc_precision;
    float drift_scale_y = 5 * loc_precision;
    float drift_scale_z = -3 * loc_precision;

    std::vector<float> drift_trajectory(3*n_timepoints);
    
    for (int i = 0; i < n_timepoints; i++)
    {
        drift_trajectory[i] = (i*drift_scale_x) / (n_timepoints - 1);
        drift_trajectory[i + n_timepoints] = (i*drift_scale_y) / (n_timepoints - 1);
        drift_trajectory[i + (2 * n_timepoints)] = (i*drift_scale_z) / (n_timepoints - 1);
    }

    std::vector<float> initial_localizations_x(n_locs);
    std::vector<float> initial_localizations_y(n_locs);
    std::vector<float> initial_localizations_z(n_locs);

    int mol_id;

    for (int i = 0; i < n_locs; i++)
    {
        mol_id = (int)floorf( ((float) i) / n_loc_per_molecule );
        initial_localizations_x[i] = mol_x_coords[mol_id] + drift_trajectory[loc_times[i]] + loc_errors_x[i];
        initial_localizations_y[i] = mol_y_coords[mol_id] + drift_trajectory[loc_times[i]+n_timepoints] + loc_errors_y[i];
        initial_localizations_z[i] = mol_z_coords[mol_id] + drift_trajectory[loc_times[i]+(2*n_timepoints)] + loc_errors_z[i];
    }

    size_t n_pairs = ((size_t) n_locs) * (n_locs-1) / 2;

    std::vector<int> pair_indices_i(n_pairs);
    std::vector<int> pair_indices_j(n_pairs);

    size_t tmp_count = 0;

    for (int i = 0; i < (n_locs-1); i++)
    {
        for (int j = i+1; j < n_locs; j++)
        {
            pair_indices_i[tmp_count] = i;
            pair_indices_j[tmp_count] = j;
            tmp_count++;
        }
    }

    // initialize the optimization

    int return_value;

    return_value = gpu_opt_storm_drift_initialize_3d(
        n_locs,
        n_timepoints,
        initial_localizations_x.data(),
        initial_localizations_y.data(),
        initial_localizations_z.data(),
        loc_times.data(),
        tmp_count,
        pair_indices_i.data(),
        pair_indices_j.data());

    float gaussian_scale = 3.0f * loc_precision;

    float output_cost_function;

    int flag_calculate_derivatives = 1;

    std::vector<float> output_derivatives(3*n_timepoints);

    return_value = gpu_opt_storm_drift_compute_3d(
        n_locs,
        n_timepoints,
        tmp_count,
        gaussian_scale,
        drift_trajectory.data(),
        &output_cost_function,
        flag_calculate_derivatives,
        output_derivatives.data());

    return_value = gpu_opt_storm_drift_free_3d();

}


int main(int argc, char const *argv[]) {
	
    int aa;

    aa = 0;

    test_3d_example();

    std::cout << std::endl << "Test completed!" << std::endl;
    std::cout << "Press ENTER to exit" << std::endl;
    std::getchar();

    return 0;

}

