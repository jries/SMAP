#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include <thrust/copy.h>
#include <algorithm>
#include <cuda_runtime.h>
//#include <nvtoolsext.h>
#include "gpu_optimize_storm_drift.h"


/* Define global variables */

int static_n_coordinates;
int static_n_timepoints;
size_t static_n_coordinate_pairs;
__device__ float * static_d_coords_x;
__device__ float * static_d_coords_y;
__device__ float * static_d_coords_z;
__device__ int * static_d_coords_time;
__device__ int * static_d_pair_indices_i;
__device__ int * static_d_pair_indices_j;


__global__ void calculate_osd_cost_function_2d(
    int const n_coordinates,
    int const n_timepoints,
    size_t const n_coordinate_pairs,
    size_t const process_start_index,
    float gaussian_scale,
    float * d_drift_trajectory,
    float * d_wa_function_values,
    int flag_calculate_derivatives,
    float * d_derivatives)
{

    int const n_threads_per_block = blockDim.x;
    int const block_index = blockIdx.x;
    int const thread_index = threadIdx.x;

    int const proc_id = process_start_index + (block_index * n_threads_per_block) + thread_index;

    // is this a valid process ID?
    bool const process_valid = (proc_id < n_coordinate_pairs);

    int const src_index = (process_valid) ? proc_id : 0;

    int const coord_index_i = static_d_pair_indices_i[src_index];
    int const coord_index_j = static_d_pair_indices_j[src_index];

    int const coord_t_i = static_d_coords_time[coord_index_i];
    float const coord_x_i = static_d_coords_x[coord_index_i] + d_drift_trajectory[coord_t_i];
    float const coord_y_i = static_d_coords_y[coord_index_i] + d_drift_trajectory[coord_t_i + n_timepoints];

    int const coord_t_j = static_d_coords_time[coord_index_j];
    float const coord_x_j = static_d_coords_x[coord_index_j] + d_drift_trajectory[coord_t_j];
    float const coord_y_j = static_d_coords_y[coord_index_j] + d_drift_trajectory[coord_t_j + n_timepoints];

    float const delta_x = coord_x_i - coord_x_j;
    float const delta_y = coord_y_i - coord_y_j;

    float const dist_sq = delta_x * delta_x + delta_y * delta_y;

    float cost_fn_value = exp(-(dist_sq) / gaussian_scale);

    // if this is a valid process, store the results
    if (process_valid)
    {
        // store the cost function
        d_wa_function_values[proc_id - process_start_index] = -cost_fn_value;

        // store the derivatives
        float old_value;

        if (flag_calculate_derivatives == 1)
        {
            if (coord_t_i != coord_t_j)
            {
                cost_fn_value = cost_fn_value * (2.0 / gaussian_scale);

                old_value = atomicAdd(d_derivatives + coord_t_i, (cost_fn_value * delta_x));
                old_value = atomicAdd(d_derivatives + coord_t_j, -(cost_fn_value * delta_x));

                old_value = atomicAdd(d_derivatives + coord_t_i + n_timepoints, (cost_fn_value * delta_y));
                old_value = atomicAdd(d_derivatives + coord_t_j + n_timepoints, -(cost_fn_value * delta_y));
            }
        }
    }

    __syncthreads();
}


__global__ void calculate_osd_cost_function_3d(
    int const n_coordinates,
    int const n_timepoints,
    size_t const n_coordinate_pairs,
    size_t const process_start_index,
    float gaussian_scale, 
    float * d_drift_trajectory, 
    float * d_wa_function_values,
    int flag_calculate_derivatives,
    float * d_derivatives)
{
    
    int const n_threads_per_block = blockDim.x;
    int const block_index = blockIdx.x;
    int const thread_index = threadIdx.x;

    int const proc_id = process_start_index + (block_index * n_threads_per_block) + thread_index;

    // is this a valid process ID?
    bool const process_valid = (proc_id < n_coordinate_pairs);

    int const src_index = (process_valid) ? proc_id : 0;

    int const coord_index_i = static_d_pair_indices_i[src_index];
    int const coord_index_j = static_d_pair_indices_j[src_index];

    int const coord_t_i = static_d_coords_time[coord_index_i];
    float const coord_x_i = static_d_coords_x[coord_index_i] + d_drift_trajectory[coord_t_i];
    float const coord_y_i = static_d_coords_y[coord_index_i] + d_drift_trajectory[coord_t_i + n_timepoints];
    float const coord_z_i = static_d_coords_z[coord_index_i] + d_drift_trajectory[coord_t_i + 2 * n_timepoints];

    int const coord_t_j = static_d_coords_time[coord_index_j];
    float const coord_x_j = static_d_coords_x[coord_index_j] + d_drift_trajectory[coord_t_j];
    float const coord_y_j = static_d_coords_y[coord_index_j] + d_drift_trajectory[coord_t_j + n_timepoints];
    float const coord_z_j = static_d_coords_z[coord_index_j] + d_drift_trajectory[coord_t_j + 2 * n_timepoints];

    float const delta_x = coord_x_i - coord_x_j;
    float const delta_y = coord_y_i - coord_y_j;
    float const delta_z = coord_z_i - coord_z_j;

    float const dist_sq = delta_x * delta_x + delta_y * delta_y + delta_z * delta_z;

    float cost_fn_value = exp(-(dist_sq)/gaussian_scale);

    // if this is a valid process, store the results
    if (process_valid)
    {
        // store the cost function
        d_wa_function_values[proc_id - process_start_index] = -cost_fn_value;

        // store the derivatives
        float old_value;

        if (flag_calculate_derivatives == 1)
        {
            if (coord_t_i != coord_t_j)
            {
                cost_fn_value = cost_fn_value * (2.0 / gaussian_scale);

                old_value = atomicAdd(d_derivatives + coord_t_i, (cost_fn_value * delta_x));
                old_value = atomicAdd(d_derivatives + coord_t_j, -(cost_fn_value * delta_x));

                old_value = atomicAdd(d_derivatives + coord_t_i + n_timepoints, (cost_fn_value * delta_y));
                old_value = atomicAdd(d_derivatives + coord_t_j + n_timepoints, -(cost_fn_value * delta_y));

                old_value = atomicAdd(d_derivatives + coord_t_i + (2 * n_timepoints), (cost_fn_value * delta_z));
                old_value = atomicAdd(d_derivatives + coord_t_j + (2 * n_timepoints), -(cost_fn_value * delta_z));
            }
        }
    }
 
    __syncthreads();
}


int gpu_opt_storm_drift_compute_2d(
    int n_coordinates,
    int n_timepoints,
    size_t n_coordinate_pairs,
    float gaussian_scale,
    float * drift_trajectory,
    float * output_cost_function,
    int flag_calculate_derivatives,
    float * output_derivatives)
{

    cudaError_t cuda_status;


    //nvtxMarkA("Start of gpu_opt_storm_drift_compute");

    if (n_coordinates != static_n_coordinates)
    {
        throw std::runtime_error("Invalid number of coordinates");
    }

    if (n_timepoints != static_n_timepoints)
    {
        throw std::runtime_error("Invalid number of time points");
    }

    if (n_coordinate_pairs != static_n_coordinate_pairs)
    {
        throw std::runtime_error("Invalid number of coordinate pairs");
    }


    // Copy the drift trajectory to the GPU

    float * d_drift_trajectory{ nullptr };

    cuda_status = cudaMalloc(&d_drift_trajectory, 2 * n_timepoints * sizeof(float));
    if (cuda_status != cudaSuccess)
    {
        throw std::runtime_error(cudaGetErrorString(cuda_status));
    }

    cuda_status = cudaMemcpy(d_drift_trajectory, drift_trajectory, 2 * n_timepoints * sizeof(float), cudaMemcpyHostToDevice);
    if (cuda_status != cudaSuccess)
    {
        throw std::runtime_error(cudaGetErrorString(cuda_status));
    }


    // Initialize an array in which to store the derivatives (this will always be a float array)

    float * d_derivatives{ nullptr };

    if (flag_calculate_derivatives == 1)
    {
        cuda_status = cudaMalloc(&d_derivatives, 2 * n_timepoints * sizeof(float));
        if (cuda_status != cudaSuccess)
        {
            throw std::runtime_error(cudaGetErrorString(cuda_status));
        }

        cuda_status = cudaMemset(d_derivatives, 0, 2 * n_timepoints * sizeof(float));
        if (cuda_status != cudaSuccess)
        {
            throw std::runtime_error(cudaGetErrorString(cuda_status));
        }
    }


    // Initialize the cost function value (stored internally as double precision)
    double tmp_cost_function = 0.0;

    // Divide the work into chunks

    int const proc_chunk_size = 16777216;
    int const n_chunks = (int)std::ceil((double)n_coordinate_pairs / (double)proc_chunk_size);

    // allocate thrust device vectors as working arrays
    thrust::device_vector<float> dev_vec_wa_function_values(proc_chunk_size);
    thrust::device_vector<double> dev_vec_wa_function_values_dbl(proc_chunk_size);

    float * d_wa_function_values = thrust::raw_pointer_cast(dev_vec_wa_function_values.data());

    for (int i = 0; i < n_chunks; i++)
    {

        size_t start_index = i*((size_t)proc_chunk_size);
        size_t end_index = std::min(start_index + proc_chunk_size - 1, n_coordinate_pairs - 1);
        int cur_chunk_n_pairs = (int)(end_index - start_index + 1);

        // clear the working array
        thrust::fill(dev_vec_wa_function_values.begin(), dev_vec_wa_function_values.end(), 0.0);

        // initialize the number of blocks and threads
        int const n_threads_per_block = 128;
        int const n_blocks = (int)std::ceil(double(cur_chunk_n_pairs) / double(n_threads_per_block));

        // calculate the gaussian contribution from each process
        calculate_osd_cost_function_2d << < n_blocks, n_threads_per_block >> > (
            n_coordinates,
            n_timepoints,
            n_coordinate_pairs,
            start_index,
            gaussian_scale,
            d_drift_trajectory,
            d_wa_function_values,
            flag_calculate_derivatives,
            d_derivatives);
        
        // copy the working array to a double precision array and sum as double precision
        thrust::copy(dev_vec_wa_function_values.begin(), dev_vec_wa_function_values.end(), dev_vec_wa_function_values_dbl.begin());

        // sum the cost function values in the double precision array and add the result to the accumulated cost function
        tmp_cost_function += thrust::reduce(dev_vec_wa_function_values_dbl.begin(), dev_vec_wa_function_values_dbl.end(), 0.0, thrust::plus<double>());

    }

    *output_cost_function = (float) tmp_cost_function;

    // copy the derivatives to host memory

    if (flag_calculate_derivatives == 1)
    {
        cuda_status = cudaMemcpy(output_derivatives, d_derivatives, 2 * n_timepoints * sizeof(float), cudaMemcpyDeviceToHost);
        if (cuda_status != cudaSuccess)
        {
            throw std::runtime_error(cudaGetErrorString(cuda_status));
        }

        cuda_status = cudaFree(d_derivatives);
        if (cuda_status != cudaSuccess)
        {
            throw std::runtime_error(cudaGetErrorString(cuda_status));
        }

    }


    cuda_status = cudaFree(d_drift_trajectory);
    if (cuda_status != cudaSuccess)
    {
        throw std::runtime_error(cudaGetErrorString(cuda_status));
    }

    // deallocate memory associated with the device vectors
    dev_vec_wa_function_values.clear();
    dev_vec_wa_function_values.shrink_to_fit();

    dev_vec_wa_function_values_dbl.clear();
    dev_vec_wa_function_values_dbl.shrink_to_fit();

    return 0;

}


int gpu_opt_storm_drift_compute_3d(
    int n_coordinates,
    int n_timepoints,
    size_t n_coordinate_pairs,
    float gaussian_scale,
    float * drift_trajectory,
    float * output_cost_function, 
    int flag_calculate_derivatives,
    float * output_derivatives)
{

    cudaError_t cuda_status;


    //nvtxMarkA("Start of gpu_opt_storm_drift_compute");

    if (n_coordinates != static_n_coordinates)
    {
        throw std::runtime_error("Invalid number of coordinates");
    }

    if (n_timepoints != static_n_timepoints)
    {
        throw std::runtime_error("Invalid number of time points");
    }

    if (n_coordinate_pairs != static_n_coordinate_pairs)
    {
        throw std::runtime_error("Invalid number of coordinate pairs");
    }


    // copy the drift trajectory to the GPU
    float * d_drift_trajectory{ nullptr };

    cuda_status = cudaMalloc(&d_drift_trajectory, 3 * n_timepoints * sizeof(float));
    if (cuda_status != cudaSuccess)
    {
        throw std::runtime_error(cudaGetErrorString(cuda_status));
    }

    cuda_status = cudaMemcpy(d_drift_trajectory, drift_trajectory, 3 * n_timepoints * sizeof(float), cudaMemcpyHostToDevice);
    if (cuda_status != cudaSuccess)
    {
        throw std::runtime_error(cudaGetErrorString(cuda_status));
    }


    // initialize an array in which to store the derivatives
    float * d_derivatives{ nullptr };

    if (flag_calculate_derivatives == 1)
    {
        cuda_status = cudaMalloc(&d_derivatives, 3 * n_timepoints * sizeof(float));
        if (cuda_status != cudaSuccess)
        {
            throw std::runtime_error(cudaGetErrorString(cuda_status));
        }

        cuda_status = cudaMemset(d_derivatives, 0, 3 * n_timepoints * sizeof(float));
        if (cuda_status != cudaSuccess)
        {
            throw std::runtime_error(cudaGetErrorString(cuda_status));
        }
    }


    // initialize the cost function value
    double tmp_cost_function = 0.0;

    // divide the work into chunks
    int const proc_chunk_size = 16777216;
    int const n_chunks = (int)std::ceil((double)n_coordinate_pairs / (double)proc_chunk_size);


    // allocate thrust device vectors as working arrays
    thrust::device_vector<float> dev_vec_wa_function_values(proc_chunk_size);
    thrust::device_vector<double> dev_vec_wa_function_values_dbl(proc_chunk_size);

    float * d_wa_function_values = thrust::raw_pointer_cast(dev_vec_wa_function_values.data());

    for (int i = 0; i < n_chunks; i++)
    {

        size_t start_index = i*((size_t)proc_chunk_size);
        size_t end_index = std::min(start_index + proc_chunk_size - 1, n_coordinate_pairs - 1);
        int cur_chunk_n_pairs = (int)(end_index - start_index + 1);

        // clear the working array
        thrust::fill(dev_vec_wa_function_values.begin(), dev_vec_wa_function_values.end(), 0.0);

        // initialize the number of blocks and threads
        int const n_threads_per_block = 128;
        int const n_blocks = (int)std::ceil(double(cur_chunk_n_pairs) / double(n_threads_per_block));

        //calculate the gaussian contribution from each process
        calculate_osd_cost_function_3d <<< n_blocks, n_threads_per_block >>> (
            n_coordinates,
            n_timepoints,
            n_coordinate_pairs,
            start_index,
            gaussian_scale,
            d_drift_trajectory,
            d_wa_function_values,
            flag_calculate_derivatives,
            d_derivatives);

        // copy the working array to a double precision array and sum as double precision
        thrust::copy(dev_vec_wa_function_values.begin(), dev_vec_wa_function_values.end(), dev_vec_wa_function_values_dbl.begin());

        // sum the cost function values in the double precision array and add the result to the accumulated cost function
        tmp_cost_function += thrust::reduce(dev_vec_wa_function_values_dbl.begin(), dev_vec_wa_function_values_dbl.end(), 0.0, thrust::plus<double>());

    }


    *output_cost_function = (float) tmp_cost_function;

    if (flag_calculate_derivatives == 1)
    {
        // copy the derivatives to host memory
        cuda_status = cudaMemcpy(output_derivatives, d_derivatives, (size_t)3 * n_timepoints * sizeof(float), cudaMemcpyDeviceToHost);
        if (cuda_status != cudaSuccess)
        {
            throw std::runtime_error(cudaGetErrorString(cuda_status));
        }

        cuda_status = cudaFree(d_derivatives);
        if (cuda_status != cudaSuccess)
        {
            throw std::runtime_error(cudaGetErrorString(cuda_status));
        }

    }


    cuda_status = cudaFree(d_drift_trajectory);
    if (cuda_status != cudaSuccess)
    {
        throw std::runtime_error(cudaGetErrorString(cuda_status));
    }

    // deallocate memory associated with the device vectors
    dev_vec_wa_function_values.clear();
    dev_vec_wa_function_values.shrink_to_fit();

    dev_vec_wa_function_values_dbl.clear();
    dev_vec_wa_function_values_dbl.shrink_to_fit();

    return 0;

}


int gpu_opt_storm_drift_initialize_2d(
    int n_coordinates,
    int n_timepoints,
    float * coordinates_x,
    float * coordinates_y,
    int * coordinates_time,
    size_t n_coordinate_pairs,
    int * pair_indices_i,
    int * pair_indices_j)
{

    //nvtxMarkA("Start of gpu_opt_storm_drift_initialize");

    cudaError_t cuda_status;


    // allocate space for storate arrays
    float * d_coordinates_x{ nullptr };
    float * d_coordinates_y{ nullptr };
    int * d_coordinates_time{ nullptr };
    int * d_pair_indices_i{ nullptr };
    int * d_pair_indices_j{ nullptr };


    cuda_status = cudaMalloc(&d_coordinates_x, n_coordinates * sizeof(float));
    if (cuda_status != cudaSuccess)
    {
        throw std::runtime_error(cudaGetErrorString(cuda_status));
    }

    cuda_status = cudaMalloc(&d_coordinates_y, n_coordinates * sizeof(float));
    if (cuda_status != cudaSuccess)
    {
        throw std::runtime_error(cudaGetErrorString(cuda_status));
    }

    cuda_status = cudaMalloc(&d_coordinates_time, n_coordinates * sizeof(int));
    if (cuda_status != cudaSuccess)
    {
        throw std::runtime_error(cudaGetErrorString(cuda_status));
    }

    cuda_status = cudaMalloc(&d_pair_indices_i, n_coordinate_pairs * sizeof(int));
    if (cuda_status != cudaSuccess)
    {
        throw std::runtime_error(cudaGetErrorString(cuda_status));
    }

    cuda_status = cudaMalloc(&d_pair_indices_j, n_coordinate_pairs * sizeof(int));
    if (cuda_status != cudaSuccess)
    {
        throw std::runtime_error(cudaGetErrorString(cuda_status));
    }


    // copy the data to the GPU 
    cuda_status = cudaMemcpy(d_coordinates_x, coordinates_x, n_coordinates * sizeof(float), cudaMemcpyHostToDevice);
    if (cuda_status != cudaSuccess)
    {
        throw std::runtime_error(cudaGetErrorString(cuda_status));
    }

    cuda_status = cudaMemcpy(d_coordinates_y, coordinates_y, n_coordinates * sizeof(float), cudaMemcpyHostToDevice);
    if (cuda_status != cudaSuccess)
    {
        throw std::runtime_error(cudaGetErrorString(cuda_status));
    }

    cuda_status = cudaMemcpy(d_coordinates_time, coordinates_time, n_coordinates * sizeof(int), cudaMemcpyHostToDevice);
    if (cuda_status != cudaSuccess)
    {
        throw std::runtime_error(cudaGetErrorString(cuda_status));
    }

    cuda_status = cudaMemcpy(d_pair_indices_i, pair_indices_i, n_coordinate_pairs * sizeof(int), cudaMemcpyHostToDevice);
    if (cuda_status != cudaSuccess)
    {
        throw std::runtime_error(cudaGetErrorString(cuda_status));
    }

    cuda_status = cudaMemcpy(d_pair_indices_j, pair_indices_j, n_coordinate_pairs * sizeof(int), cudaMemcpyHostToDevice);
    if (cuda_status != cudaSuccess)
    {
        throw std::runtime_error(cudaGetErrorString(cuda_status));
    }


    // store the device pointer addresses in global variables on the device
    cuda_status = cudaMemcpyToSymbol(static_d_coords_x, &d_coordinates_x, sizeof(d_coordinates_x), 0, cudaMemcpyHostToDevice);
    if (cuda_status != cudaSuccess)
    {
        throw std::runtime_error(cudaGetErrorString(cuda_status));
    }

    cuda_status = cudaMemcpyToSymbol(static_d_coords_y, &d_coordinates_y, sizeof(d_coordinates_y), 0, cudaMemcpyHostToDevice);
    if (cuda_status != cudaSuccess)
    {
        throw std::runtime_error(cudaGetErrorString(cuda_status));
    }

    cuda_status = cudaMemcpyToSymbol(static_d_coords_time, &d_coordinates_time, sizeof(d_coordinates_time), 0, cudaMemcpyHostToDevice);
    if (cuda_status != cudaSuccess)
    {
        throw std::runtime_error(cudaGetErrorString(cuda_status));
    }

    cuda_status = cudaMemcpyToSymbol(static_d_pair_indices_i, &d_pair_indices_i, sizeof(d_pair_indices_i), 0, cudaMemcpyHostToDevice);
    if (cuda_status != cudaSuccess)
    {
        throw std::runtime_error(cudaGetErrorString(cuda_status));
    }

    cuda_status = cudaMemcpyToSymbol(static_d_pair_indices_j, &d_pair_indices_j, sizeof(d_pair_indices_j), 0, cudaMemcpyHostToDevice);
    if (cuda_status != cudaSuccess)
    {
        throw std::runtime_error(cudaGetErrorString(cuda_status));
    }

    // store global scalar variables
    static_n_coordinates = n_coordinates;
    static_n_timepoints = n_timepoints;
    static_n_coordinate_pairs = n_coordinate_pairs;


    //nvtxMarkA("End of gpu_opt_storm_drift_initialize");

    return 0;
}


int gpu_opt_storm_drift_initialize_3d(
    int n_coordinates,
    int n_timepoints,
    float * coordinates_x,
    float * coordinates_y,
    float * coordinates_z, 
    int * coordinates_time, 
    size_t n_coordinate_pairs, 
    int * pair_indices_i, 
    int * pair_indices_j)
{

    //nvtxMarkA("Start of gpu_opt_storm_drift_initialize");

    cudaError_t cuda_status;


    // allocate space for storate arrays
    float * d_coordinates_x{ nullptr };
    float * d_coordinates_y{ nullptr };
    float * d_coordinates_z{ nullptr };
    int * d_coordinates_time{ nullptr };
    int * d_pair_indices_i{ nullptr };
    int * d_pair_indices_j{ nullptr };


    cuda_status = cudaMalloc(&d_coordinates_x, n_coordinates * sizeof(float));
    if (cuda_status != cudaSuccess)
    {
        throw std::runtime_error(cudaGetErrorString(cuda_status));
    }

    cuda_status = cudaMalloc(&d_coordinates_y, n_coordinates * sizeof(float));
    if (cuda_status != cudaSuccess)
    {
        throw std::runtime_error(cudaGetErrorString(cuda_status));
    }

    cuda_status = cudaMalloc(&d_coordinates_z, n_coordinates * sizeof(float));
    if (cuda_status != cudaSuccess)
    {
        throw std::runtime_error(cudaGetErrorString(cuda_status));
    }

    cuda_status = cudaMalloc(&d_coordinates_time, n_coordinates * sizeof(int));
    if (cuda_status != cudaSuccess)
    {
        throw std::runtime_error(cudaGetErrorString(cuda_status));
    }

    cuda_status = cudaMalloc(&d_pair_indices_i, n_coordinate_pairs * sizeof(int));
    if (cuda_status != cudaSuccess)
    {
        throw std::runtime_error(cudaGetErrorString(cuda_status));
    }

    cuda_status = cudaMalloc(&d_pair_indices_j, n_coordinate_pairs * sizeof(int));
    if (cuda_status != cudaSuccess)
    {
        throw std::runtime_error(cudaGetErrorString(cuda_status));
    }


    // copy the data to the GPU 
    cuda_status = cudaMemcpy(d_coordinates_x, coordinates_x, n_coordinates * sizeof(float), cudaMemcpyHostToDevice);
    if (cuda_status != cudaSuccess)
    {
        throw std::runtime_error(cudaGetErrorString(cuda_status));
    }

    cuda_status = cudaMemcpy(d_coordinates_y, coordinates_y, n_coordinates * sizeof(float), cudaMemcpyHostToDevice);
    if (cuda_status != cudaSuccess)
    {
        throw std::runtime_error(cudaGetErrorString(cuda_status));
    }

    cuda_status = cudaMemcpy(d_coordinates_z, coordinates_z, n_coordinates * sizeof(float), cudaMemcpyHostToDevice);
    if (cuda_status != cudaSuccess)
    {
        throw std::runtime_error(cudaGetErrorString(cuda_status));
    }

    cuda_status = cudaMemcpy(d_coordinates_time, coordinates_time, n_coordinates * sizeof(int), cudaMemcpyHostToDevice);
    if (cuda_status != cudaSuccess)
    {
        throw std::runtime_error(cudaGetErrorString(cuda_status));
    }

    cuda_status = cudaMemcpy(d_pair_indices_i, pair_indices_i, n_coordinate_pairs * sizeof(int), cudaMemcpyHostToDevice);
    if (cuda_status != cudaSuccess)
    {
        throw std::runtime_error(cudaGetErrorString(cuda_status));
    }

    cuda_status = cudaMemcpy(d_pair_indices_j, pair_indices_j, n_coordinate_pairs * sizeof(int), cudaMemcpyHostToDevice);
    if (cuda_status != cudaSuccess)
    {
        throw std::runtime_error(cudaGetErrorString(cuda_status));
    }


    // store the device pointer addresses in global variables on the device
    cuda_status = cudaMemcpyToSymbol(static_d_coords_x, &d_coordinates_x, sizeof(d_coordinates_x), 0, cudaMemcpyHostToDevice);
    if (cuda_status != cudaSuccess)
    {
        throw std::runtime_error(cudaGetErrorString(cuda_status));
    }

    cuda_status = cudaMemcpyToSymbol(static_d_coords_y, &d_coordinates_y, sizeof(d_coordinates_y), 0, cudaMemcpyHostToDevice);
    if (cuda_status != cudaSuccess)
    {
        throw std::runtime_error(cudaGetErrorString(cuda_status));
    }

    cuda_status = cudaMemcpyToSymbol(static_d_coords_z, &d_coordinates_z, sizeof(d_coordinates_z), 0, cudaMemcpyHostToDevice);
    if (cuda_status != cudaSuccess)
    {
        throw std::runtime_error(cudaGetErrorString(cuda_status));
    }

    cuda_status = cudaMemcpyToSymbol(static_d_coords_time, &d_coordinates_time, sizeof(d_coordinates_time), 0, cudaMemcpyHostToDevice);
    if (cuda_status != cudaSuccess)
    {
        throw std::runtime_error(cudaGetErrorString(cuda_status));
    }

    cuda_status = cudaMemcpyToSymbol(static_d_pair_indices_i, &d_pair_indices_i, sizeof(d_pair_indices_i), 0, cudaMemcpyHostToDevice);
    if (cuda_status != cudaSuccess)
    {
        throw std::runtime_error(cudaGetErrorString(cuda_status));
    }

    cuda_status = cudaMemcpyToSymbol(static_d_pair_indices_j, &d_pair_indices_j, sizeof(d_pair_indices_j), 0, cudaMemcpyHostToDevice);
    if (cuda_status != cudaSuccess)
    {
        throw std::runtime_error(cudaGetErrorString(cuda_status));
    }


    // store global variables
    static_n_coordinates = n_coordinates;
    static_n_timepoints = n_timepoints;
    static_n_coordinate_pairs = n_coordinate_pairs;


    //nvtxMarkA("End of gpu_opt_storm_drift_initialize");

    return 0;
}



int gpu_opt_storm_drift_free_2d()
{

    cudaError_t cuda_status;

    // free device memory
    static_n_coordinates = 0;
    static_n_timepoints = 0;
    static_n_coordinate_pairs = 0;

    float * d_coordinates_x{ nullptr };
    float * d_coordinates_y{ nullptr };
    int * d_coordinates_time{ nullptr };
    int * d_pair_indices_i{ nullptr };
    int * d_pair_indices_j{ nullptr };

    cuda_status = cudaMemcpyFromSymbol(&d_coordinates_x, static_d_coords_x, sizeof(d_coordinates_x), 0, cudaMemcpyDeviceToHost);
    cuda_status = cudaMemcpyFromSymbol(&d_coordinates_y, static_d_coords_y, sizeof(d_coordinates_y), 0, cudaMemcpyDeviceToHost);
    cuda_status = cudaMemcpyFromSymbol(&d_coordinates_time, static_d_coords_time, sizeof(d_coordinates_time), 0, cudaMemcpyDeviceToHost);
    cuda_status = cudaMemcpyFromSymbol(&d_pair_indices_i, static_d_pair_indices_i, sizeof(d_pair_indices_i), 0, cudaMemcpyDeviceToHost);
    cuda_status = cudaMemcpyFromSymbol(&d_pair_indices_j, static_d_pair_indices_j, sizeof(d_pair_indices_j), 0, cudaMemcpyDeviceToHost);

    cuda_status = cudaFree(d_coordinates_x);
    cuda_status = cudaFree(d_coordinates_y);
    cuda_status = cudaFree(d_coordinates_time);
    cuda_status = cudaFree(d_pair_indices_i);
    cuda_status = cudaFree(d_pair_indices_j);

    //cuda_status = cudaDeviceReset();

    return 0;
}




int gpu_opt_storm_drift_free_3d()
{

    cudaError_t cuda_status;

    // free device memory
    static_n_coordinates = 0;
    static_n_timepoints = 0;
    static_n_coordinate_pairs = 0;
    
    float * d_coordinates_x{ nullptr };
    float * d_coordinates_y{ nullptr };
    float * d_coordinates_z{ nullptr };
    int * d_coordinates_time{ nullptr };
    int * d_pair_indices_i{ nullptr };
    int * d_pair_indices_j{ nullptr };

    cuda_status = cudaMemcpyFromSymbol(&d_coordinates_x, static_d_coords_x, sizeof(d_coordinates_x), 0, cudaMemcpyDeviceToHost);
    cuda_status = cudaMemcpyFromSymbol(&d_coordinates_y, static_d_coords_y, sizeof(d_coordinates_y), 0, cudaMemcpyDeviceToHost);
    cuda_status = cudaMemcpyFromSymbol(&d_coordinates_z, static_d_coords_z, sizeof(d_coordinates_z), 0, cudaMemcpyDeviceToHost);
    cuda_status = cudaMemcpyFromSymbol(&d_coordinates_time, static_d_coords_time, sizeof(d_coordinates_time), 0, cudaMemcpyDeviceToHost);
    cuda_status = cudaMemcpyFromSymbol(&d_pair_indices_i, static_d_pair_indices_i, sizeof(d_pair_indices_i), 0, cudaMemcpyDeviceToHost);
    cuda_status = cudaMemcpyFromSymbol(&d_pair_indices_j, static_d_pair_indices_j, sizeof(d_pair_indices_j), 0, cudaMemcpyDeviceToHost);

    cuda_status = cudaFree(d_coordinates_x);
    cuda_status = cudaFree(d_coordinates_y);
    cuda_status = cudaFree(d_coordinates_z);
    cuda_status = cudaFree(d_coordinates_time);
    cuda_status = cudaFree(d_pair_indices_i);
    cuda_status = cudaFree(d_pair_indices_j);

    //cuda_status = cudaDeviceReset();

    return 0;
}


