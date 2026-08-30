#include "gpu_optimize_storm_drift.h"


int gpu_opt_storm_drift_initialize_2d_portable(int argc, void *argv[])
{

    return gpu_opt_storm_drift_initialize_2d(
        *(int *)argv[0],
        *(int *)argv[1],
        (float *)argv[2],
        (float *)argv[3],
        (int *)argv[5],
        *(size_t *)argv[6],
        (int *)argv[7],
        (int *)argv[8]);

}

int gpu_opt_storm_drift_compute_2d_portable(int argc, void *argv[])
{

    return gpu_opt_storm_drift_compute_2d(
        *(int *)argv[0],
        *(int *)argv[1],
        *(size_t *)argv[2],
        *(float *)argv[3],
        (float *)argv[4],
        (float *)argv[5],
        *(int *)argv[6],
        (float *)argv[7]);

}

int gpu_opt_storm_drift_free_2d_portable(int argc, void *argv[])
{

    return gpu_opt_storm_drift_free_2d();

}

int gpu_opt_storm_drift_initialize_3d_portable(int argc, void *argv[])
{

	return gpu_opt_storm_drift_initialize_3d(
		*(int *)argv[0],
        *(int *)argv[1],
        (float *)argv[2],
        (float *)argv[3],
        (float *)argv[4],
        (int *)argv[5],
        *(size_t *)argv[6],
        (int *)argv[7],
        (int *)argv[8]);

}

int gpu_opt_storm_drift_compute_3d_portable(int argc, void *argv[])
{

    return gpu_opt_storm_drift_compute_3d(
        *(int *)argv[0],
        *(int *)argv[1],
        *(size_t *)argv[2],
        *(float *)argv[3],
        (float *)argv[4],
        (float *)argv[5],
        *(int *)argv[6],
        (float *)argv[7]);

}

int gpu_opt_storm_drift_free_3d_portable(int argc, void *argv[])
{

    return gpu_opt_storm_drift_free_3d();

}
