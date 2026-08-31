
;cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
;
;   dv_gpu_osd_initialize_3d
;   
;cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

function dv_gpu_osd_initialize_3d, n_timepoints, $
                                   coordinates_x, $
                                   coordinates_y, $
                                   coordinates_z, $
                                   coordinates_t, $
                                   pair_indices_i, $
                                   pair_indices_j, $
                                   binary_dir = binary_dir, $
                                   execution_time = execution_time
                      

    compile_opt idl2, strictarrsubs
    @dv_func_err_handler

    if N_elements(binary_dir) eq 0 then $
        message, 'Binary file directory not specified'

    tmp_timer = tic()

    library_name  = binary_dir + 'gpu_optimize_storm_drift.dll'
    function_name = 'gpu_opt_storm_drift_initialize_3d_portable'
        
    n_coordinates = N_elements(coordinates_x)
    n_pairs = N_elements(pair_indices_i)

    input_coordinates_x = float(coordinates_x)
    input_coordinates_y = float(coordinates_y)
    input_coordinates_z = float(coordinates_z)
    input_coordinates_t = long(coordinates_t)
    input_pair_indices_i = long(pair_indices_i)
    input_pair_indices_j = long(pair_indices_j)
    
    input_n_coordinates = long(n_coordinates)
    input_n_timepoints = long(n_timepoints)
    input_n_pairs = ulong64(n_pairs)
   
             
    ; call the dll
            
    tmp =  call_external(library_name, $
                         function_name, $
                         input_n_coordinates, $
                         input_n_timepoints, $
                         input_coordinates_x, $
                         input_coordinates_y, $
                         input_coordinates_z, $
                         input_coordinates_t, $
                         input_n_pairs, $
                         input_pair_indices_i, $
                         input_pair_indices_j, $
                         RETURN_TYPE = 3, $
                         /VERBOSE)
    

    if tmp ne 0 then message, 'Error code ' + strtrim(string(tmp),2)

    execution_time = toc(tmp_timer)

    return, tmp
    
end

