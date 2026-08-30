;cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
;
;   dv_gpu_osd_compute_3d
;   
;cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

function dv_gpu_osd_compute_3d, n_coordinates, $
                                n_timepoints, $
                                n_coordinate_pairs, $
                                gaussian_scale, $
                                drift_trajectory, $
                                output_cost_function, $
                                flag_calculate_derivatives, $
                                output_derivatives, $
                                binary_dir = binary_dir, $
                                execution_time = execution_time
                      

    compile_opt idl2, strictarrsubs
    @dv_func_err_handler

    if N_elements(binary_dir) eq 0 then $
        message, 'Binary file directory not specified'

    tmp_timer = tic()

    library_name  = binary_dir + 'gpu_optimize_storm_drift.dll'
    function_name = 'gpu_opt_storm_drift_compute_3d_portable'
        

    if size(drift_trajectory,/TYPE) ne 4 then message, 'Invalid drift trajectory type'

    input_n_coordinates = long(n_coordinates)
    input_n_timepoints = long(n_timepoints)
    input_n_pairs = ulong64(n_coordinate_pairs)
    input_gauss_scale = float(gaussian_scale)
    input_flag_calculate_derivatives = long(flag_calculate_derivatives)
    

    output_cost_function = float(0.0)
    
    if flag_calculate_derivatives eq 1 then begin
    
        output_derivatives = fltarr(3*n_timepoints,/NOZERO)
             
    endif else begin
        
        output_derivatives = 0.0D
        
    endelse
             
    ; call the dll
            
    tmp =  call_external(library_name, $
                         function_name, $
                         input_n_coordinates, $
                         input_n_timepoints, $
                         input_n_pairs, $
                         input_gauss_scale, $
                         drift_trajectory, $
                         output_cost_function, $
                         input_flag_calculate_derivatives, $
                         output_derivatives, $
                         RETURN_TYPE = 3, $
                         /VERBOSE)
    

    if tmp ne 0 then message, 'Error code ' + strtrim(string(tmp),2)

    if flag_calculate_derivatives eq 0 then output_derivatives = []

    execution_time = toc(tmp_timer)

    return, tmp
    
end

