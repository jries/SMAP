;cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
;
;   dv_gpu_osd_free_3d
;   
;cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

function dv_gpu_osd_free_3d, binary_dir = binary_dir, $
                             execution_time = execution_time
                      

    compile_opt idl2, strictarrsubs
    @dv_func_err_handler

    if N_elements(binary_dir) eq 0 then $
        message, 'Binary file directory not specified'

    tmp_timer = tic()

    library_name  = binary_dir + 'gpu_optimize_storm_drift.dll'
    function_name = 'gpu_opt_storm_drift_free_3d_portable'
        
             
    ; call the dll
            
    tmp =  call_external(library_name, $
                         function_name, $
                         RETURN_TYPE = 3, $
                         /VERBOSE)
    

    if tmp ne 0 then message, 'Error code ' + strtrim(string(tmp),2)

    execution_time = toc(tmp_timer)

    return, tmp
    
end

