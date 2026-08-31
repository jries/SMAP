;ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
;
;   dv_an_proc_optimize_drift_minimal
;                    
;   time_seg_def: 0: define according to the number of time segments
;                 1: define according to the number of localizations
;                 2: define according to specified equal time windows
;                 
;   n_segments_in: number of time segments into which the localization data 
;                  is divided (time_seg_def = 0)
;   
;   n_loc_per_segment_in: number of localizations per time segment (used to 
;                         define the number of segments when time_seg_def = 1)
;   
;   n_frm_per_segment_in: number of frames per time segment (used to 
;                         define the number of segments when time_seg_def = 2)
;
;ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


function dv_an_proc_optimize_drift_minimal_opt_func, drift_correction

    compile_opt idl2, strictarrsubs

    COMMON dv_an_proc_optimize_drift_minimal_cb, cb_n_loc, $
                                                   cb_n_timepoints, $
                                                   cb_n_pairs, $
                                                   cb_gauss_scale, $
                                                   cb_binary_dir, $
                                                   cb_progress_bar

    result = dv_gpu_osd_compute_3d(cb_n_loc, $
                                   cb_n_timepoints, $
                                   cb_n_pairs, $
                                   cb_gauss_scale, $
                                   drift_correction, $
                                   output_cost_function, $
                                   0, $
                                   output_derivatives, $
                                   binary_dir = cb_binary_dir)

    if obj_valid(cb_progress_bar) then begin
        
        labeltxt = 'Drift correction (optimization), cost function: ' $
            + string(output_cost_function)
        
        cb_progress_bar.label_text = labeltxt
        cb_progress_bar.Update_text, labeltxt
    
    endif
    
    print, output_cost_function
    
    return, output_cost_function

end


function dv_an_proc_optimize_drift_minimal_opt_func_derivative_gpu, drift_correction

    compile_opt idl2, strictarrsubs

    COMMON dv_an_proc_optimize_drift_minimal_cb, cb_n_loc, $
                                                   cb_n_timepoints, $
                                                   cb_n_pairs, $
                                                   cb_gauss_scale, $
                                                   cb_binary_dir, $
                                                   cb_progress_bar

    result = dv_gpu_osd_compute_3d(cb_n_loc, $
                                   cb_n_timepoints, $
                                   cb_n_pairs, $
                                   cb_gauss_scale, $
                                   drift_correction, $
                                   output_cost_function, $
                                   1, $
                                   output_derivatives, $
                                   binary_dir = cb_binary_dir)
    
    return, output_derivatives

end




function dv_an_proc_optimize_drift_minimal_LBFGS_opt_func, drift_correction, $
                                                           gradient_values = gradient_values

    compile_opt idl2, strictarrsubs

    COMMON dv_an_proc_optimize_drift_minimal_cb, cb_n_loc, $
                                                   cb_n_timepoints, $
                                                   cb_n_pairs, $
                                                   cb_gauss_scale, $
                                                   cb_binary_dir, $
                                                   cb_progress_bar

    result = dv_gpu_osd_compute_3d(cb_n_loc, $
                                   cb_n_timepoints, $
                                   cb_n_pairs, $
                                   cb_gauss_scale, $
                                   drift_correction, $
                                   output_cost_function, $
                                   1, $
                                   output_derivatives, $
                                   binary_dir = cb_binary_dir)

    if obj_valid(cb_progress_bar) then begin
        
        labeltxt = 'Drift correction (optimization), cost function: ' $
            + string(output_cost_function)
        
        cb_progress_bar.label_text = labeltxt
        cb_progress_bar.Update_text, labeltxt
    
    endif
    
    ;print, output_cost_function
    
    gradient_values = temporary(output_derivatives)
    return, output_cost_function

end



function dv_an_proc_optimize_drift_minimal, loc_x, $
                                            loc_y, $
                                            loc_z, $
                                            loc_t, $
                                            n_segments_in, $
                                            pixel_size_um, $
                                            initial_gaussian_scale_um, $
                                            max_drift_um, $                 
                                            optimize_gauss_scale = optimize_gauss_scale, $
                                            smoothing_width = smoothing_width, $
                                            linear_interpolation = linear_interpolation, $
                                            spline_interpolation = spline_interpolation, $
                                            lbfgsb_algorithm = lbfgsb_algorithm, $
                                            start_frame_out = start_frame_out, $
                                            center_frame_out = center_frame_out, $
                                            end_frame_out = end_frame_out, $
                                            displacements_out = displacements_out, $
                                            displacements_smoothed_out = displacements_smoothed_out, $
                                            displacements_unsmoothed_out = displacements_unsmoothed_out, $
                                            gaussian_scale_um_out = gaussian_scale_um_out, $
                                            skip_user_query = skip_user_query, $
                                            optimization_cancelled = optimization_cancelled, $
                                            binary_dir = binary_dir


    compile_opt idl2, strictarrsubs
    @dv_func_err_handler


    COMMON dv_an_proc_optimize_drift_minimal_cb, cb_n_loc, $
                                                 cb_n_timepoints, $
                                                 cb_n_pairs, $
                                                 cb_gauss_scale, $
                                                 cb_binary_dir, $
                                                 cb_progress_bar


    optimization_cancelled = 0
    scan_data_present = 0

    if N_elements(optimize_gauss_scale) eq 0 then optimize_gauss_scale = 0
    if N_elements(smoothing_width) eq 0 then smoothing_width = 1
    if N_elements(linear_interpolation) eq 0 then linear_interpolation = 0
    if N_elements(spline_interpolation) eq 0 then spline_interpolation = 0
    if N_elements(lbfgsb_algorithm) eq 0 then lbfgsb_algorithm = 0
    if N_elements(skip_user_query) eq 0 then skip_user_query = 0
    if N_elements(binary_dir) eq 0 then binary_dir = ''
      
      

    n_localizations = N_elements(loc_x)

    n_loc_per_segment = floor(n_localizations / float(n_segments_in))

    n_segments = n_segments_in
    
    nloc_str = strtrim(string(n_loc_per_segment),2)
    nseg_str = strtrim(string(n_segments),2)
    
    msgtxt = 'Drift correct with ' + nseg_str + ' time windows, ' + $
             nloc_str + ' localizations per window.  Proceed?'
             
             
    if skip_user_query eq 0 then begin
             
        result = DIALOG_MESSAGE(msgtxt, /QUESTION)
            
        if result eq 'No' then begin
            
            optimization_cancelled = 1
            return, 0

        endif
    
    endif
    
    
    ; sort the localizations by time
    
    sort_index = sort(loc_t)
    
    loc_x_sorted = loc_x[sort_index]
    loc_y_sorted = loc_y[sort_index]
    loc_z_sorted = loc_z[sort_index]
    loc_t_sorted = loc_t[sort_index]
    
    
    ; define segments
    
    loc_segment = loc_t_sorted
    loc_segment_mean_time = fltarr(n_segments)
    loc_segment_start_time = fltarr(n_segments)
    loc_segment_end_time = fltarr(n_segments)
    
    for i = 0, n_segments-1 do begin
        
        si = (i * n_loc_per_segment) < (n_localizations-1)
        ei = (si + n_loc_per_segment - 1) < (n_localizations-1)
        
        if i eq (n_segments-1) then ei = (n_localizations-1)
        
        loc_segment[si:ei] = i
        
        loc_segment_mean_time[i] = mean(loc_t_sorted[si:ei])
        
        loc_segment_start_time[i] = loc_t_sorted[si]
        loc_segment_end_time[i] = loc_t_sorted[ei]
        
    endfor

    localizations_x = loc_x_sorted
    localizations_y = loc_y_sorted
    localizations_z = loc_z_sorted
    localizations_segment = loc_segment

    rawdata_xy_pixsize = pixel_size_um
    
    rawdata_nframes = ceil(max(loc_t_sorted)) + 1
    
    tmp_sfrm = loc_segment_start_time
    tmp_center_frm = loc_segment_mean_time
    tmp_efrm = loc_segment_end_time
    
    segment_per_loc_array = loc_segment
    
    segment_per_frm_array = lonarr(rawdata_nframes)
    
    tmp_si = 0
    
    for i = 0, n_segments-1 do begin
        
        si = tmp_si
        ei = ceil(loc_segment_end_time[i])
        
        segment_per_frm_array[si:ei] = i
        
        tmp_si = ei + 1
        
    endfor
    
    
    gaussian_scale = initial_gaussian_scale_um / rawdata_xy_pixsize
    max_drift_scale = max_drift_um / rawdata_xy_pixsize
    
    hist_binsize = max_drift_scale
    
    x_loc_min = min(localizations_x, max=x_loc_max)
    y_loc_min = min(localizations_y, max=y_loc_max)
    z_loc_min = min(localizations_z, max=z_loc_max)
    
    nbins_x = ceil((x_loc_max-x_loc_min)/hist_binsize) > 1
    nbins_y = ceil((y_loc_max-y_loc_min)/hist_binsize) > 1
    nbins_z = ceil((z_loc_max-z_loc_min)/hist_binsize) > 1
    
    hist_x_max = x_loc_min + nbins_x * hist_binsize
    hist_y_max = y_loc_min + nbins_y * hist_binsize
    hist_z_max = z_loc_min + nbins_z * hist_binsize
    
    
    n_loc = N_elements(localizations_x)
    tmpdat = fltarr(3,n_loc,/NOZERO)
    tmpdat[0,*] = localizations_x
    tmpdat[1,*] = localizations_y
    tmpdat[2,*] = localizations_z
    
    
    ; first pair list
    
    hist_a_min = [x_loc_min, y_loc_min, z_loc_min]
    hist_a_max = [hist_x_max, hist_y_max, hist_z_max]
    hist_a_nbins = [nbins_x, nbins_y, nbins_z]
    
    h = hist_nd(tmpdat, MIN=hist_a_min, MAX=hist_a_max, NBINS=hist_a_nbins, REVERSE_INDICES=ri)
    
    pair_indices_i_ptr_arr = ptrarr(nbins_x,nbins_y,nbins_z)
    pair_indices_j_ptr_arr = ptrarr(nbins_x,nbins_y,nbins_z)
    bin_n_pairs_arr = lonarr(nbins_x,nbins_y,nbins_z)
    
    for i = 0, nbins_x-1 do begin
        for j = 0, nbins_y-1 do begin
            for k = 0, nbins_z-1 do begin
                
                ind=[i+nbins_x*(j+nbins_y*k)]
                si = ri[ind]
                ei = ri[ind+1]-1
                
                if ei gt si then begin
                
                    ; there are at least two localizations
                
                    tmp_indices = ri[si:ei]

                    bin_n_loc = N_elements(tmp_indices)

                    bin_n_pairs = (bin_n_loc * (bin_n_loc-1)) / 2

                    ; calculate the pairs within the current bin

                    bin_pair_indices_i = lonarr(bin_n_pairs,/NOZERO)
                    bin_pair_indices_j = lonarr(bin_n_pairs,/NOZERO)
    
                    tmp_count = 0
                    for m = 0, bin_n_loc-2 do begin
                        
                        tmp_n_entries = (bin_n_loc - m) - 1                        
                        bin_pair_indices_i[tmp_count] = replicate(tmp_indices[m], tmp_n_entries)
                        bin_pair_indices_j[tmp_count] = tmp_indices[lindgen(tmp_n_entries) + (m+1)]
                        tmp_count += tmp_n_entries

                    endfor
                        
                    pair_indices_i_ptr_arr[i,j,k] = ptr_new(bin_pair_indices_i,/NO_COPY)
                    pair_indices_j_ptr_arr[i,j,k] = ptr_new(bin_pair_indices_j,/NO_COPY)
                    bin_n_pairs_arr[i,j,k] = bin_n_pairs

                endif

            endfor
        endfor
        

        ;print, float(i)/(nbins_x-1)

    endfor
    
    
    n_pairs = total(bin_n_pairs_arr,/INTEGER)
    
    pair_indices_i = lonarr(n_pairs,/NOZERO)
    pair_indices_j = lonarr(n_pairs,/NOZERO)
    
    tmp_count = 0
    
    for i = 0, nbins_x-1 do begin
        for j = 0, nbins_y-1 do begin
            for k = 0, nbins_z-1 do begin
                
                tmp_n_pairs = bin_n_pairs_arr[i,j,k]
                
                if tmp_n_pairs gt 0 then begin
                    
                    pair_indices_i[tmp_count] = *pair_indices_i_ptr_arr[i,j,k]
                    pair_indices_j[tmp_count] = *pair_indices_j_ptr_arr[i,j,k]
                    tmp_count += tmp_n_pairs
                    
                endif

            endfor
        endfor
    endfor
    
    ptr_free, pair_indices_i_ptr_arr
    ptr_free, pair_indices_j_ptr_arr
    

    cb_n_loc = n_loc
    cb_n_timepoints = n_segments
    cb_n_pairs = n_pairs
    cb_gauss_scale = gaussian_scale
    cb_binary_dir = binary_dir
    cb_progress_bar = obj_new()

    
    ; initialize the GPU functions
    result = dv_gpu_osd_initialize_3d(n_segments, $
                                      localizations_x, $
                                      localizations_y, $
                                      localizations_z, $
                                      localizations_segment, $
                                      pair_indices_i, $
                                      pair_indices_j, $
                                      binary_dir = binary_dir)
    
    
    ; do the optimization
    
    gtol = 1.0e-7
    
    current_drift_correction = fltarr(3*n_segments)
    stored_gauss_scale = cb_gauss_scale
    
    
    ; catch convergence errors during the optimization
    
    convergence_error_count = 0
    optimization_aborted = 0
    
    CATCH, error_status
    
    if error_status ne 0 then begin
        
        PRINT, 'Error index: ', error_status
        PRINT, 'Error message: ', !ERROR_STATE.MSG
        convergence_error_count += 1
        
        if convergence_error_count eq 2 then begin
            
            ; after two errors, restart the optimization with 
            ; a larger gaussian scale value
            current_drift_correction = fltarr(3*n_segments)
            cb_gauss_scale = cb_gauss_scale * 2
            stored_gauss_scale = cb_gauss_scale
       
        endif
       
        if convergence_error_count gt 5 then begin
            
            ; after five errors, cancel the opimization
            CATCH, /CANCEL
            PRINT, 'Aborting optimization'
            optimization_aborted = 1
            GOTO, ABORT_OPTIMIZATION
        
        endif
        
    endif
    
    
    ; set parameters for the LBFGS algorithm
    
    optimization_function = 'dv_an_proc_optimize_drift_minimal_LBFGS_opt_func'

    bound_types = lonarr(3*n_segments)
    lower_bounds = dblarr(3*n_segments)
    upper_bounds = dblarr(3*n_segments)
    
    bound_types[*] = 2
    lower_bounds[*] = -abs(2*max_drift_scale)
    upper_bounds[*] = abs(2*max_drift_scale)
    
    input_factr = double(1e3)
    
    
    
    if lbfgsb_algorithm eq 0 then begin
    
        ; dfpmin (Numerical Recipes, section 10.9)
    
        dfpmin, current_drift_correction, $
                gtol, $
                fmin, $
                'dv_an_proc_optimize_drift_minimal_opt_func', $
                'dv_an_proc_optimize_drift_minimal_opt_func_derivative_gpu', $
                ITER = n_iterations_out, $
                ITMAX = 100
                
        print, 'Optimized cost function value: ', fmin
                
    endif else begin

        ; LBFGSB

        verbose_output = 1

        dv_lbfgsb_min, optimization_function, $
                       current_drift_correction, $
                       bound_types = bound_types, $
                       lower_bounds = lower_bounds, $
                       upper_bounds = upper_bounds, $
                       input_factr = input_factr, $
                       output_converged = output_converged, $
                       output_stopped = output_stopped, $
                       output_error = output_error, $
                       output_message = output_message, $
                       output_state = output_state, $
                       output_n_iter = n_iterations_out, $
                       output_fcn_value = output_fcn_value, $
                       verbose = verbose_output, $
                       binary_dir = binary_dir
                       
        print, 'Optimized cost function value: ', output_fcn_value
                       
    endelse

    print, 'initial gauss_scale ', cb_gauss_scale
    print, 'initial n_iterations ', n_iterations_out


    if optimize_gauss_scale eq 1 then begin

        ; divide the gaussian scale by a constant factor and repeat the process
        ; until the drift trajectory is no longer changing
        
        quita = 0
        
        last_drift_rms_change = !values.f_infinity
        
        while (quita eq 0) do begin
        
            last_drift_correction = current_drift_correction
            stored_gauss_scale = cb_gauss_scale
        
            cb_gauss_scale = cb_gauss_scale / 1.5
            
                    
            if lbfgsb_algorithm eq 0 then begin
            
                ; dfpmin
            
                dfpmin, current_drift_correction, $
                        gtol, $
                        fmin, $
                        'dv_an_proc_optimize_drift_minimal_opt_func', $
                        'dv_an_proc_optimize_drift_minimal_opt_func_derivative_gpu', $
                        ITER = n_iterations_out
                        
                print, 'Optimized cost function value: ', fmin
                        
            endif else begin
                    
                ; LBFGSB
                    
                dv_lbfgsb_min, optimization_function, $
                               current_drift_correction, $
                               bound_types = bound_types, $
                               lower_bounds = lower_bounds, $
                               upper_bounds = upper_bounds, $
                               input_factr = input_factr, $
                               output_converged = output_converged, $
                               output_stopped = output_stopped, $
                               output_error = output_error, $
                               output_message = output_message, $
                               output_state = output_state, $
                               output_n_iter = n_iterations_out, $
                               output_fcn_value = output_fcn_value, $
                               verbose = 1, $
                               binary_dir = binary_dir
                               
                print, 'Optimized cost function value: ', output_fcn_value
                               
            endelse
                    
                    
            drift_rms_change = sqrt(mean((last_drift_correction - current_drift_correction)^2))
            
            if drift_rms_change gt last_drift_rms_change then quita = 1
    
            last_drift_rms_change = drift_rms_change
    
            print, 'gauss_scale (um)', cb_gauss_scale * rawdata_xy_pixsize
            print, 'n_iterations ', n_iterations_out
            print, 'drift_rms_change ', drift_rms_change * rawdata_xy_pixsize
    
        endwhile
        
        ; discard the last iteration of the drift correction optimization
        current_drift_correction = last_drift_correction

    endif

ABORT_OPTIMIZATION:

    ; free the GPU functions
    result = dv_gpu_osd_free_3d(binary_dir = binary_dir)

    ; reset the common block variables
    cb_n_loc = 0
    cb_n_timepoints = 0
    cb_n_pairs = 0
    cb_gauss_scale = 0.0
    cb_binary_dir = ''
    cb_progress_bar = obj_new()


    x_drift_opt = -current_drift_correction[0:(n_segments-1)]
    y_drift_opt = -current_drift_correction[n_segments:((2*n_segments)-1)]
    z_drift_opt = -current_drift_correction[(2*n_segments):((3*n_segments)-1)]

;    ; plot the drift data
;    ydata = fltarr(n_segments,3)
;    ydata[0,0] = x_drift_opt
;    ydata[0,1] = y_drift_opt
;    ydata[0,2] = z_drift_opt
;    xdata = center_frame_out[lindgen(n_segments)]
;    result = daxview_plot_xy_data(xdata,ydata)


    raw_drift_x = x_drift_opt[segment_per_frm_array]
    raw_drift_y = y_drift_opt[segment_per_frm_array]
    raw_drift_z = z_drift_opt[segment_per_frm_array]


    ; smooth the data with a boxcar averaging filter
    
    smoothed_x = smooth(x_drift_opt, smoothing_width, /EDGE_TRUNCATE)
    smoothed_y = smooth(y_drift_opt, smoothing_width, /EDGE_TRUNCATE)
    smoothed_z = smooth(z_drift_opt, smoothing_width, /EDGE_TRUNCATE)

    if linear_interpolation eq 1 or spline_interpolation eq 1 then begin
    
        spl_interp_mod = (spline_interpolation eq 1) and (N_elements(smoothed_x) ge 4)

        smoothed_x = interpol(smoothed_x, tmp_center_frm, lindgen(rawdata_nframes), SPLINE=spl_interp_mod)
        smoothed_y = interpol(smoothed_y, tmp_center_frm, lindgen(rawdata_nframes), SPLINE=spl_interp_mod)
        smoothed_z = interpol(smoothed_z, tmp_center_frm, lindgen(rawdata_nframes), SPLINE=spl_interp_mod)
 
    endif else begin
        
        smoothed_x = smoothed_x[segment_per_frm_array]
        smoothed_y = smoothed_y[segment_per_frm_array]
        smoothed_z = smoothed_z[segment_per_frm_array]
        
    endelse
 
                  
    displacements_smoothed_out = fltarr(3,rawdata_nframes,/NOZERO)

    displacements_smoothed_out[0,*] = smoothed_x
    displacements_smoothed_out[1,*] = smoothed_y
    displacements_smoothed_out[2,*] = smoothed_z
    
    start_frame_out = tmp_sfrm
    center_frame_out = tmp_center_frm
    end_frame_out = tmp_efrm

    displacements_out = fltarr(3,n_segments,/NOZERO)    
    displacements_out[0,*] = x_drift_opt
    displacements_out[1,*] = y_drift_opt
    displacements_out[2,*] = z_drift_opt

    displacements_unsmoothed_out = fltarr(3,rawdata_nframes,/NOZERO)
    displacements_unsmoothed_out[0,*] = raw_drift_x
    displacements_unsmoothed_out[1,*] = raw_drift_y
    displacements_unsmoothed_out[2,*] = raw_drift_z

    gaussian_scale_um_out = stored_gauss_scale * rawdata_xy_pixsize

    ; evaluate the result
    
    xd_min = min(x_drift_opt, max=xd_max)
    yd_min = min(y_drift_opt, max=yd_max)
    zd_min = min(z_drift_opt, max=zd_max)

    xd_scale_um = (xd_max-xd_min)*rawdata_xy_pixsize
    yd_scale_um = (yd_max-yd_min)*rawdata_xy_pixsize
    zd_scale_um = (zd_max-zd_min)*rawdata_xy_pixsize

    optimization_succeeded = optimization_aborted eq 0 and $
                             xd_scale_um le max_drift_um and $
                             yd_scale_um le max_drift_um and $
                             zd_scale_um le max_drift_um

    return, optimization_succeeded
    
end
    
    