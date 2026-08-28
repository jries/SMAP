if isfield(results, 'Single_tracking_analysis')
    stvar = results.Single_tracking_analysis;
    track_summary=[];
    track_details=[];
    for k=1:length(stvar)
        track_summary=vertcat(track_summary,stvar{k}.summarytable);
        track_details=vertcat(track_details,stvar{k}.trackstable);
    end
end
if isfield(results, 'Dual_tracking_analysis')
    dtvar = results.Dual_tracking_analysis;
    dual_track_summary=[];
    for k=1:length(dtvar)
        dual_track_summary=vertcat(dual_track_summary,dtvar{k});
    end
end
if isfield(results, 'Locstatistics')
    dtvar = results.Locstatistics;
    statistics=[];
    for k=1:length(dtvar)
        statistics=vertcat(statistics,dtvar{k});
    end
end