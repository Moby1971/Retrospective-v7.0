function footer = makemrdfooter(inputfooter,par)


parameters = {':NO_SAMPLES no_samples, ',':NO_VIEWS no_views, ',':NO_VIEWS_2 no_views_2, ', ... 
                ':NO_ECHOES no_echoes, ',':EXPERIMENT_ARRAY no_experiments, ',':NO_AVERAGES no_averages, ', ...
                ':VAR pe1_order, ',':VAR slice_nav, ',':VAR radial_on, ', ... 
                ':VAR frame_loop_on, ',':VAR tr, ',':VAR te, ', ...
                ':BATCH_SLICES batch_slices, ',':NO_SLICES no_slices, '
                };

replacepars = {par.NoSamples,par.NoViews,par.NoViews2, ... 
                par.NoEchoes,par.NoExperiments,par.NoAverages, ... 
                par.peorder,par.slicenav,par.radialon, ... 
                par.frameloopon,par.tr,par.te, ...
                par.batchslices,par.NoSlices
                };


for i = 1:length(parameters)
    
    txt = parameters{i};
    var = replacepars{i};
    
    pos = strfind(inputfooter,txt);
    
    if ~isempty(pos)
        oldtxtlength = strfind(inputfooter(pos+length(txt):pos+length(txt)+6),char(13))-1;
        newtext = [num2str(var),'     '];
        newtext = newtext(1:6);
        inputfooter = replaceBetween(inputfooter,pos+length(txt),pos+length(txt)+oldtxtlength-1,newtext);
    end
    
end

footer = inputfooter;


end