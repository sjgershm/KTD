function x = KTD_make_stimulus(onset,dur,trial_length)
    
    x = zeros(trial_length,1);
    if ~isempty(onset)
        x(onset:(onset+dur-1)) = 1;
    end