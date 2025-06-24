function analog_data = to_analog(time_data, N, total_time)

    time_res = total_time/N;
    
    analog_data = zeros(1, N);
        
    for ti = 1: length(time_data)
         idx = floor(mod(time_data(ti)/time_res, N)) + 1;
         analog_data(idx) = 1;
    end
end