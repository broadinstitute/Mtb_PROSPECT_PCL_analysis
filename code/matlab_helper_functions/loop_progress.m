function loop_progress(ii,tot,interval)

% Fuction to print loop progress

if ii==1 || mod(ii,interval)==0
        disp(sprintf('%d/%d',ii,tot))
end