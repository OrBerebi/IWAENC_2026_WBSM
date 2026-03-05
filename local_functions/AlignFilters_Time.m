function alighned_signal = AlignFilters_Time(ref,signal)
%TIME_DOMAIN_ALIGHMENT Summary of this function goes here
    N = length(ref(:,1));
    tau_min = -1000;
    tau_max = 999;
    tau_vec = tau_min:tau_max;
    cov_val = zeros(1,length(tau_vec));
    ear_ref = ref';
    ear_signal = signal';

    % Find Corellation of delay
    for tau=tau_vec
        tmp_signal = circshift(ear_signal,tau,2);
        tmp_signal = tmp_signal(:,1:N);
        cov_val(tau-tau_min+1) = trace(ear_ref*tmp_signal');
    end
    
    % Align to the max correlation delay index
    [~,index] = max(cov_val);
    tau_final = index+tau_min-1;
    ear_signal = circshift(ear_signal,tau_final,2);
    alighned_signal = ear_signal(:,1:N)';
    a = (trace(alighned_signal'*ref)/trace(alighned_signal'*alighned_signal));
    alighned_signal = a*alighned_signal;
end

