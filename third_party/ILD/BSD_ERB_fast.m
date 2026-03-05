function BSD = BSD_ERB_fast(x_ref, x_test, ERBstruct)

    fs = ERBstruct.fs;
    f_lim = ERBstruct.f_lim;
    f_c = ERBstruct.f_c;
    C = ERBstruct.C;
    n = ERBstruct.n;

    % match IR length
    if size(x_ref,1)~=size(x_test,1)
        if size(x_ref,1) > size(x_test,1)
            x_test(end+1:size(x_ref,1),:) = 0;
        else
            x_ref(end+1:size(x_test,1),:) = 0;
        end
    end

    % match IR channels
    if size(x_ref,2) > size(x_test,2)
        x_test = repmat(x_test, [1 size(x_ref,2)]);
    end

    N = size(x_ref,1);
    c = size(x_ref,2);

    % allocate
    X_out_ref = zeros(n,c);
    X_out_test = zeros(n,c);

    % filterbank magnitude responses (in dB) for cutoff detection
    C_db_tmp = db(C);

    % Loop over filters
    for k = 1:n
        % determine needed length
        C_db = C_db_tmp(:,k);
        Nmin = find(flipud(C_db) > max(C_db)-60, 1, 'first');
        Nmin = 2^nextpow2(size(C,1) - Nmin + 1);

        Ncur = 2^nextpow2(1/f_c(k) * 4 * fs);
        Ncur = max([Ncur Nmin N]);

        % FFT frequency indices
        f_lim_id(1) = max([round(f_lim(1)/(fs/Ncur))+1 2]);
        f_lim_id(2) = round(f_lim(2)/(fs/Ncur))+1;

        % filter signals
        Hk = abs(fft(C(:,k), Ncur));
        X_ref = repmat(Hk, [1, c]) .* (abs(fft(x_ref, Ncur)).^2);
        X_test = repmat(Hk, [1, c]) .* (abs(fft(x_test, Ncur)).^2);

        % band energy
        X_out_ref(k,:)  = sum(X_ref(f_lim_id(1):f_lim_id(2),:));
        X_out_test(k,:) = sum(X_test(f_lim_id(1):f_lim_id(2),:));
    end

    % ERB error
    BSD = abs(10*log10(abs(X_out_test)./abs(X_out_ref+eps)));

end
