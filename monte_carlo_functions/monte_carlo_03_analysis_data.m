function monte_carlo_03_analysis_data(base_path, results_dir, results_txt_file, condition_name, plot_flag)

    conf.fs_proc = 16e3;
    conf.f_low_limit  = 20;
    conf.f_mid_limit  = 1500; 
    conf.f_high_limit = conf.fs_proc/2;

    f_band = [conf.f_mid_limit,conf.f_high_limit];
    ERBstruct = prepare_ERB(f_band, conf.fs_proc);
    ERBstruct.fs = conf.fs_proc;
    ERBstruct.f_band = f_band;

    if ~exist(results_dir, 'dir'), error('Results directory not found: %s', results_dir); end
    file_list = dir(fullfile(results_dir, 'eval_scene_*.mat'));
    N_files = length(file_list);

    if N_files == 0
        warning('No result files found in %s', results_dir);
        return;
    end

    results_mix  = [];
    freq_data = struct();

    for i = 1:N_files
        fname = fullfile(file_list(i).folder, file_list(i).name);
        
        % --- TRY-CATCH BLOCK FOR THE ENTIRE FILE ---
        try
            tmp = load(fname, 'eval_data');
            data = tmp.eval_data;
            
            fs = data.fs;
            mix_ref = data.SpeechMixed.ref;
            algos = fieldnames(data.SpeechMixed);
            algos(strcmp(algos, 'ref')) = [];
            
            [bms_ref, ild_ref, ic_ref, ~] = compute_binaural_metrics(mix_ref, fs, conf);
            new_row = struct('ID', data.ID);
            
            for a = 1:length(algos)
                algo_name = algos{a};
                sig_test = data.SpeechMixed.(algo_name);
                
                % Check for finiteness before entering sub-functions
                if ~all(isfinite(sig_test(:))) || ~all(isfinite(mix_ref(:)))
                    warning('Skipping %s in file %s: Non-finite values detected.', algo_name, file_list(i).name);
                    continue; 
                end

                % Metrics calculation
                [m_algo, m_algo_f] = calculate_all_metrics(mix_ref, sig_test, fs, conf, ERBstruct);
                [m_algo_stft, m_algo_stft_f] = calculate_all_metrics_stft(mix_ref, sig_test, fs, conf, ERBstruct);
                [bms_test, ild_test, ic_test, f_erb] = compute_binaural_metrics(sig_test, fs, conf);
                bin_errs = compute_binaural_errors(bms_ref, ild_ref, ic_ref, bms_test, ild_test, ic_test);
                
                % --- ACCUMULATE FREQUENCY DATA ---
                % (Keep your existing accumulation logic here...)
                if ~isfield(freq_data, algo_name)
                    freq_data.(algo_name).nmse = []; freq_data.(algo_name).stft_nmse = [];
                    freq_data.(algo_name).mag = [];  freq_data.(algo_name).stft_mag = [];
                    freq_data.(algo_name).lsd = [];  freq_data.(algo_name).stft_lsd = [];
                    freq_data.(algo_name).bsd = [];  freq_data.(algo_name).stft_bsd = [];
                    freq_data.(algo_name).bms_sq = []; freq_data.(algo_name).ild_sq = []; freq_data.(algo_name).ic_sq = [];
                end
                freq_data.(algo_name).nmse(:, end+1) = m_algo_f.nmse;
                freq_data.(algo_name).stft_nmse(:, end+1) = m_algo_stft_f.nmse;
                freq_data.(algo_name).mag(:, end+1) = m_algo_f.mag;
                freq_data.(algo_name).stft_mag(:, end+1) = m_algo_stft_f.mag;
                freq_data.(algo_name).lsd(:, end+1) = m_algo_f.lsd;
                freq_data.(algo_name).stft_lsd(:, end+1) = m_algo_stft_f.lsd;
                freq_data.(algo_name).bsd(:, end+1) = m_algo_f.bsd;
                freq_data.(algo_name).stft_bsd(:, end+1) = m_algo_stft_f.bsd;
                freq_data.(algo_name).bms_sq(:, end+1) = (bms_ref - bms_test).^2;
                freq_data.(algo_name).ild_sq(:, end+1) = (ild_ref - ild_test).^2;
                freq_data.(algo_name).ic_sq(:, end+1)  = (ic_ref - ic_test).^2;

                if ~isfield(freq_data, 'axes')
                    freq_data.axes.nmse = m_algo_f.f_nmse; freq_data.axes.stft_nmse = m_algo_stft_f.f_nmse;
                    freq_data.axes.mag  = m_algo_f.f_mag;  freq_data.axes.stft_mag  = m_algo_stft_f.f_mag;
                    freq_data.axes.lsd  = m_algo_f.f_lsd;  freq_data.axes.stft_lsd  = m_algo_stft_f.f_lsd;
                    freq_data.axes.bsd  = ERBstruct.f_c;   freq_data.axes.erb       = f_erb;
                end

                new_row.(sprintf('stft_NMSE_%s', algo_name)) = m_algo_stft.nmse;
                new_row.(sprintf('stft_MagErr_%s', algo_name)) = m_algo_stft.mag;
                new_row.(sprintf('stft_LSD_%s', algo_name))   = m_algo_stft.lsd;
                new_row.(sprintf('stft_BSD_%s', algo_name))   = m_algo_stft.bsd;
                new_row.(sprintf('RMSE_BMS_%s', algo_name)) = bin_errs.bms_rmse;
                new_row.(sprintf('RMSE_ILD_%s', algo_name)) = bin_errs.ild_rmse;
                new_row.(sprintf('RMSE_IC_%s', algo_name))  = bin_errs.ic_rmse;
            end
            
            if isempty(results_mix), results_mix = new_row; else, results_mix(end+1) = new_row; end

        catch ME
            % If ANY error occurs (pwelch, memory, dimensions), print a warning and skip
            warning('CRITICAL ERROR in file %s. Skipping this example. Error: %s', ...
                file_list(i).name, ME.message);
            continue; 
        end
    end

    T_mix  = struct2table(results_mix);

    % Open text file for appending
    fid = fopen(results_txt_file, 'a');
    fprintf(fid, '\n======================================================\n');
    fprintf(fid, 'CONDITION: %s\n', condition_name);
    fprintf(fid, '======================================================\n');
    
    print_comparison_table_to_file(fid, 'MIXED SIGNALS', T_mix);
    fclose(fid);
    
    if plot_flag
        freq_data.conf = conf;
        plot_frequency_metrics(freq_data, algos);
    end

end

function print_comparison_table_to_file(fid, title_str, T)
    fprintf(fid, '\n--- %s ---\n', title_str);
    var_names = T.Properties.VariableNames;
    algo_cols = var_names(startsWith(var_names, 'stft_NMSE_'));
    algos = strrep(algo_cols, 'stft_NMSE_', '');
    
    metrics = {'stft_NMSE', 'stft_MagErr', 'stft_LSD','stft_BSD', ...
               'RMSE_BMS', 'RMSE_ILD', 'RMSE_IC'};
    
    header_str = sprintf('%-15s | ', 'Metric');
    for a = 1:length(algos)
        header_str = [header_str, sprintf('%-20s | ', algos{a})];
    end
    fprintf(fid, '%s\n', header_str);
    fprintf(fid, '%s\n', repmat('-', 1, length(header_str)));
    
    for i = 1:length(metrics)
        m = metrics{i};
        row_str = sprintf('%-15s | ', m);
        for a = 1:length(algos)
            col_name = sprintf('%s_%s', m, algos{a});
            if any(strcmp(var_names, col_name))
                vals = T.(col_name);
                val_str = sprintf('%6.2f +/- %4.2f', mean(vals), std(vals));
            else
                val_str = sprintf('%15s', 'N/A');
            end
            row_str = [row_str, sprintf('%-20s | ', val_str)];
        end
        fprintf(fid, '%s\n', row_str);
    end
    fprintf(fid, '%s\n\n', repmat('-', 1, length(header_str)));



    
end

% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================

function sig_out = MatchEnergy(sig_ref, sig_in, fs, f_cutoff)
    % Matches energy based on LOW FREQUENCIES ONLY (< f_cutoff)
    f_cutoff_low = f_cutoff(1);
    f_cutoff_high = f_cutoff(2);
    % 1. Compute FFT
    nfft = 2^nextpow2(size(sig_ref,1));
    X_ref = fft(sig_ref, nfft);
    X_in  = fft(sig_in, nfft);
    
    % 2. Define Frequency Axis
    f_axis = linspace(0, fs, nfft)';
    
    % 3. Select Low Freq Bins
    [~,idx_low] = min(abs(f_axis - f_cutoff_low));
    [~,idx_high] = min(abs(f_axis - f_cutoff_high));
    
    % 4. Compute Energy in Low Band
    E_ref_low = sum(abs(X_ref(idx_low:idx_high, :)).^2, 1);
    E_in_low  = sum(abs(X_in(idx_low:idx_high, :)).^2, 1);
    
    % 5. Calculate Gain Factor
    gains = mean(sqrt(E_ref_low ./ (E_in_low + 1e-12)));
    
    % 6. Apply Gain to Full Time Signal
    sig_out = sig_in .* gains;
end

function [err_spec, mag_err_spec, ref_spec] = get_spectral_energy(sig_ref, sig_test, nfft)
    if size(sig_ref,1) < size(sig_ref,2), sig_ref = sig_ref.'; end
    if size(sig_test,1) < size(sig_test,2), sig_test = sig_test.'; end
    
    win = hann(size(sig_ref,1));
    S_ref = fft(sig_ref .* win, nfft);
    S_test = fft(sig_test .* win, nfft);
    
    S_ref = S_ref(1:end/2+1, :);
    S_test = S_test(1:end/2+1, :);
    
    % 1. Complex Error (NMSE) -> sum(|Ref - Test|^2)
    err_spec = sum(abs(S_ref - S_test).^2, 2);
    
    % 2. Magnitude Error -> sum(| |Ref| - |Test| |^2)
    mag_err_spec = sum(abs(abs(S_ref) - abs(S_test)).^2, 2);
    
    % 3. Reference Power
    ref_spec = sum(abs(S_ref).^2, 2);
end

function [out, out_f] = calculate_all_metrics(sig_ref, sig_test, fs, conf, ERBstruct)
    nfft = 2^nextpow2(size(sig_ref, 1));
    X_ref  = fft(sig_ref, nfft);
    X_test = fft(sig_test, nfft);
    
    f_axis = linspace(0, fs, nfft)';
    idx_pos = 1:floor(nfft/2)+1;
    X_ref  = X_ref(idx_pos, :);
    X_test = X_test(idx_pos, :);
    f_axis = f_axis(idx_pos);
    
    idx_nmse = find(f_axis >= conf.f_low_limit & f_axis <= conf.f_mid_limit);
    idx_mag  = find(f_axis > conf.f_mid_limit  & f_axis <= conf.f_high_limit);
    
    [out.nmse, raw_nmse] = calc_nmse_band(X_ref(idx_nmse, :), X_test(idx_nmse, :));
    [out.mag,  raw_mag]  = calc_mag_band(X_ref(idx_mag, :), X_test(idx_mag, :));
    [out.lsd,  raw_lsd]  = calc_lsd(X_ref(idx_mag, :), X_test(idx_mag, :));
    [out.bsd,  out_f.bsd]  = calc_bsd_wrapper(sig_ref, sig_test, ERBstruct); % ERB is inherently fixed-size
    
    % --- FIX FOR VARYING SIGNAL LENGTHS ---
    % Interpolate the output vectors onto a fixed standardized frequency grid (2000 points)
    N_plot_bins = 2000; 
    out_f.f_nmse = linspace(conf.f_low_limit, conf.f_mid_limit, N_plot_bins)';
    out_f.f_mag  = linspace(conf.f_mid_limit, conf.f_high_limit, N_plot_bins)';
    out_f.f_lsd  = out_f.f_mag;
    out_f.f_bsd  = ERBstruct.f_c;
    
    % Interpolate the raw data onto the fixed grid
    out_f.nmse = interp1(f_axis(idx_nmse), raw_nmse(:), out_f.f_nmse, 'linear', 'extrap');
    out_f.mag  = interp1(f_axis(idx_mag),  raw_mag(:),  out_f.f_mag,  'linear', 'extrap');
    out_f.lsd  = interp1(f_axis(idx_mag),  raw_lsd(:),  out_f.f_lsd,  'linear', 'extrap');
end

function [score, score_f] = calc_nmse_band(X_ref, X_test)
    err = abs(X_ref - X_test).^2;
    pow = abs(X_ref).^2;
    score = 10 * log10(mean(sum(err, 1) ./ (sum(pow, 1) + 1e-12)));
    score_f = 10 * log10(mean(err ./ (pow + 1e-12), 2) + 1e-12);
end

function [score, score_f] = calc_mag_band(X_ref, X_test)
    err = abs(abs(X_ref) - abs(X_test)).^2;
    pow = abs(X_ref).^2;
    score = 10 * log10(mean(sum(err, 1) ./ (sum(pow, 1) + 1e-12)));
    score_f = 10 * log10(mean(err ./ (pow + 1e-12), 2) + 1e-12);
end

function [score, score_f] = calc_lsd(X_ref, X_test)
    P_ref  = abs(X_ref).^2 + 1e-12;
    P_test = abs(X_test).^2 + 1e-12;
    log_diff_sq = (10*log10(P_test ./ P_ref)).^2;
    score = mean(sqrt(mean(log_diff_sq, 1))); 
    score_f = mean(sqrt(log_diff_sq), 2);
end

function [score, score_f] = calc_bsd_wrapper(t_ref, t_test, ERBstruct)
    bsd_val = BSD_ERB_fast(t_ref, t_test, ERBstruct);
    score = mean(mean(bsd_val)); 
    score_f = mean(bsd_val, 2);
end

function [out, out_f] = calculate_all_metrics_stft(sig_ref, sig_test, fs, conf, ERBstruct)


    win_size = 512; 
    overlap = win_size / 2; % Recommended: 50% overlap for audio metrics
    nfft = 512;
    window = hann(win_size); % Recommended: Hann window to prevent spectral leakage
    
    num_channels = size(sig_ref, 2);
    
    % Get Freq Axis
    [~, f_axis, ~] = spectrogram(sig_ref(:, 1), window, overlap, nfft, fs);
    idx_nmse = find(f_axis >= conf.f_low_limit & f_axis <= conf.f_mid_limit);
    idx_mag  = find(f_axis > conf.f_mid_limit  & f_axis <= conf.f_high_limit);
    
    out_f.f_nmse = f_axis(idx_nmse);
    out_f.f_mag  = f_axis(idx_mag);
    out_f.f_lsd  = f_axis(idx_mag);
    
    nmse_ears = zeros(1, num_channels); 
    mag_ears = zeros(1, num_channels); 
    lsd_ears = zeros(1, num_channels);
    
    nmse_f_ears = zeros(length(idx_nmse), num_channels);
    mag_f_ears  = zeros(length(idx_mag), num_channels);
    lsd_f_ears  = zeros(length(idx_mag), num_channels);
    
    for ch = 1:num_channels
        [S_ref, ~, ~] = spectrogram(sig_ref(:, ch), window, overlap, nfft, fs);
        [S_test, ~, ~] = spectrogram(sig_test(:, ch), window, overlap, nfft, fs);
        ref_sq = abs(S_ref).^2;
        
        % 1. NMSE STFT (Fixed ratio of means)
        err_sq = abs(S_ref - S_test).^2;
        mse_t = sum(err_sq(idx_nmse, :), 1);
        ref_t = sum(ref_sq(idx_nmse, :), 1);
        act = ref_t > 1e-12;
        if any(act)
            %nmse_ears(ch) = mean(mse_t(act) ./ ref_t(act));
            nmse_f_ears(:, ch) = mean(err_sq(idx_nmse, act), 2) ./ (mean(ref_sq(idx_nmse, act), 2) + 1e-12);
            
            nmse_ears(ch) = squeeze(mean(10 * log10(nmse_f_ears(:, ch)+1e-12),1));
        end
        
        % 2. Mag STFT (Fixed ratio of means)
        err_mag = abs(abs(S_ref) - abs(S_test)).^2;
        mse_mag_t = sum(err_mag(idx_mag, :), 1);
        ref_mag_t = sum(ref_sq(idx_mag, :), 1);
        act_m = ref_mag_t > 1e-12;
        if any(act_m)
            %mag_ears(ch) = mean(mse_mag_t(act_m) ./ ref_mag_t(act_m));
            mag_f_ears(:, ch) = mean(err_mag(idx_mag, act_m), 2) ./ (mean(ref_sq(idx_mag, act_m), 2) + 1e-12);
            mag_ears(ch) = squeeze(mean(10 * log10(mag_f_ears(:, ch)+1e-12),1));
        end
        
        % 3. LSD STFT (Valid)
        P_ref  = abs(S_ref(idx_mag, :)).^2 + 1e-12;
        P_test = abs(S_test(idx_mag, :)).^2 + 1e-12;
        log_diff_sq = (10 * log10(P_test ./ P_ref)).^2;
        %lsd_ears(ch) = mean(sqrt(mean(log_diff_sq(:, act_m), 1)));
        lsd_f_ears(:, ch) = sqrt(mean(log_diff_sq(:, act_m), 2));
        lsd_ears(ch) = squeeze(mean(lsd_f_ears(:, ch),1));
    end
    
    out.nmse = mean(nmse_ears) ;
    out_f.nmse = 10 * log10(mean(nmse_f_ears, 2) + 1e-12);
    
    out.mag = mean(mag_ears);
    out_f.mag = 10 * log10(mean(mag_f_ears, 2) + 1e-12);
    
    out.lsd = mean(lsd_ears);
    out_f.lsd = mean(lsd_f_ears, 2);
    
    % 4. BSD STFT (Fixed active frame division)
    % Note: Keeping this zero-overlap to match original explicit block handling.
    num_frames = floor(size(sig_ref, 1) / win_size);
    pad_len = max(win_size, size(ERBstruct.C, 1)); 
    bsd_frame_scores = zeros(num_frames, 1);
    bsd_f_frames = zeros(length(ERBstruct.f_c), num_frames);
    v_frames = 0;
    
    for fr = 1:num_frames
        idx = (fr-1)*win_size + 1 : fr*win_size;
        frame_ref = sig_ref(idx, :);
        
        if sum(frame_ref(:).^2) > 1e-8
            v_frames = v_frames + 1;
            frame_ref_pad  = [frame_ref; zeros(pad_len - win_size, size(frame_ref,2))];
            frame_test_pad = [sig_test(idx, :); zeros(pad_len - win_size, size(sig_test,2))];
            bsd_val = BSD_ERB_fast(frame_ref_pad, frame_test_pad, ERBstruct);
            
            bsd_frame_scores(v_frames) = mean(bsd_val(:)); % Pack tightly using v_frames index
            bsd_f_frames(:, v_frames)  = mean(bsd_val, 2);
        end
    end
    
    if v_frames > 0
        % Calculate using exact v_frames denominator to prevent perfect-score exclusion
        %out.bsd = sum(bsd_frame_scores(1:v_frames)) / v_frames; 
        out_f.bsd = sum(bsd_f_frames(:, 1:v_frames), 2) / v_frames;
        out.bsd   = mean(out_f.bsd,1);
    else
        out.bsd = 0; 
        out_f.bsd = zeros(length(ERBstruct.f_c), 1);
    end
end



function BSD = BSD_ERB_fast(x_ref, x_test, ERBstruct)
    fs = ERBstruct.fs;
    C = ERBstruct.C;
    n = ERBstruct.n;
    len = max(size(x_ref,1), size(x_test,1));
    x_ref(end+1:len,:) = 0;
    x_test(end+1:len,:) = 0;
    N = len;
    c_ears = size(x_ref,2);
    X_out_ref = zeros(n,c_ears);
    X_out_test = zeros(n,c_ears);
    for k = 1:n
        Ncur = size(x_ref,1); 
        Hk = abs(fft(C(:,k), Ncur));
        X_ref_band  = (repmat(Hk, [1, c_ears]) .* abs(fft(x_ref, Ncur))).^2;
        X_test_band = (repmat(Hk, [1, c_ears]) .* abs(fft(x_test, Ncur))).^2;
        X_out_ref(k,:)  = sum(X_ref_band, 1);
        X_out_test(k,:) = sum(X_test_band, 1);
    end
    BSD = abs(10*log10(X_out_test ./ (X_out_ref + 1e-12)));
end



function print_comparison_table(title_str, T)
    fprintf('\n--- %s ---\n', title_str);
    
    % Find all evaluated algorithms by scanning table headers
    var_names = T.Properties.VariableNames;
    algo_cols = var_names(startsWith(var_names, 'NMSE_'));
    algos = strrep(algo_cols, 'NMSE_', '');
    
    % List of all metrics calculated
    % metrics = {'NMSE', 'MagErr', 'LSD', 'BSD', ...
    %            'stft_NMSE', 'stft_MagErr', 'stft_LSD', 'stft_BSD', ...
    %            'RMSE_BMS', 'RMSE_ILD', 'RMSE_IC', 'Bias_IC'};

    metrics = {'stft_NMSE', 'stft_MagErr', 'stft_LSD','stft_BSD', ...
               'RMSE_BMS', 'RMSE_ILD', 'RMSE_IC'};
    
    
    % Create dynamic header
    header_str = sprintf('%-15s | ', 'Metric');
    for a = 1:length(algos)
        header_str = [header_str, sprintf('%-20s | ', algos{a})];
    end
    fprintf('%s\n', header_str);
    fprintf('%s\n', repmat('-', 1, length(header_str)));
    
    % Print data rows dynamically
    for i = 1:length(metrics)
        m = metrics{i};
        row_str = sprintf('%-15s | ', m);
        for a = 1:length(algos)
            col_name = sprintf('%s_%s', m, algos{a});
            if any(strcmp(var_names, col_name))
                vals = T.(col_name);
                val_str = sprintf('%6.2f +/- %4.2f', mean(vals), std(vals));
            else
                val_str = sprintf('%15s', 'N/A');
            end
            row_str = [row_str, sprintf('%-20s | ', val_str)];
        end
        fprintf('%s\n', row_str);
    end
    fprintf('%s\n\n', repmat('-', 1, length(header_str)));
end

function ERBstruct = create_erb_struct(fs, f_lim, n_bands)
    ERBstruct.fs = fs;
    ERBstruct.f_lim = f_lim;
    ERBstruct.n = n_bands;
    low_erb = 21.4 * log10(4.37 * f_lim(1) / 1000 + 1);
    high_erb = 21.4 * log10(4.37 * f_lim(2) / 1000 + 1);
    erb_val = linspace(low_erb, high_erb, n_bands)';
    f_c = (10.^(erb_val / 21.4) - 1) / 4.37 * 1000;
    ERBstruct.f_c = f_c;
    T = 1/fs;
    len_ir = round(0.1 * fs); 
    t = (0:len_ir-1)' * T;
    C = zeros(len_ir, n_bands);
    for k = 1:n_bands
        fc = f_c(k);
        bw = 24.7 * (4.37 * fc / 1000 + 1);
        b  = 1.019 * bw;
        env = t.^(4-1) .* exp(-2*pi*b*t);
        ir  = env .* cos(2*pi*fc*t);
        C(:,k) = ir / max(abs(fft(ir)));
    end
    ERBstruct.C = C;
end


function [BMS_erb, ILD_erb, IC_erb, f_erb] = compute_binaural_metrics(p_t, fs, conf)
    % Extracts BMS, ILD, and IC from a binaural signal [N_samples x 2]
    % using 512 samples with no overlap.
    
    win_size = 512;
    overlap = 0;
    nfft = 512;
    
    % Left and Right channel spectrograms (using Rectangular window for no overlap)
    [S_L, F, ~] = spectrogram(p_t(:, 1), rectwin(win_size), overlap, nfft, fs);
    [S_R, ~, ~] = spectrogram(p_t(:, 2), rectwin(win_size), overlap, nfft, fs);
    
    % Expected values (mean across time frames) -> Eq 33
    c_ll = mean(abs(S_L).^2, 2);
    c_rr = mean(abs(S_R).^2, 2);
    c_lr = mean(S_L .* conj(S_R), 2);
    
    % Compute metrics per frequency bin -> Eq 34, 35, 36
    BMS_f = 10 * log10(c_ll + c_rr + 1e-12);
    ILD_f = 10 * log10((c_ll + 1e-12) ./ (c_rr + 1e-12));
    IC_f  = real(c_lr) ./ (sqrt(c_ll .* c_rr) + 1e-12);
    
    % Average over ERB scale 
    erb_scale = 21.4 * log10(4.37 * F / 1000 + 1);
    
    % Define ERB band edges based on conf limits
    num_erb_bands = 32; 
    f_min_erb = 21.4 * log10(4.37 * conf.f_low_limit / 1000 + 1);
    f_max_erb = 21.4 * log10(4.37 * conf.f_high_limit / 1000 + 1);
    erb_edges = linspace(f_min_erb, f_max_erb, num_erb_bands + 1);
    
    BMS_erb = zeros(num_erb_bands, 1);
    ILD_erb = zeros(num_erb_bands, 1);
    IC_erb  = zeros(num_erb_bands, 1);
    f_erb   = zeros(num_erb_bands, 1); 
    
    for b = 1:num_erb_bands
        % Find STFT bins that fall into the current ERB band
        idx = find(erb_scale >= erb_edges(b) & erb_scale < erb_edges(b+1));
        
        if isempty(idx)
            continue; 
        end
        
        % Average the metric values within this ERB band
        BMS_erb(b) = mean(BMS_f(idx));
        ILD_erb(b) = mean(ILD_f(idx));
        IC_erb(b)  = mean(IC_f(idx));
        f_erb(b)   = (10^((erb_edges(b) + erb_edges(b+1))/2 / 21.4) - 1) / 4.37 * 1000;
    end
end

function metrics = compute_binaural_errors(BMS_ref, ILD_ref, IC_ref, BMS_test, ILD_test, IC_test)
    % Computes RMSE and biases across frequency bands
    valid_idx = BMS_ref ~= 0; 
    N_f = sum(valid_idx);
    
    if N_f == 0
        metrics.bms_rmse = NaN;
        metrics.ild_rmse = NaN;
        metrics.ic_rmse  = NaN;
        metrics.ic_bias  = NaN;
        return;
    end
    
    metrics.bms_rmse = sqrt((1/N_f) * sum((BMS_ref(valid_idx) - BMS_test(valid_idx)).^2));
    metrics.ild_rmse = sqrt((1/N_f) * sum((ILD_ref(valid_idx) - ILD_test(valid_idx)).^2));
    metrics.ic_rmse  = sqrt((1/N_f) * sum((IC_ref(valid_idx) - IC_test(valid_idx)).^2));
    
    % Delta IC Bias: IC_test - IC_ref
    % Positive bias = test is too coherent (lacking envelopment)
    metrics.ic_bias = (1/N_f) * sum(IC_test(valid_idx) - IC_ref(valid_idx));
end


function plot_frequency_metrics(freq_data, algos)
    % FIGURE 1: Standard Objective Metrics
    figure('Name', 'Objective Error Metrics vs Frequency');
    colors = lines(length(algos));
    
    % Subplot 1: NMSE
    subplot(2,2,1); hold on; grid on;
    for a = 1:length(algos)
        %plot(freq_data.axes.nmse, mean(freq_data.(algos{a}).nmse, 2), 'Color', colors(a,:), 'LineWidth', 1.5, 'DisplayName', [algos{a} ' (FFT)']);
        plot(freq_data.axes.stft_nmse, mean(freq_data.(algos{a}).stft_nmse, 2), 'Color', colors(a,:), 'LineWidth', 1.5, 'DisplayName', [algos{a} ]);
    end
    set(gca, 'XScale', 'log'); xlim([freq_data.conf.f_low_limit freq_data.conf.f_mid_limit]);
    xlabel('Frequency (Hz)'); ylabel('NMSE (dB)'); title('NMSE'); legend('Location', 'best');
    
    % Subplot 2: MagErr
    subplot(2,2,2); hold on; grid on;
    for a = 1:length(algos)
        %plot(freq_data.axes.mag, mean(freq_data.(algos{a}).mag, 2), 'Color', colors(a,:), 'LineWidth', 1.5, 'DisplayName', [algos{a} ' (FFT)']);
        plot(freq_data.axes.stft_mag, mean(freq_data.(algos{a}).stft_mag, 2),  'Color', colors(a,:), 'LineWidth', 1.5, 'DisplayName', [algos{a} ]);
    end
    set(gca, 'XScale', 'log'); xlim([freq_data.conf.f_mid_limit freq_data.conf.f_high_limit]);
    xlabel('Frequency (Hz)'); ylabel('Mag Error (dB)'); title('Magnitude Error'); legend('Location', 'best');
    
    % Subplot 3: LSD
    subplot(2,2,3); hold on; grid on;
    for a = 1:length(algos)
        %plot(freq_data.axes.lsd, mean(freq_data.(algos{a}).lsd, 2), 'Color', colors(a,:), 'LineWidth', 1.5, 'DisplayName', [algos{a} ' (FFT)']);
        plot(freq_data.axes.stft_lsd, mean(freq_data.(algos{a}).stft_lsd, 2),'Color', colors(a,:), 'LineWidth', 1.5, 'DisplayName', [algos{a} ]);
    end
    set(gca, 'XScale', 'log'); xlim([freq_data.conf.f_mid_limit freq_data.conf.f_high_limit]);
    xlabel('Frequency (Hz)'); ylabel('LSD (dB)'); title('Log-Spectral Distance'); legend('Location', 'best');
    
    % Subplot 4: BSD
    subplot(2,2,4); hold on; grid on;
    for a = 1:length(algos)
        %plot(freq_data.axes.bsd, mean(freq_data.(algos{a}).bsd, 2), 'Color', colors(a,:), 'LineWidth', 1.5, 'DisplayName', [algos{a} ' (FFT)']);
        plot(freq_data.axes.bsd, mean(freq_data.(algos{a}).stft_bsd, 2),  'Color', colors(a,:), 'LineWidth', 1.5, 'DisplayName', [algos{a} ]);
    end
    set(gca, 'XScale', 'log'); xlim([freq_data.conf.f_mid_limit freq_data.conf.f_high_limit]);
    xlabel('Frequency (Hz)'); ylabel('BSD'); title('Binaural Spectral Distance'); legend('Location', 'best');
    
    
    % FIGURE 2: Binaural Errors
    figure('Name', 'Binaural RMSE vs Frequency (ERB)');
    
    subplot(1,3,1); hold on; grid on;
    for a = 1:length(algos), plot(freq_data.axes.erb, sqrt(mean(freq_data.(algos{a}).bms_sq, 2)), 'Color', colors(a,:), 'LineWidth', 2, 'DisplayName', algos{a}); end
    set(gca, 'XScale', 'log'); xlim([freq_data.conf.f_low_limit freq_data.conf.f_high_limit]);
    xlabel('Frequency (Hz)'); ylabel('RMSE (dB)'); title('BMS RMSE'); legend('Location', 'best');
    
    subplot(1,3,2); hold on; grid on;
    for a = 1:length(algos), plot(freq_data.axes.erb, sqrt(mean(freq_data.(algos{a}).ild_sq, 2)), 'Color', colors(a,:), 'LineWidth', 2, 'DisplayName', algos{a}); end
    set(gca, 'XScale', 'log'); xlim([freq_data.conf.f_low_limit freq_data.conf.f_high_limit]);
    xlabel('Frequency (Hz)'); ylabel('RMSE (dB)'); title('ILD RMSE'); legend('Location', 'best');
    
    subplot(1,3,3); hold on; grid on;
    for a = 1:length(algos), plot(freq_data.axes.erb, sqrt(mean(freq_data.(algos{a}).ic_sq, 2)), 'Color', colors(a,:), 'LineWidth', 2, 'DisplayName', algos{a}); end
    set(gca, 'XScale', 'log'); xlim([freq_data.conf.f_low_limit freq_data.conf.f_high_limit]);
    xlabel('Frequency (Hz)'); ylabel('RMSE'); title('IC RMSE'); legend('Location', 'best');
end



