function verify_universal_resampling(ATF_direct, HRTF_direct)
% VERIFY_UNIVERSAL_RESAMPLING
%
% Visualizes how the "Centered Time-Padding" logic correctly handles
% BOTH the Split Impulse (ATF) and the Centered Impulse (HRTF).

    % Extract representative impulses
    % 1. ATF (Split)
    H_atf_f = squeeze(ATF_direct(1, :)); 
    % Force to Time Domain (N=1024 approx)
    nIn = length(H_atf_f) * 2 - 2; % Approx reconstruction
    H_full = [H_atf_f, conj(H_atf_f(end-1:-1:2))];
    h_atf_t = real(ifft(H_full));
    
    % 2. HRTF (Centered)
    H_hrtf_f = squeeze(HRTF_direct(1, :));
    H_full_h = [H_hrtf_f, conj(H_hrtf_f(end-1:-1:2))];
    h_hrtf_t = real(ifft(H_full_h));
    
    % Target Size
    nOut = 4096;
    
    % --- Apply Universal Fix ---
    h_atf_new  = resample_universal(h_atf_t, nOut);
    h_hrtf_new = resample_universal(h_hrtf_t, nOut);

    % --- Plotting ---
    figure('Color','w', 'Name', 'Universal Resampling Verification');
    
    % Row 1: ATF (The Problem Child)
    subplot(2,2,1);
    plot(h_atf_t, 'r', 'LineWidth', 1.5);
    title('Original ATF (Split Impulse)'); grid on; xlim([1 length(h_atf_t)]);
    
    subplot(2,2,2);
    plot(h_atf_new, 'b', 'LineWidth', 1.5);
    title('Resampled ATF (Fixed)'); 
    subtitle('Tail wraps correctly to the NEW end (4096)');
    grid on; xlim([1 nOut]);
    
    % Row 2: HRTF (The Standard)
    subplot(2,2,3);
    plot(h_hrtf_t, 'r', 'LineWidth', 1.5);
    title('Original HRTF (Centered)'); grid on; xlim([1 length(h_hrtf_t)]);
    
    subplot(2,2,4);
    plot(h_hrtf_new, 'b', 'LineWidth', 1.5);
    title('Resampled HRTF (Fixed)'); 
    subtitle('Pulse stays centered and intact');
    grid on; xlim([1 nOut]);
    
    fprintf('VERIFICATION:\n');
    fprintf('1. ATF Plot (Top Right): Should have spikes at 0 and 4096. Middle should be 0.\n');
    fprintf('2. HRTF Plot (Bottom Right): Should look exactly like Bottom Left, just with more space.\n');
end

function h_out = resample_universal(h_in, N_new)
    N_old = length(h_in);
    
    % 1. Shift Center to 0
    h_centered = fftshift(h_in);
    
    % 2. Pad Edges
    pad_total = N_new - N_old;
    pad_L = floor(pad_total/2);
    pad_R = ceil(pad_total/2);
    h_pad = [zeros(1, pad_L), h_centered, zeros(1, pad_R)]; % Row vector assumpt
    if iscolumn(h_in), h_pad = [zeros(pad_L,1); h_centered; zeros(pad_R,1)]; end
    
    % 3. Unshift
    h_out = ifftshift(h_pad);
end