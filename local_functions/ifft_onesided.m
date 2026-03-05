function h_t = ifft_onesided(H, nfft, dim)
% IFFT_ONESIDED Reconstructs time-domain signal from one-sided FFT.
%
%   h_t = ifft_onesided(H)
%       Assumes H is from an EVEN length FFT (N = 2*(size(H,dim)-1)).
%       Operates along dimension 2 by default.
%
%   h_t = ifft_onesided(H, nfft)
%       Reconstructs based on specific nfft (handles both ODD and EVEN).
%       If nfft is Odd, it uses odd-reconstruction logic.
%       If nfft is Even, it uses even-reconstruction logic.
%
%   h_t = ifft_onesided(H, nfft, dim)
%       Operates along specified dimension 'dim'.

    % 1. Default Dimension: If not specified, assume dim 2
    if nargin < 3 || isempty(dim)
        dim = 2;
    end
    
    % Get the number of frequency bins
    n_bins = size(H, dim);

    % 2. Default NFFT: If not specified, assume standard Even case
    if nargin < 2 || isempty(nfft)
        nfft = 2 * (n_bins - 1);
    end

    % 3. Reconstruction Logic
    if mod(nfft, 2) == 0
        % --- EVEN CASE (e.g., 257 bins -> 512 points) ---
        % Structure: [DC, PosFreqs, Nyquist]
        % We need to reconstruct: [DC, Pos, Nyq, Conj(Flip(Pos))]
        
        % Check if input size matches expected even bins
        if n_bins ~= nfft/2 + 1
            error('Input size %d does not match expected size %d for Even NFFT %d', n_bins, nfft/2+1, nfft);
        end
        
        % Extract components using generic slicing
        idx_dc  = 1;
        idx_pos = 2 : (n_bins - 1);
        idx_nyq = n_bins;
        
        % Create slices
        dc  = extract_slice(H, idx_dc, dim);
        pos = extract_slice(H, idx_pos, dim);
        nyq = extract_slice(H, idx_nyq, dim);
        
        % Reconstruct full spectrum
        % cat(dim, DC, Pos, Nyq, ConjFlipPos)
        H_full = cat(dim, dc, pos, nyq, conj(flip(pos, dim)));
        
    else
        % --- ODD CASE (e.g., 257 bins -> 513 points) ---
        % Structure: [DC, PosFreqs] (No Nyquist bin in Odd FFTs)
        % We need to reconstruct: [DC, Pos, Conj(Flip(Pos))]
        
        % Check if input size matches expected odd bins
        if n_bins ~= (nfft + 1) / 2
             error('Input size %d does not match expected size %d for Odd NFFT %d', n_bins, (nfft+1)/2, nfft);
        end
        
        % Extract components
        idx_dc  = 1;
        idx_pos = 2 : n_bins;
        
        % Create slices
        dc  = extract_slice(H, idx_dc, dim);
        pos = extract_slice(H, idx_pos, dim);
        
        % Reconstruct full spectrum
        % cat(dim, DC, Pos, ConjFlipPos)
        H_full = cat(dim, dc, pos, conj(flip(pos, dim)));
    end

    % 4. Compute Inverse FFT
    % 'symmetric' ensures output is real (discards floating point imaginary noise)
    h_t = ifft(H_full, nfft, dim, 'symmetric');

end

% Helper function to keep the main code clean 
function slice = extract_slice(arr, indices, dim)
    % Creates a dynamic index for any dimension
    S.type = '()';
    S.subs = repmat({':'}, 1, ndims(arr));
    S.subs{dim} = indices;
    slice = subsref(arr, S);
end