function ERBstruct = prepare_ERB(f_lim, fs)

    if nargin < 1, f_lim = [50 20000]; end
    if nargin < 2, fs = 48e3; end

    % number of filters and center frequencies
    n = round(21.4*log10(0.004367*f_lim(2)+1) ...
             -21.4*log10(0.004367*f_lim(1)+1));
    f_c = ERBSpace(f_lim(1), fs/2, n);
    fcoefs = MakeERBFilters(fs, n, f_lim(1));

    % sort ascending frequency
    fcoefs = flipud(fcoefs);
    f_c    = flipud(f_c);

    % get rid of filters above upper f_lim
    n      = find(f_c <= f_lim(2), 1, 'last');
    f_c    = f_c(1:n);
    fcoefs = fcoefs(1:n,:);

    % estimate needed filter length
    Nmax = 2^nextpow2( 1/f_c(1) * 4 * fs );

    % get filter impulse responses
    C = zeros(Nmax,1);
    C(1) = 1;
    C = ERBFilterBank(C', fcoefs)';

    % make minimum phase
    C = AKphaseManipulation(C, fs, 'min', 4, 0);

    % precompute filter FFT magnitudes for different lengths
    ERBstruct.fs    = fs;
    ERBstruct.f_lim = f_lim;
    ERBstruct.f_c   = f_c;
    ERBstruct.C     = C;       % filterbank IRs
    ERBstruct.n     = n;

end
