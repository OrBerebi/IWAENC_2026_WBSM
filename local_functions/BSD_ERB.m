
function [BSD,f_c,n] = BSD_ERB(x_ref, x_test, f_lim, fs)

    % default parameter
    if ~exist('f_lim', 'var')
        f_lim = [50 20000];
    end
    if ~exist('fs', 'var')
        fs = 44100;
    end

    % throw a warning if input appears to be zero phase
    % (all checks if all angles for one spectrum are zero,
    %  any checks if this is true for any channel)
    xPhase      = any(all(abs(angle(fft(x_ref)))  < 1e-10));
    targetPhase = any(all(abs(angle(fft(x_test))) < 1e-10));

    % don't throw the warning if it is zero phase but a dirac impulse
    xDirac      = all( all( x_ref(2:end,:) < 1e-10 ) );
    targetDirac = all( all( x_test(2:end,:) < 1e-10 ) );

    if xPhase && ~xDirac
        warning('AKerbError:Input', 'At least one channel of ''x'' appears to be zero phase. This can produce errors, and can be avoided by circshift(x, [round(size(x,1)/2) 0])')
    end
    if targetPhase  && ~targetDirac
        warning('AKerbError:Input', 'At least one channel of ''target'' appears to be zero phase. This can produce errors, and can be avoided by circshift(x, [round(size(target,1)/2) 0])')
    end


    % match IR length
    if size(x_ref,1)~=size(x_test,1)
        if size(x_ref,1)>size(x_test,1)
            x_test(end+1:size(x_ref,1),:) = 0;
        else
            x_ref(end+1:size(x_test,1),:) = 0;
        end
    end

    % match IR channels
    if size(x_ref,2) > size(x_test,2)
        x_test = repmat(x_test, [1 size(x_ref,2)]);
    end

    % get length and channels
    N = size(x_ref,1);  % ir length
    c = size(x_ref,2);  % ir num. of channels

    % ---------------------------------------- 1. get FIR filter in freq domain
    % get ERB filter coefficients
    n      = round(21.4*log10(0.004367*f_lim(2)+1)-21.4*log10(0.004367*f_lim(1)+1));    % number of filters
    f_c    = ERBSpace(f_lim(1), fs/2, n);                                               % center frequencies
    fcoefs = MakeERBFilters(fs,n,f_lim(1));                                             % filter coefficients

    % sort with ascending frequency
    fcoefs = flipud(fcoefs);
    f_c    = flipud(f_c);

    % get rid of filters above upper f_lim
    n      = find(f_c <= f_lim(2), 1, 'last');
    f_c    = f_c(1:n);
    fcoefs = fcoefs(1:n,:);

    % estimate needed filter length (4 times the largest cycle)
    Nmax = 2^nextpow2( 1/f_c(1) * 4 * fs );

    % get filter impulse responses
    C = zeros(Nmax,1);
    C(1) = 1;
    C = ERBFilterBank(C',fcoefs)';

    % make minimum phase
    C = AKphaseManipulation(C, fs, 'min', 4, 0); %O.B

    % filter x with filerbank and calculate error in each band
    x_out = zeros(n,c,2);

    X_out_ref = zeros(n,c);
    X_out_test = zeros(n,c);

    C_db_tmp = db(C);
    C_db = zeros(size(C_db_tmp));
    % ------------------------------------- 2. filter signals and get ERB error
    for k = 1:n
        
        % find point where filter decayed for 60 dB
        C_db = C_db_tmp(:,k);
        %C_db = 10*log10(C(:,k)); %O.B
        Nmin = find(flipud(C_db) > max(C_db)-60, 1, 'first');
        Nmin = 2^nextpow2(Nmax - Nmin + 1);
        
        % needed filter length at current center frequency
        Ncur = 2^nextpow2( 1/f_c(k) * 4 * fs );
        Ncur = max([Ncur Nmin N]);
        
        % get indeces of upper and lower frequency limit
        f_lim_id(1) = round(f_lim(1)/(fs/Ncur))+1;
        f_lim_id(1) = max([f_lim_id(1) 2]);         % make sure we dont use the 0 Hz bin
        f_lim_id(2) = round(f_lim(2)/(fs/Ncur))+1;
        
        % filter input signals
        %t1 = repmat(abs(fft(C(:,k), Ncur)), [1, c]);
        %t2 = (abs(fft(x_l, Ncur)).^2);
        X_ref = repmat(abs(fft(C(:,k), Ncur)), [1, c]) .* (abs(fft(x_ref, Ncur)).^2);
        X_test = repmat(abs(fft(C(:,k), Ncur)), [1, c]) .* (abs(fft(x_test, Ncur)).^2);
        
        % get energy
        X_ref = sum(X_ref(f_lim_id(1):f_lim_id(2),:));
        X_test = sum(X_test(f_lim_id(1):f_lim_id(2),:));
        
        
        % get ERB error
        X_out_ref(k,:) = X_ref;
        X_out_test(k,:) = X_test;

    end

    BSD = abs(10*log10(abs(X_out_test)./abs(X_out_ref+eps)));

end

