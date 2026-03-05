function hobj_rot_BSM= process_hrtf(BSMobj,HRTFpath,WignerDpath)
    hobj = sofaToEaro(convertStringsToChars(HRTFpath));   % hobj is HRIR earo object
    if ~(BSMobj.head_rot_az == 0)
        load(WignerDpath);
        %N_HRTF_rot = BSMobj.N_PW;
        N_HRTF_rot = 32;
        DN = (N_HRTF_rot + 1)^2; % size of the wignerD matrix
        D_allAngles = D(:, 1 : DN);
    else
        N_HRTF_rot = BSMobj.N_PW;
        %N_HRTF_rot = 42;
        %N_HRTF_rot = 35;
        %N_HRTF_rot = 12; % FOR HATS
        
        %DN = (N_HRTF_rot + 1)^2; % size of the wignerD matrix
        %D_allAngles = eye(DN);
    end

    hobj.shutUp = true;
    %%Interpolate HRTF to frequencies
    hobj_freq_grid = hobj;
    if strcmp(hobj_freq_grid.dataDomain{1},'FREQ'), hobj_freq_grid=hobj_freq_grid.toTime(); end
    % resample HRTF to fs
    hobj_freq_grid = hobj_freq_grid.resampleData(BSMobj.fs);
    hobj_freq_grid = hobj_freq_grid.toFreq(BSMobj.filt_samp);
    hobj_freq_grid.taps = BSMobj.filt_samp;
    % Trim negative frequencies
    hobj_freq_grid.data = hobj_freq_grid.data(:, 1:ceil(BSMobj.filt_samp/2)+1, :);

    %hobj_rot = BSM_toolbox.RotateHRTF(hobj_freq_grid, N_HRTF_rot, D_allAngles, BSMobj.head_rot_az);
    
    % Interpolate HRTF to BSM grid
    hobj_rot_BSM = hobj_freq_grid.toSH(N_HRTF_rot,'SRC');
    hobj_rot_BSM = hobj_rot_BSM.toSpace('SRC', BSMobj.th_BSMgrid_vec, BSMobj.ph_BSMgrid_vec);
    hobj_rot_BSM.sourceGrid.quadWeight = []; % I dont have the Sampling wights for the BSMgrid
    hobj_rot_BSM.nData = size(hobj_rot_BSM.data,1);
    hobj_rot_BSM.sourceGrid.r = hobj.sourceGrid.r(1:length(BSMobj.ph_BSMgrid_vec));

    % tmp_hobj = hobj_rot_BSM.toTime(BSMobj.filt_samp,true);
    % tmp_hobj_2 = hobj_freq_grid.toTime(BSMobj.filt_samp,true);
    % tmp_hobj = hobj_rot_BSM.toTime(BSMobj.filt_samp,true);
    % 
    % [dataID1, azimuth1, elevation1]= tmp_hobj.closestIrSrc(-90,90);
    % [dataID2, azimuth2, elevation2]= tmp_hobj_2.closestIrSrc(-90,90);


    % figure
    % plot(squeeze(tmp_hobj_2.data(dataID2,:,:)))
    % ylim([-1,1])
    % xlim([0,400])
    % figure
    % plot(squeeze(tmp_hobj.data(dataID1,:,:)))
    % ylim([-1,1])
    % xlim([0,400])
    % figure
    % plot(squeeze(hobj.data(dataID2,:,:)))
    % ylim([-1,1])
    % xlim([0,400])



end