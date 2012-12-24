function Z_ScanMain
    % Z-scan fitting program
    clear all;
    % init & GO!
    silicaOpts = initSilicaOpts();
    ZScan_fit(silicaOpts);
end

function opt = initSilicaOpts
    opt = Options;
    opt.dataType = 'silica';
    opt.lambda = 800; %nm
    opt.filename = 'silica2.txt';
end