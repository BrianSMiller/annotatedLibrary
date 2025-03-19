function sp = spectroParams(type)
varNames = {'nrow','ncol','pre','post','sampleRate','nfft','noverlap','lowFreq','highFreq'};
varTypes = {'int32' ,'int32', 'double','double','double','int32','int32','double','double'};
sp = table('Size',[1 9],'VariableTypes',varTypes,'VariableNames',varNames);

switch(lower(type))
    case {'blue','abw','pbw','bmd','bm-d'}
        sp.nrow = 4; % Number of rows per page
        sp.ncol = 4; % Number of columns per page

        sp.pre = 5;
        sp.post = 5;

        sp.sampleRate = 250;
        sp.nfft = 128;
        sp.noverlap = 112;
        sp.lowFreq = 30;
        sp.highFreq = 125;

    case {'srw','srwUp','upcall','srw_up','narw','narw-up','narwup'}
        sp.nrow = 4; % Number of rows per page
        sp.ncol = 4; % Number of columns per page

        sp.pre = 1;
        sp.post = 1;

        sp.sampleRate = 1000;
        sp.nfft = 256;
        sp.noverlap = 224;
        sp.lowFreq = 30;
        sp.highFreq = 500;
end
