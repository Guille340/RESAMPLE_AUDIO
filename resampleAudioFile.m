%  xr = RESAMPLEAUDIOFILE(filePath,channel,samples,P,Q,varargin)
%
%  DESCRIPTION
%  Reads a segment of the specified CHANNEL in the audio file located in 
%  FILEPATH and resamples it into frequency FS*P/Q. The segment starts at 
%  sample SAMPLES(1) and contains at total number of samples equal to 
%  SAMPLES(2). Input parameters P and Q are the upsampling and downsampling 
%  factors. The function supports audio files in WAV and RAW formats.
%
%  INPUT VARIABLES
%  - filePath [string]: absolute path of audio file
%  - channel [integer number]: selected audio channel. Use ¦channel¦ = []
%    to read all channels.
%  - samples [firstSample nSamples]: section of the audio file to be
%    read. The first element is the starting audio sample (firstSample =
%    1,2,...) and the second element the number of samples to be read.
%  - P [integer number]: upsampling factor (fr = fs*P/Q)
%  - Q [integer number]: downsampling factor (fr = fs*P/Q)
%  - sampleRate (varargin{1}) [integer number]: sampling rate (RAW).
%  - numChannels (varargin{2}) [integer number]: number of channels (RAW).
%  - bitsPerSample (varargin{3}) [integer number]: bit resolution (RAW).
%  - byteOrder (varargin{4}) [char]: endianess (RAW).
%    ¬ 'l': little endian
%    ¬ 'b': big endian
%
%  OUTPUT VARIABLES
%  - xr [numeric vector]: resampled audio. The data is returned as 'single'
%    class.
%
%  INTERNALLY CALLED FUNCTIONS
%  - readwavHeader
%  - readAudioFile
%
%  CONSIDERATIONS & LIMITATIONS
%  - Use #resamplingFactors to estimate the optimal upsampling and
%    downsampling factors based on memory and processing speed limitations.
%    Whenever possible, choose low values for P and Q that also result
%    in an integer resampling frequency (fr = fs*P/Q). The maximum of P
%    and Q determines the order of the low-pass recovery/anti-alias filter,
%    and P has also an important impact on the memory usage.
%
%  FUNCTION CALLS
%  1. xr = resampleAudioFile(filePath,channel,samples,P,Q)
%     - For 'WAV' format. Variable input arguments FS, NCH, NBIT
%       and EN can also be included in the call but they'll be ignored
%       if FILEPATH has .wav extension.
%
%  2. xr = resampleAudioFile(filePath,channel,samples,P,Q,fs,nch,nbit,en)
%     - Only for 'RAW' format.
%
%  See also READAUDIOFILE, READWAVHEADER

%  VERSION 1.1
%  Updated: 14 Aug 2021
%  Author: Guillermo Jimenez Arranz
%  - Fixed bug with overlap-save method that created a discontinuity
%    between processed windows.
%  - Cleaned code.
%  
%  VERSION 1.0
%  Guillermo Jimenez Arranz
%  email: gjarranz@gmail.com
%  06 Apr 2018

function xr = resampleAudioFile(filePath,channel,samples,P,Q,varargin)

[~,~,fext] = fileparts(filePath);

switch fext
    case {'.3gp','.aa','.aac','.aax','.act','.aiff','.amr','.ape',...
            '.au','.awb','.dct','.dss','.dvf','.flac','.gsm','.iklax',...
            '.ivs','.m4a','.m4b','.m4p','.mmf','.mp3','.mpc','.msv',...
            '.ogg','.oga','.mogg','.opus','.ra','.rm','.sln','.tta',...
            '.vox','.wma','.wv','.webm','.8svx'}
        error('Audio format not supported')

    case '.wav' % Import WAV file parameters
        
        % Error control: number of input arguments
        if nargin < 5, error('Not enough input arguments'); end
        if nargin > 5, error('Too many input arguments'); end

        % Read WAV header
        wavHeader = readwavHeader(filePath);

        % Size of audio data [bytes]
        fileStats = dir(filePath);
        if strcmp(wavHeader.subchunk2ID,'data')
            dataSize = wavHeader.subchunk2size;
            % subchunk2size = 0 when WAV file does not close properly (corrupt)
            if ~wavHeader.subchunk2size 
                dataSize = floor((fileStats.bytes-44)/wavHeader.blockAlign) ...
                    *wavHeader.blockAlign; % size of audio data [bytes]
                warning(['Size of audio data disagrees with header '...
                    'information: file is potentially corrupted. The '...
                    'data will be imported; note that this might result '...
                    'in a few initial samples being incorrectly read']);
            end
        elseif strcmp(wavHeader.subchunk2ID,'junk')
            dataSize = floor((fileStats.bytes-44)/wavHeader.blockAlign)...
                *wavHeader.blockAlign; % size of audio data [bytes]
        else
            error('Not recognised string for SUBCHUNK2ID')
        end

        % Audio parameters
        fs = wavHeader.sampleRate;
        nch = wavHeader.numChannels;
        nbit = wavHeader.bitsPerSample;
        en = wavHeader.byteOrder;

    otherwise % Import RAW file parameters
        % Error control: number of input arguments
        if nargin < 9, error('Not enough input arguments'); end
        if nargin > 9, error('Too many input arguments'); end

        % Size of audio data [bytes]
        fileStats = dir(filePath);
        dataSize = fileStats.bytes;

        % Audio parameters
        fs = varargin{1};
        nch = varargin{2};
        nbit = varargin{3};
        en = varargin{4};
end

% RESAMPLE AUDIO FILE
if ~isequal(P,Q)
    % Determine the Section of Audio to be Read
    allSamples = dataSize*8/(nbit*nch);
    if ~isempty(samples)
        firstSample = min(max(samples(1),1),allSamples); % first audio sample to read (1 to ALLSAMPLES)
        nSamples = min(samples(2),allSamples - firstSample + 1); % number of audio samples to read (1 to ALLSAMPLES)
           
        % Error Control: Start Sample and Number of Samples
        if firstSample ~= samples(1)
            warning(['The starting sample SAMPLE(1) is out of bounds. '...
                'SAMPLE(1) = %d will be used'],firstSample)
        end        
        if nSamples ~= samples(2)
            warning(['The number of samples SAMPLE(2) is out of bounds. '...
                'SAMPLE(2) = %d will be used'],nSamples)
        end
    else
        firstSample = 1;
        nSamples = allSamples;
    end

    % Design of Combined Recovery/Anti-Alias Filter
    fnpass = 1/max([P Q])-eps; % normalised -6dB passband frequency (fnpass = fpass/(0.5*fs))
    fnstop = 1.25*fnpass; % normalised stopband frequency (fnstop = fstop/(0.5*fs))
    atten = 60; % stopband attenuation
    N = 2*atten/(22*(fnstop-fnpass)); % number of taps for the FIR filter (Fred Harris "rule of thumb" method)
    N = ceil((N - 1)/P)*P + 1; % make the overlap multiple of the interpolation factor (N-1 ~ P)
    b = fir1(N-1,fnpass)'; % coefficients of low-pass FIR filter (order N - 1)

    % Calculate Optimal Length of Overlapping FFT Window (L)
    nBitsSingle = 32; % number of bits for single-precision data (audio)
    windowBytes = 52428800; % 50 MB window
    nOverlap = N - 1; % window overlap (nOverlap ~ P)
    Lt = nSamples*P + nOverlap; % length of upsampled audio file (with overlap)
    M = round(windowBytes*8/nBitsSingle * 1/P)*P; % length of upsampled audio chunk (without overlap)
    L = M + nOverlap; % length of upsampled audio chunk (with overlap)
    
    % Limit the Maximum Length of the Upsampled Window to Lt
    if L > Lt 
        L = Lt;
        M = L - nOverlap;
    end
    
    % Make M >= nOverlap to ensure only the first window is zero-padded
    if M < nOverlap 
        error(['The filter order exceeds the maximum accepted for '...
            'efficient processing']);
    end

    % Filtering, Overlap-Save and Resampling Parameters
    nWindows = ceil(nSamples*P/M); % number of windows
    zeropad = nWindows*M - nSamples*P; % length of zero-padding
    phase = 0; % downsampling phase (variable)
    sr1 = 1; % starting sample for the resampled audio vector xr

    % Calculate the Fourier Transform of the Filter
    H = conj(fft(b,L)); % conjugate FFT of FIR filter b padded with M-1 zeros

    % NOTE: in theory you should use H = FFT(b,L) when convolving, since
    % the complex conjugate is meant for correlation. However, using the
    % conjugate here is a simple trick to shift the distortion effects of
    % the FIR filter to the end of the recovered segment for them to be
    % discarded. In this way, the distorted portion of the segment is
    % confined to the last samples of the filtered signal.
    
    % Initialise overlap section with mean value of first samples
    xu_overlap = zeros(nOverlap,1); % overlapping section (initialised with the average value)

    % RESAMPLE SELECTED AUDIO (RAW and WAV)
    hfig = waitbar(0,'','Name','Resampling'); % initialise waitbar
    for m = 1:nWindows-1
        % Read audio chunk
        iChunkSample1 = (m-1)*M/P + firstSample; % index of first sample to be read
        iChunkSample2 = m*M/P + firstSample - 1; % index of last sample to be read
        nChunkSamples = iChunkSample2 - iChunkSample1 + 1; % number of samples to read
        x =  readAudioFile(filePath,channel,[iChunkSample1 nChunkSamples],...
            'float',fs,nch,nbit,en);
        
        % Upsample
        xu = zeros(length(x)*P,1);
        xu(1:P:end) = x * P;
        xu = [xu_overlap ; xu];  %#ok<AGROW>
        xu_overlap = xu(M+1:L); % overlap to be used in next audio section
        clear x
        
        % Filter
        Xu = fft(xu);
        xf = ifft(Xu.*H);
        clear xu Xu
        
        % Remove Initial Overlap (zero-pad) from First Window
        if m == 1
            xf(1:round(nOverlap/2)) = [];
        end
        
        % Downsample
        xd = downsample(xf(1:end - nOverlap),Q,phase);
        sr2 = sr1 + length(xd) - 1;
        xr(sr1:sr2) = xd;
        clear xf xd
        
        % Parameters for next iteration
        phase = ceil((M - phase)/Q)*Q - (M - phase);
        sr1 = sr2 + 1;
        
        % Display processing status
        waitbarString = sprintf('Window %d/%d',m,nWindows);
        waitbar(m/nWindows,hfig,waitbarString); % show processing progress
    end

    % Read last audio chunk
    iChunkSample1 = (nWindows-1)*M/P + firstSample; % index of first sample to be read
    iChunkSample2 = firstSample + nSamples - 1; % index of last sample to be read
    nChunkSamples = iChunkSample2 - iChunkSample1 + 1; % number of samples to read
    x =  readAudioFile(filePath,channel,[iChunkSample1 nChunkSamples],...
        'float',fs,nch,nbit,en);
    
    % Upsample
    xu = zeros(length(x)*P,1);
    xu(1:P:end) = x * P;
    xu = [xu_overlap ; xu ; zeros(zeropad,1)];
    clear x
    
    % Filter
    Xu = fft(xu);
    xf = ifft(Xu.*H);
    clear xu Xu
    
    % Downsample
    xd = downsample(xf(1:end - zeropad),Q,phase);
    sr2 = sr1 + length(xd) - 1;
    xr(sr1:sr2) = xd;
    clear xf xd
    
    % Display processing status
    waitbarString = sprintf('Window %d/%d',nWindows,nWindows);
    waitbar(1,hfig,waitbarString); % show processing progress
    delete(hfig) % delete waitbar
    
    % Limit Signal Amplitude
    xrmax = 2^(nbit-1) - 1;
    xrmin = -2^(nbit-1);
    xr(xr > xrmax) = xrmax;
    xr(xr < xrmin) = xrmin;
    
else % No resampling (read audio file directly)
    xr =  readAudioFile(filePath,channel,samples,'float',fs,nch,nbit,en);
end
