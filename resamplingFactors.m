%  [P,Q] = RESAMPLINGFACTORS(fs,fr,tolerance,dataSize)
%
%  DESCRIPTION: Calculates the upsampling and downsampling factors P and
%  Q required to convert a signal from its original sampling rate FS
%  to a new sampling rate FR = FS*P/Q. 
%
%  Sometimes, an exact sampling rate conversion is not possible and FS*P/Q 
%  will differ from FR. A non-exact conversion may is neccesary due to 
%  limitations in the available memory. Higher memory usage is associated
%  with larger values of P (upsampling factor) and size of the signal
%  to be resampled.
% 
%  The input parameters DATASIZE and TOLERANCE are used to control the 
%  conversion accuracy.
%
%  RESAMPLINGFACTORS can be used to estimate the optimal values for the
%  input upsampling and downsampling factors (P,Q) in RESAMPLEAUDIOFILE.
% 
%  INPUT ARGUMENTS
%  - fs [integer number]: sampling frequency [Hz]
%  - fr [integer number]: resampling frequency [Hz]
%  - tolerance [number]: error tolerance for the resampling frequency FR
%    [%]. P and Q will be selected to keep FS*P/Q within the tolerance limits
%    (i.e. (1-TOLERANCE/100)*FR < FS*P/Q < (1 + TOLERANCE/100)*FR).
%  - dataSize [number]: size of the signal to be resampled [bytes]. It is
%    recommended to resample the signal in small segments to limit the
%    memory usage and increase the processing speed. For example, for a
%    segment of a single-precision signal of TAU = 1 second duration and 
%    sampled at FS = 96000 Hz, DATASIZE = FS*4 = 384000 Bytes.
%
%  OUTPUT ARGUMENTS
%  - P [number]: upsampling factor
%  - Q [number]: downsampling factor
%
%  INTERNALLY CALLED FUNCTIONS
%  - readwavHeader
%  - readAudioFile
%
%  FUNCTION CALLS
%  [P,Q] = RESAMPLINGFACTORS(fs,fr,tolerance,dataSize)
%
%  See also RESAMPLEAUDIOFILE

function [P,Q] = resamplingFactors(fs,fr,tolerance,dataSize)

% Estimate of maximum upsampling factor PMAX
[~,sysmem] = memory;
freeRam = sysmem.PhysicalMemory.Available / 2; % RAM free for processing (taken as 1/2 * available RAM) [bytes]
Pmax = round(freeRam/dataSize - 1); % maximum upsampling factor

% Calculate Initial Upsampling and Downsampling Factors (P,Q)
k = gcd(fs,fr); % greatest common divisor of FS and FR
P0 = fr/k; % upsampling factor (temporal)
Q0 = fs/k; % downsampling factor (temporal)
Q_factors = factor(Q0);

% Initialise Upsampling and Downsampling Factors (P,Q)
P = P0; % upsampling factor
Q = Q0; % downsampling factor

% Compute Conditions for Loop Continuity (Memory, Tolerance, and Factor)
memoryFlag = P0 > Pmax; % RAM memory limit
toleranceFlag = (fs*P0/Q0 >= fr*(1-tolerance/100)) ...
    && (fs*P0/Q0 <= fr*(1+tolerance/100)); % resampling rate tolerance
factorFlag = length(Q_factors) > 1; % factorisation limit

% Adjust Upsampling and Downsampling Factors (P,Q)
if (toleranceFlag || memoryFlag) && factorFlag
    % Divide(P,Q) by a Factor that Reasonably Keeps the Ratio P/Q
    [~,ind] = min(abs(round(P0./Q_factors).*Q_factors - P0));
    P0 = round(P0./Q_factors(ind));
    Q0 = Q0/Q_factors(ind);
    
    % Calculate Upsampling and Downsampling Factors (P,Q)
    k = gcd(P0,Q0);
    P0 = P0/k;
    Q0 = Q0/k;
    Q_factors = factor(Q0);
    
    % Compute Conditions for Loop Continuity (Memory, Tolerance, and Factor)
    memoryFlag = P0 > Pmax; % RAM memory limit
    toleranceFlag = (fs*P0/Q0 >= fr*(1-tolerance/100)) ...
        || (fs*P0/Q0 <= fr*(1+tolerance/100)); % resampling rate tolerance
    factorFlag = length(Q_factors) > 1; % factorisation limit

    while (toleranceFlag || memoryFlag) && factorFlag
        % Update Upsampling and Downsampling Factors (P,Q)
        P = P0;
        Q = Q0;
        
        % Divide(P,Q) by a Factor that Reasonably Keeps the Ratio P/Q
        [~,ind] = min(abs(round(P0./Q_factors).*Q_factors - P0));
        P0 = round(P0./Q_factors(ind));
        Q0 = Q0/Q_factors(ind);
        
        % Calculate Upsampling and Downsampling Factors (P,Q)
        k = gcd(P0,Q0);
        P0 = P0/k;
        Q0 = Q0/k;
        Q_factors = factor(Q0);
        
        % Compute Conditions for Loop Continuity (Memory, Tolerance, and Factor)
        memoryFlag = P0 > Pmax; % RAM memory limit
        toleranceFlag = (fs*P0/Q0 > fr*(1-tolerance/100)) ...
            && (fs*P0/Q0 < fr*(1+tolerance/100)); % resampling rate tolerance
        factorFlag = length(Q_factors) > 1; % factorisation limit
    end
end
