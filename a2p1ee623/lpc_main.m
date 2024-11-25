% Main Script
% Parameters
Fs = 16000;            % Sampling rate (16 kHz for wideband speech)
frame_duration = 0.03; % Frame size (30 ms)
frame_shift = 0.02;    % Frame shift (20 ms overlap)
order = 12;            % LPC order
dct_order = 50;        % Number of DCT coefficients

% Load speech signal
[x, Fs] = audioread('ht.wav');%recorded sample audiofile
x = x / max(abs(x)); % Normalize signal

% Frame processing parameters
frame_size = round(frame_duration * Fs);
frame_step = round(frame_shift * Fs);
num_frames = floor((length(x) - frame_size) / frame_step) + 1;

% Pre-compute global energy threshold
global_threshold = mean(x.^2) * 0.1; % Dynamic energy threshold

% Initialize synthesized signals
synthesized_plain = zeros(size(x));
synthesized_voice_excited = zeros(size(x));

% Hamming window
window = hamming(frame_size);

for i = 1:num_frames
    % Extract current frame
    frame_start = (i-1) * frame_step + 1;
    frame_end = frame_start + frame_size - 1;
    frame = x(frame_start:frame_end) .* window;
    
    % LPC Analysis
    [a, G, residual] = lpc_analysis(frame, order);
    
    % Voiced/Unvoiced Detection
    if sum(frame.^2) > global_threshold
        % Voiced: Generate impulse train
        pitch_period = Fs / 100; % Approximate pitch period (10 ms)
        excitation = zeros(size(frame));
        excitation(1:round(pitch_period):end) = 1;
    else
        % Unvoiced: Generate white noise
        excitation = randn(size(frame));
    end
    
    % Plain LPC Synthesis
    synthesized_frame_plain = lpc_synthesis(excitation, a, G);
    synthesized_plain(frame_start:frame_end) = synthesized_plain(frame_start:frame_end) + synthesized_frame_plain;

    % Compress residual using DCT
    dct_residual = dct(residual);
    dct_residual_truncated = dct_residual(1:dct_order);

    % Decompress residual
    decompressed_residual = idct([dct_residual_truncated; zeros(length(residual) - dct_order, 1)]);

    % Voice-Excited LPC Synthesis
    synthesized_frame_voice = lpc_synthesis(decompressed_residual, a, G);
    synthesized_voice_excited(frame_start:frame_end) = synthesized_voice_excited(frame_start:frame_end) + synthesized_frame_voice;
end

% Normalize synthesized signals
synthesized_plain = synthesized_plain / max(abs(synthesized_plain));
synthesized_voice_excited = synthesized_voice_excited / max(abs(synthesized_voice_excited));

% Save synthesized signals
audiowrite('synthesized_plain.wav', synthesized_plain, Fs);
audiowrite('synthesized_voice_excited.wav', synthesized_voice_excited, Fs);

% Load Original and Synthesized Signals
original = x; % Use the original loaded speech signal
synth_plain = synthesized_plain;
synth_voice_excited = synthesized_voice_excited;

% Calculate Segmental SNR
segSNR_plain = segsnr(original, synth_plain, Fs, 'FrameDuration', 0.02);
segSNR_voice_excited = segsnr(original, synth_voice_excited, Fs, 'FrameDuration', 0.02);

% Display Results
fprintf('Segmental SNR (Plain LPC): %.2f dB\n', segSNR_plain);
fprintf('Segmental SNR (Voice-excited LPC): %.2f dB\n', segSNR_voice_excited);

% Measure runtime for Plain LPC Synthesis
tic;
for i = 1:num_frames
    frame_start = (i-1) * frame_step + 1;
    frame_end = frame_start + frame_size - 1;
    frame = x(frame_start:frame_end) .* window;
    
    % LPC Analysis
    [a, G, ~] = lpc_analysis(frame, order);
    
    % Voiced/Unvoiced Detection
    if sum(frame.^2) > global_threshold
        % Voiced: Generate impulse train
        pitch_period = Fs / 100; % Approximate pitch period (10 ms)
        excitation = zeros(size(frame));
        excitation(1:round(pitch_period):end) = 1;
    else
        % Unvoiced: Generate white noise
        excitation = randn(size(frame));
    end
    
    % Plain LPC Synthesis
    synthesized_frame_plain = lpc_synthesis(excitation, a, G);
    synthesized_plain(frame_start:frame_end) = synthesized_plain(frame_start:frame_end) + synthesized_frame_plain;
end
runtime_plain = toc;

% Measure runtime for Voice-Excited LPC Synthesis
tic;
for i = 1:num_frames
    frame_start = (i-1) * frame_step + 1;
    frame_end = frame_start + frame_size - 1;
    frame = x(frame_start:frame_end) .* window;
    
    % LPC Analysis
    [a, G, residual] = lpc_analysis(frame, order);
    
    % Compress residual using DCT
    dct_residual = dct(residual);
    dct_residual_truncated = dct_residual(1:dct_order);

    % Decompress residual
    decompressed_residual = idct([dct_residual_truncated; zeros(length(residual) - dct_order, 1)]);

    % Voice-Excited LPC Synthesis
    synthesized_frame_voice = lpc_synthesis(decompressed_residual, a, G);
    synthesized_voice_excited(frame_start:frame_end) = synthesized_voice_excited(frame_start:frame_end) + synthesized_frame_voice;
end
runtime_voice_excited = toc;

% Display runtimes
fprintf('Runtime (Plain LPC): %.3f seconds\n', runtime_plain);
fprintf('Runtime (Voice-Excited LPC): %.3f seconds\n', runtime_voice_excited);


% Local Function Definitions
function [a, G, residual] = lpc_analysis(frame, order)
    a = lpc(frame, order);
    residual = filter(a, 1, frame);
    G = sqrt(sum(residual.^2) / length(residual));
end

function synthesized = lpc_synthesis(residual, a, G)
    synthesized = filter(1, a, residual) * G;
end

function segSNR = segsnr(original, synthesized, Fs, varargin)
    p = inputParser;
    addParameter(p, 'FrameDuration', 0.02);
    parse(p, varargin{:});
    frame_duration = p.Results.FrameDuration;

    frame_size = round(frame_duration * Fs);
    num_frames = floor(length(original) / frame_size);

    segSNR_frames = zeros(num_frames, 1);
    for i = 1:num_frames
        start_idx = (i-1)*frame_size + 1;
        end_idx = start_idx + frame_size - 1;
        if end_idx > length(original)
            break;
        end
        frame_orig = original(start_idx:end_idx);
        frame_synth = synthesized(start_idx:end_idx);
        if sum(frame_orig.^2) > 1e-6
            segSNR_frames(i) = 10 * log10(sum(frame_orig.^2) / sum((frame_orig - frame_synth).^2));
        end
    end
    segSNR = mean(segSNR_frames, 'omitnan');
end
