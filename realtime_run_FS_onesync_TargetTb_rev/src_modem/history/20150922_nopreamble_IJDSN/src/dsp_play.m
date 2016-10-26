function [obj]=dsp_play(obj, fs, len, stat, frame)

if strcmp(stat, 'start')
    %Initialization of recording obj.
    obj = dsp.AudioPlayer('SampleRate', fs, ...
        'DeviceDataType', '16-bit integer', ...
        'BufferSizeSource', 'Property', ...
        'BufferSize', len, ...
        'QueueDuration', 0);
    disp('Speaker playing initialized...');
    
    frame = 0;
end

if strcmp(stat, 'stop')
    release(obj);
    disp('Speaker playing stopped...');
end

if strcmp(stat, 'record')
    frame_int = cast(frame', 'int16');
    step(obj, frame_int);
end