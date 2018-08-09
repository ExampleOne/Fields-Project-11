function [data, labels] = readTAC(path)
%READTAC Reads a TAC file, and takes the input into a "data" variable.
%Labels labels the different brain regions.
%   Here data(:, 1) is the starttime of each frame, data(:, 2) is the end
%   time of each frame, and data(ii + 2) is the TAC of the iith brain
%   region (ie. the brain region labelled as labels{ii}). Here labels is a
%   cell array.
    fullData = importdata(path, '\t');
    data = [fullData.data];
    labels = fullData.textdata(1:end);
end

