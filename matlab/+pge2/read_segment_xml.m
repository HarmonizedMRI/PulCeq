function seq = read_pulse_sequence_xml(filename)
% READ_PULSE_SEQUENCE_XML Reads an ICE Pulse Sequence XML file.
%
%   seq = READ_PULSE_SEQUENCE_XML(filename)
%
%   Returns a struct array where each element corresponds to a <sequencer>
%   block in the XML file, with fields:
%       - id
%       - title
%       - xtitle
%       - ytitle
%       - waveform_name
%       - time  (vector)
%       - value (vector)
%
%   Example:
%       seq = read_pulse_sequence_xml('scan.xml.0001');

    if nargin < 1
        [file, path] = uigetfile('*.xml.*', 'Select an XML file');
        if isequal(file, 0)
            error('No file selected.');
        end
        filename = fullfile(path, file);
    end

    % Read the XML file into a DOM
    try
        xDoc = xmlread(filename);
    catch ME
        error('Failed to read XML file %s.\nError: %s', filename, ME.message);
    end

    % Get all <sequencer> elements
    seqNodes = xDoc.getElementsByTagName('sequencer');
    nSeq = seqNodes.getLength;
    seq = repmat(struct(), 1, nSeq);

    for i = 1:nSeq
        thisSeq = seqNodes.item(i-1);

        % Get sequencer attributes
        seq(i).id     = char(thisSeq.getAttribute('id'));
        seq(i).title  = char(thisSeq.getAttribute('title'));
        seq(i).xtitle = char(thisSeq.getAttribute('xtitle'));
        seq(i).ytitle = char(thisSeq.getAttribute('ytitle'));

        % Get the <data> child node
        dataNodes = thisSeq.getElementsByTagName('data');
        if dataNodes.getLength > 0
            dataNode = dataNodes.item(0);
            seq(i).waveform_name = char(dataNode.getAttribute('waveform'));

            % Extract text content and parse into numeric arrays
            txt = char(dataNode.getTextContent);
            lines = strtrim(splitlines(txt));
            lines = lines(~cellfun('isempty', lines));

            % Convert each line to numeric [time, value]
            data = cellfun(@(L) sscanf(L, '%f'), lines, 'UniformOutput', false);
            for ii = 1:length(data)
                data{ii} = data{ii}';
            end
            data = vertcat(data{:});

            seq(i).time = data(:,1);
            seq(i).value = data(:,2);
        else
            seq(i).waveform_name = '';
            seq(i).time = [];
            seq(i).value = [];
        end
    end
end

