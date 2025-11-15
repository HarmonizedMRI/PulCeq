function thetaVals = readThetaRegisters(filename)
% readThetaRegisters - Extract THETA register values (OPC 406)
%   from a text dump containing multiple r_*_frq sections.
%
% Usage:
%   thetaVals = readThetaRegisters('regs.txt')

    % Read file as text
    fid = fopen(filename, 'r');
    if fid < 0
        error('Cannot open file: %s', filename);
    end
    txt = textscan(fid, '%s', 'Delimiter', '\n', 'Whitespace', '');
    fclose(fid);
    lines = txt{1};
    
    thetaVals = struct('name', {}, 'theta', {});
    
    i = 1;
    currentName = '';
    while i <= numel(lines)
        line = strtrim(lines{i});
        
        % Detect NAME line
        if startsWith(line, 'NAME:')
            currentName = strtrim(extractAfter(line, 'NAME:'));
        end
        
        % Detect THETA register
        if contains(line, 'OPC') && contains(line, '406')
            % Next two lines are the data bytes
            if i+2 <= numel(lines)
                dat1 = extractByteValue(lines{i+1});
                dat0 = extractByteValue(lines{i+2});
                % Combine bytes: assume 16-bit value (high byte first)
                theta = bitor(bitshift(dat1, 8), dat0);
                theta = bitor(bitshift(dat1, 8), dat0);
                
                thetaVals(end+1).name = currentName; %#ok<AGROW>
                thetaVals(end).theta = theta;
            end
        end
        i = i + 1;
    end
end

function val = extractByteValue(line)
% Extract integer value from a line like:
%   DAT       798      DATA BYTE 1 (0xCC)  (=52224)
%   or        DAT 600 DATA BYTE 0 (0x 0) (=0)
    expr = '=\s*(\d+)';
    tokens = regexp(line, expr, 'tokens');
    if ~isempty(tokens)
        val = str2double(tokens{1}{1});
    else
        val = NaN;
    end
end

