function DAx = loadPoly(fname, n, m, type)
% loadPoly  Load a polynomial map file written by the COSY/DA propagator.
%
%   DAx = loadPoly(fname, n, m, type)
%
%   Inputs:
%     fname - path to the polynomial map file (string)
%     n     - number of DA variables (polynomial dimension)
%     m     - number of polynomials per batch (e.g. 6 for a state vector)
%     type  - 0: pure DA polynomials
%             1: Taylor models (each polynomial also carries remainder bounds)
%
% Outputs: %     DAx   - column array of structs of length (total_batches * m), each with:
%               .C    - coefficient column vector (one entry per monomial)
%               .E    - exponent matrix (K x n, one row per monomial)
%             additionally, when type == 1:
%               .Rinf - lower remainder bound
%               .Rsup - upper remainder bound
%
%   File format:
%     Each polynomial block starts with a line beginning with 'I' or 'A'.
%     Each coefficient line contains: [index coeff ... exponents...].
%     A line beginning with '-' marks the end of a block.
%     When n == 16, coefficient data spans two consecutive lines.
%     When n > 16, same two-line layout as n == 16.

file = char(textread(fname, '%s', 'delimiter', '\n')); %#ok<DTXTRD>

% Find the start positions of polynomial blocks
pos = find((file(:,1) == 'I') | (file(:,1) == 'A'));
pos = pos(1:m);

DAx = [];

if n < 16
    % --- Standard single-line format ---
    for i = 1:m:length(pos)
        for j = 0:(m-1)
            k = pos(i+j) + 1;
            DAx_temp(j+1).C = [];  %#ok<AGROW>
            DAx_temp(j+1).E = [];  %#ok<AGROW>

            % Read monomial entries until block-end marker '-'
            while file(k,1) ~= '-'
                temp = str2num(file(k,:)); %#ok<ST2NM>
                DAx_temp(j+1).C = [DAx_temp(j+1).C; temp(2)];
                DAx_temp(j+1).E = [DAx_temp(j+1).E; temp(4:(4+n-1))];
                k = k + 1;
            end

            % Represent the zero polynomial explicitly
            if isempty(DAx_temp(j+1).C)
                DAx_temp(j+1).C = 0;
                DAx_temp(j+1).E = zeros(1, n);
            end

            % Taylor model: read remainder interval [ Rinf, Rsup ]
            if type
                temp  = file(k+2, :);
                pos1  = find(temp == '[');
                pos2  = find(temp == ',');
                pos3  = find(temp == ']');
                DAx_temp(j+1).Rinf = str2double(temp(pos1+1:pos2-1));
                DAx_temp(j+1).Rsup = str2double(temp(pos2+1:pos3-1));
            end
        end
        DAx = [DAx; DAx_temp]; %#ok<AGROW>
    end

elseif n == 16
    % --- Two-line format for n = 16 ---
    for i = 1:m:length(pos)
        for j = 0:(m-1)
            k = pos(i+j) + 1;
            DAx_temp(j+1).C = [];  %#ok<AGROW>
            DAx_temp(j+1).E = [];  %#ok<AGROW>

            while file(k,1) ~= '-'
                temp = str2num(file(k,:)); %#ok<ST2NM>
                DAx_temp(j+1).C = [DAx_temp(j+1).C; temp(2)];
                DAx_temp(j+1).E = [DAx_temp(j+1).E; temp(4:(4+n-1))];
                k = k + 2;  % advance by 2 lines per monomial
            end

            if isempty(DAx_temp(j+1).C)
                DAx_temp(j+1).C = 0;
                DAx_temp(j+1).E = zeros(1, n);
            end

            if type
                temp  = file(k+2, :);
                pos1  = find(temp == '[');
                pos2  = find(temp == ',');
                pos3  = find(temp == ']');
                DAx_temp(j+1).Rinf = str2double(temp(pos1+1:pos2-1));
                DAx_temp(j+1).Rsup = str2double(temp(pos2+1:pos3-1));
            end
        end
        DAx = [DAx; DAx_temp]; %#ok<AGROW>
    end

else
    % --- Two-line format for n > 16: exponents split across two lines ---
    for i = 1:m:length(pos)
        for j = 0:(m-1)
            k = pos(i+j) + 1;
            DAx_temp(j+1).C = [];  %#ok<AGROW>
            DAx_temp(j+1).E = [];  %#ok<AGROW>

            while file(k,1) ~= '-'
                temp1 = str2num(file(k,:));   %#ok<ST2NM>
                temp2 = str2num(file(k+1,:)); %#ok<ST2NM>
                DAx_temp(j+1).C = [DAx_temp(j+1).C; temp1(2)];
                % Exponents span the tail of line 1 and the head of line 2
                DAx_temp(j+1).E = [DAx_temp(j+1).E; temp1(4:end-1) temp2(1:end-1)];
                k = k + 2;
            end

            if isempty(DAx_temp(j+1).C)
                DAx_temp(j+1).C = 0;
                DAx_temp(j+1).E = zeros(1, n);
            end

            if type
                temp  = file(k+2, :);
                pos1  = find(temp == '[');
                pos2  = find(temp == ',');
                pos3  = find(temp == ']');
                DAx_temp(j+1).Rinf = str2double(temp(pos1+1:pos2-1));
                DAx_temp(j+1).Rsup = str2double(temp(pos2+1:pos3-1));
            end
        end
        DAx = [DAx; DAx_temp]; %#ok<AGROW>
    end
end
