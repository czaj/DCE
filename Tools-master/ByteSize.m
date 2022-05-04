function b = ByteSize(theVariable, returnType, fid)
% getByteSize returns the mem.usage of the provided variable(theVariable) to the given file identifier. 
% returnType is assigned meaningfully according to the byte size if not stated
% Output is written to screen if fid is 1, empty or not provided.
s = whos('theVariable');
b = s.bytes;
if nargin == 1 || isempty(returnType)
    scale = floor(log(b)/log(1024));
    switch scale
        case 0
            returnType = 'byte';
        case 1
            returnType = 'kb';
        case 2
            returnType = 'mb';
        case 3
            returnType = 'gb';
        case 4
            returnType = 'tb';
        case -inf
            % Size occasionally returned as zero (eg some Java objects).
            returnType = 'byte';
            warning('Size occasionally returned as zero (eg some Java objects). Bytes assumed');
        otherwise
            returnType = 'petabytes';
            warning('Over 1024 petabyte. petabytes assumed');
    end
end
switch returnType
    case {'b','byte','bytes'}
        b = s.bytes;
    case {'kb','kbs','kilobyte','kilobytes'}
        b = b / 1024;
    case {'mb','mbs','megabyte','megabytes'}
        b = b / 1024^2;
    case {'gb','gbs','gigabyte','gigabytes'}
        b = b / 1024^3;
    case {'tb','tbs','terabyte','terabytes'}
        b = b / 1024^4;
    case {'pb','pbs','petabyte','petabytes'}
        b = b / 1024^5;
    otherwise
        returnType = 'bytes';
end
if nargin <= 2 || isempty(fid) || fid == 1
    fprintf(1,[num2str(b) ' ' returnType '\n']);
elseif nargin > 2 && ~isempty(fid) && fid > 2
    try
        fprintf(fid,[num2str(b) ' ' returnType '\n']);
    catch
        warning(['fid(' num2str(fid) ') could not be edited. Hence the output will be written on the screen.']);
        fprintf(1,[num2str(b) ' ' returnType '\n']);
    end
end
end