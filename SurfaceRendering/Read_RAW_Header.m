% Read_RAW_Header
% Returns the image header information and some extra information
% for any raw image file.
% Usage: stHeader = Read_RAW_Header(filename)
% Input
%   filename - RAW format file name
% Output
%   stHeader - header structure containing information on the data
%
% Example usage:
% fullFileName_RAW = fullfile(cd, 'SimulatedData_Distance.VOL')
% stHeader = Read_RAW_Header(fullFileName_RAW)
%
%------------------------------------------------------------------------------
%example: Read_RAW_3DArray('C:/Users/Zahra/Dropbox/Shuang/SurfaceRendering_matlab/gabardine_od_supertile.vol')
function stHeader = Read_RAW_Header(fullFilename)
	% Check syntax.  Must have at least one input argument, the full filename. 
	if (nargin ~= 1)
		error('Usage: stHeader = Read_RAW_Header(fullFilename)');
	end 

	if (ischar(fullFilename)~=1)
		error('Requires a string filename as an argument.');
    end 
    
	stHeader.Endian = 'Little-endian';
	stHeader.EndianArg = 'ieee-le';
    
	% Open the file for reading.
	fileHandleID = fopen(fullFilename, 'rb', stHeader.EndianArg);
	if (fileHandleID == -1)
		error(['Error opening ', fullFilename, ' for input.']);
	end

	% Go to the beginning of the file.
	% Shouldn't be necessary, but can't hurt.
	fseek(fileHandleID, 0, 'bof');

	% Read the total number of voxels in the image.
	% Read bytes 1-3. ASCII Bytes ’V’, ’O’, and ’L’
	V = fread(fileHandleID, 1, 'char*1');
    O = fread(fileHandleID, 1, 'char*1');
    L = fread(fileHandleID, 1, 'char*1');
    stHeader.Format = 'VOL';
    
    %file format version number
    stHeader.version = fread(fileHandleID, 1, '*ubit8');
       
  	% The next 4 bytes are the Encoding identifier 
	% Read bytes 5-8.: 1-Dense float32-based
	stHeader.representation = fread(fileHandleID, 1, '*int32');

	% Read in the dimensions for the different directions.
	% They'll be in bytes 9-20.
	stHeader.x_size = fread(fileHandleID, 1, '*int32');
	stHeader.y_size = fread(fileHandleID, 1, '*int32');
	stHeader.z_size = fread(fileHandleID, 1, '*int32');
 	stHeader.TotalVoxels = stHeader.x_size * stHeader.y_size * stHeader.z_size;
    
    %Read number of channels:
    stHeader.channels = fread(fileHandleID, 1, '*ubit32');
    
    %Read the bounding box
    stHeader.xmin = fread(fileHandleID, 1, '*single');
    stHeader.ymin = fread(fileHandleID, 1, '*single');
    stHeader.zmin = fread(fileHandleID, 1, '*single');
    stHeader.xmax = fread(fileHandleID, 1, '*single');
    stHeader.ymax = fread(fileHandleID, 1, '*single');
    stHeader.zmax = fread(fileHandleID, 1, '*single');

	% Assign some other useful information.
	fileInfo = dir(fullFilename); 
	% MATLAB returns the information in a structure with fields:
	%    name
	%    date
	%    bytes
	%    isdir
	%    datenum
	stHeader.HeaderSize = 48;
	fileSizeInBytes = fileInfo.bytes;
	dataSizeInBytes = double(fileSizeInBytes) - double(stHeader.HeaderSize);
	stHeader.BytesPerVoxel = int32(round(dataSizeInBytes / double(stHeader.TotalVoxels)));
	stHeader.FileSizeInBytes = fileSizeInBytes;
	stHeader.DataSizeInBytes = dataSizeInBytes;

	% Close the file.
	fclose(fileHandleID);	% Close the file.