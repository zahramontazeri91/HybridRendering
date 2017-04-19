%  Read_RAW_3DArray
%  Reads a RAW file from disk into a 3D array with control over cropping and subsampling.
%
%  Usage:
%	[data3D, stHeader] = Read_RAW_3DArray(fullFileName, stInputParameters);
%
%  Inputs:
%    fullFileName - full path and filename of the RAW file.
%    stInputParameters - a MATLAB structure that is dynamically created.
%                        defaults get applied for the missing structure members.
%                        See example below.
%  Output:
%    data3D - 3-D array of voxel values, either uint8 or uint16.
%
%  Example usage:
%	Create a MATLAB structure dynamically to hold the parameters.
%	stInputParameters.XStart = 1;
%	stInputParameters.XEnd = 1024;
%	stInputParameters.YStart = 1;
%	stInputParameters.YEnd = 1024;
%	stInputParameters.ZStart = 10;
%	stInputParameters.ZEnd = 50;
%	stInputParameters.Subsample = 4;
%	fullFileName_RAW = 'SimulatedData_Distance8Bit.R
%	[data3D, stHeader] = Read_RAW_3DArray(fullFileName_RAW, stInputParameters);
%
%------------------------------------------------------------------------------
 
%example: Read_RAW_3DArray('C:/Users/Zahra/Dropbox/Shuang/SurfaceRendering_matlab/gabardine_od_supertile.vol')

function [data3D, stHeader] = Read_RAW_3DArray(fullFileName)
    z_min = 1;
    z_max = 170;
	% Check syntax.  Must have at least one input argument, the full filename. 
	if (nargin ~= 1)
		error('Usage: stHeader = Read_RAW_Header(fullFilename)');
	end 

	if (ischar(fullFileName)~=1)
		error('Requires a string filename as an argument.');
    end 

	% Read in header and get 3D image array dimensions.
	stHeader = Read_RAW_Header(fullFileName);
	% Extract out sizes to more conveniently-named local variables.
	% Note: fread() requires that x_size and y_size be doubles.
	x_size = double(stHeader.x_size);
	y_size = double(stHeader.y_size);
	z_size = double(stHeader.z_size);

	% Open the image data file for reading.
	fileHandle = fopen(fullFileName, 'rb', stHeader.EndianArg);
	if (fileHandle == -1)
		error(['Read_RAW_3DArray() reports error opening ', fullFileName, ' for input.']);
	end

	% Skip past header of stHeader.HeaderSize bytes.
	bytesToSkip = int32(stHeader.HeaderSize);
	fseek(fileHandle, bytesToSkip, 'bof');

    %Read in data in one 1D Vector
    vol_1D = fread(fileHandle,'*float32');

    vol_3D = zeros(x_size, y_size, z_size);
    hdr = zeros(x_size,y_size);
    
    for z = [z_min:z_max]%according to the cropped z value obtained using cropZ.m
            for y = [1:y_size]
                for x = [1:x_size]
                    vol_3D(x,y,z) = vol_1D((z-1)*(x_size*y_size) + (y-1)*x_size + x);
                end            
            end
    end
    
%     %display the slices from the side view:
%     temp = zeros(x_size,z_size);
%     for i = 1:y_size
%         for j = 1:z_size
%             temp(:,j) = vol_3D(:,i,j);
%         end
%         imwrite(temp,['Vertical Slices/slice_',num2str(i),'.png']);
%     end


%%%Capture the bottom surface (the dataset captured upsidedown, so this
%is the top side of the sample -> having diagonal patern)
    for x = [1:x_size]
            for y = [1:y_size]
                %for z = [z_size:-1:1] 
                for z = [z_max:-1:z_min] %according to the cropped z value obtained using cropZ.m
                    %find i s.t. SigmaVj(from i to n) > 0.95 Sigma Vj(from 1 to n)
                    if sum(vol_3D(x,y,1:z)) > 0.95*(sum(vol_3D(x,y,:)))
                       hdr(x,y)= z;
                       %disp(z);                 
                    end
                end            
            end
    end


%%%Capture the top surface
%     for x = [1:x_size]
%             for y = [1:y_size]
%                 for z = [1:z_size]
%                     %find i s.t. SigmaVj(from i to n) > 0.95 Sigma Vj(from 1 to n)
%                     if sum(vol_3D(x,y,z:z_size)) > 0.95*(sum(vol_3D(x,y,:)))
%                        hdr(x,y)= z_size-z;
%                        %disp(z);                 
%                     end
% 
%                 end            
%             end
%     end    


      exrwrite( hdr, 'input/height.exr' );
      imagesc(hdr);
      colorbar;
      axis on;
      grid on;
    
	% Close the file.
	fclose(fileHandle);