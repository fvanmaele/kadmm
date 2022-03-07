function [A, b, u_orig] = deblurringLoadData( dataset, instance, color, filterType, filterSize )
%[A, b, u_orig] = deblurringLoadData( dataset, instance, color, filterType, filterSize )
%
%   INPUT:
%   -------------------------------------------------------------
%   dataset           {'BSDS500'}
%   instance          the instance number in the dataset
%   color             true | false = grayscale
%   filterType        {'Gaussian'}, {'Box'}, {'Motion'}
%   filterSize        must be an odd number
%
%   OUTPUT:
%   --------------------------------------------------------------
%   A                 linear operator   NxN
%   b                 observed data     NxC
%   u_orig            original image    WxHxC

    % Default value for filterSize
    if nargin < 5
        filterSize = 5;
    end
    
    % Default value for filterType
    if nargin < 4
        filterType = 'Gaussian';
    end
    
    % Default value for color
    if nargin < 3
        color = false;
    end 
    
    % Default value for instance
    if nargin < 2
        instance = 1;
    end
    
    % Default value for dataset
    if nargin < 1
        dataset = 'BSDS500';
    end
    
    if filterSize < 1 || mod(filterSize,2) == 0
        error('Error: filterSize must be an positive odd number' )
    end

    %% Load Data
    mfn  = mfilename;
    mffn = mfilename('fullpath');
    datapath = [mffn(1:end-numel(mfn)),'/data/',dataset,'/'];
    files    = dir([datapath,'*.jpg*']);
    if size(files,1) < 1
        error(['Error: no image data found for dataset ' dataset])
    end
    if instance < 1 || instance > size(files,1)
        error(['Error: instance must be in the range of [1,' num2str(size(files,1)) '] for dataset ' dataset])
    end
    
    u_orig   = im2double(imread([datapath,files(instance).name]));
    
    % convert to grayscale
    if color == false
        u_orig = rgb2gray(u_orig);
    end 
    
    height = size(u_orig, 1);
    width = size(u_orig, 2);
    pixels = height*width;    
    channels = size(u_orig,3);    
    %% Generate filter mask
    switch filterType
      case 'Gaussian'
        var = 10;
        filter = fspecial('gaussian', [filterSize filterSize], var);
      case 'Box'
        filter = ones(filterSize, filterSize);
        filter = filter./sum(filter(:));
      case 'Motion'
        filter = fliplr( diag(ones(filterSize-1,1), -1) + diag(ones(filterSize,1), 0) + diag(ones(filterSize-1,1), 1) );
        filter = filter./sum(filter(:));
      otherwise
        error('Error: unsupported type of filter')
    end
    
    %% Transform filter mask into corresponding matrix A    
    A = MOPfilter2( height, width, filter, 'mirror', ones(2,2)*floor(size(filter,1)/2) );

    b = A*reshape(u_orig, pixels, channels);

end