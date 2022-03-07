%%
function u_orig = loadData(instance)
    dataset = 'BSDS500';
    path = ['data/', dataset,'/'];
    files = dir([path, '*.jpg*']);
    
    if size(files, 1) < 1
        error(['Error: no image data found for dataset ' dataset])
    end
    if instance < 1 || instance > size(files,1)
        error(['Error: instance must be in the range of [1,' num2str(size(files,1)) '] for dataset ' dataset])
    end
    
    u_orig = im2double(imread([path, files(instance).name]));
    u_orig = rgb2gray(u_orig);
end