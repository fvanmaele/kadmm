%%
res_var = 1.3935e-04;
res0 = loadData(103);
res1 = imresize(res0, 0.1);
res2 = imresize(res0, 0.25);
res3 = imresize(res0, 0.5);
res4 = imresize(res0, 0.75);

%%
t0 = benchmark(res0, 100);
t1 = benchmark(res1, 100);
scale01 = length(res0(:)) / length(res1(:));
t2 = benchmark(res2, 100);
scale02 = length(res0(:)) / length(res2(:));
t3 = benchmark(res3, 100);
scale03 = length(res0(:)) / length(res3(:));
t4 = benchmark(res4, 100);
scale04 = length(res0(:)) / length(res4(:));

%%
function elapsedTime = benchmark(resI, iterations)
    % setup
    [resIP, resIA] = prepareData(resI);
    resIV = unifrnd(0, 1, length(resIP(:)), 3);
    resID = unifrnd(0, 1, length(resIP(:)), 3);
    % timed function
    tic
    kadmmTV(resIA, resIP(:), resIV, resID, 16, 1, 0.01, 1e-3, iterations);
    elapsedTime = toc;
end

%%
function [res, A] = prepareData(res0)
    height = size(res0, 1);
    width = size(res0, 2);
    A = forward_Sv2d(height, width, 21, 4, 'periodic', 'corners');
    res = reshape(A * res0(:), height, width);
    res = imnoise(res, 'gaussian', 0, 1.3935e-04);
end