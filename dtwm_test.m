function [] = dtwm_test()

dbstop if error; 

%test signals
x = [1 2 3 4 4 5]; 
y = [1 2 2 2 3 4 5 6];


try
    mex dtwm_costmatrix.c;
catch
    fprintf('mex file for costmatrix d/n compile successfully.\nmust call dtwm with ''use_mex'',''false''\n');
    return;
end

tic
map = dtwm(x,y);
toc

subplot(2,2,[1 3]);
plot_distance(x,y,map,'grid',true);

subplot(2,2,2);
plot_warp(x,y,map);

subplot(2,2,4);
plot_warp(x,y,map,'mode','warped');



end
