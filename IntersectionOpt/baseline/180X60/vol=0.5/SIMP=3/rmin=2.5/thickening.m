%% Thickening of the part
clear all
close all

load('baseline.mat','xPhys')
nely = size(xPhys,1);
nelx = size(xPhys,2);
thick = 20;
xPhys_thick = repmat(xPhys,1,1,thick);

figure(1)
colormap(gray); imagesc(1-xPhys_thick(:,:,1)); caxis([0 1]); axis equal;axis tight manual; axis off; drawnow;daspect([1 1 1])
figure(2)
colormap(gray); imagesc(1-xPhys_thick(:,:,floor(thick/2))); caxis([0 1]); axis equal;axis tight manual; axis off; drawnow;daspect([1 1 1])
figure(3)
colormap(gray); imagesc(1-xPhys_thick(:,:,thick)); caxis([0 1]); axis equal;axis tight manual; axis off; drawnow;daspect([1 1 1])

save('comp.mat','xPhys_thick')