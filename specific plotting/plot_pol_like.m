[x,y] = ndgrid(bmvec,avec);
balike = squeeze(sum(like,1));
% for ii = 1:6
% figure(ii+2); clf; surf(x,y,balike(:,:,ii*2));
% end
%%
bamar = 1;
for jj = 1:1:jmax
bamar = bamar.*balike(:,:,jj)/sum(sum(balike(:,:,jj)));
bamar = 100*bamar/sum(sum(bamar));
end
figure(9);clf; surf(x,y,bamar);
%%
bamar = 0;
for jj = 1:1:jmax
bamar = bamar+balike(:,:,jj)/sum(sum(balike(:,:,jj)));
end
figure(10);clf; surf(x,y,bamar);