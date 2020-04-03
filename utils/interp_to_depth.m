map2dnowsh = squeeze(map3dnow(:,:,didxsh));
map2dnowdp = squeeze(map3dnow(:,:,didxdp));
mapcatnow = cat(1, map2dnowsh(:)', map2dnowdp(:)');
map2dnow = interp1(depth(didxsh:didxdp), mapcatnow, newdepth);
map2dnow = reshape(map2dnow, size(map2dnowsh));
