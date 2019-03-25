function NaturalScenesToContrast(localPath)
addpath('nsutils');

if nargin < 1
    localPath = pwd;
end

% load images
cd(localPath)
files = dir('imageData\*.mat');
scenes = [];
for ff = 1:size(files,1)
    fpath = fullfile(localPath,'imageData', files(ff).name);
    data = load(fpath);
    scenes = cat(3, scenes, data.projection);
    fprintf('Loaded file: %s\n', fpath);
end

%     scenes = data.projection;
xRes = 360/size(scenes,2);

x = (0:xRes:360-xRes)';
y = (0:xRes:xRes*size(scenes,1)-xRes)';

% individual photo receptors
filtStd = 5/(2*sqrt(2*log(2)));

xFilt = ifftshift(normpdf(x,180,filtStd));
yFilt = ifftshift(normpdf(y,(y(end)+xRes)/2,filtStd));

xyFiltMat = yFilt*xFilt';

xyFiltTensor = repmat(xyFiltMat,[1 1 size(scenes,3)]);

% for mean estimation
filtStdContrast = 20;

xFiltContrast = ifftshift(normpdf(x,180,filtStdContrast));
yFiltContrast = ifftshift(normpdf(y,(y(end)+xRes)/2,filtStdContrast));

xyFiltMatContrast = yFiltContrast*xFiltContrast';

xyFiltTensorContrast = repmat(xyFiltMatContrast,[1 1 size(scenes,3)]);

% perform the convolution
fftScenes = fft2(scenes);
fftFilt = fft2(xyFiltTensor);
fftFiltContrast = fft2(xyFiltTensorContrast);

filteredScenes = real(ifft2(fftScenes.*fftFilt));
filteredScenesContrast = real(ifft2(fftScenes.*fftFiltContrast));

finalContrast = (filteredScenes-filteredScenesContrast)./filteredScenesContrast;

MakeFigure;
subplot(1,2,1);
imagesc(log(scenes(:,:,99)))
subplot(1,2,2);
imagesc(finalContrast(:,:,99));
colormap(gray);

if ~isfolder('savedData'), mkdir('savedData'); end
spath = fullfile(localPath, '\savedData\combinedFiltered2D.mat');
save(spath, 'scenes', 'finalContrast', '-v7.3');
fprintf('Saved contrast images to %s\n', spath);

end