function [] = display3dImage(image, colorMap, resizeFactor, titleStrings)
% titleStrings must be a string array to display titles for multiple
% concatenated base images

% plot initial image
slice0 = round(size(image, 3)/2);
picH = imshow(absNorm(myImresizeNearest(image(:,:,slice0))),...
                                    [0, 1], 'border', 'loose');
             
colormap(colorMap)

% set figure position and borders
set(gca, 'Units', 'Pixels')
axPos = get(gca, 'position');
borderWidth = [40, 50, 40, 70];     % left, top, right, bottom]
set(gcf, 'position', [100, 100, axPos(3) + sum(borderWidth([1, 3])), axPos(4) +  sum(borderWidth([2, 4]))])
set(gca, 'position', [borderWidth(1), borderWidth(4), axPos(3), axPos(4)])

% create GUI
axPos = get(gca, 'position');
sliderPos(2) = axPos(2) - 30;
sliderPos(3) = min(axPos(3), 512);
sliderPos(4) = 20;
sliderPos(1) = axPos(1) + round((axPos(3) - sliderPos(3))/2);

sliderH = uicontrol('style','slider','position', sliderPos,...
                    'min', 1, 'max', size(image, 3),...
                    'value', slice0,...
                    'sliderstep', [1, 1]/max((size(image, 3) - 1), 1));

textH = uicontrol('style','text', 'position', sliderPos + [0, -30, 0, 5],...
                  'string', num2str(slice0, 'slice %d'), ...
                  'FontSize', 17);
                
                
addlistener(sliderH, 'Value', 'PostSet', @callbackfn);

% add titles
titleWidth = axPos(3)/numel(titleStrings);

titlePos = zeros(1, 4);
titlePos(2) = axPos(2) + axPos(4) + 12;
titlePos(3) = titleWidth;
titlePos(4) = 28;

for ii = 1:numel(titleStrings)
    
    titlePos(1) = axPos(1) + (ii - 1)*titleWidth;
    
    % to place a text box using pixel units
    uicontrol(gcf, 'style', 'text', 'position', titlePos,...
              'string', titleStrings(ii), 'FontSize', 17)
    
end


function out = myImresizeNearest(in)
% to help presentation of low resolution data
% equivalent to image processing toolbox function
% imresize with 'method', 'nearest' option
dimIn = size(in);
assert(length(size(dimIn)) == 2, 'only accepts 2D images')

out = repmat(in, [1, 1, resizeFactor, resizeFactor]);
out = permute(out, [3, 1, 4, 2]);
out = reshape(out, dimIn*resizeFactor);

end

function callbackfn(source, eventdata)
% exectutes GUI update
slice = max(round(get(eventdata.AffectedObject, 'Value')), 1);
set(sliderH, 'Value', slice)
picH.CData  = absNorm(myImresizeNearest(image(:,:,slice)));
textH.String = num2str(slice, 'slice %d');

end

end