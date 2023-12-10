%% Rings Segmentation Image:
% Segments the rings around the chosen area in the image.
%
% Syntax:
% seg_img = PanGUI_RingSegImg(img,img_L)
%
% Input:
% img - Original image, uint8[0,255].
% img_L - The segmented areas, values of [0,1,2,3,4].
%
% Output:
% seg_img - The segmented image, double[0,255].

function seg_img = PanGUI_RingSegImg(img,img_L)
seg_img = img;
seg_img(:,:,1) = ...
    (img_L == 1).*0.5.*(255 + double(img(:,:,1))) + ...
    (img_L == 2).*0.5.*(0 + double(img(:,:,1))) + ...
    (img_L == 3).*0.5.*(0 + double(img(:,:,1))) + ...
    (img_L == 4).*0.5.*(255 + double(img(:,:,1))) + ...
    (img_L ~= 1).*(img_L ~= 2).*(img_L ~= 3).*(img_L ~= 4).*double(img(:,:,1));
seg_img(:,:,2) = ...
    (img_L == 1).*0.5.*(0 + double(img(:,:,2))) + ...
    (img_L == 2).*0.5.*(255 + double(img(:,:,2))) + ...
    (img_L == 3).*0.5.*(0 + double(img(:,:,2))) + ...
    (img_L == 4).*0.5.*(255 + double(img(:,:,2))) + ...
    (img_L ~= 1).*(img_L ~= 2).*(img_L ~= 3).*(img_L ~= 4).*double(img(:,:,2));
seg_img(:,:,3) = ...
    (img_L == 1).*0.5.*(0 + double(img(:,:,3))) + ...
    (img_L == 2).*0.5.*(0 + double(img(:,:,3))) + ...
    (img_L == 3).*0.5.*(255 + double(img(:,:,3))) + ...
    (img_L == 4).*0.5.*(0 + double(img(:,:,3))) + ...
    (img_L ~= 1).*(img_L ~= 2).*(img_L ~= 3).*(img_L ~= 4).*double(img(:,:,3));
end