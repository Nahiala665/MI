%% Intensity
Im_max = max(max(Im)); % to have the max value of intensity -> do you use that later ? 

%% contour
figure
subplot(1,3,1), imshow(Im_bin,[]), title('Binarized')
subplot(1,3,2), imcontour(Im_int,1,'m')
subplot(1,3,3), imshow(Im_bin,[]), hold on, imcontour(Im_int,1,'m'), title('Contours')
% do we use the contours ? 