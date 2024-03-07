function [X,Y,T,AUC] = ROC_EM(t_detection,t_detection_wt)

pixel_num = numel(t_detection);
% Plot the ROC curve
t_detection = reshape(t_detection,pixel_num,1);
t_detection_wt = reshape(t_detection_wt,pixel_num,1);
output_values = [t_detection; t_detection_wt];
labels_GT = [0*t_detection; ones(size(t_detection_wt))];
[X,Y,T,AUC] = perfcurve(labels_GT,output_values,1);


% interpulate data
X = interp(X,4);
Y = interp(Y,4);

end