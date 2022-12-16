function [Sensitivity, Specificity, Similarity] = SegmentationPerformance(A,B)
    % A is the ground truth, B is the segmented shape.
    % A and B need to be binary images
    sumindex = A + B;
    % True Positive
    TP = length(find(sumindex == 2));
    % True Negative
    TN = length(find(sumindex == 0));
    substractindex = A - B;
    % False Positive
    FP = length(find(substractindex == -1));
    %False Negative
    FN = length(find(substractindex == 1));

    % Evaluation of the segmentation
    Sensitivity = TP/(TP+FN);
    Specificity = TN/(TN+FP);
    Similarity = dice(A,B);

end
   