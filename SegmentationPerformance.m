function [Sensitivity, Specificity, Similarity] = SegmentationPerformance(A,B)
    % A is the ground truth, B is the segmented result.
    % A and B need to be binary images
    sumindex = A + B;
    TP = length(find(sumindex == 2));
    TN = length(find(sumindex == 0));
    substractindex = A - B;
    FP = length(find(substractindex == -1));
    FN = length(find(substractindex == 1));
    Sensitivity = TP/(TP+FN);
    Specificity = TN/(TN+FP);
    Similarity = dice(A,B);
end
   