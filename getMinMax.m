function [min_right, max_left, min_top]= getMinMax(right, left, top)

R_right=[]; C_right=[];
for i = 1:size(right,1)
    for j = 1:size(right,2)
        if right(i,j) == 1
           R_right(end+1) = i;
           C_right(end+1) = j;
        end
    end
end
min_right = min(C_right);

R_left=[]; C_left=[];
for i = 1:size(left,1)
    for j = 1:size(left,2)
        if left(i,j) == 1
           R_left(end+1) = i;
           C_left(end+1) = j;
        end
    end
end
max_left = max(C_left);

R_top=[]; C_top=[];
for i = 1:size(top,1)
    for j = 1:size(top,2)
        if top(i,j) == 1
           R_top(end+1) = i;
           C_top(end+1) = j;
        end
    end
end
min_top= min(R_top);

end