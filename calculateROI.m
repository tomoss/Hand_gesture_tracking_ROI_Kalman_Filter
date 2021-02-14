function [min_col, min_row, max_col, max_row, w, h]= calculateROI(SM)

R=[]; C=[];
for i = 1:size(SM,1)
    for j = 1:size(SM,2)
        if SM(i,j) == 1
           R(end+1) = i;
           C(end+1) = j;
        end
    end
end

max_col = max(C);
min_col = min(C);
max_row = max(R);
min_row = min(R);

w = max_col - min_col;
h = max_row - min_row;

end