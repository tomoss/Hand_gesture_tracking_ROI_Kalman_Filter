function [min_col_left, min_col_right, min_row_top, right_widthbbox, left_widthbbox, top_highbbox] = getRegions(min_col, min_row, max_col, max_row, w, h)

min_col_left  = max_col;
min_col_right = min_col - 0.2*w; 
min_row_top   = min_row - 0.2*h;

right_widthbbox = 0.2*w;
left_widthbbox  = 0.2*w;
top_highbbox    = 0.4*h;

end