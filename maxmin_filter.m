function [filtered] = maxmin_filter(im, mode, win_size)
% maxmin_filter - 最大最小值滤波
%
% input:
%   - im: h*w, 2D待滤波数据
%   - mode: str, 滤波方式, max/min
%   - win_size: int, 滤波窗口半径
% output:
%   - filtered: h*w, 滤波后数据
%
% docs:
%   - 对数据进行镜像填充, 避免滤波后数据大小与原始数据不一致
%

mode = lower(mode);
modes = struct('max', @max, 'min', @min);
if isfield(modes, mode)
    func = getfield(modes, mode);
else
    error('%s is not correct parameter, must be [max, min]!', mode);
end

[h, w, c] = size(im);
if c > 1
    error('%s must be 2D!', im);
end
filtered = zeros(h, w);

% 填充
im_padded = padarray(im, [win_size, win_size], Inf);

% filter
for row = win_size+1:win_size+h
    for col = win_size+1:win_size+w
        win = im_padded(row-win_size:row+win_size, col-win_size:col+win_size);
        win = win(:);
        filtered(row-win_size, col-win_size) = func(win);
    end
end

end