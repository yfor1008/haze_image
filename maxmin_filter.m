function [filtered] = maxmin_filter(im, mode, win_size)
% maxmin_filter - �����Сֵ�˲�
%
% input:
%   - im: h*w, 2D���˲�����
%   - mode: str, �˲���ʽ, max/min
%   - win_size: int, �˲����ڰ뾶
% output:
%   - filtered: h*w, �˲�������
%
% docs:
%   - �����ݽ��о������, �����˲������ݴ�С��ԭʼ���ݲ�һ��
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

% ���
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