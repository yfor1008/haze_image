function [dehazed_im] = dehaze(varargin)
% dehaze - 去雾
%
% usage: 
%   [dehazed_im] = dehaze(im, win_size, ratio, w, thres, t0)
%
% input:
%   - im: h*w*3, rgb图像
%   - win_size: int, 最小值滤波窗口半径
%   - ratio: float, 取有雾像素比例, 文章中为0.001
%   - w: float, 透射率限制系数, 文章中为0.95, 值越小去雾效果越不明显
%   - thres: int, 对大气光值限制阈值, 默认为220
%   - t0: float, 透射率限制阈值, 默认为0.1
% output:
%   - dehazed_im: h*w*3, 去雾后图像
%
% docs:
%   - 实现何凯明的《Single Image Haze Removal Using Dark Channel Prior》
%   - step1: 计算暗原色, 利用公式5
%   - step2: 估计大气光值A
%   - step3: 估计透射率, 利用公式12
%   - step4: 对透射率进行调整, 得到更精细透射率图像, 原文中使用soft matting, 这里使用导向滤波替代
%   - step5: 恢复得到辐射图像, 利用公式22
%

% 解析参数
[im, win_size, ratio, w, thres, t0] = parse_varargin(varargin{:});

im = double(im);

% step1
dc = get_dark_channel(im, win_size);

% step2
A = estimate_A(im, dc, ratio, thres);

% step3
transmission = estimate_t(im, A, win_size, w);

% step4
% transmission = matte(im, transmission); % 非常慢
filterRadius = ceil(min(size(im,1), size(im,2)) / 50);
% transmission = guide_filter(transmission, im, filterRadius, 10^-4); % 滤波窗口半径至少是最小值滤波的4倍
transmission = imguidedfilter(transmission, im, [filterRadius*2+1, filterRadius*2+1], 10^-4, 4); % 滤波窗口半径至少是最小值滤波的4倍

% step5
dehazed_im = get_radiance(im, A, transmission, t0);
dehazed_im = uint8(dehazed_im);

end

function [dark_channel] = get_dark_channel(im, win_size)
% get_dark_channel - 获取暗原色
%
% input:
%   - im: h*w*3, rgb图像
%   - win_size: int, 滤波窗口半径
% output:
%   - dark_channel: h*w, 暗原色
%
% docs:
%   - step1: 暗原色先验, rgb通道中最小的数
%   - step2: 暗原色, 最小滤波
%

[h, w, ~] = size(im);
dark_channel = zeros(h, w);

% step1
for row = 1:h
    for col = 1:w
        dark_channel(row, col) = min(im(row, col, :));
    end
end

% step2
dark_channel = maxmin_filter(dark_channel, 'min', win_size);

end

function [A] = estimate_A(im, dc, ratio, thres)
% estimate_A - 估计大气光值
%
% input：
%   - im：h*w*c, 原始图像, c=3时为RGB图像, c=1时为gray图像
%   - dc: h*w, 暗原色
%   - ratio: float, 取有雾像素比例
%   - thres: int, A的阈值
% output:
%   - A: c*1, 大气光值
%
% docs:
%   - step1: 从暗原色中按亮度取前0.1%的像素, 这些像素大部分是有雾的
%   - step2: 这些像素中, 在原始图像寻找这些像素对应的具有最高亮度的值作为A值, 这里取这些像素的均值
%   - step3: 这里增加对A的限制
%

[h, w, c] = size(im);

% step1
im_size = h * w;
bright_dc = reshape(dc, im_size, 1); % 变成列, 方便处理
[~, idx] = sort(bright_dc, 'descend'); % 对dc进行排序
topNum = floor(im_size * ratio);

% step2
bright_im = reshape(im, im_size, c);
A = zeros(c, 1);
for ind = 1:topNum % 取前ratio数据
    A = A + bright_im(idx(ind));
end
A = A / topNum;

% step3
A = min(A, thres);

end

function [transmission] = estimate_t(im, A, win_size, w)
% estimate_t - 估计透射率
%
% input:
%   - im: h*w*3, rgb图像
%   - A: c*1, 大气光值, c=3时为rgb图像, c=1时为gray图像
%   - win_size: int, 滤波窗口半径
%   - w: float, [0,1]限制参数
% output:
%   - transmission: h*w, 透射率
%
% docs:
%   - 公式12

c = length(A);
im3 = zeros(size(im));
if c == 3
    for idx = 1:c
        im3(:,:,idx) = im(:,:,idx) / A(idx);
    end
else
    im3 = im3 / A;
end

transmission = 1 - w * get_dark_channel(im3, win_size);

end

function [radiance] = get_radiance(im, A, t, t0)
% get_radiance - 得到辐射图
%
% input:
%   - im: h*w*3, rgb图像
%   - A: c*1, 大气光值, c=3时为rgb图像, c=1时为gray图像
%   - t: h*w, 透射率
%   - t0: float, 透射率阈值
% output:
%   - radiance: h*w*3, 恢复得到的rgb图像
%
% docs:
%   - 利用公式22
%

c = length(A);
radiance = zeros(size(im));

t = max(t, t0);

if c == 3
    for idx = 1:3
        radiance(:,:,idx) = (im(:,:,idx) - A(idx)) ./ t + A(idx);
    end
else
    for idx = 1:3
        radiance(:,:,idx) = (im(:,:,idx) - A) ./ t + A;
    end
end

end

function [im, win_size, ratio, w, thres, t0] = parse_varargin(varargin)
% parse_varargin - 解析参数
%

num = length(varargin);
switch num
case 1
    im = varargin{1};
    win_size = ceil(min(size(im,1), size(im,2)) / 400 * 15);
    ratio = 0.001;
    w = 0.95;
    thres = 220;
    t0 = 0.1;
case 2
    im = varargin{1};
    win_size = varargin{2};
    ratio = 0.001;
    w = 0.95;
    thres = 220;
    t0 = 0.1;
case 3
    im = varargin{1};
    win_size = varargin{2};
    ratio = varargin{3};
    w = 0.95;
    thres = 220;
    t0 = 0.1;
case 4
    im = varargin{1};
    win_size = varargin{2};
    ratio = varargin{3};
    w = varargin{4};
    thres = 220;
    t0 = 0.1;
case 5
    im = varargin{1};
    win_size = varargin{2};
    ratio = varargin{3};
    w = varargin{4};
    thres = varargin{5};
    t0 = 0.1;
case 6
    im = varargin{1};
    win_size = varargin{2};
    ratio = varargin{3};
    w = varargin{4};
    thres = varargin{5};
    t0 = varargin{6};
otherwise
    error('parameters nor correct!');
end
test = 0;


end    
    