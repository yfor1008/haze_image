function [dehazed_im] = dehaze(varargin)
% dehaze - ȥ��
%
% usage: 
%   [dehazed_im] = dehaze(im, win_size, ratio, w, thres, t0)
%
% input:
%   - im: h*w*3, rgbͼ��
%   - win_size: int, ��Сֵ�˲����ڰ뾶
%   - ratio: float, ȡ�������ر���, ������Ϊ0.001
%   - w: float, ͸��������ϵ��, ������Ϊ0.95, ֵԽСȥ��Ч��Խ������
%   - thres: int, �Դ�����ֵ������ֵ, Ĭ��Ϊ220
%   - t0: float, ͸����������ֵ, Ĭ��Ϊ0.1
% output:
%   - dehazed_im: h*w*3, ȥ���ͼ��
%
% docs:
%   - ʵ�ֺο����ġ�Single Image Haze Removal Using Dark Channel Prior��
%   - step1: ���㰵ԭɫ, ���ù�ʽ5
%   - step2: ���ƴ�����ֵA
%   - step3: ����͸����, ���ù�ʽ12
%   - step4: ��͸���ʽ��е���, �õ�����ϸ͸����ͼ��, ԭ����ʹ��soft matting, ����ʹ�õ����˲����
%   - step5: �ָ��õ�����ͼ��, ���ù�ʽ22
%

% ��������
[im, win_size, ratio, w, thres, t0] = parse_varargin(varargin{:});

im = double(im);

% step1
dc = get_dark_channel(im, win_size);

% step2
A = estimate_A(im, dc, ratio, thres);

% step3
transmission = estimate_t(im, A, win_size, w);

% step4
% transmission = matte(im, transmission); % �ǳ���
filterRadius = ceil(min(size(im,1), size(im,2)) / 50);
% transmission = guide_filter(transmission, im, filterRadius, 10^-4); % �˲����ڰ뾶��������Сֵ�˲���4��
transmission = imguidedfilter(transmission, im, [filterRadius*2+1, filterRadius*2+1], 10^-4, 4); % �˲����ڰ뾶��������Сֵ�˲���4��

% step5
dehazed_im = get_radiance(im, A, transmission, t0);
dehazed_im = uint8(dehazed_im);

end

function [dark_channel] = get_dark_channel(im, win_size)
% get_dark_channel - ��ȡ��ԭɫ
%
% input:
%   - im: h*w*3, rgbͼ��
%   - win_size: int, �˲����ڰ뾶
% output:
%   - dark_channel: h*w, ��ԭɫ
%
% docs:
%   - step1: ��ԭɫ����, rgbͨ������С����
%   - step2: ��ԭɫ, ��С�˲�
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
% estimate_A - ���ƴ�����ֵ
%
% input��
%   - im��h*w*c, ԭʼͼ��, c=3ʱΪRGBͼ��, c=1ʱΪgrayͼ��
%   - dc: h*w, ��ԭɫ
%   - ratio: float, ȡ�������ر���
%   - thres: int, A����ֵ
% output:
%   - A: c*1, ������ֵ
%
% docs:
%   - step1: �Ӱ�ԭɫ�а�����ȡǰ0.1%������, ��Щ���ش󲿷��������
%   - step2: ��Щ������, ��ԭʼͼ��Ѱ����Щ���ض�Ӧ�ľ���������ȵ�ֵ��ΪAֵ, ����ȡ��Щ���صľ�ֵ
%   - step3: �������Ӷ�A������
%

[h, w, c] = size(im);

% step1
im_size = h * w;
bright_dc = reshape(dc, im_size, 1); % �����, ���㴦��
[~, idx] = sort(bright_dc, 'descend'); % ��dc��������
topNum = floor(im_size * ratio);

% step2
bright_im = reshape(im, im_size, c);
A = zeros(c, 1);
for ind = 1:topNum % ȡǰratio����
    A = A + bright_im(idx(ind));
end
A = A / topNum;

% step3
A = min(A, thres);

end

function [transmission] = estimate_t(im, A, win_size, w)
% estimate_t - ����͸����
%
% input:
%   - im: h*w*3, rgbͼ��
%   - A: c*1, ������ֵ, c=3ʱΪrgbͼ��, c=1ʱΪgrayͼ��
%   - win_size: int, �˲����ڰ뾶
%   - w: float, [0,1]���Ʋ���
% output:
%   - transmission: h*w, ͸����
%
% docs:
%   - ��ʽ12

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
% get_radiance - �õ�����ͼ
%
% input:
%   - im: h*w*3, rgbͼ��
%   - A: c*1, ������ֵ, c=3ʱΪrgbͼ��, c=1ʱΪgrayͼ��
%   - t: h*w, ͸����
%   - t0: float, ͸������ֵ
% output:
%   - radiance: h*w*3, �ָ��õ���rgbͼ��
%
% docs:
%   - ���ù�ʽ22
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
% parse_varargin - ��������
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
    