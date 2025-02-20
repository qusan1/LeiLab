%%%%%%%%%%%%%%%%Data read%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------
% OldRawData=textread('STEAM0185.txt');
%-----------------------------------------------------------------
%    OldRawData=load('TS_cell_5_10G.txt');
%-----------------------------------------------------------------  
%   OldRawData=importdata('TS_cell_5_10G.mat');
%-----------------------------------------------------------------
fid=fopen('Data/data_r436_c1364.bin','rb+'); %rb+ 读写打开一个二进制文件，只允许读写数据
RawData = fread(fid,'int16');
fclose ('all');
%-----------------------------------------------------------------

%-------------整体插值部分---------------------------------
RawData=InsertPoint(RawData,1);
%--------------------------------------------------------
PulPol = 1;             %"1" means positive pulses
if size(RawData, 1) > 1   
    RawData = RawData * PulPol;
else
    RawData = RawData' * PulPol;
end

RawData_f = fft(RawData);
RawData_f(1)=0; %将直流分量直接变0， 省去上面正负对称操作;
RepRate = find(RawData_f == max(RawData_f(1 : round(end / 2))), 1, 'first');%找出RawData_f最大值的位置（横坐标代表频率）
Width = round(0.4* length(RawData) / RepRate) * 2;

PassBand = round(1.5 * RepRate);
filter = [0; ones(PassBand, 1); zeros(length(RawData) - 1 - 2 * PassBand, 1); ones(PassBand, 1)];%前一个PassBand，后一个PassBand
FilterData = real(ifft(RawData_f .* filter));      %ifft 为离散快速傅里叶逆变换

SubData = FilterData(1 : Width * 2.5);
[Peaks, Locs] = findpeaks(SubData);                %Peaks对应峰值，Locs对应峰位数
if Locs(1) - Width / 2 > 0
    FirstPulsePos = Locs(1) - Width / 2;
else
    FirstPulsePos = Locs(2) - Width / 2;
end
FirstPulse = RawData(FirstPulsePos : FirstPulsePos + Width);

CrossCor = [];
for i = length(RawData)- Width * 2.5 : length(RawData)- Width
    CrossCor = [CrossCor, sum(FirstPulse .* RawData(i : i + Width))];
end
LastPulsePos = find(CrossCor == max(CrossCor), 1, 'last') + length(RawData)- Width * 2.5 - 1;

FilterData = FilterData(FirstPulsePos : LastPulsePos - 1);
ColNum = length(findpeaks(FilterData)) + 1;
Duration = (LastPulsePos - FirstPulsePos) / (ColNum - 1);
StartPoint = round((0 : ColNum - 1) * Duration) + FirstPulsePos - 1;%四舍五入取每一个脉冲开始的点
ImgIndex = reshape(repmat((0 : Width)', 1, ColNum) + repmat(StartPoint, Width + 1, 1), 1, []);%用repmat函数构造位置矩阵，再用reshape函数变为向量

Img = RawData(ImgIndex);   %按位置取
Img = reshape(Img, Width + 1, []);
figure(1);
imagesc(Img);

%-----------------剪切图像y部分----------------------------------------------------------------------------------
Pulse_threhold_img_ratio=0.19;
Img_first_column=Img(:,1);
Column_Start_position=find( Img_first_column>(min(Img_first_column)+(max(Img_first_column)-min(Img_first_column))*Pulse_threhold_img_ratio) );

%----
if ~isempty(Column_Start_position)   %防止出BUG 不为空才剪切
Img=Img(Column_Start_position(1):min(Column_Start_position(end),Width), : );
end

Background1 = mean(Img(:, 5 : 25), 2);
Background2 = mean(Img(:, end - 25 : end-5), 2);
if mean(Background1) >= mean(Background2)
    Background = Background1;
else
    Background = Background2;
end
Img = Img - repmat(Background, 1, ColNum);

%-------------归一化在 0-255 -----------------------------------------
waveMax = max(max(Img)); waveMin = min(min(Img));
Img = (Img - waveMin)*(2^8-1)/(waveMax - waveMin);
Img = uint8(Img);
%--------------------------------------------------------------------

Img = imresize(Img,[400,400]);
figure(2);
imshow(Img);
imwrite(Img,'Img/image.png')
  
