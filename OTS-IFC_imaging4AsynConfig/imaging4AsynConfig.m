%%%%%%%%%%%%%%%%Data read%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------
% OldRawData=textread('STEAM0185.txt');
%-----------------------------------------------------------------
%    OldRawData=load('TS_cell_5_10G.txt');
%-----------------------------------------------------------------  
%   OldRawData=importdata('TS_cell_5_10G.mat');
%-----------------------------------------------------------------
fid=fopen('Data/data_r436_c1364.bin','rb+'); %rb+ ��д��һ���������ļ���ֻ�����д����
RawData = fread(fid,'int16');
fclose ('all');
%-----------------------------------------------------------------

%-------------�����ֵ����---------------------------------
RawData=InsertPoint(RawData,1);
%--------------------------------------------------------
PulPol = 1;             %"1" means positive pulses
if size(RawData, 1) > 1   
    RawData = RawData * PulPol;
else
    RawData = RawData' * PulPol;
end

RawData_f = fft(RawData);
RawData_f(1)=0; %��ֱ������ֱ�ӱ�0�� ʡȥ���������ԳƲ���;
RepRate = find(RawData_f == max(RawData_f(1 : round(end / 2))), 1, 'first');%�ҳ�RawData_f���ֵ��λ�ã����������Ƶ�ʣ�
Width = round(0.4* length(RawData) / RepRate) * 2;

PassBand = round(1.5 * RepRate);
filter = [0; ones(PassBand, 1); zeros(length(RawData) - 1 - 2 * PassBand, 1); ones(PassBand, 1)];%ǰһ��PassBand����һ��PassBand
FilterData = real(ifft(RawData_f .* filter));      %ifft Ϊ��ɢ���ٸ���Ҷ��任

SubData = FilterData(1 : Width * 2.5);
[Peaks, Locs] = findpeaks(SubData);                %Peaks��Ӧ��ֵ��Locs��Ӧ��λ��
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
StartPoint = round((0 : ColNum - 1) * Duration) + FirstPulsePos - 1;%��������ȡÿһ�����忪ʼ�ĵ�
ImgIndex = reshape(repmat((0 : Width)', 1, ColNum) + repmat(StartPoint, Width + 1, 1), 1, []);%��repmat��������λ�þ�������reshape������Ϊ����

Img = RawData(ImgIndex);   %��λ��ȡ
Img = reshape(Img, Width + 1, []);
figure(1);
imagesc(Img);

%-----------------����ͼ��y����----------------------------------------------------------------------------------
Pulse_threhold_img_ratio=0.19;
Img_first_column=Img(:,1);
Column_Start_position=find( Img_first_column>(min(Img_first_column)+(max(Img_first_column)-min(Img_first_column))*Pulse_threhold_img_ratio) );

%----
if ~isempty(Column_Start_position)   %��ֹ��BUG ��Ϊ�ղż���
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

%-------------��һ���� 0-255 -----------------------------------------
waveMax = max(max(Img)); waveMin = min(min(Img));
Img = (Img - waveMin)*(2^8-1)/(waveMax - waveMin);
Img = uint8(Img);
%--------------------------------------------------------------------

Img = imresize(Img,[400,400]);
figure(2);
imshow(Img);
imwrite(Img,'Img/image.png')
  
