%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description:This example code provides a calibration algorithm for reconstructing OTS-IFC cell images from raw data obtained via asynchronous sampling.
% Authors: Lei Cheng, Zhou Jiehua(zhoujiehua@whu.edu.cn)
% Affiliation: LeiLab, Institute of Industrial Science, Wuhan University
% Date: 25.02.20
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%----Data read-------------------------------------------------------------
fid=fopen('Data/data_r327_c1099.bin','rb+'); %rb+ 
RawData = fread(fid,'int16');
fclose ('all');
%--------------------------------------------------------------------------

%----Interpolation---------------------------------------------------------
RawData=InsertPoint(RawData,0);%"0" means no interpolation
%--------------------------------------------------------------------------

%----Make sure the waveform is a positive pulse----------------------------
PulPol = 1;             %"1" means positive pulses
if size(RawData, 1) > 1   
    RawData = RawData * PulPol;
else
    RawData = RawData' * PulPol;
end
%--------------------------------------------------------------------------

%----Step 1----------------------------------------------------------------
%fast Fourier transform (FFT) is performed on the original array to obtain the approximate integer pulse repetition period and the pulse numbers of this segment of data
RawData_f = fft(RawData);
RawData_f(1)=0; 
RepRate = find(RawData_f == max(RawData_f(1 : round(end / 2))), 1, 'first');
Width = round(0.5* length(RawData) / RepRate) * 2;
%--------------------------------------------------------------------------

%----Step 2----------------------------------------------------------------
%the fundamental frequency signal of the pulse can be extracted by filtering out the high-frequency components from the original signal, followed by performing the inverse fast Fourier transform (IFFT)
PassBand = round(1.5 * RepRate);
filter = [0; ones(PassBand, 1); zeros(length(RawData) - 1 - 2 * PassBand, 1); ones(PassBand, 1)];
FilterData = real(ifft(RawData_f .* filter));      
SubData = FilterData(1 : Width * 2.5);
[Peaks, Locs] = findpeaks(SubData);              
if Locs(1) - Width / 2 > 0
    FirstPulsePos = Locs(1) - Width / 2;
else
    FirstPulsePos = Locs(2) - Width / 2;
end
%--------------------------------------------------------------------------

%----Step 3----------------------------------------------------------------
%based on the approximate integer pulse repetition period, select out the first one pulse data array and the last two pulses data array from the original array, and calculate the correlation of the two data arrays
FirstPulse = RawData(FirstPulsePos : FirstPulsePos + Width);
CrossCor = [];
for i = length(RawData)- Width * 2.5 : length(RawData)- Width
    CrossCor = [CrossCor, sum(FirstPulse .* RawData(i : i + Width))];
end
%--------------------------------------------------------------------------

%----Step 4----------------------------------------------------------------
%intercept the original array as a new array at the location of the maximum correlation coefficient and the start position to ensure that the new array contains an integer multiple of the number of pulses
LastPulsePos = find(CrossCor == max(CrossCor), 1, 'last') + length(RawData)- Width * 2.5 - 1;
FilterData = FilterData(FirstPulsePos : LastPulsePos - 1);
ColNum = length(findpeaks(FilterData)) + 1;
%--------------------------------------------------------------------------

%----Step 5----------------------------------------------------------------
%select every pulse from the new array according to the start position index to stack the image matrix
Duration = (LastPulsePos - FirstPulsePos) / (ColNum - 1);
StartPoint = round((0 : ColNum - 1) * Duration) + FirstPulsePos - 1;%
ImgIndex = reshape(repmat((0 : Width)', 1, ColNum) + repmat(StartPoint, Width + 1, 1), 1, []);%
Img = RawData(ImgIndex);   %
Img = reshape(Img, Width + 1, []);
%--------------------------------------------------------------------------
figure(1);
imagesc(Img);

%----image cropping--------------------------------------------------------
Pulse_threhold_img_ratio=0.19;
Img_first_column=Img(:,1);
Column_Start_position=find( Img_first_column>(min(Img_first_column)+(max(Img_first_column)-min(Img_first_column))*Pulse_threhold_img_ratio) );
if ~isempty(Column_Start_position)   %
Img=Img(Column_Start_position(1):min(Column_Start_position(end),Width), : );
end
%--------------------------------------------------------------------------

%----Remove background-----------------------------------------------------
Background1 = mean(Img(:, 5 : 25), 2);
Background2 = mean(Img(:, end - 25 : end-5), 2);
if mean(Background1) >= mean(Background2)
    Background = Background1;
else
    Background = Background2;
end
Img = Img - repmat(Background, 1, ColNum);
%--------------------------------------------------------------------------

%----normalization to 8-bit------------------------------------------------
waveMax = max(max(Img)); waveMin = min(min(Img));
Img = (Img - waveMin)*(2^8-1)/(waveMax - waveMin);
Img = uint8(Img);
%--------------------------------------------------------------------------
figure(2);
imshow(Img);

%----Save image------------------------------------------------------------
Img = imresize(Img,[200,200]);
imwrite(Img,'Img/image.png');
%--------------------------------------------------------------------------
figure(3);
imshow(Img);
