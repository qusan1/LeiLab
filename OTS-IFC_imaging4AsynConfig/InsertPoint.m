%Change at 220904 Make sure the out put is a column vector
function [Waveform_back] = InsertPoint( Waveform,insertpoint )

x = 1:1:length(Waveform); 
% insertpoint=0;
v =double( Waveform);
 xq=1:(1/(insertpoint+1)):length(Waveform);
Waveform_back= interp1(x,v,xq);  %Default is 'linear'
% Waveform_back= interp1(x,v,xq,'spline');

if size(Waveform_back, 2) > 1  %Make sure it's a column vector 
Waveform_back = Waveform_back';
end
end