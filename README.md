# Light-field-refocusing
Matlab code for refocusing images for Light field images

%R. Sharma, S. Perry, and E. Cheng, "Noise-Resilient Depth Estimation for Light Field Images Using Focal Stack and FFT Analysis," Sensors, vol. 22, no. 5, p. 1993, 2022.

% Light field refocusing algorithm that can generated focal stack with
% slope difference of 0.01

%INPUT : LF = Light field image [9,9,H,W,color channel]
%        Slope = depth value to refocus. (synthetic image = -4 to +4)
%                                        (Lytro image = -2 to +2)
%
%OUTPUT : ImgOut_med = refocused image  [H,W,color channel]
