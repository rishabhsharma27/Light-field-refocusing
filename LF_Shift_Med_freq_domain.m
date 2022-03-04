%R. Sharma, S. Perry, and E. Cheng, "Noise-Resilient Depth Estimation for Light Field Images Using Focal Stack and FFT Analysis," Sensors, vol. 22, no. 5, p. 1993, 2022.

% Light field refocusing algorithm that can generated focal stack with
% slope difference of 0.01

%INPUT : LF = Light field image [9,9,H,W,color channel]
%        Slope = depth value to refocus. (synthetic image = -4 to +4)
%                                        (Lytro image = -2 to +2)
%
%OUTPUT : ImgOut_med = refocused image  [H,W,color channel]
                                                    
            

function [ImgOut_med] = LF_Shift_Med_freq_domain( LF, Slope )


    if(size(LF,3)==size(LF,4))

        w = hann(size(LF,3));
        w=w.*w';
        for i=1:size(LF,1)
            for j=1:size(LF,2)

                a=im2single(squeeze(LF(i,j,:,:,:)));

                LF_sep{i,j}=(ifftshift(fftshift(fft2(a)).*w)); 

            end
        end
        w_count = [1,size(LF,3)]; 
    elseif(size(LF,3)>size(LF,4))
        
        w = hann(size(LF,3));
        w=w.*w';
        for i=1:size(LF,1)
            for j=1:size(LF,2)
                
                a=im2single(squeeze(LF(i,j,:,:,:)));
                a=padarray(a,[0 (size(LF,3)-size(LF,4))],'pre' );
                LF1(i,j,:,:,:) = a; 
                LF_sep{i,j}=(ifftshift(fftshift(fft2(a)).*w)); 

            end
        end
        w_count = [2,size(LF,4)];
        LF=LF1;
        
    elseif(size(LF,4)>size(LF,3))
        w = hann(size(LF,4));
        w=w.*w';
        for i=1:size(LF,1)
            for j=1:size(LF,2)
                
                a=im2single(squeeze(LF(i,j,:,:,:)));
                a=padarray(a,[(size(LF,4)-size(LF,3)) 0],'pre' );
                LF1(i,j,:,:,:) = a; 
                LF_sep{i,j}=(ifftshift(fftshift(fft2(a)).*w)); 

            end
        end
        w_count = [3,size(LF,3)]; 
        LF=LF1;
    end
        
        

FiltOptions = struct('UpsampRate', 1);  % Change 

% if (FiltOptions.UpsampRate~=1)
%     LFSize = size(LF);
%     a1=(squeeze(LF(1,1,:,:,:)));
%     a1 = imresize(a1,FiltOptions.UpsampRate);
%     LFSize(1,3) = size(a1,1);
%     LFSize(1,4) = size(a1,2);
% else
    LFSize = size(LF);
% end


NColChans = size(LF,5);


TVSlope = Slope ;
SUSlope = Slope ;

%%
v = linspace(1,LFSize(3), round(LFSize(3)*FiltOptions.UpsampRate));
u = linspace(1,LFSize(4), round(LFSize(4)*FiltOptions.UpsampRate));


%%


NewLFSize = LFSize;
NewLFSize(3:4) = [length(v), length(u)];



%%

VOffsetVec=zeros(1,LFSize(1));

UOffsetVec=zeros(1,LFSize(2));

for i=1:LFSize(1)
    VOffsetVec(i)=(i-1-floor(LFSize(1)/2)).*TVSlope;
    UOffsetVec(i)=(i-1-floor(LFSize(2)/2)).*SUSlope;
end


%%
    LFOut = zeros(NewLFSize, 'like', LF);


    [a,b,~]=size(LF_sep{1,1});
    Nr = ifftshift([-fix(a/2):ceil(a/2)-1]);
    Nc = ifftshift([-fix(b/2):ceil(b/2)-1]);
    [xF,yF] = meshgrid(Nr,Nc);
    
   
          
for( TIdx = 1:LFSize(1) )
	VOffset = VOffsetVec(TIdx)*FiltOptions.UpsampRate;
    
    for( SIdx = 1:LFSize(2) )
		UOffset = UOffsetVec(SIdx)*FiltOptions.UpsampRate;
        

          in=LF_sep{TIdx,SIdx}; %// Define input signal
          

            x0=-UOffset; %// Define shifts
            y0=-VOffset;
            

            H=exp(1i*2*pi*(x0*xF/a+y0*yF/b));
            IF_image = real(ifft2(in.*H));
            
            x0=round(x0);
            y0=round(y0);
            if (x0>0)
               IF_image(:,end-abs(x0)+1:end,:)=0;
            elseif(x0<0)
                IF_image(:,1:abs(x0),:)=0;
            end
            
            if (y0>0)
               IF_image(end-abs(y0)+1:end,:,:)=0;
            elseif(y0<0)
                IF_image(1:abs(y0),:,:)=0;
            end

                LFOut(TIdx,SIdx, :,:, :)=IF_image;

        
    end
end

        LF = LFOut;

		t = reshape(LF(:,:,:,:,1:NColChans), [prod(LFSize(1:2)), NewLFSize(3:4), NColChans]);
		ImgOut_med = squeeze(nanmedian(t));
        if( w_count(1,1)==1)
            ImgOut_med = ImgOut_med;
        elseif( w_count(1,1)==2)
            ImgOut_med = ImgOut_med(:,size(ImgOut_med,2)-w_count(1,2)+1:end,:);
        elseif( w_count(1,1)==3)
            ImgOut_med = ImgOut_med(size(ImgOut_med,2)-w_count(1,2)+1:end,:,:);
        end


end

