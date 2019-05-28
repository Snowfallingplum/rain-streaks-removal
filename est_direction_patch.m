function direction=est_direction_patch(RainMap)
%% Rainmap is a 3-d data,in this algorithm get initial RainMap got by Dt(Rainy)

Size=size(RainMap);
window_size=floor(0.25*Size(1));
strides=floor(0.1*Size(1));
number=floor((Size(1)-window_size)/strides)+1;
direct=zeros(number*number,Size(3));
for i=1:Size(3)
    for j=1:number*number
        x_top_left=mod((j-1),number)*strides+1;
        y_top_left=floor((j-1)/number)*strides+1;
        patch_img=RainMap(x_top_left:x_top_left+window_size,y_top_left:y_top_left+window_size,i);
        direct(j,i)=est_theta_Fourier(patch_img);
    end
end
direct=direct(:);
direct=sort(direct);
%direct=direct(0.2*size(direct):size(direct));%some patch has no droprain,so discard 20% data
direction=median(direct(0.5*size(direct):size(direct)));

