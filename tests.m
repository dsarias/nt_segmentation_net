%trash code
tic
clc
image= image_backup;
%applies initial diffusion to make area around nunchuck also affected by
%turbulance
image=diffusion(image,1200,0,8);
%creates a turbulance mask
turbulance_value = round(3*rand)+1;
turbulance_weights = turbulance(1200,turbulance_value);
%Applies turbulance mask
image =double(image) .* turbulance_weights/max(max(turbulance_weights));
%Fixes values of the image to be within 8-bit ranges with a randomized
%maximum brightness
brightness_factor = 45 + rand*210;
image=image/max(max(image))*brightness_factor;
image=uint8(image); %converts image back to 

%to make nanotubes mesh with the noise better - stand out less
image=diffusion(image,1200,1,4); %applies diffusion

%applies noise to the whole image
noise_mean = rand*0.25';
image = imnoise(image,'gaussian',noise_mean,noise_mean*rand*.01);

imshow(image)


toc

function image = turbulance(resolution,zoom_initial)
    image = zeros(resolution, resolution);
    zoom = zoom_initial;
  while zoom >= 1
    noise_image = zeros(resolution/zoom, resolution/zoom);
    noise_image = imnoise(noise_image,'gaussian',0.5,0.5);
    noise_image = imresize(noise_image,[resolution, resolution]);
    image = image + noise_image * zoom;
        %imshow(image)
    pause(1)  

    zoom = zoom/2.0;
  end

  image = image/max(max(image)) *255;
  
  

end
  
function out=diffusion(P,resolution,nt_min, nt_variable)
%solves the 2D diffusion equation to spread out the nunchuck
Nt=round(nt_min+rand*nt_variable);%round(15*rand+20); %number of time steps (random with min of 15)
dt=0.1; %time step size-arbitrary units
D=0.8; %diffusion coefficient
dx=1; %step size

P=im2double(P); %converts imnage to double in order to make the calculations

PNew=P; %douplicates image

for k=1:Nt %time loop
    for i=2:resolution-1 %space loops
        for j=2:resolution-1
            %if(ismember([j i],dots ,'rows')==1)
             %   continue
            %end
            PNew(i,j)=double(P(i,j))+(D*dt)/(dx^2)*(P(i+1,j)+P(i-1,j)+P(i,j+1)+P(i,j-1)-4*P(i,j));
            %finite difference method for diffusion equation
        end
    end
    P=PNew; %passes new values of image for next loop
end

P=im2uint8(P); %converts image back to 8bit

out=P;
end 