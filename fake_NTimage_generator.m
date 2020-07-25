clc
close

%{ 
TODO: 
Add ramdom variations in brightness along nanotubes
    Make the shape more variable - by adding a small diffusion step before
    noise is added to arms - to increase the size and randomize shape a bit

Change the function that write muiltiple images so that it prints the pixel
data too

Create function/find a way to add the nanotube fragments to the image
%}

%{ 
Pixel Label ID labels: 
1 - nanotube
2 - fragment (nanotube fragment)
3 - background
%}

tic
resolution = 2600; %not the actual last resolution, only used for scaling



[image,pixel_labels]=generate_image(resolution);

imshow(image)
imwrite(image,'/home/mmbl/Documents/dsa/Segmentation_NN/fakedata/fakedatatest.tif'); %writes image to right folder
imwrite(pixel_labels,'/home/mmbl/Documents/dsa/Segmentation_NN/fakedata/fakedatatest_labels.tif'); %writes image to right folder

%}


%{
N=100 %number of images that will be generated
%binSize=5 %angle bin size
%numOfBins=360/binSize %number of angles that 180 will be split into

save_path = '../fakedata/training_images';
nt_image_path = strcat(save_path,'/nanotube_images'); %path where generated nanotube images will be saved
pld_path = strcat(save_path,'/pixel_label_data'); % folder where the pixel label data will be saved

makeFolders(save_path,nt_image_path,pld_path) %will delete old training data and prepare new folders


parfor j=1:N
	image=generate_image(resolution); %generates the artificial nanotube image
	savePath=strcat(nt_image_path,'/',num2str(j),'.tif'); %path to where the image will be saved
	imwrite(image,savePath); %writes image to right folder
    
    %{
    %progress bar display    
    clc
    bi=round((j/N)*50);
    bar=repmat('|',1,bi);
    empbar=repmat('_',1,50-bi);
    disp(strcat(bar,empbar,num2str(round(((j/N)*100),3)),"%"))
    pause(1)
    %}
end
    
runTime=toc;

strcat(num2str(N),' images in: ',num2str(toc/60),' mins')


%}

%________________________FUNCTIONS_______________________________________________
%________________________________________________________________________________
%________________________________________________________________________________
function [imageout,pixel_dataout]=generate_image(resolution)
image=zeros(resolution,'uint8'); %blank canavas for 8bit image
pixel_data=zeros(resolution,'uint8'); %will contain the pixel label data


[image,pixel_data1]=nanoTubes(image,pixel_data,resolution);
[image,pixel_data2]=nanoTubes_small(image,pixel_data,resolution);
%image=nanotube_fragments(image,resolution);

image=diffusion(image,resolution); %applies diffusion

%mask=zeros(resolution,'uint8');
%mask(image>0)=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%mask=double(image)./double(max(max(image)));

%mask(mask>0) = mask(mask>0) + rand*0.5;

%mask(mask>1)=1;
%image = imnoise(image,'gaussian',0.5,0.5);
%image=image.*mask
%image=uint8(double(image).*mask);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
final_resolution = 1200;
image=imresize(image,[final_resolution,final_resolution]);

pixel_data=pixel_data1+pixel_data2;

noise_mean = rand*0.25';
image = imnoise(image,'gaussian',noise_mean,noise_mean*rand*.01);

imageout=image;
pixel_dataout=pixel_data;

end

function [imageout,pixel_dataout]=nanoTubes(image,pixel_data,resolution) 
%This function will generate the well formed nanotubes

final_resolution = 1200;
limit_initial = 0.99; %limit that determines whether another tube wil be drawn in the loop
limit_factor = 0.999;  %factor by which the limit is decreased every iteration
brightness_min = 40; %minimum brightness of tubes drawn
brightness_variable = 205; %variability of tube (brigtness_min + rand*brightness_variable)
size_min = 32; %minimum size
size_variable = 43; %variability in size (same as brightness)
label = 1; %label for viable nanotubes

[image,pixel_data]=nanoTubes_general(image,pixel_data,label,resolution,limit_initial,limit_factor,brightness_min,brightness_variable,size_min,size_variable);

pixel_data =imresize(pixel_data,[final_resolution,final_resolution]); %resizing pixel labels to final resolution
pixel_data(pixel_data>0)=label;

imageout=image;
pixel_dataout=pixel_data;

end 

function [imageout,pixel_dataout]=nanoTubes_small(image,pixel_data,resolution) 
%This function will generate the well formed nanotubes

final_resolution = 1200;
limit_initial = 0.99; %limit that determines whether another tube wil be drawn in the loop
limit_factor = 0.9;  %factor by which the limit is decreased every iteration
brightness_min = 40; %minimum brightness of tubes drawn
brightness_variable = 205; %variability of tube (brigtness_min + rand*brightness_variable)
size_min = 8; %minimum size
size_variable = 22; %variability in size (same as brightness)
label = 2; %label for nanotubes that are too small for analysis

[image,pixel_data]=nanoTubes_general(image,pixel_data,label,resolution,limit_initial,limit_factor,brightness_min,brightness_variable,size_min,size_variable);

pixel_data =imresize(pixel_data,[final_resolution,final_resolution]); %resizing pixel labels to final resolution
pixel_data(pixel_data>0)=label;

imageout=image;
pixel_dataout=pixel_data;
end 


function [imageout,pixel_dataout]=nanoTubes_general(image,pixel_data,label,resolution,limit_initial,limit_factor,brightness_min,brightness_variable,size_min,size_variable)
%This function will generate the well formed nanotubes

luck=1; %so that the odds are rolled at least once
num=0; %will count the number of nanotubes generated
limit=limit_initial; %initial threashold to determine whether arm will be drawn

%x & y positions of the start of the nanotube random number sampled from a gaussian
positionx_mean=resolution/2*(1/2 +rand); 
positiony_mean=resolution/2*(1/2 +rand); 
position_std=resolution/8;

while luck==1

    if rand<limit %chance that a tube will be drawn
        num=num+1; %adds to count
        center=[round(normrnd(positionx_mean,position_std)),round(normrnd(positiony_mean,position_std))]; %randomized starting position
        angleArm=360*rand; %random direction
        brightness=brightness_variable*rand+brightness_min; %brightness of the nanotube
        size=rand*size_variable+size_min; % size of the nanotube (in pixels)
        [image,pixel_data]=arm(image,pixel_data,label,angleArm,brightness,size,center,brightness,resolution);
    else
        luck=0;

    end
    limit=limit*limit_factor; %lowers the limit every time a new tube is drawn
end

num;
imageout=image;
pixel_dataout=pixel_data;

end 


function out=nanotube_fragments(image,resolution) 
%This function will generate nanotube fragments that appear in the frame

luck=1; %so that the odds are rolled at least once
num=0; %will count the number of nanotube fragements generated
limit=0.999; %initial threashold to determine whether fragement will be drawn
position_mean=resolution/2;
position_std=resolution/8;

while luck==1

    if rand<limit %chance that a tube will be drawn
        num=num+1; %adds to count
        center=[round(normrnd(position_mean,position_std)),round(normrnd(position_mean,position_std))]; %randomized starting position
        angleArm=360*rand; %random direction
        brightness=25*100*rand; %brightness of the fragment
        size=3 +rand*2; % size of the nanotube (in pixels)
        image=arm(image,angleArm,brightness,size,center,brightness,resolution);
    else
        luck=0;

    end
    limit=limit*0.999; %lowers the limit every time a new tube is drawn
end

num;
out=image;

end 

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



function [imageout, pixel_dataout]=arm(image,pixel_data,label,angleInitial,brightness,size,center,vertexBrightness,resolution)

point=[0,0];%because we must start random walk from the center before we perform the rotation
Length=size; %size of arm-depending on input and random variable
angle=angleInitial;%initial angle of the arm-input

%from Amber's filament simulation program
n_steps=round(Length)-1;%number of pixels
sigma_i=3.7;%the sigma of the Gaussian distribution from which we draw bend angle of each step

theta_i=normrnd(0,sigma_i,[1,n_steps]);
theta_i(1)=0;

running_sum=cumsum(theta_i);%the cumulative sum of angles for each filament
x=zeros(1,n_steps);%make empty arrays for storing x and y coordinates for random walks
y=zeros(1,n_steps);

%in this for-loop, fill x and y arrays using the simulated angles
for i=1:n_steps
    if i==1
        x=point(1)+cosd(running_sum(i));
        y=point(2)+sind(running_sum(i));
    else
        x(i)=x(i-1)+cosd(running_sum(i));
        y(i)=y(i-1)+sind(running_sum(i));
    end
end

%smooth the x-y curve.
for counter=1:30
    x=smooth(x);y=smooth(y);
end

if length(x)<8%for "backtubes", Length could be very short, so we don't need to do this first rotation to "correct"
    R = [cosd(angle) -sind(angle); sind(angle) cosd(angle)];%rotation matrix
    rotated = (R*[x,y]')';

    %translate the filament so it starts at the desired point.
    x=rotated(:,1)+center(1);
    y=rotated(:,2)+center(2);
else
    %now measure the inital angle
    for i=1:3%look at three vectors each of Length 2 (or 3, between 3 neighboring dots)
        vector=[x(i+4)-x(i+2),y(i+4)-y(i+2)];
        initial_angle(i)=atand(vector(2)/vector(1));
    end
    off_angle=mean(initial_angle);%calculate the inital angle that the filament is off by
    R = [cosd(-off_angle) -sind(-off_angle); sind(-off_angle) cosd(-off_angle)];%rotation matrix so the filament is at 0 degrees
    rotated = (R*[x,y]')';

    %rotate to achieve the desired orientation. Note the handedness.
    R = [cosd(angle) -sind(angle); sind(angle) cosd(angle)];%rotation matrix
    rotated_again = (R*rotated')';

    %translate the filament so it starts at the desired point.
    x=rotated_again(:,1)+center(1);
    y=rotated_again(:,2)+center(2);
end

evalue=100; %1.8*rand+0.2; %value used in the exponetial that determines the arms initial brightness
for i=1:n_steps 
    if i==1 %does not mark first point to create the characteristic void in the middle of the nunchuck
        continue
    end
    if round(x(i))-1<1 || round(x(i))+1>resolution || round(y(i))-1<1 || round(y(i))+1>resolution
        break
    end

    %if brightness2==-1
    %    image(round(y(i))-1:round(y(i))+1,round(x(i))-1:round(x(i))+1)=uint8(brightness);%marks a 2x2 area at next point with input brightness
    %else
        
        brightnessLevel=(brightness-vertexBrightness)*(1-exp(-(evalue)*(i-1)))+vertexBrightness;
        image(round(y(i))-1:round(y(i))+1,round(x(i))-1:round(x(i))+1)=uint8(brightnessLevel);%marks a 2x2 area at next point with input brightness
        pixel_data(round(y(i))-2:round(y(i))+2,round(x(i))-2:round(x(i))+2)=uint8(250); %marks the same pixels of the arm with the label on the pixel data image

    %end
    
end

imageout=image;
pixel_dataout=pixel_data;
end

function makeFolders(save_path,nt_image_path,pld_path)
	%This function will make the folders where the images and pixel data will be stored


if exist(save_path,'dir')
    disp("True")
    rmdir(save_path)
end

mkdir(save_path) %makes inital folder

mkdir(nt_image_path)
mkdir(pld_path)

end


%{


for i=1:100
luck=1; %so that the odds are rolled at least once

limit=0.99; %initial threashold to determine whether arm will be drawn
num=0; %will count the number of nanotubes generated
while luck==1

    if rand<limit %chance that a tube will be drawn
        num=num+1; %adds to count

    else
        luck=0;

    end
    limit=limit*0.90; %lowers the limit every time a new tube is drawn
end
list2(i)=num;
end

function out=nanoTubes_small(image,resolution) 
%This function will generate nanoteubes that are too small for analysis

luck=1; %so that the odds are rolled at least once
num=0; %will count the number of nanotubes generated
limit=0.99; %initial threashold to determine whether arm will be drawn

position_mean=resolution/2;
position_std=resolution/8;

while luck==1

    if rand<limit %chance that a tube will be drawn
        num=num+1; %adds to count
        center=[round(normrnd(position_mean,position_std)),round(normrnd(position_mean,position_std))]; %randomized starting position
        angleArm=360*rand; %random direction
        brightness=205*rand+40; %brightness of the nanotube
        size=rand*22+8; % size of the nanotube (in pixels)
        image=arm(image,angleArm,brightness,size,center,brightness,resolution);
    else
        luck=0;

    end
    limit=limit*0.90; %lowers the limit every time a new tube is drawn
end

num;
out=image;

end 


function out=nanoTubes(image,resolution) 
%This function will generate the well formed nanotubes

luck=1; %so that the odds are rolled at least once
num=0; %will count the number of nanotubes generated
limit=0.99; %initial threashold to determine whether arm will be drawn

position_mean=resolution/2;
position_std=resolution/8;

while luck==1

    if rand<limit %chance that a tube will be drawn
        num=num+1; %adds to count
        center=[round(normrnd(position_mean,position_std)),round(normrnd(position_mean,position_std))]; %randomized starting position
        angleArm=360*rand; %random direction
        brightness=205*rand+40; %brightness of the nanotube
        size=rand*43+32; % size of the nanotube (in pixels)
        image=arm(image,angleArm,brightness,size,center,brightness,resolution);
    else
        luck=0;

    end
    limit=limit*0.90; %lowers the limit every time a new tube is drawn
end

num;
out=image;

end 

function makeFolders(binSize)
numOfFolders=360/binSize;

if exist('FakeNunchuckImages')==7
    %disp("True")
    rmdir FakeNunchuckImages s
end
mkdir FakeNunchuckImages %makes inital folder

for i=1:numOfFolders
    angle=(i-1)*binSize+(binSize/2)-180; %so that folder name is in the middle of the bin
    %disp(angle)
    folderData=['FakeNunchuckImages/' num2str(angle)]; %path of new folder with name
    mkdir(folderData)
end
end

function out=reformatImages(img)
    img(201:227,201:227)=uint8(0);
    img = cat(3, img, img, img);
    out=img;
end

function out=noise(image,resolution)
resolution=1200
image=imresize(image,[1200,1200]);
%next loop adds a random value (20-30) to each pixel to simulate noise
noiseLevel=randi(20,'uint8')+20;
noiseBaseline=randi(noiseLevel/2,'uint8')+noiseLevel/3;
%noiseBaseline=noiseLevel/3
for i=1:resolution
    for j=1:resolution
        %image(i,j)=image(i,j)+randi(40,'uint8')+uint8(40);
        image(i,j)=image(i,j)+randi(noiseLevel,'uint8')+uint8(noiseBaseline);
    end 
end

%next loop will increase each pixel to simulate differences in image
%brightness
brightnessIncrease=uint8(40*rand); %random amount that each pixel will be increased by
for i=1:resolution
    for j=1:resolution
        image(i,j)=image(i,j)+brightnessIncrease;
    end 
end
out=image;
end


function [out,center]=nunchuck(image,angleNun,resolution) %creates a nunchuch in the image/


center=[600,600]; %randomized center (40x40 square)
angleArm=360*rand; %intial angle at which first arm will come out

brightness1=randi(195,'uint8')+uint8(60);%brightness for the first(brighter arm) arm
brightness=double(brightness1);
brightness=normrnd(brightness/2,10); %picks second brightness usign a normal distribution
brightness2=uint8(brightness);
basebrightness=brightness2-0.8*rand*brightness2;

image1=arm(image,angleArm,brightness1,20,center,basebrightness,resolution); %draws first arm (originally 35)
angleNun=180-angleNun; %switched how the angle is defined 
angleArm=angleArm+angleNun; %angle of second arm according to nunchuck angle desired
%brightness=brightness/2; %brightness for the second arm should be half of brighter arm
image2=arm(image,angleArm,brightness2,20,center,basebrightness,resolution); %creates the new arm at desired angle (originally 25)
    %dimmer and shorter(most likely) -last two numbers-to recreate single labeled arm
    %image1(center(2)-20:center(2)+20,center(1)-20:center(1)+20)=uint8(50);
out=image1+image2; %adds two images to create a superposition of the arms


end 

%}

%{
Things of note:
-Includes randomized arm sizes (taking into account tendency for dimmer arm
to be shorter)
-Randomized arms curvature
-Difference in arm bnrightness
-Simulated noise
-Randomized image brightness
-Randomized nunchuck center position
%}