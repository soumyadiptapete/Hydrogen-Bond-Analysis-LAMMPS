%% Code for counting the number of hydrogen bonds existing in between water and hydronium molecules in a LAMMPS molecular dynamics trajectory
%% variables DP and lambda are obtained from the address of the file. They can also be input directly
clc;close all; clear all;
n_cores=32;
parpool(n_cores);
%getting DP and lambda from path address automatically
path=cd;

DP_str='ab';lambda_str='cd';
k=strfind(path,'Nafion_');
x=length('Nafion_');
DP_str(1)=path(k(2)+x);
if path(k(2)+x+1)~='/'
    DP_str(2)=path(k(2)+x+1);
else
    DP_str(2)=[];
end
DP=str2num(DP_str)

k=strfind(path,'lambda_');
x=length('lambda_');
lambda_str(1)=path(k+x);
if path(k+x+1)~='/'
    lambda_str(2)=path(k+x+1);
else
    lambda_str(2)=[];
end
lambda=str2num(lambda_str)

if DP==10
    path=strcat('/home/sengupt/Nafion_ph/Nafion_',num2str(DP),'/lambda_',num2str(lambda))
else
    path=strcat('/home/sengupt/Nafion_ph/Nafion_',num2str(DP),'/lambda_',num2str(lambda),'/restart')
end
cd(path);

fid=fopen('wrap.lammpstrj','r');

path=strcat('/home/sengupt/Nafion_ph/Nafion_',num2str(DP),'/lambda_',num2str(lambda),'/HB_cnt')
cd(path);


dist_hb=[1.59 2.27];%distance-criteria for h-bonding
angle_hb=140;%angle-criteria for h-bonding
max_lags=1;
max_ts=3000;
chains=20;

atoms_per_chain=692-DP;

sulfur_atoms=(10-DP)*chains;


hydronium_mols=DP*chains;
water_mols=lambda*chains*10-hydronium_mols;

water_start=atoms_per_chain*chains+1;
water_end=water_start+water_mols*3-1;


hydronium_start=water_end+1;
hydronium_end=hydronium_start+hydronium_mols*4-1;

hydrogen_atoms=hydronium_mols*3;
oxygen_atoms=hydronium_mols;



water=water_start:water_end; %water atoms index
oxygen=water_start+1:3:water_end-1; %oxygen of water atoms index

hydronium=hydronium_start:hydronium_end;

water_len=length(water)/3;
hydronium_len=length(hydronium)/4;
 
hb_ww_cnt=zeros(max_ts,1);
if DP~=0
hb_wh_cnt=zeros(max_ts,1);
hb_hh_cnt=zeros(max_ts,1);
end

max_hb_ww_cnt=nchoosek(water_len,2)*4;
if DP~=0
max_hb_wh_cnt=water_len*hydronium_len*5;
max_hb_hh_cnt=nchoosek(hydronium_len,2)*6;
end

t=0; cnt_steps=0;
while ~feof(fid) && cnt_steps<max_ts


    id=fgetl(fid);
    if strcmp(id,'ITEM: TIMESTEP')==1
        id=fgetl(fid);
        t=t+1;
        Timestep=str2num(id);
    end
    
    if strcmp(id,'ITEM: NUMBER OF ATOMS')==1
        
        id=fgetl(fid);
        
    end
    
    if strcmp(id,'ITEM: BOX BOUNDS pp pp pp')==1 || strcmp(id,'ITEM: BOX BOUNDS pp pp ff')
        
        id=fgetl(fid);id=fgetl(fid);id=fgetl(fid);
        
    end
     
     if strcmp(id,'ITEM: ATOMS id type x y z q ')
             
        cnt_steps=cnt_steps+1;
        cnt_steps
        %skip lines for all atoms before the first water atom
        for i=1:water(1)-1
            id=fgetl(fid);
        end
      
        water_count=0;
        
        %read the data from a particular timestep for waters
        for i=1:length(water)/3 %loop for the number of water molecules
         
            water_count=water_count+1; % count of water molecules
            
            for j=1:3
                id=fgetl(fid);
                water_mol(water_count).atom_data(j,:)=str2num(id);
            end

        end
 
 if DP~=0
        hydronium_count=0;
        
        %read the data from a particular timestep for hydronium
        for i=1:length(hydronium)/4 %loop for the number of hydronium molecules
         
            hydronium_count=hydronium_count+1; % count of hydronium molecules
            
            for j=1:4
                id=fgetl(fid);
                hydronium_mol(hydronium_count).atom_data(j,:)=str2num(id);
            end

        end
 end
        
        
   
%water-water hydrogen bonds

parfor (i=1:water_len,n_cores)
 
x(i)=compute_hb_ww(i,water_mol,water_len,dist_hb,angle_hb);

end
 
hb_ww_cnt(cnt_steps,1)=sum(x);
clear x;

if DP~=0    
% %water hydronium hydrogen bonds

parfor (i=1:water_len,n_cores)
   
  x(i)=compute_hb_wh(i,water_mol,hydronium_mol,hydronium_len,dist_hb,angle_hb);  

end

hb_wh_cnt(cnt_steps,1)=sum(x);
clear x;

%  
% %hydronium-hydronium hydrogen bonds

parfor (i=1:hydronium_len,n_cores)
    
x(i)=compute_hb_hh(i,hydronium_mol,hydronium_len,dist_hb,angle_hb);

end

hb_hh_cnt(cnt_steps,1)=sum(x);
clear x;
end
 
clear water_mol 

if DP~=0
clear hydronium_mol
end
     end
   
     
end

hb_ww_cnt_norm=hb_ww_cnt/max_hb_ww_cnt;
if DP~=0
hb_wh_cnt_norm=hb_wh_cnt/max_hb_wh_cnt;
hb_hh_cnt_norm=hb_hh_cnt/max_hb_hh_cnt;
end
if DP~=0
save HB_cnt.mat -v7.3 hb_ww_cnt hb_wh_cnt hb_hh_cnt hb_ww_cnt_norm hb_wh_cnt_norm hb_hh_cnt_norm
else
save HB_cnt.mat -v7.3 hb_ww_cnt hb_ww_cnt_norm 
end

