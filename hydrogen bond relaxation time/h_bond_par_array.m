clc;close all; clear all;
warning('off','all');
n_cores=32;parpool(n_cores);
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

fid=fopen('wrap_hb.lammpstrj','r');

path=strcat('/home/sengupt/Nafion_ph/Nafion_',num2str(DP),'/lambda_',num2str(lambda),'/HB_acf')
cd(path);

chains=20;
dist_hb=[1.59 2.27];%distance-criteria for h-bonding
angle_hb=140;%angle-criteria for h-bonding

max_ts=1000;%Maximum time steps 1000*0.005ps=5 ps
n_skip_step=2;%skip very other step

atoms_per_chain=692-DP;


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

%number of hydrogen bonds
n_hb_ww=nchoosek(water_len,2)*4;
if DP~=0
n_hb_wh=water_len*hydronium_len*5;
end

%%preallocate the hydrogen bond arrays 
hb_ww=zeros(max_ts,n_hb_ww,'uint8');
if DP~=0
hb_wh=zeros(max_ts,n_hb_wh,'uint8');
end

%preallocate water residence time arrays

rsqr_ww_exp=zeros(n_hb_ww,1,'single');
tau_ww_exp=zeros(n_hb_ww,1,'single');


%%start reading data and computing hydrogen bond existence data for all
%%possible hydrogen bonds
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
     
    
    if strcmp(id,'ITEM: ATOMS id type x y z q ') && mod(t,n_skip_step)==0
        
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
        x=zeros(water_len,(water_len-1)*2);
        parfor (i=1:water_len,n_cores)
            
            x(i,:)=compute_hb_ww(i,water_mol,water_len,dist_hb,angle_hb);
            
        end
        
        hb_ww(cnt_steps,:)=x(:);
        clear x;
        
        if DP~=0
            % %water hydronium hydrogen bonds
            x=zeros(water_len,hydronium_len*5);
            parfor (i=1:water_len,n_cores)
                
                x(i,:)=compute_hb_wh(i,water_mol,hydronium_mol,hydronium_len,dist_hb,angle_hb);
                
            end
            
            hb_wh(cnt_steps,:)=x(:);
            clear x;
            
            
        end
        
        clear water_mol
        
        if DP~=0
            clear hydronium_mol
        end
    end
    
    
end

%ACF_computation
max_lags=max_ts-5;
time_const=(5/1000)*n_skip_step; %sampling interval in ps
disp('water_water_acf_starting..........');
%water-water

for i=1:n_hb_ww
    [rsqr_ww_exp(i,1),tau_ww_exp(i,1)]= compute_acf_ww(hb_ww,max_lags,time_const,i);
end

clear hb_ww;

save rtf_ww_exp.mat -v7.3 rsqr_ww_exp tau_ww_exp 

if DP~=0
    %water-hydronium
    n_hb_wh=size(hb_wh,2);
    rsqr_wh_exp=zeros(n_hb_wh,1);
    tau_wh_exp=zeros(n_hb_wh,1);
 
    disp('water_hydronium_acf_starting..........');
    
    for i=1:n_hb_wh
        [rsqr_wh_exp(i,1),tau_wh_exp(i,1)]= compute_acf_wh(hb_wh,max_lags,time_const,i);
    end
    
    clear hb_wh;
  
    save rtf_wh_exp.mat -v7.3 rsqr_wh_exp tau_wh_exp
    
    
end




