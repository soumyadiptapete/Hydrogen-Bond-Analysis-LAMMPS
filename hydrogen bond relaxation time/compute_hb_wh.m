function x=compute_hb_wh(i,water_mol,hydronium_mol,hydronium_len,dist_hb,angle_hb)
x=zeros(1,hydronium_len*5);
hb_wh_count=0;
    for j=1:hydronium_len
       
        donor_h_index=[1 3];%water is donor
        for k=donor_h_index
            hb_wh_count=hb_wh_count+1;
            %donor is water (index i) and acceptor is hydronium
            %(index j)
            di=i;ai=j;
            dist=sum((water_mol(di).atom_data(k,3:5)-hydronium_mol(ai).atom_data(2,3:5)).^2);
            if dist>=dist_hb(1)^2 && dist<=dist_hb(2)^2 %distance between donor hydrogen and acceptor oxygen
                v1=water_mol(di).atom_data(k,3:5)-water_mol(di).atom_data(2,3:5);%vector connecting donor hydrogen-donor oxygen
                v2=water_mol(di).atom_data(k,3:5)-hydronium_mol(ai).atom_data(2,3:5);%vector connecting donor hydrogen-acceptor oxygen
                angle=acos(dot(v1, v2) / (norm(v1) * norm(v2)));
                angle=rad2deg(angle);
                if angle>=angle_hb
                   x(hb_wh_count)=1;
                  
                end
            else
                x(hb_wh_count)=0;
            end
        end
        
        donor_h_index=[1 3 4];%hydronium is donor
        
        for k=donor_h_index
            hb_wh_count=hb_wh_count+1;
            %acceptor is water and donor is hydronium
            di=j;ai=i;
            dist=sum((hydronium_mol(di).atom_data(k,3:5)-water_mol(ai).atom_data(2,3:5)).^2);
            if dist>=dist_hb(1)^2 && dist<=dist_hb(2)^2
                v1=hydronium_mol(di).atom_data(k,3:5)-hydronium_mol(di).atom_data(2,3:5);
                v2=hydronium_mol(di).atom_data(k,3:5)-water_mol(ai).atom_data(2,3:5);
                angle=acos(dot(v1, v2) / (norm(v1) * norm(v2)));
                angle=rad2deg(angle);
                if angle>=angle_hb
                 x(hb_wh_count)=1;
                end
            else
               x(hb_wh_count)=0;
            end
        end
    end
end