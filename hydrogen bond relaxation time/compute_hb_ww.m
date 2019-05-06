function x=compute_hb_ww(i,water_mol,water_len,dist_hb,angle_hb)
% computes hydrogen bonds between water molecules with water as both donor or acceptor
hb_ww_count=0;
x=zeros(1,(water_len-1)*2);
for j=1:water_len

        if i~=j

            di=i; ai=j; %di-donor index, ai-acceptor index, donor is the molecule which gives the h-atom
            donor_h_index=[1 3];% water atom data 1st and last atoms are hydrogen atoms
            for k=donor_h_index
                hb_ww_count=hb_ww_count+1;
                       
                    dist=sum((water_mol(di).atom_data(k,3:5)-water_mol(ai).atom_data(2,3:5)).^2);
                if dist>=dist_hb(1)^2 && dist<=dist_hb(2)^2 %distance between donor hydrogen and acceptor oxygen
                    v1=water_mol(di).atom_data(k,3:5)-water_mol(di).atom_data(2,3:5);%vector connecting donor hydrogen-donor oxygen
                    v2=water_mol(di).atom_data(k,3:5)-water_mol(ai).atom_data(2,3:5);%vector connecting donor hydrogen-acceptor oxygen
                    angle=acos(dot(v1, v2) / (norm(v1) * norm(v2)));
                    angle=rad2deg(angle);
                    if angle>=angle_hb
                        x(hb_ww_count)=1;
                       
                    end
                else
                    x(hb_ww_count)=0;
                end
            end
        end

end
end
