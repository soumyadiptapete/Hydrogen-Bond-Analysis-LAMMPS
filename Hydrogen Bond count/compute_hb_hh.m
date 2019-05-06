function x=compute_hb_hh(i,hydronium_mol,hydronium_len,dist_hb,angle_hb)
x=0;
    for j=1:hydronium_len
        hb_hh_count=0;
        if i~=j
            
            donor_h_index=[1 3 4];
            for k=donor_h_index
                hb_hh_count=hb_hh_count+1;
                
                di=i; ai=j;%di=donor index,ai=acceptor index
                dist=sum((hydronium_mol(di).atom_data(k,3:5)-hydronium_mol(ai).atom_data(2,3:5)).^2);
             
                if dist>=dist_hb(1)^2 && dist<=dist_hb(2)^2
                    v1=hydronium_mol(di).atom_data(k,3:5)-hydronium_mol(di).atom_data(2,3:5);
                    v2=hydronium_mol(di).atom_data(k,3:5)-hydronium_mol(ai).atom_data(2,3:5);
                    angle=acos(dot(v1, v2) / (norm(v1) * norm(v2)));
                    angle=rad2deg(angle);
                    if angle>=angle_hb
                        x=x+1;
                      
                    end
        
                end
            end
        end
    end
end