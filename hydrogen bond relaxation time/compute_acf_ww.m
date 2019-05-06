function [rsqr_ww_exp,tau_ww_exp]= compute_acf_ww(hb_ww,max_lags,time_const,i)    
warning('off','all');
start_point=hb_ww(1:20:100,i);
if sum(start_point)>0
        i
        ACF_ww=zeros(max_lags+1,1);
        ACF_ww=single(ACF_ww);
        temp=xcorr(hb_ww(:,i),max_lags,'unbiased');
        ACF_ww=temp(max_lags+1:2*max_lags+1);
        clear temp
        ACF_ww=(ACF_ww)/(sum(hb_ww(:,i).^2))/length(hb_ww(:,i));
       
        x=((1:length(ACF_ww))-1)*time_const;
 
        %exponential fit
        [f1,gof1]=fit(x',ACF_ww,'exp1','MaxIter',5000,'MaxFunEvals',5000);
        rsqr_ww_exp=gof1.adjrsquare;
        tau_ww_exp=-1/f1.b;
        
      clear ACF_ww
else
       

        rsqr_ww_exp=0;
        tau_ww_exp=0;
        
end

end