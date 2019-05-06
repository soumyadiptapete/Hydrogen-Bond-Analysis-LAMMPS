function [rsqr_wh_exp,tau_wh_exp]= compute_acf_wh(hb_wh,max_lags,time_const,i) 
warning('off','all');   

start_point=hb_wh(1:20:100,i);
if sum(start_point)>0
        i
        ACF_wh=zeros(max_lags+1,1);
        ACF_wh=single(ACF_wh);
        temp=xcorr(hb_wh(:,i),max_lags,'unbiased');
        
        ACF_wh=temp(max_lags+1:2*max_lags+1);
        clear temp
        ACF_wh=(ACF_wh)/(sum(hb_wh(:,i).^2))/length(hb_wh(:,i));
      
        x=((1:length(ACF_wh))-1)*time_const;
        

        %exponential fit
        [f1,gof1]=fit(x',ACF_wh,'exp1','MaxIter',5000,'MaxFunEvals',5000);
        rsqr_wh_exp=gof1.adjrsquare;
        tau_wh_exp=-1/f1.b;
        clear ACF_wh
else
        

        rsqr_wh_exp=0;
        tau_wh_exp=0;
        
end

end