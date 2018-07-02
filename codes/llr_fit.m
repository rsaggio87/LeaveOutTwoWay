function sigma_predict= llr_fit(Lambda_P,Lambda_B,y,eta_h,subsample_llr_fit,K,movers,T,KGrid)

if nargin<=8
    KGrid=1000;
end


%% Set up dimensions
NT=size(eta_h,1);
maxT=max(T);

%% Run LLR
tic
sigma_i=y.*eta_h;
Pii=diag(Lambda_P);
Bii=diag(Lambda_B);
hBest=1/NT^(1/3);
if subsample_llr_fit == 0
    f = fit([Pii Bii],sigma_i,'lowess','Normalize','on','Span',hBest);
    sigma_predict = feval(f,[Pii, Bii]);
    selnan=isnan(sigma_predict);
    sigma_predict(selnan)=sigma_i(selnan);
end

if subsample_llr_fit == 1 && K==0
    sigma_predict=zeros(NT,1);
    for ti=1:maxT
            sel=(T==ti & ~movers);
            sigma_predict(sel)=mean(sigma_i(sel));
    end
    sigma_use=sigma_i(movers);
    Xuse=[Bii(movers) Pii(movers)];
    cellS=size(Xuse,1);
    hBest=1/cellS^(1/3);
    f = fit(Xuse,sigma_use,'lowess','Normalize','on','Span',hBest);
    sigma_predict(movers) = feval(f,[Bii(movers), Pii(movers)]);
    selnan=isnan(sigma_predict);
    sigma_predict(selnan)=sigma_i(selnan);
end

if subsample_llr_fit == 1 && K>0
    sigma_predict=zeros(NT,1);
    %LLR on Stayers.
    sigma_use=sigma_i(~movers);
    Xuse=[Bii(~movers) Pii(~movers)];
    cellS=size(Xuse,1);
    hBest=1/cellS^(1/3);
    f = fit(Xuse,sigma_use,'lowess','Normalize','on','Span',hBest);
    sigma_predict(~movers) = feval(f,[Bii(~movers), Pii(~movers)]);
    %LLR on Movers.
    Xuse=[Bii(~movers) Pii(~movers)];
    cellS=size(Xuse,1);
    hBest=1/cellS^(1/3);
    f = fit(Xuse,sigma_use,'lowess','Normalize','on','Span',hBest);
    sigma_predict(movers) = feval(f,[Bii(movers), Pii(movers)]);
    selnan=isnan(sigma_predict);
    sigma_predict(selnan)=sigma_i(selnan);
end

if subsample_llr_fit == 2 && K==0
    sigma_predict=zeros(NT,1);
    for ti=1:maxT
            sel=(T==ti & ~movers);
            sigma_predict(sel)=mean(sigma_i(sel));
    end
    Bii=Bii(movers);
    Pii=Pii(movers);
    g_B = group_equally(Bii, KGrid);
    g_P = group_equally(Pii, KGrid);
    [~,~,qq]=unique([g_B,g_P],'rows');
    weight=accumarray(qq,1);
    Bii_bin=accumarray(qq,Bii,[],@mean);
    Pii_bin=accumarray(qq,Pii,[],@mean);
    sigma_use=accumarray(qq,sigma_i(movers),[],@mean);
    Xuse=[Bii_bin Pii_bin];
    cellS=size(Bii_bin,1);
    hBest=1/cellS^(1/3);
    f = fit(Xuse,sigma_use,'lowess','Normalize','on','Span',hBest,'Weights',weight);
    sigma_predict(movers) = feval(f,[Bii, Pii]);
    selnan=isnan(sigma_predict);
    sigma_predict(selnan)=sigma_i(selnan);
end

if subsample_llr_fit == 2 && K>0
    sigma_predict=zeros(NT,1);
    Bii_old=Bii;
    Pii_old=Pii;
    %Binned LLR on stayers.
    Bii=Bii_old(~movers);
    Pii=Pii_old(~movers);
    sigma_use=sigma_i(~movers);
    g_B = group_equally(Bii, KGrid);
    g_P = group_equally(Pii, KGrid);
    [~,~,qq]=unique([g_B,g_P],'rows');
    weight=accumarray(qq,1);
    Bii_bin=accumarray(qq,Bii,[],@mean);
    Pii_bin=accumarray(qq,Pii,[],@mean);
    sigma_use=accumarray(qq,sigma_use,[],@mean);
    Xuse=[Bii_bin Pii_bin];
    cellS=size(Bii_bin,1);
    hBest=1/cellS^(1/3);
    f = fit(Xuse,sigma_use,'lowess','Normalize','on','Span',hBest,'Weights',weight);
    sigma_predict(~movers) = feval(f,[Bii, Pii]);
    %Binned LLR on movers.
    Bii=Bii_old(movers);
    Pii=Pii_old(movers);
    sigma_use=sigma_i(movers);
    g_B = group_equally(Bii, KGrid);
    g_P = group_equally(Pii, KGrid);
    [~,~,qq]=unique([g_B,g_P],'rows');
    weight=accumarray(qq,1);
    Bii_bin=accumarray(qq,Bii,[],@mean);
    Pii_bin=accumarray(qq,Pii,[],@mean);
    sigma_use=accumarray(qq,sigma_use,[],@mean);
    Xuse=[Bii_bin Pii_bin];
    cellS=size(Bii_bin,1);
    hBest=1/cellS^(1/3);
    f = fit(Xuse,sigma_use,'lowess','Normalize','on','Span',hBest,'Weights',weight);
    sigma_predict(movers) = feval(f,[Bii, Pii]);
    selnan=isnan(sigma_predict);
    sigma_predict(selnan)=sigma_i(selnan);
end

if subsample_llr_fit == 3
    g_B = group_equally(Bii, KGrid);
    g_P = group_equally(Pii, KGrid);
    [~,~,qq]=unique([g_B,g_P],'rows');
    weight=accumarray(qq,1);
    Bii_bin=accumarray(qq,Bii,[],@mean);
    Pii_bin=accumarray(qq,Pii,[],@mean);
    sigma_use=accumarray(qq,sigma_i,[],@mean);
    Xuse=[Bii_bin Pii_bin];
    cellS=size(Bii_bin,1);
    hBest=1/cellS^(1/3);
    f = fit(Xuse,sigma_use,'lowess','Normalize','on','Span',hBest,'Weights',weight);
    sigma_predict= feval(f,[Bii, Pii]);
    selnan=isnan(sigma_predict);
    sigma_predict(selnan)=sigma_i(selnan);
end


disp('Time to Perform non-parametric fit')
toc
disp('Fraction of negative predicted values for tilde sigma_i^2')
mean(sigma_predict<0)
end    


