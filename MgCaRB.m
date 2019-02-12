%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MgCaRB: Compute salinity and pH adjusted planktonic Mg/Ca SST
%
% age in ka
% Mg/Ca in mmol/mol. If one column then analytical uncertainty prescribed to be +/-2%
% species ID (1=G. ruber, 2=T. sacculifer, 3=G. bulloides, 4=O.universa, 5='generic')
% relative - absolute (enter 0) or delta temperature (enter start and end age)
% TAlk - modern site alkalinity
% DpCO2 - modern site CO2 disequilibrium if CO2 is to be used (or 0)
% d11B - d11B pH record if available as d11BOH4 (enter 0 if pCO2 is to be
%    used). Can be two column with uncertainties.
% sal - modern site salinity or a vector of salinity values
% unc - Monte-Carlo uncertainty simulation?
% plot - plot data? (1/0)
%
function MgTemp = MgCaRB(age,MgCa,species,relative,TAlk,DpCO2,d11B,sal,unc,fig)

load('lookup_data','seaLevel','pCO2')

if numel(sal)==1    % if salinity data not prescribed
    salinity = interp1(seaLevel(:,1),seaLevel(:,2:4),age); % interpolate sea level from input
    salinity = sal + salinity./min(seaLevel(:,2)).*0.7; % salinity scaled from sea level
else
    salinity = sal;
end
atmosCO2 = interp1(pCO2(:,1),pCO2(:,2:2),age) + DpCO2;    % CO2 interpolated from input

d11Bsw = 39.61;

if species==1   % calibration choice - G.ruber (w)
    a = 0.035; b = 0.064; c = -0.87; d = -0.03;
    a1 = 0.013; b1 = 0.0; c1 = 0.12; d1 = 0;
    intCov = @(n) 1.563359 - 25.016445.*n;   % covariance of intercept & T-sens. for unc. prop. (see SI)
elseif species==2   % T.sacculifer
    a = 0.052; b = 0.065; c = 0.03; d = -0.28;
    a1 = 0.012; b1 = 0.0; c1 = 0.13; d1 = 0;
    intCov = @(n) 1.390751 - 25.769725.*n;
elseif species==3   % G. bulloides
    a = 0.033; b = 0.064; c = -0.87; d = 0.14;
    a1 = 0.007; b1 = 0.0; c1 = 0.13; d1 = 0;
    intCov = @(n) 1.543317 - 21.816096.*n;
elseif species==4   % O.universa
    a = 0.040; b = 0.075; c = -0.5; d = 0.48;
    a1 = 0.014; b1 = 0.0; c1 = 0.12; d1 = 0;
    intCov = @(n) 2.091911 - 21.679394.*n;
else                % generic
    a = 0.034; b = 0.062; c = -0.73; d = 0;
    a1 = 0.007; b1 = 0.0; c1 = 0.07; d1 = 0;
    intCov = @(n) 1.647334 - 23.57042.*n;
end

if relative(1)~=0 %&& numel(relative)>1
    if max(size(relative))==1
        relative = [0 relative];
    end
    HolAge = age(age>relative(1) & age<relative(2));
    if size(HolAge,1)==0
        error('No Holocene data! Cannot compute relative T.')
    end
    [~,idx]=ismember(HolAge,age);
end

if unc==0
    pH = NaN(size(age,1),1);
    MgTemp = NaN(size(age,1),1);
    T = repmat(25,100,1);
    if d11B==0      % if no d11B record, use atmospheric CO2
        for i = 1:size(age,1)
            t = 10;
            z = 1;
            while t>0.001   %iteratively solve for temp+pH
                z = z+1;
                carbSys = CO2SYS(TAlk,atmosCO2(i,1),1,4,salinity(i,1),T(z-1),T(z-1),0,0,15,1,1,4,1);
                pH(i) = carbSys(:,3);       % pH not output
                T(z) = MgCa(i)./exp(a.*(salinity(i,1)-35) + (pH(i)-8).*c + d);
                T(z) = log(T(z))./b;
                t = abs(T(z)-T(z-1));
            end
            MgTemp(i) = T(z);
        end
    else        % use d11B
        for i = 1:size(age,1)
            t = 10;
            z = 1;
            while t>0.001
                z = z+1;
                KB = (-8966.9 -2890.53*salinity(i,1)^0.5 - 77.942*salinity(i,1) +...
                    1.728*salinity(i,1)^1.5 - 0.0996*salinity(i,1)^2)/(T(z-1)+273.15) + ...
                    (148.0248 + 137.1942*salinity(i,1)^0.5 + 1.62142*salinity(i,1)) + ...
                    (-24.4344 - 25.085*salinity(i,1)^0.5 - 0.2474*salinity(i,1))*log((T(z-1)+273.15)) + ...
                    0.053105*salinity(i,1)^0.5*(T(z-1)+273.15);
                pKB = -log10(exp(KB));
                pH(i) = pKB - log10(-((d11Bsw-d11B(i))/(d11Bsw-(1.0272*d11B(i)) - 1000*(1.0272 - 1))));
                T(z) = MgCa(i)./exp(a.*(salinity(i,1)-35) + (pH(i)-8).*c + d);
                T(z) = log(T(z))./b;
                t = abs(T(z)-T(z-1));
            end
            MgTemp(i) = T(z);
        end
    end
    if relative(1)~=0 && numel(relative)>1
        HolTemp = nanmean(MgTemp(idx(1):idx(size(idx,1)),j));
        MgTempRel(:,j) = MgTemp(:,j) - HolTemp;
    end

    
else    % perform uncertianty analysis
    mCN = 1000;
    pH = NaN(size(age,1),mCN);
    MgTemp = NaN(size(age,1),mCN);
    MgTempRel = NaN(size(age,1),mCN); 
    T = repmat(25,100,mCN);
    errInt = [normrnd(a,a1/2,mCN,1) normrnd(b,b1/2,mCN,1) normrnd(c,c1/2,mCN,1) ...
        normrnd(d,d1/2,mCN,1)];
    errInt(:,4) = intCov(errInt(:,2));  % replace intercept unc. with covariance
    if d11B==0  % use pCO2
        tic;
        f = waitbar(0,'running...','Name','Uncertainty simulation');
        errEx = NaN(4,size(age,1),mCN);
        for i = 1:size(age,1)   % normal random external errors
            if size(MgCa,2)==1  % Mg/Ca uncertainty not given,use 2% of value
                errEx(1,i,:) = normrnd(MgCa(i,1),MgCa(i,1).*0.01,1,mCN);
            else        % use value
                errEx(1,i,:) = MgCa(i,1)-MgCa(i,2) + 2.*MgCa(i,2).*rand(1,mCN);
            end
                errEx(2,i,:) = normrnd(0,20,1,mCN); % pCO2 unc.
                errEx(3,i,:) = normrnd(salinity(i,1),1,1,mCN); % sal. unc.
                errEx(4,i,:) = normrnd(TAlk,50,1,mCN); % TAlk unc. 
        end
        for j = 1:mCN
            for i = 1:size(age,1)
                t = 10;
                z = 1;
                while t>0.001   %iteratively solve for temp+pH
                    z = z+1;
                    carbSys = CO2SYS(errEx(4,i,j),atmosCO2(i,1)+errEx(2,i,j),1,4,...
                        errEx(3,i,j),T(z-1,j),T(z-1,j),0,0,15,1,1,4,1);
                    pH(i,j) = carbSys(:,3);       % pH not output
                    T(z,j) = (errEx(1,i,j))./exp(errInt(j,1).*(errEx(3,i,j)-35) + ...
                        (pH(i,j)-8).*errInt(j,3) + d);% (d+err(1,5)*d1));
                    T(z,j) = log(T(z,j))./errInt(j,2);
                    t = abs(T(z,j)-T(z-1,j));
                end
                MgTemp(i,j) = T(z,j);
            end
            waitbar(j/mCN,f)
        end
        close(f)
        toc
    else
        errEx = NaN(5,size(age,1),mCN);
        for i = 1:size(age,1)   % normal random external errors
            if size(MgCa,2)==1  % Mg/Ca uncertainty not given,use 2% of value
                errEx(1,i,:) = normrnd(MgCa(i,1),MgCa(i,1).*0.01,1,mCN);
            else        % use value
                errEx(1,i,:) = MgCa(i,1)-MgCa(i,2) + 2.*MgCa(i,2).*rand(1,mCN);
            end
                errEx(2,i,:) = normrnd(0,20,1,mCN); % pCO2 unc.
                errEx(3,i,:) = normrnd(salinity(i,1),1,1,mCN); % sal. unc.
                errEx(4,i,:) = normrnd(TAlk,50,1,mCN); % TAlk unc.  
            if size(d11B,2)==1  % d11B uncertainty not given, no. unc. prop
                errEx(5,i,:) = normrnd(d11B(i,1),0,1,mCN);
            else        % use value
                errEx(5,i,:) = d11B(i,1)-d11B(i,2) + 2.*d11B(i,2).*rand(1,mCN);
            end
        end
        for j = 1:mCN
            for i = 1:size(age,1)
            t = 10;
            z = 1;
            while t>0.001
                z = z+1;
                KB = (-8966.9 -2890.53*errEx(3,i,j)^0.5 - 77.942*errEx(3,i,j) +...
                    1.728*errEx(3,i,j)^1.5 - 0.0996*errEx(3,i,j)^2)/(T(z-1,j)+273.15) + ...
                    (148.0248 + 137.1942*errEx(3,i,j)^0.5 + 1.62142*errEx(3,i,j)) + ...
                    (-24.4344 - 25.085*errEx(3,i,j)^0.5 - 0.2474*errEx(3,i,j))*log((T(z-1,j)+273.15)) + ...
                    0.053105*errEx(3,i,j)^0.5*(T(z-1,j)+273.15);
                pKB = -log10(exp(KB));
                pH(i,j) = pKB - log10(-((d11Bsw-errEx(5,i,j))/(d11Bsw-(1.0272*errEx(5,i,j)) - 1000*(1.0272 - 1))));
                T(z,j) = (errEx(1,i,j))./exp(errInt(j,1).*(errEx(3,i,j)-35) + ...
                        (pH(i,j)-8).*errInt(j,3) + d);% (d+err(1,5)*d1));
                T(z,j) = log(T(z,j))./errInt(j,2);
                t = abs(T(z,j)-T(z-1,j));
            end
            MgTemp(i,j) = T(z,j);
            end
        end
    end
    if relative(1)~=0 && numel(relative)>1 || numel(relative)>1
        for j = 1:mCN
            HolTemp = nanmean(MgTemp(idx(1):idx(size(idx,1)),j));
            MgTempRel(:,j) = MgTemp(:,j) - HolTemp;
        end
    end
end


if fig==1  % draw figure?
    close(figure(1))   % 1SD 2SD envelopes
    H = figure(1);
    if relative(1)~=0 || numel(relative)>1
        fill([age ; flipud(age)],...
            [prctile(MgTempRel,98,2) ; flipud(prctile(MgTempRel,2,2))],...
            'b','facealpha',0.1,'edgecolor','w')
        hold on
        fill([age ; flipud(age)],...
            [prctile(MgTempRel,84,2) ; flipud(prctile(MgTempRel,16,2))],...
            'b','facealpha',0.1,'edgecolor','w')
        plot(age,prctile(MgTempRel,50,2),'-b','linewidth',0.5)
        scatter(age,prctile(MgTempRel,50,2),20,'white','filled',...
            'markeredgecolor','b','linewidth',0.5)
        set(gcf,'color','white')
        xlabel('age')
        ylabel('\DeltaT (\circC)')
    else
        fill([age ; flipud(age)],...
            [prctile(MgTemp,98,2) ; flipud(prctile(MgTemp,2,2))],...
            'b','facealpha',0.1,'edgecolor','w')
        hold on
        fill([age ; flipud(age)],...
            [prctile(MgTemp,84,2) ; flipud(prctile(MgTemp,16,2))],...
            'b','facealpha',0.1,'edgecolor','w')
        plot(age,prctile(MgTemp,50,2),'-b','linewidth',0.5)
        scatter(age,prctile(MgTemp,50,2),20,'white','filled',...
            'markeredgecolor','b','linewidth',0.5)
        set(gcf,'color','white')
        xlabel('age')
        ylabel('T (\circC)')
    end
end

MgOut = [prctile(MgTempRel,98,2) prctile(MgTempRel,84,2) prctile(MgTempRel,50,2) prctile(MgTempRel,16,2) prctile(MgTempRel,2,2)];

MgOut = [prctile(MgTemp,98,2) prctile(MgTemp,84,2) prctile(MgTemp,50,2) prctile(MgTemp,16,2) prctile(MgTemp,2,2)];

pHOut = [prctile(pH,98,2) prctile(pH,84,2) prctile(pH,50,2) prctile(pH,16,2) prctile(pH,2,2)];

close(figure(2))   % 1SD 2SD envelopes
    H = figure(2);
        fill([age ; flipud(age)],...
            [prctile(pH,98,2) ; flipud(prctile(pH,2,2))],...
            'b','facealpha',0.1,'edgecolor','w')
        hold on
        fill([age ; flipud(age)],...
            [prctile(pH,84,2) ; flipud(prctile(pH,16,2))],...
            'b','facealpha',0.1,'edgecolor','w')
        plot(age,prctile(pH,50,2),'-b','linewidth',0.5)
        scatter(age,prctile(pH,50,2),20,'white','filled',...
            'markeredgecolor','b','linewidth',0.5)
        set(gcf,'color','white')
        xlabel('age')
        ylabel('pH (total)')