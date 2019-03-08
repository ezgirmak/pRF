%function getRSP(directory, subjID, datatype)
clear all close all
directory = '/project_space/VISCOG_FOLDERS';
subjID = 197348;
dd=3;

atlas  = load(fullfile(directory, 'atlas.mat'));

%data directory
dataDr = fullfile(directory,num2str(subjID),'pRF/');
%result directory
resDir = fullfile(directory, 'RSM');
i_roi = [find([atlas.wang2015 == 2]);  find([atlas.wang2015 == 3])]; %select only V1
disp(['roi size = ', num2str(length(i_roi))]);

cd(dataDr);

%% load data
%find the file for original
file = dir('*_2mm_*_angle_eccentricity.mat');
file = file.name;
load(file);
orig = collated;
i_orig = find([orig.pRF(:).didFit]);

%find the file for filtered
file = dir('*_2mm_filtered_*_angle_eccentricity.mat')
file = file.name;
load(file);
filt = collated;
i_filt = find([filt.pRF(:).didFit]);

for dd=1:3
    if dd==1
        datatype = 'orig';
        data = orig;
        i_fit = i_orig;
    elseif dd==2
        datatype = 'filt';
        data = filt;
        i_fit = i_filt;
    else
        datatype = 'both';
        data = orig;
        i_fit = intersect(find([data.pRF(:).didFit]) , find([filt.pRF(:).didFit]));
    end
    
    
    %% find values in the ROI
    
    disp(['fitted data = ', num2str(length(i_fit))]);
    i_data = [];
    for ii=1:length(i_fit)
        if ismember(data.pRF(i_fit(ii)).id, i_roi)
            i_data = cat(1, i_data, i_fit(ii));
        end
    end
    
    disp(['data 2 analyse = ', num2str(length(i_data))]);
    
    %% calculate mean time courses in polar bins
    
    alist = linspace(-pi,pi, 20);
    rlist = linspace(0,12, 16);
    %rlist = [linspace(log10(.001), log10(12), 25)].^10;
    
    for r=1:length(rlist)-1
        disp(['finding ', num2str(r)])
        for a=1:length(alist)-1
            seed_ind = find([data.pRF.angle]>alist(a) & [data.pRF.angle]<=alist(a)+1 ...
                & [data.pRF.radius]>rlist(r) & [data.pRF.radius]<=rlist(r)+1 & [data.pRF(:).didFit] ==1);
            for s=1:length(data.scan)
                if datatype == 'both'
                    ave_tc( r, a, s).tc =  mean([data.scan(s).vtc(:, seed_ind)-filt.scan(s).vtc(:, seed_ind)], 2);
                else
                    ave_tc( r, a, s).tc =  mean(data.scan(s).vtc(:, seed_ind), 2);
                end
                ave_tc( r, a, s).n =  length(seed_ind);
            end
        end
    end
    %% calculate cross correlations across every bin
    
    
    
    
    %%
    
    nr = length(rlist)-1;
    na = length(alist)-1;
    ns = size(ave_tc,3);
    
    midr = floor(nr/2);  % 'center' radius
    mida = floor(na/2);
    %
    c = zeros(nr,na);
    n = zeros(nr,na);
    for r=1:nr
        disp(['finding ', num2str(r)])
        for a=1:na
            
            for s=1:ns
                tmpcor = corr([ave_tc(r,a,s).tc],[ave_tc(:,:,s).tc]);
                tmpcor = reshape(tmpcor,nr,na);
                da=(1:na)+a-mida; % difference in angles
                da = mod(da-1,na)+1;
                
                dr = (1:nr)+r-midr; %difference in radius
                gooddr = dr>=1 & dr<=nr;
                shiftcor = NaN*ones(nr,na);
                shiftcor(gooddr,:) = tmpcor(dr(gooddr),da);
                
                
                id = ~isnan(shiftcor);
                
                c(id) = c(id)+shiftcor(id);
                n(id) = n(id)+1;
            end
        end
    end
    pol = c./n;
    
    
    
    
    %% plot it
    
    f1 = figure(dd);
    
    pol(1) =0;
    imagesc(rlist(1:end-1), alist(1:end-1), pol');
    colormap(hot)
    colorbar
    ylabel('angle');
    xlabel('radius')
    title([num2str(subjID) ' ' datatype]);
    saveName = [resDir '/' num2str(subjID) '_' datatype  '_pol.mat'];
    save(saveName, 'pol', '-v7.3');
    %%
    naa = 10;
    aa = linspace(0,(alist(2)-alist(1)),naa);  %arc
    wons = ones(1,naa);
    
    colList = hot(256);
    
    mincor = -.3;
    maxcor = 1;
    figure(dd)
    clf
    subplot('Position',[0.05,.2,.9,.8])
    hold on
    
    for ai=1:na
        for ri=1:nr
            if ~isnan(pol(ri,ai))
                colid = floor(255*(pol(ri,ai)-mincor)/(maxcor-mincor));
                col = colList(colid,:);
            else
                col = [.5,.5,.5];
            end
            r=[rlist(ri),rlist(ri+1)*wons,rlist(ri)*wons];
            a=[alist(ai),alist(ai)+aa,alist(ai)+fliplr(aa)];
            patch(r.*cos(a),r.*sin(a),col,'EdgeColor','none')
        end
    end
    axis equal
    set(gca,'Color',[.5,.5,.5]);
    set(gca,'XTick',[]);
    set(gca,'YTick',[]);
    
    %     subplot('Position',[0.05,.05,.9,.1]);
    %
    %     imagesc(linspace(mincor,maxcor),ones(1,256),linspace(mincor,maxcor))
    %     colormap(colList)
    %     set(gca,'YTick',[]);
    %     axis equal
    %     axis tight
    %     set(gca,'YLim',[.98,1.02]);
    saveas(f1, [saveName '_RSM.jpg'], 'jpg');
    
end

