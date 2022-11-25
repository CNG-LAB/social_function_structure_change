 
% Code for https://www.biorxiv.org/content/10.1101/2020.11.11.377895v4 
% Questions valk@cbs.mpg.de



for load_data_and_dependencies  = 1
    %BrainSpace: brainspace.readthedocs.io
    %Colorbrewer colormaps: https://de.mathworks.com/matlabcentral/fileexchange/58350-cbrewer2
    %Surfstat / brainstat: brainstat.readthedocs.io
    %RaincloudPlot: https://github.com/RainCloudPlots/RainCloudPlots
    %fdr correction code
end

for load_surface = 1
    % load:  mesh
    % -----------------
    cd(path to '/fsaverage5/surf/')
    SPHERE      = SurfStatReadSurf({'lh.sphere','rh.sphere'});
    S           = SurfStatReadSurf({'lh.inflated','rh.inflated'});
    SW          = SurfStatReadSurf({'lh.white','rh.white'});
    SP          = SurfStatReadSurf({'lh.pial','rh.pial'});
    SM.tri      = SW.tri;
    SM.coord    = (SW.coord + SP.coord)./2;
    SInf        = SW;
    SInf.coord  = 0.2 *SW.coord + 0.8* S.coord;
    
end

for load_aux_data  = 1
%     Schaefer 400 parcellation
%     Behavioral data (attention, compassion, ToM)
%     Sex, IDs, Timepoint, quality control
end

for load_brain_data_func = 1
   
    for cortex = 1
%         Ingredients:
%         - HCP mean functional connectome
%         - functional connectomes of study
%         - cortical thickness data of the sample, as published in Valk,
%         2017

        rowskeep      = (1:400);
        isthere_fc    = find(cortex400(:,1,2)~=0);
        alignedrs     = zeros(992,400,10);
        noalignedrs   = zeros(992,400,10);
        for i = 1:992
            tic
            try
                other_rs = squeeze(cortex400(i,:,:));
                other_rs(eye(size(other_rs))==1) = 0;
                if(mean(mean(other_rs)) ~= 0 | ~isnan(mean(mean(other_rs))))
                    Gp = GradientMaps('kernel','na','approach','dm','alignment','pa');
                    Gp = Gp.fit({fc400m, other_rs(1:400,1:400)});
                    alignedrs(i,:,:) = Gp.aligned{:,2};
                    noalignedrs(i,:,:) = Gp.gradients{:,2};
                else
                end
            catch
            end
            toc
        end
        
        %flip directions if needed
        [r] = corr(squeeze(alignedrs(1,:,1))',squeeze(alignedrs(:,:,1))');
        G1400 = alignedrs(:,:,1).*sign(r)';
        [r] = corr(squeeze(-alignedrs(1,:,2))',squeeze(alignedrs(:,:,2))');
        G2400 = alignedrs(:,:,2).*sign(r)';
        [r] = corr(squeeze(-alignedrs(1,:,3))',squeeze(alignedrs(:,:,3))');
        G3400 = alignedrs(:,:,3).*sign(r)';
        
        for i = 1:992
            GGG400(i,:) = sqrt((G1400(i,:).^2)+(G2400(i,:).^2)+(G3400(i,:).^2));
        end
                        
    end
        
    
    for ctx_control = 1
        GGG400corr = zeros(992,400);
        for j = 1
            measure_used = GGG400;
            keeproi3 =  zeros(992,1);
            keeproi2 = ~isnan(mean(measure_used,2)).* (nanmean(measure_used,2)~=0);
            keeproi3 = ~isnan(mean(ct400,2)).* (nanmean(ct400,2)~=0);
            keeptmp = keeproi2.*keeproi3;
            for i = 1:400
                seed = ct400(find(keeptmp),i);
                M    = 1 + term(seed);
                slm  = SurfStatLinMod(measure_used(find(keeptmp),i),M);
                GGG400corr(find(keeptmp),i) = measure_used(find(keeptmp),i) - slm.X*slm.coef;
            end
        end
    end
  
    
    for gsr_cortex = 1
        rsgsr = load('/Users/sofievalk/Documents/Data/Resource_stuff/yeo400gsr.mat')
        
        alignedrs = zeros(992,400,10);
        noalignedrs = zeros(992,400,10);
        for i = 1:992
            try
                other_rs = squeeze(rsgsr.yeo400conn(i,:,:));
                other_rs(eye(size(other_rs))==1) = 0;
                if(mean(mean(other_rs)) ~= 0 | ~isnan(mean(mean(other_rs))))
                    Gp = GradientMaps('kernel','na','approach','dm','alignment','pa');
                    Gp = Gp.fit({fc400m, other_rs});
                    
                    alignedrs(i,:,:) = Gp.aligned{:,2};
                    noalignedrs(i,:,:) = Gp.gradients{:,2};
                else
                end
            catch
            end
        end
        
        
        [r] = corr(squeeze(alignedrs(1,:,1))',squeeze(alignedrs(:,:,1))');
        G1allgsr = alignedrs(:,:,1).*sign(r)';
        [r] = corr(squeeze(-alignedrs(1,:,2))',squeeze(alignedrs(:,:,2))');
        G2allgsr = alignedrs(:,:,2).*sign(r)';
        [r] = corr(squeeze(-alignedrs(1,:,3))',squeeze(alignedrs(:,:,3))');
        G3allgsr = alignedrs(:,:,3).*sign(r)';
        
        for i = 1:992
            GGGgsr(i,:) = sqrt((G1allgsr(i,:).^2)+(G2allgsr(i,:).^2)+(G3allgsr(i,:).^2));
        end
    end
end

for load_brain_data_struct = 1
    % Ingredients
    % 1. microstructural profiles, can be computed using
    % micapipe.readthedocs.io or
    % https://github.com/CNG-LAB/CNG_cingulate_gradients for example.
    
    % More ingredients: 
    % Freesurfer data, see also Valk, 2017
     
    
    
    for CTX_control = 1
        MPC_layers_ctx = zeros(992,12,400);
        for j = 1:12
            measure_used = squeeze(MPC_layers(:,j,:));
            keeproi3 =  zeros(992,1);
            keeproi2 = ~isnan(mean(measure_used,2)).* (nanmean(measure_used,2)~=0);
            keeproi3 = ~isnan(mean(ct400,2)).* (nanmean(ct400,2)~=0);
            keeptmp = keeproi2.*keeproi3;
            for i = 1:400
                seed = ct400(find(keeptmp),i);
                M    = 1 + term(seed);
                slm  = SurfStatLinMod(measure_used(find(keeptmp),i),M);
                MPC_layers_ctx(find(keeptmp),j,i) = measure_used(find(keeptmp),i) - slm.X*slm.coef;
            end
        end
    end
    
    
end

% ggg
for  make_change_func = 1
    KZ       = zeros(992,400);
    keeproi2 = mintersect(isthere_fc,find(~isnan(mean(GGG400,2))),find(fmri_qc==1));
    keeproi  = zeros(992,1);
    keeproi(keeproi2) = 1; 
    
    KZ = GGG400;  
    minit = zeros(size(id));
    tp_t0_there = minit;
    tp_t1_there = minit;
    tp_t2_there = minit;
    tp_t3_there = minit;
    daysfrombl  = minit;
    daysfromt1  = minit;
    daysfromt2  = minit;
    DCChange_bl  = zeros(size(KZ))-666;
    DCChange_t1  = zeros(size(KZ))-666;
    DCChange_t2  = zeros(size(KZ))-666;
    
    %     who have data there?
    for i = 1:length(id)
        if keeproi(i,:) == 1
            i
            inter0 		  = intersect(find(tpnum==0),find(idnum==idnum(i)));
            inter1        = intersect(find(tpnum==1),find(idnum==idnum(i)));
            inter2        = intersect(find(tpnum==2),find(idnum==idnum(i)));
            inter3        = intersect(find(tpnum==3),find(idnum==idnum(i)));
            
            tp_t0_there(i,1)    = ~isempty(inter0) .* keeproi(i)==1;
            if ~isempty(inter0)
                tp_t0_there(i,1) =  tp_t0_there(i,1) .* keeproi(inter0);
            end
            tp_t1_there(i,1)    = ~isempty(inter1) .* keeproi(i)==1;
            if ~isempty(inter1)
                tp_t1_there(i,1) =  tp_t1_there(i,1) .* keeproi(inter1);
            end
            tp_t2_there(i,1)    = ~isempty(inter2) .* keeproi(i)==1;
            if ~isempty(inter2)
                tp_t2_there(i,1) =  tp_t2_there(i,1) .* keeproi(inter2);
            end
            tp_t3_there(i,1)    = ~isempty(inter3) .* keeproi(i)==1;
            if ~isempty(inter3)
                tp_t3_there(i,1) =  tp_t3_there(i,1) .* keeproi(inter3);
            end
            if keeproi(inter0)~=0
                daysfrombl(i,1)    = days(i)-days(inter0);
                DCChange_bl(i,:,:)    = (KZ(i,:,:) - KZ(inter0,:,:));
            end
            if keeproi(inter1)~=0
                daysfromt1(i,1)    = days(i)-days(inter1);
                DCChange_t1(i,:,:)    = (KZ(i,:,:) - KZ(inter1,:,:));
            end
            if keeproi(inter2)~=0
                daysfromt2(i,1)    = days(i)-days(inter2);
                DCChange_t2(i,:,:)    = (KZ(i,:,:) - KZ(inter2,:,:));
            end
        else
        end
    end
    
    ZChange_last = zeros(size(DCChange_bl));
    ZChange_last(strcmp(group4,'Group_1') & strcmp(tp,'T1'),:,:) = DCChange_bl(strcmp(group4,'Group_1') & strcmp(tp,'T1'),:,:);
    ZChange_last(strcmp(group4,'Group_1') & strcmp(tp,'T2'),:,:) = DCChange_t1(strcmp(group4,'Group_1') & strcmp(tp,'T2'),:,:);
    ZChange_last(strcmp(group4,'Group_1') & strcmp(tp,'T3'),:,:) = DCChange_t2(strcmp(group4,'Group_1') & strcmp(tp,'T3'),:,:);
    
    ZChange_last(strcmp(group4,'Group_2') & strcmp(tp,'T1'),:,:) = DCChange_bl(strcmp(group4,'Group_2') & strcmp(tp,'T1'),:,:);
    ZChange_last(strcmp(group4,'Group_2') & strcmp(tp,'T2'),:,:) = DCChange_t1(strcmp(group4,'Group_2') & strcmp(tp,'T2'),:,:);
    ZChange_last(strcmp(group4,'Group_2') & strcmp(tp,'T3'),:,:) = DCChange_t2(strcmp(group4,'Group_2') & strcmp(tp,'T3'),:,:);
    
    ZChange_last(strcmp(group4,'Group3') & strcmp(tp,'T1'),:,:)  = DCChange_bl(strcmp(group4,'Group3')  & strcmp(tp,'T1'),:,:);
    ZChange_last(strcmp(group4,'Group3') & strcmp(tp,'T2'),:,:)  = DCChange_t1(strcmp(group4,'Group3') & strcmp(tp,'T2'),:,:);
    ZChange_last(strcmp(group4,'Group3') & strcmp(tp,'T3'),:,:)  = DCChange_t2(strcmp(group4,'Group3') & strcmp(tp,'T3'),:,:);
    
    ZChange_last(strcmp(group4,'Control') & strcmp(tp,'T1'),:,:) = DCChange_bl(strcmp(group4,'Control') & strcmp(tp,'T1'),:,:);
    ZChange_last(strcmp(group4,'Control') & strcmp(tp,'T2'),:,:) = DCChange_t1(strcmp(group4,'Control') & strcmp(tp,'T2'),:,:);
    ZChange_last(strcmp(group4,'Control') & strcmp(tp,'T3'),:,:) = DCChange_t2(strcmp(group4,'Control') & strcmp(tp,'T3'),:,:);
    
    ggg_400_bl = DCChange_bl;
    ggg_400_last = ZChange_last;
end

% ggg gsr corrected
for  make_change_func = 1
    KZ       = zeros(992,400);
    keeproi2 = mintersect(isthere_fc,find(~isnan(mean(GGG400,2))),find(fmri_qc==1));
    keeproi  = zeros(992,1);
    keeproi(keeproi2) = 1; 
    
    KZ = GGGgsr;  
    minit = zeros(size(id));
    tp_t0_there = minit;
    tp_t1_there = minit;
    tp_t2_there = minit;
    tp_t3_there = minit;
    daysfrombl  = minit;
    daysfromt1  = minit;
    daysfromt2  = minit;
    DCChange_bl  = zeros(size(KZ))-666;
    DCChange_t1  = zeros(size(KZ))-666;
    DCChange_t2  = zeros(size(KZ))-666;
    
    %     who have data there?
    for i = 1:length(id)
        if keeproi(i,:) == 1
            i
            inter0 		  = intersect(find(tpnum==0),find(idnum==idnum(i)));
            inter1        = intersect(find(tpnum==1),find(idnum==idnum(i)));
            inter2        = intersect(find(tpnum==2),find(idnum==idnum(i)));
            inter3        = intersect(find(tpnum==3),find(idnum==idnum(i)));
            
            tp_t0_there(i,1)    = ~isempty(inter0) .* keeproi(i)==1;
            if ~isempty(inter0)
                tp_t0_there(i,1) =  tp_t0_there(i,1) .* keeproi(inter0);
            end
            tp_t1_there(i,1)    = ~isempty(inter1) .* keeproi(i)==1;
            if ~isempty(inter1)
                tp_t1_there(i,1) =  tp_t1_there(i,1) .* keeproi(inter1);
            end
            tp_t2_there(i,1)    = ~isempty(inter2) .* keeproi(i)==1;
            if ~isempty(inter2)
                tp_t2_there(i,1) =  tp_t2_there(i,1) .* keeproi(inter2);
            end
            tp_t3_there(i,1)    = ~isempty(inter3) .* keeproi(i)==1;
            if ~isempty(inter3)
                tp_t3_there(i,1) =  tp_t3_there(i,1) .* keeproi(inter3);
            end
            if keeproi(inter0)~=0
                daysfrombl(i,1)    = days(i)-days(inter0);
                DCChange_bl(i,:,:)    = (KZ(i,:,:) - KZ(inter0,:,:));
            end
            if keeproi(inter1)~=0
                daysfromt1(i,1)    = days(i)-days(inter1);
                DCChange_t1(i,:,:)    = (KZ(i,:,:) - KZ(inter1,:,:));
            end
            if keeproi(inter2)~=0
                daysfromt2(i,1)    = days(i)-days(inter2);
                DCChange_t2(i,:,:)    = (KZ(i,:,:) - KZ(inter2,:,:));
            end
        else
        end
    end
    
    ZChange_last = zeros(size(DCChange_bl));
    ZChange_last(strcmp(group4,'Group_1') & strcmp(tp,'T1'),:,:) = DCChange_bl(strcmp(group4,'Group_1') & strcmp(tp,'T1'),:,:);
    ZChange_last(strcmp(group4,'Group_1') & strcmp(tp,'T2'),:,:) = DCChange_t1(strcmp(group4,'Group_1') & strcmp(tp,'T2'),:,:);
    ZChange_last(strcmp(group4,'Group_1') & strcmp(tp,'T3'),:,:) = DCChange_t2(strcmp(group4,'Group_1') & strcmp(tp,'T3'),:,:);
    
    ZChange_last(strcmp(group4,'Group_2') & strcmp(tp,'T1'),:,:) = DCChange_bl(strcmp(group4,'Group_2') & strcmp(tp,'T1'),:,:);
    ZChange_last(strcmp(group4,'Group_2') & strcmp(tp,'T2'),:,:) = DCChange_t1(strcmp(group4,'Group_2') & strcmp(tp,'T2'),:,:);
    ZChange_last(strcmp(group4,'Group_2') & strcmp(tp,'T3'),:,:) = DCChange_t2(strcmp(group4,'Group_2') & strcmp(tp,'T3'),:,:);
    
    ZChange_last(strcmp(group4,'Group3') & strcmp(tp,'T1'),:,:)  = DCChange_bl(strcmp(group4,'Group3')  & strcmp(tp,'T1'),:,:);
    ZChange_last(strcmp(group4,'Group3') & strcmp(tp,'T2'),:,:)  = DCChange_t1(strcmp(group4,'Group3') & strcmp(tp,'T2'),:,:);
    ZChange_last(strcmp(group4,'Group3') & strcmp(tp,'T3'),:,:)  = DCChange_t2(strcmp(group4,'Group3') & strcmp(tp,'T3'),:,:);
    
    ZChange_last(strcmp(group4,'Control') & strcmp(tp,'T1'),:,:) = DCChange_bl(strcmp(group4,'Control') & strcmp(tp,'T1'),:,:);
    ZChange_last(strcmp(group4,'Control') & strcmp(tp,'T2'),:,:) = DCChange_t1(strcmp(group4,'Control') & strcmp(tp,'T2'),:,:);
    ZChange_last(strcmp(group4,'Control') & strcmp(tp,'T3'),:,:) = DCChange_t2(strcmp(group4,'Control') & strcmp(tp,'T3'),:,:);
    
    ggg_400_gsr = ZChange_last;
end

% ggg cortical thickness corrected
for  make_change_func = 1
    KZ       = zeros(992,400);
    keeproi2 = mintersect(isthere_fc,find(~isnan(mean(GGG400corr,2))),find(fmri_qc==1));
    keeproi  = zeros(992,1);
    keeproi(keeproi2) = 1; 
    
    KZ = GGG400corr;  
    minit = zeros(size(id));
    tp_t0_there = minit;
    tp_t1_there = minit;
    tp_t2_there = minit;
    tp_t3_there = minit;
    daysfrombl  = minit;
    daysfromt1  = minit;
    daysfromt2  = minit;
    DCChange_bl  = zeros(size(KZ))-666;
    DCChange_t1  = zeros(size(KZ))-666;
    DCChange_t2  = zeros(size(KZ))-666;
    
    %     who have data there?
    for i = 1:length(id)
        if keeproi(i,:) == 1
            i
            inter0 		  = intersect(find(tpnum==0),find(idnum==idnum(i)));
            inter1        = intersect(find(tpnum==1),find(idnum==idnum(i)));
            inter2        = intersect(find(tpnum==2),find(idnum==idnum(i)));
            inter3        = intersect(find(tpnum==3),find(idnum==idnum(i)));
            
            tp_t0_there(i,1)    = ~isempty(inter0) .* keeproi(i)==1;
            if ~isempty(inter0)
                tp_t0_there(i,1) =  tp_t0_there(i,1) .* keeproi(inter0);
            end
            tp_t1_there(i,1)    = ~isempty(inter1) .* keeproi(i)==1;
            if ~isempty(inter1)
                tp_t1_there(i,1) =  tp_t1_there(i,1) .* keeproi(inter1);
            end
            tp_t2_there(i,1)    = ~isempty(inter2) .* keeproi(i)==1;
            if ~isempty(inter2)
                tp_t2_there(i,1) =  tp_t2_there(i,1) .* keeproi(inter2);
            end
            tp_t3_there(i,1)    = ~isempty(inter3) .* keeproi(i)==1;
            if ~isempty(inter3)
                tp_t3_there(i,1) =  tp_t3_there(i,1) .* keeproi(inter3);
            end
            if keeproi(inter0)~=0
                daysfrombl(i,1)    = days(i)-days(inter0);
                DCChange_bl(i,:,:)    = (KZ(i,:,:) - KZ(inter0,:,:));
            end
            if keeproi(inter1)~=0
                daysfromt1(i,1)    = days(i)-days(inter1);
                DCChange_t1(i,:,:)    = (KZ(i,:,:) - KZ(inter1,:,:));
            end
            if keeproi(inter2)~=0
                daysfromt2(i,1)    = days(i)-days(inter2);
                DCChange_t2(i,:,:)    = (KZ(i,:,:) - KZ(inter2,:,:));
            end
        else
        end
    end
    
    ZChange_last = zeros(size(DCChange_bl));
    ZChange_last(strcmp(group4,'Group_1') & strcmp(tp,'T1'),:,:) = DCChange_bl(strcmp(group4,'Group_1') & strcmp(tp,'T1'),:,:);
    ZChange_last(strcmp(group4,'Group_1') & strcmp(tp,'T2'),:,:) = DCChange_t1(strcmp(group4,'Group_1') & strcmp(tp,'T2'),:,:);
    ZChange_last(strcmp(group4,'Group_1') & strcmp(tp,'T3'),:,:) = DCChange_t2(strcmp(group4,'Group_1') & strcmp(tp,'T3'),:,:);
    
    ZChange_last(strcmp(group4,'Group_2') & strcmp(tp,'T1'),:,:) = DCChange_bl(strcmp(group4,'Group_2') & strcmp(tp,'T1'),:,:);
    ZChange_last(strcmp(group4,'Group_2') & strcmp(tp,'T2'),:,:) = DCChange_t1(strcmp(group4,'Group_2') & strcmp(tp,'T2'),:,:);
    ZChange_last(strcmp(group4,'Group_2') & strcmp(tp,'T3'),:,:) = DCChange_t2(strcmp(group4,'Group_2') & strcmp(tp,'T3'),:,:);
    
    ZChange_last(strcmp(group4,'Group3') & strcmp(tp,'T1'),:,:)  = DCChange_bl(strcmp(group4,'Group3')  & strcmp(tp,'T1'),:,:);
    ZChange_last(strcmp(group4,'Group3') & strcmp(tp,'T2'),:,:)  = DCChange_t1(strcmp(group4,'Group3') & strcmp(tp,'T2'),:,:);
    ZChange_last(strcmp(group4,'Group3') & strcmp(tp,'T3'),:,:)  = DCChange_t2(strcmp(group4,'Group3') & strcmp(tp,'T3'),:,:);
    
    ZChange_last(strcmp(group4,'Control') & strcmp(tp,'T1'),:,:) = DCChange_bl(strcmp(group4,'Control') & strcmp(tp,'T1'),:,:);
    ZChange_last(strcmp(group4,'Control') & strcmp(tp,'T2'),:,:) = DCChange_t1(strcmp(group4,'Control') & strcmp(tp,'T2'),:,:);
    ZChange_last(strcmp(group4,'Control') & strcmp(tp,'T3'),:,:) = DCChange_t2(strcmp(group4,'Control') & strcmp(tp,'T3'),:,:);
    
    ggg_400_last_ctx = ZChange_last;
end

% g1
for  make_change_func_g1 = 1
    KZ       = zeros(992,400);
    keeproi2 = mintersect(isthere_fc,find(~isnan(mean(GGG400,2))),find(fmri_qc==1));
    keeproi  = zeros(992,1);
    keeproi(keeproi2) = 1; 
    KZ = G1400;  
    minit = zeros(size(id));
    tp_t0_there = minit;
    tp_t1_there = minit;
    tp_t2_there = minit;
    tp_t3_there = minit;
    daysfrombl  = minit;
    daysfromt1  = minit;
    daysfromt2  = minit;
    DCChange_bl  = zeros(size(KZ))-666;
    DCChange_t1  = zeros(size(KZ))-666;
    DCChange_t2  = zeros(size(KZ))-666;
    
    %     who have data there?
    for i = 1:length(id)
        if keeproi(i,:) == 1
            i
            inter0 		  = intersect(find(tpnum==0),find(idnum==idnum(i)));
            inter1        = intersect(find(tpnum==1),find(idnum==idnum(i)));
            inter2        = intersect(find(tpnum==2),find(idnum==idnum(i)));
            inter3        = intersect(find(tpnum==3),find(idnum==idnum(i)));
            
            tp_t0_there(i,1)    = ~isempty(inter0) .* keeproi(i)==1;
            if ~isempty(inter0)
                tp_t0_there(i,1) =  tp_t0_there(i,1) .* keeproi(inter0);
            end
            tp_t1_there(i,1)    = ~isempty(inter1) .* keeproi(i)==1;
            if ~isempty(inter1)
                tp_t1_there(i,1) =  tp_t1_there(i,1) .* keeproi(inter1);
            end
            tp_t2_there(i,1)    = ~isempty(inter2) .* keeproi(i)==1;
            if ~isempty(inter2)
                tp_t2_there(i,1) =  tp_t2_there(i,1) .* keeproi(inter2);
            end
            tp_t3_there(i,1)    = ~isempty(inter3) .* keeproi(i)==1;
            if ~isempty(inter3)
                tp_t3_there(i,1) =  tp_t3_there(i,1) .* keeproi(inter3);
            end
            if keeproi(inter0)~=0
                daysfrombl(i,1)    = days(i)-days(inter0);
                DCChange_bl(i,:,:)    = (KZ(i,:,:) - KZ(inter0,:,:));
            end
            if keeproi(inter1)~=0
                daysfromt1(i,1)    = days(i)-days(inter1);
                DCChange_t1(i,:,:)    = (KZ(i,:,:) - KZ(inter1,:,:));
            end
            if keeproi(inter2)~=0
                daysfromt2(i,1)    = days(i)-days(inter2);
                DCChange_t2(i,:,:)    = (KZ(i,:,:) - KZ(inter2,:,:));
            end
        else
        end
    end
    
    ZChange_last = zeros(size(DCChange_bl));
    ZChange_last(strcmp(group4,'Group_1') & strcmp(tp,'T1'),:,:) = DCChange_bl(strcmp(group4,'Group_1') & strcmp(tp,'T1'),:,:);
    ZChange_last(strcmp(group4,'Group_1') & strcmp(tp,'T2'),:,:) = DCChange_t1(strcmp(group4,'Group_1') & strcmp(tp,'T2'),:,:);
    ZChange_last(strcmp(group4,'Group_1') & strcmp(tp,'T3'),:,:) = DCChange_t2(strcmp(group4,'Group_1') & strcmp(tp,'T3'),:,:);
    
    ZChange_last(strcmp(group4,'Group_2') & strcmp(tp,'T1'),:,:) = DCChange_bl(strcmp(group4,'Group_2') & strcmp(tp,'T1'),:,:);
    ZChange_last(strcmp(group4,'Group_2') & strcmp(tp,'T2'),:,:) = DCChange_t1(strcmp(group4,'Group_2') & strcmp(tp,'T2'),:,:);
    ZChange_last(strcmp(group4,'Group_2') & strcmp(tp,'T3'),:,:) = DCChange_t2(strcmp(group4,'Group_2') & strcmp(tp,'T3'),:,:);
    
    ZChange_last(strcmp(group4,'Group3') & strcmp(tp,'T1'),:,:)  = DCChange_bl(strcmp(group4,'Group3')  & strcmp(tp,'T1'),:,:);
    ZChange_last(strcmp(group4,'Group3') & strcmp(tp,'T2'),:,:)  = DCChange_t1(strcmp(group4,'Group3') & strcmp(tp,'T2'),:,:);
    ZChange_last(strcmp(group4,'Group3') & strcmp(tp,'T3'),:,:)  = DCChange_t2(strcmp(group4,'Group3') & strcmp(tp,'T3'),:,:);
    
    ZChange_last(strcmp(group4,'Control') & strcmp(tp,'T1'),:,:) = DCChange_bl(strcmp(group4,'Control') & strcmp(tp,'T1'),:,:);
    ZChange_last(strcmp(group4,'Control') & strcmp(tp,'T2'),:,:) = DCChange_t1(strcmp(group4,'Control') & strcmp(tp,'T2'),:,:);
    ZChange_last(strcmp(group4,'Control') & strcmp(tp,'T3'),:,:) = DCChange_t2(strcmp(group4,'Control') & strcmp(tp,'T3'),:,:);
    
    g1_400_bl = DCChange_bl;
    g1_400_last = ZChange_last;
end

% g2
for  make_change_func_g2 = 1
    KZ       = zeros(992,400);
    keeproi2 = mintersect(isthere_fc,find(~isnan(mean(GGG400,2))),find(fmri_qc==1));
    keeproi  = zeros(992,1);
    keeproi(keeproi2) = 1; 
    KZ = G2400;  
    minit = zeros(size(id));
    tp_t0_there = minit;
    tp_t1_there = minit;
    tp_t2_there = minit;
    tp_t3_there = minit;
    daysfrombl  = minit;
    daysfromt1  = minit;
    daysfromt2  = minit;
    DCChange_bl  = zeros(size(KZ))-666;
    DCChange_t1  = zeros(size(KZ))-666;
    DCChange_t2  = zeros(size(KZ))-666;
    
    %     who have data there?
    for i = 1:length(id)
        if keeproi(i,:) == 1
            i
            inter0 		  = intersect(find(tpnum==0),find(idnum==idnum(i)));
            inter1        = intersect(find(tpnum==1),find(idnum==idnum(i)));
            inter2        = intersect(find(tpnum==2),find(idnum==idnum(i)));
            inter3        = intersect(find(tpnum==3),find(idnum==idnum(i)));
            
            tp_t0_there(i,1)    = ~isempty(inter0) .* keeproi(i)==1;
            if ~isempty(inter0)
                tp_t0_there(i,1) =  tp_t0_there(i,1) .* keeproi(inter0);
            end
            tp_t1_there(i,1)    = ~isempty(inter1) .* keeproi(i)==1;
            if ~isempty(inter1)
                tp_t1_there(i,1) =  tp_t1_there(i,1) .* keeproi(inter1);
            end
            tp_t2_there(i,1)    = ~isempty(inter2) .* keeproi(i)==1;
            if ~isempty(inter2)
                tp_t2_there(i,1) =  tp_t2_there(i,1) .* keeproi(inter2);
            end
            tp_t3_there(i,1)    = ~isempty(inter3) .* keeproi(i)==1;
            if ~isempty(inter3)
                tp_t3_there(i,1) =  tp_t3_there(i,1) .* keeproi(inter3);
            end
            if keeproi(inter0)~=0
                daysfrombl(i,1)    = days(i)-days(inter0);
                DCChange_bl(i,:,:)    = (KZ(i,:,:) - KZ(inter0,:,:));
            end
            if keeproi(inter1)~=0
                daysfromt1(i,1)    = days(i)-days(inter1);
                DCChange_t1(i,:,:)    = (KZ(i,:,:) - KZ(inter1,:,:));
            end
            if keeproi(inter2)~=0
                daysfromt2(i,1)    = days(i)-days(inter2);
                DCChange_t2(i,:,:)    = (KZ(i,:,:) - KZ(inter2,:,:));
            end
        else
        end
    end
    
    ZChange_last = zeros(size(DCChange_bl));
    ZChange_last(strcmp(group4,'Group_1') & strcmp(tp,'T1'),:,:) = DCChange_bl(strcmp(group4,'Group_1') & strcmp(tp,'T1'),:,:);
    ZChange_last(strcmp(group4,'Group_1') & strcmp(tp,'T2'),:,:) = DCChange_t1(strcmp(group4,'Group_1') & strcmp(tp,'T2'),:,:);
    ZChange_last(strcmp(group4,'Group_1') & strcmp(tp,'T3'),:,:) = DCChange_t2(strcmp(group4,'Group_1') & strcmp(tp,'T3'),:,:);
    
    ZChange_last(strcmp(group4,'Group_2') & strcmp(tp,'T1'),:,:) = DCChange_bl(strcmp(group4,'Group_2') & strcmp(tp,'T1'),:,:);
    ZChange_last(strcmp(group4,'Group_2') & strcmp(tp,'T2'),:,:) = DCChange_t1(strcmp(group4,'Group_2') & strcmp(tp,'T2'),:,:);
    ZChange_last(strcmp(group4,'Group_2') & strcmp(tp,'T3'),:,:) = DCChange_t2(strcmp(group4,'Group_2') & strcmp(tp,'T3'),:,:);
    
    ZChange_last(strcmp(group4,'Group3') & strcmp(tp,'T1'),:,:)  = DCChange_bl(strcmp(group4,'Group3')  & strcmp(tp,'T1'),:,:);
    ZChange_last(strcmp(group4,'Group3') & strcmp(tp,'T2'),:,:)  = DCChange_t1(strcmp(group4,'Group3') & strcmp(tp,'T2'),:,:);
    ZChange_last(strcmp(group4,'Group3') & strcmp(tp,'T3'),:,:)  = DCChange_t2(strcmp(group4,'Group3') & strcmp(tp,'T3'),:,:);
    
    ZChange_last(strcmp(group4,'Control') & strcmp(tp,'T1'),:,:) = DCChange_bl(strcmp(group4,'Control') & strcmp(tp,'T1'),:,:);
    ZChange_last(strcmp(group4,'Control') & strcmp(tp,'T2'),:,:) = DCChange_t1(strcmp(group4,'Control') & strcmp(tp,'T2'),:,:);
    ZChange_last(strcmp(group4,'Control') & strcmp(tp,'T3'),:,:) = DCChange_t2(strcmp(group4,'Control') & strcmp(tp,'T3'),:,:);
    
    g2_400_bl = DCChange_bl;
    g2_400_last = ZChange_last;
end

% g3
for  make_change_func_g3 = 1
    KZ       = zeros(992,400);
    keeproi2 = mintersect(isthere_fc,find(~isnan(mean(GGG400,2))),find(fmri_qc==1));
    keeproi  = zeros(992,1);
    keeproi(keeproi2) = 1; 
    KZ = G3400;  
    minit = zeros(size(id));
    tp_t0_there = minit;
    tp_t1_there = minit;
    tp_t2_there = minit;
    tp_t3_there = minit;
    daysfrombl  = minit;
    daysfromt1  = minit;
    daysfromt2  = minit;
    DCChange_bl  = zeros(size(KZ))-666;
    DCChange_t1  = zeros(size(KZ))-666;
    DCChange_t2  = zeros(size(KZ))-666;
    
    %     who have data there?
    for i = 1:length(id)
        if keeproi(i,:) == 1
            i
            inter0 		  = intersect(find(tpnum==0),find(idnum==idnum(i)));
            inter1        = intersect(find(tpnum==1),find(idnum==idnum(i)));
            inter2        = intersect(find(tpnum==2),find(idnum==idnum(i)));
            inter3        = intersect(find(tpnum==3),find(idnum==idnum(i)));
            
            tp_t0_there(i,1)    = ~isempty(inter0) .* keeproi(i)==1;
            if ~isempty(inter0)
                tp_t0_there(i,1) =  tp_t0_there(i,1) .* keeproi(inter0);
            end
            tp_t1_there(i,1)    = ~isempty(inter1) .* keeproi(i)==1;
            if ~isempty(inter1)
                tp_t1_there(i,1) =  tp_t1_there(i,1) .* keeproi(inter1);
            end
            tp_t2_there(i,1)    = ~isempty(inter2) .* keeproi(i)==1;
            if ~isempty(inter2)
                tp_t2_there(i,1) =  tp_t2_there(i,1) .* keeproi(inter2);
            end
            tp_t3_there(i,1)    = ~isempty(inter3) .* keeproi(i)==1;
            if ~isempty(inter3)
                tp_t3_there(i,1) =  tp_t3_there(i,1) .* keeproi(inter3);
            end
            if keeproi(inter0)~=0
                daysfrombl(i,1)    = days(i)-days(inter0);
                DCChange_bl(i,:,:)    = (KZ(i,:,:) - KZ(inter0,:,:));
            end
            if keeproi(inter1)~=0
                daysfromt1(i,1)    = days(i)-days(inter1);
                DCChange_t1(i,:,:)    = (KZ(i,:,:) - KZ(inter1,:,:));
            end
            if keeproi(inter2)~=0
                daysfromt2(i,1)    = days(i)-days(inter2);
                DCChange_t2(i,:,:)    = (KZ(i,:,:) - KZ(inter2,:,:));
            end
        else
        end
    end
    
    ZChange_last = zeros(size(DCChange_bl));
    ZChange_last(strcmp(group4,'Group_1') & strcmp(tp,'T1'),:,:) = DCChange_bl(strcmp(group4,'Group_1') & strcmp(tp,'T1'),:,:);
    ZChange_last(strcmp(group4,'Group_1') & strcmp(tp,'T2'),:,:) = DCChange_t1(strcmp(group4,'Group_1') & strcmp(tp,'T2'),:,:);
    ZChange_last(strcmp(group4,'Group_1') & strcmp(tp,'T3'),:,:) = DCChange_t2(strcmp(group4,'Group_1') & strcmp(tp,'T3'),:,:);
    
    ZChange_last(strcmp(group4,'Group_2') & strcmp(tp,'T1'),:,:) = DCChange_bl(strcmp(group4,'Group_2') & strcmp(tp,'T1'),:,:);
    ZChange_last(strcmp(group4,'Group_2') & strcmp(tp,'T2'),:,:) = DCChange_t1(strcmp(group4,'Group_2') & strcmp(tp,'T2'),:,:);
    ZChange_last(strcmp(group4,'Group_2') & strcmp(tp,'T3'),:,:) = DCChange_t2(strcmp(group4,'Group_2') & strcmp(tp,'T3'),:,:);
    
    ZChange_last(strcmp(group4,'Group3') & strcmp(tp,'T1'),:,:)  = DCChange_bl(strcmp(group4,'Group3')  & strcmp(tp,'T1'),:,:);
    ZChange_last(strcmp(group4,'Group3') & strcmp(tp,'T2'),:,:)  = DCChange_t1(strcmp(group4,'Group3') & strcmp(tp,'T2'),:,:);
    ZChange_last(strcmp(group4,'Group3') & strcmp(tp,'T3'),:,:)  = DCChange_t2(strcmp(group4,'Group3') & strcmp(tp,'T3'),:,:);
    
    ZChange_last(strcmp(group4,'Control') & strcmp(tp,'T1'),:,:) = DCChange_bl(strcmp(group4,'Control') & strcmp(tp,'T1'),:,:);
    ZChange_last(strcmp(group4,'Control') & strcmp(tp,'T2'),:,:) = DCChange_t1(strcmp(group4,'Control') & strcmp(tp,'T2'),:,:);
    ZChange_last(strcmp(group4,'Control') & strcmp(tp,'T3'),:,:) = DCChange_t2(strcmp(group4,'Control') & strcmp(tp,'T3'),:,:);
    
    g3_400_bl = DCChange_bl;
    g3_400_last = ZChange_last;
end
   
% layers
for  make_change_struct = 1
    layer_400_last = zeros(992,12,400);
    layer_400_bl = zeros(992,12,400);

    for ii = 1:12
        measure_used = squeeze(MPC_layers(:,ii,:));
        KZ = zeros(992,400);   
        keeproi3 =  zeros(992,1);
        keeproi2 = ~isnan(mean(measure_used,2)).* (nanmean(measure_used,2)~=0)
        keeproi3(mintersect(isthere_fc,find(~isnan(mean(GGG400,2))),find(fmri_qc==1))) = 1;       
        keeptmp = keeproi2.*keeproi3;
        
        r_tmp = corr(measure_used(keeptmp==1,:)',mean(measure_used(keeptmp==1,:))','type','spearman')
        keeproi = zeros(992,1);
        keeproi(keeptmp==1) = 1;
        KZ    = measure_used;      
        minit = zeros(size(id));
        tp_t0_there = minit;
        tp_t1_there = minit;
        tp_t2_there = minit;
        tp_t3_there = minit;
        daysfrombl  = minit;
        daysfromt1  = minit;
        daysfromt2  = minit;
        DCChange_bl  = zeros(size(KZ))-666;
        DCChange_t1  = zeros(size(KZ))-666;
        DCChange_t2  = zeros(size(KZ))-666;
        
        %     who have data there?
        for i = 1:length(id)
            i
            if keeproi(i) == 1
                i
                inter0 		  = intersect(find(tpnum==0),find(idnum==idnum(i)));
                inter1        = intersect(find(tpnum==1),find(idnum==idnum(i)));
                inter2        = intersect(find(tpnum==2),find(idnum==idnum(i)));
                inter3        = intersect(find(tpnum==3),find(idnum==idnum(i)));
                
                tp_t0_there(i,1)    = ~isempty(inter0) .* keeproi(i)==1;
                if ~isempty(inter0)
                    tp_t0_there(i,1) =  tp_t0_there(i,1) .* keeproi(inter0);
                end
                tp_t1_there(i,1)    = ~isempty(inter1) .* keeproi(i)==1;
                if ~isempty(inter1)
                    tp_t1_there(i,1) =  tp_t1_there(i,1) .* keeproi(inter1);
                end
                tp_t2_there(i,1)    = ~isempty(inter2) .* keeproi(i)==1;
                if ~isempty(inter2)
                    tp_t2_there(i,1) =  tp_t2_there(i,1) .* keeproi(inter2);
                end
                tp_t3_there(i,1)    = ~isempty(inter3) .* keeproi(i)==1;
                if ~isempty(inter3)
                    tp_t3_there(i,1) =  tp_t3_there(i,1) .* keeproi(inter3);
                end
                if keeproi(inter0)~=0
                    daysfrombl(i,1)    = days(i)-days(inter0);
                    DCChange_bl(i,:,:)    = (KZ(i,:,:) - KZ(inter0,:,:));
                end
                if keeproi(inter1)~=0
                    daysfromt1(i,1)    = days(i)-days(inter1);
                    DCChange_t1(i,:,:)    = (KZ(i,:,:) - KZ(inter1,:,:));
                end
                if keeproi(inter2)~=0
                    daysfromt2(i,1)    = days(i)-days(inter2);
                    DCChange_t2(i,:,:)    = (KZ(i,:,:) - KZ(inter2,:,:));
                end
            else
            end
        end
        
        ZChange_last = zeros(size(DCChange_bl));
        ZChange_last(strcmp(group4,'Group_1') & strcmp(tp,'T1'),:,:) = DCChange_bl(strcmp(group4,'Group_1') & strcmp(tp,'T1'),:,:);
        ZChange_last(strcmp(group4,'Group_1') & strcmp(tp,'T2'),:,:) = DCChange_t1(strcmp(group4,'Group_1') & strcmp(tp,'T2'),:,:);
        ZChange_last(strcmp(group4,'Group_1') & strcmp(tp,'T3'),:,:) = DCChange_t2(strcmp(group4,'Group_1') & strcmp(tp,'T3'),:,:);
        
        ZChange_last(strcmp(group4,'Group_2') & strcmp(tp,'T1'),:,:) = DCChange_bl(strcmp(group4,'Group_2') & strcmp(tp,'T1'),:,:);
        ZChange_last(strcmp(group4,'Group_2') & strcmp(tp,'T2'),:,:) = DCChange_t1(strcmp(group4,'Group_2') & strcmp(tp,'T2'),:,:);
        ZChange_last(strcmp(group4,'Group_2') & strcmp(tp,'T3'),:,:) = DCChange_t2(strcmp(group4,'Group_2') & strcmp(tp,'T3'),:,:);
        
        ZChange_last(strcmp(group4,'Group3') & strcmp(tp,'T1'),:,:)  = DCChange_bl(strcmp(group4,'Group3')  & strcmp(tp,'T1'),:,:);
        ZChange_last(strcmp(group4,'Group3') & strcmp(tp,'T2'),:,:)  = DCChange_t1(strcmp(group4,'Group3') & strcmp(tp,'T2'),:,:);
        ZChange_last(strcmp(group4,'Group3') & strcmp(tp,'T3'),:,:)  = DCChange_t2(strcmp(group4,'Group3') & strcmp(tp,'T3'),:,:);
        
        ZChange_last(strcmp(group4,'Control') & strcmp(tp,'T1'),:,:) = DCChange_bl(strcmp(group4,'Control') & strcmp(tp,'T1'),:,:);
        ZChange_last(strcmp(group4,'Control') & strcmp(tp,'T2'),:,:) = DCChange_t1(strcmp(group4,'Control') & strcmp(tp,'T2'),:,:);
        ZChange_last(strcmp(group4,'Control') & strcmp(tp,'T3'),:,:) = DCChange_t2(strcmp(group4,'Control') & strcmp(tp,'T3'),:,:);
        
        
        layer_400_last(:,ii,:) = ZChange_last;
        layer_400_bl(:,ii,:) = DCChange_bl;
    end
end

% layers cortical thickness corrected
for  make_change_struct = 1
    layer_400_last_ctx = zeros(992,12,400);

    for ii = 1:12
        measure_used = squeeze(MPC_layers_ctx(:,ii,:));
        KZ = zeros(992,400);   
        keeproi3 =  zeros(992,1);
        keeproi2 = ~isnan(mean(measure_used,2)).* (nanmean(measure_used,2)~=0)
        keeproi3(mintersect(isthere_fc,find(~isnan(mean(GGG400,2))),find(fmri_qc==1))) = 1;       
        keeptmp = keeproi2.*keeproi3;
        
        r_tmp = corr(measure_used(keeptmp==1,:)',mean(measure_used(keeptmp==1,:))','type','spearman')
        keeproi = zeros(992,1);
        keeproi(keeptmp==1) = 1;
        KZ    = measure_used;      
        minit = zeros(size(id));
        tp_t0_there = minit;
        tp_t1_there = minit;
        tp_t2_there = minit;
        tp_t3_there = minit;
        daysfrombl  = minit;
        daysfromt1  = minit;
        daysfromt2  = minit;
        DCChange_bl  = zeros(size(KZ))-666;
        DCChange_t1  = zeros(size(KZ))-666;
        DCChange_t2  = zeros(size(KZ))-666;
        
        %     who have data there?
        for i = 1:length(id)
            i
            if keeproi(i) == 1
                i
                inter0 		  = intersect(find(tpnum==0),find(idnum==idnum(i)));
                inter1        = intersect(find(tpnum==1),find(idnum==idnum(i)));
                inter2        = intersect(find(tpnum==2),find(idnum==idnum(i)));
                inter3        = intersect(find(tpnum==3),find(idnum==idnum(i)));
                
                tp_t0_there(i,1)    = ~isempty(inter0) .* keeproi(i)==1;
                if ~isempty(inter0)
                    tp_t0_there(i,1) =  tp_t0_there(i,1) .* keeproi(inter0);
                end
                tp_t1_there(i,1)    = ~isempty(inter1) .* keeproi(i)==1;
                if ~isempty(inter1)
                    tp_t1_there(i,1) =  tp_t1_there(i,1) .* keeproi(inter1);
                end
                tp_t2_there(i,1)    = ~isempty(inter2) .* keeproi(i)==1;
                if ~isempty(inter2)
                    tp_t2_there(i,1) =  tp_t2_there(i,1) .* keeproi(inter2);
                end
                tp_t3_there(i,1)    = ~isempty(inter3) .* keeproi(i)==1;
                if ~isempty(inter3)
                    tp_t3_there(i,1) =  tp_t3_there(i,1) .* keeproi(inter3);
                end
                if keeproi(inter0)~=0
                    daysfrombl(i,1)    = days(i)-days(inter0);
                    DCChange_bl(i,:,:)    = (KZ(i,:,:) - KZ(inter0,:,:));
                end
                if keeproi(inter1)~=0
                    daysfromt1(i,1)    = days(i)-days(inter1);
                    DCChange_t1(i,:,:)    = (KZ(i,:,:) - KZ(inter1,:,:));
                end
                if keeproi(inter2)~=0
                    daysfromt2(i,1)    = days(i)-days(inter2);
                    DCChange_t2(i,:,:)    = (KZ(i,:,:) - KZ(inter2,:,:));
                end
            else
            end
        end
        
        ZChange_last = zeros(size(DCChange_bl));
        ZChange_last(strcmp(group4,'Group_1') & strcmp(tp,'T1'),:,:) = DCChange_bl(strcmp(group4,'Group_1') & strcmp(tp,'T1'),:,:);
        ZChange_last(strcmp(group4,'Group_1') & strcmp(tp,'T2'),:,:) = DCChange_t1(strcmp(group4,'Group_1') & strcmp(tp,'T2'),:,:);
        ZChange_last(strcmp(group4,'Group_1') & strcmp(tp,'T3'),:,:) = DCChange_t2(strcmp(group4,'Group_1') & strcmp(tp,'T3'),:,:);
        
        ZChange_last(strcmp(group4,'Group_2') & strcmp(tp,'T1'),:,:) = DCChange_bl(strcmp(group4,'Group_2') & strcmp(tp,'T1'),:,:);
        ZChange_last(strcmp(group4,'Group_2') & strcmp(tp,'T2'),:,:) = DCChange_t1(strcmp(group4,'Group_2') & strcmp(tp,'T2'),:,:);
        ZChange_last(strcmp(group4,'Group_2') & strcmp(tp,'T3'),:,:) = DCChange_t2(strcmp(group4,'Group_2') & strcmp(tp,'T3'),:,:);
        
        ZChange_last(strcmp(group4,'Group3') & strcmp(tp,'T1'),:,:)  = DCChange_bl(strcmp(group4,'Group3')  & strcmp(tp,'T1'),:,:);
        ZChange_last(strcmp(group4,'Group3') & strcmp(tp,'T2'),:,:)  = DCChange_t1(strcmp(group4,'Group3') & strcmp(tp,'T2'),:,:);
        ZChange_last(strcmp(group4,'Group3') & strcmp(tp,'T3'),:,:)  = DCChange_t2(strcmp(group4,'Group3') & strcmp(tp,'T3'),:,:);
        
        ZChange_last(strcmp(group4,'Control') & strcmp(tp,'T1'),:,:) = DCChange_bl(strcmp(group4,'Control') & strcmp(tp,'T1'),:,:);
        ZChange_last(strcmp(group4,'Control') & strcmp(tp,'T2'),:,:) = DCChange_t1(strcmp(group4,'Control') & strcmp(tp,'T2'),:,:);
        ZChange_last(strcmp(group4,'Control') & strcmp(tp,'T3'),:,:) = DCChange_t2(strcmp(group4,'Control') & strcmp(tp,'T3'),:,:);
        
        
        layer_400_last_ctx(:,ii,:) = ZChange_last;
    end
end

% save the functional masks in folder
for functional_networks = 1
    % functional decoding
for func_decoding = 1
    % ingredients: cifti-matlab-master
    % functional data
    % schaefer 400 in 32k space
    
    addpath(path to '/cifti-matlab-master')
    
    HCP400_7 = ciftiopen(path to '/data/Schaefer2018_400Parcels_7Networks_order.dlabel.nii','/Applications/workbench/bin_macosx64/wb_command')
    
    files_lh = dir(path to '/data/*association-test_z_lh.shape.gii');
    files_rh = dir(path to '/data/*association-test_z_rh.shape.gii');
    
    maps_ns = zeros(400,5);
    for j = 1:5
        
        L = gifti([path to '/data/', files_lh(j).name]);
        R = gifti([path to '/data/', files_rh(j).name]);
        funcmap= [L.cdata;R.cdata];
        
        for i = 1:400
            maps_ns(i,j) = mean(funcmap(find(HCP400_7.cdata==i)));
        end
        maps_full(:,j) = funcmap;
    end

    for i  =  1:5
    maps_ns(find(maps_ns(:,i) < prctile(maps_ns(:,i),90)),i) = 0;
    end
    
    maps_ns2 = maps_ns;
    maps_ns(:,2) = maps_ns2(:,4);
    maps_ns(:,3) = maps_ns2(:,2);
    maps_ns(:,4) = maps_ns2(:,3);
    
    
    for j = 1:2
        tmp = zeros(20484,1);
        for i = 1:200
            tmp(find(parcels400==i+1)) = nanmean(maps_ns(i,j));
        end
        for i = 1:200
            tmp(find(parcels400==i+1001)) = nanmean(maps_ns(i+200,j));
        end
        f=figure,
        BoSurfStatViewData(tmp,SN,'')
        colormap((cbrewer('seq','YlOrBr',11)))
        %exportfigbo(f,[RPATH 'F1.attention_map', num2str(j) ,'.png'],'png', 10)
    end
    
    for j = 3:4
        tmp = zeros(20484,1);
        for i = 1:200
            tmp(find(parcels400==i+1)) = nanmean(maps_ns(i,j));
        end
        for i = 1:200
            tmp(find(parcels400==i+1001)) = nanmean(maps_ns(i+200,j));
        end
        f=figure,
        BoSurfStatViewData(tmp,SN,'')
        colormap((cbrewer('seq','Reds',11)))
        exportfigbo(f,[RPATH 'F1.affect_map', num2str(j) ,'.png'],'png', 10)
    end
    
    for j = 5
        tmp = zeros(20484,1);
        for i = 1:200
            tmp(find(parcels400==i+1)) = nanmean(maps_ns(i,j));
        end
        for i = 1:200
            tmp(find(parcels400==i+1001)) = nanmean(maps_ns(i+200,j));
        end
        f=figure,
        BoSurfStatViewData(tmp,SN,'')
        colormap(flipud(cbrewer('seq','Greens',11)))
        exportfigbo(f,[RPATH 'F1.tom_map.png'],'png', 10)
    end

    
end

end


% F1
for f1 = 1
   
    keepbl = mintersect(isthere_fc,find(~isnan(mean(GGG400,2))),find(fmri_qc==1),...
        find(tpnum==0), find(~isnan(mean(measure_used,2)).* (nanmean(measure_used,2)~=0)));
    
    
    maps_ns_bin = maps_ns;
    maps_ns_bin(maps_ns>0) = 1
    M = 1;
    M2 = 1 + term(maps_ns_bin)
    slm0 = SurfStatLinMod(mean(GGG400(keepbl,:))',M)
    slm1 = SurfStatLinMod(mean(GGG400(keepbl,:))',M2)
    slm = SurfStatF(slm0,slm1)
    SurfStatP(slm)
    
      for human_scatter = 1
                    to_brain_m = GGG400;
                    change_p = [];
                    change_p.c = nanmean(to_brain_m(keepbl,find(maps_ns(:,1))),2);
                    change_p.v = nanmean(to_brain_m(keepbl,find(maps_ns(:,2))),2);
                    change_p.sn = nanmean(to_brain_m(keepbl,find(maps_ns(:,3))),2);
                    change_p.da = nanmean(to_brain_m(keepbl,find(maps_ns(:,4))),2);
                    change_p.i = nanmean(to_brain_m(keepbl,find(maps_ns(:,5))),2);
                    cl = [0.9,    0.8,    0.1; 0.9844,    0.8828,    0.0195;     0.6471    0.0588    0.0824 ;0.9882    0.7333    0.6314; 0.1922,    0.6392,    0.3294];
                    fig_position = [200 200 600 400]; % coordinates for figures
                    d{1} = change_p.c';
                    d{2} = change_p.v';
                    d{3} = change_p.sn';
                    d{4} = change_p.da';
                    d{5} = change_p.i';
                    means = cellfun(@mean, d);
                    means_ggg = means;
                    variances = cellfun(@std, d);
                    var_ggg = variances;
                    f = figure('Position', fig_position);
                    h1 = raincloud_plot(d{1}, 'box_on', 1, 'color', cl(1,:), 'alpha', 0.5,...
                        'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15,...
                        'box_col_match', 0);
                    h2 = raincloud_plot(d{2}, 'box_on', 1, 'color', cl(2,:), 'alpha', 0.5,...
                        'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0);
                    h3 = raincloud_plot(d{3}, 'box_on', 1, 'color', cl(3,:), 'alpha', 0.5,...
                        'box_dodge', 1, 'box_dodge_amount', .55, 'dot_dodge_amount', .55, 'box_col_match', 0);
                    h4 = raincloud_plot(d{4}, 'box_on', 1, 'color', cl(4,:), 'alpha', 0.5,...
                        'box_dodge', 1, 'box_dodge_amount', .75, 'dot_dodge_amount', .75,...
                        'box_col_match', 0);
                    h5 = raincloud_plot(d{5}, 'box_on', 1, 'color', cl(5,:), 'alpha', 0.5,...
                        'box_dodge', 1, 'box_dodge_amount', .95, 'dot_dodge_amount', .95, 'box_col_match', 0);
                    set(gca,'Xlim',[0.02 0.16], 'YLim', [-45 45]);
                    box off
                    exportfigbo(f,[RPATH 'F1.eccentricity_netwoorks.png'],'png', 10)
                end

end

for f2 = 1
    
    %GGG
        Xn            = squeeze(layer_400_last(:,1,:)); % check for the presence of microstructure
        Xn(isnan(Xn)) = 0;
        Xn(isinf(Xn)) = 0;
        Xm            = ggg_400_last;
        Xm(isnan(Xm)) = 0;
        Xm(isinf(Xm)) = 0;
        keep1   = union(find(strcmp(groupN,'Affect')),    find(strcmp(groupN,'Perspective')));
        keep1   = union(keep1,find(strcmp(groupN,'Presence')))
        keep1   = union(keep1,find(strcmp(groupN,'Control1')))
        keep2   = intersect(find(abs(mean(Xn,2)) <666), find(sum(Xn,2)~=0));
        keep3   = intersect(find(~strcmp(groupN,'Group3')),find(tpnum>0));  
        keep4   = intersect(find(abs(mean(Xm,2)) <666), find(sum(Xm,2)~=0));
        keep    = mintersect(keep1, keep2,keep3,keep4);
        
        ik      = id(keep,:);
        ink     = idnum(keep);
        tk      = tpnum(keep,:);
        tnk     = tp(keep,:);
        gk      = group(keep);
        g4k     = group4(keep);
        gNk     = groupN(keep);
        ak      = age(keep);
        sk      = cellstr(sex(keep));   
        
        A       = term(ak);
        S       = term(sk);
        G       = term(gk);
        G4      = term(g4k);
        GN      = term(gNk);
        Sub     = term(ink);
        Tn      = term(tnk);
        
        M        = 1 + A + S + GN + random(Sub) + I;
        
        slm      = SurfStatLinMod(ggg_400_last(keep,:),M);
        slm      = SurfStatT(slm,(GN.Presence))
        p        = 2*(1 - tcdf(abs(slm.t),slm.df))

        tmp = zeros(20484,1);
        for i = 1:200
            tmp(find(parcels400==i+1)) = slm.t(i);
        end
        for i = 1:200
            tmp(find(parcels400==i+1001)) =slm.t(i+200);
        end
        
        f=figure,
        BoSurfStatViewData(tmp,SN,'')
        BoSurfStatColLim([-4 4])
        colormap(flipud(cbrewer('div','RdBu',11)))
       % exportfigbo(f,[RPATH 'F1.Presence.t.png'],'png', 10)       
        
        slm      = SurfStatLinMod(ggg_400_last(keep,:),M);
        slm      = SurfStatT(slm,(GN.Perspective))
        p        = 2*(1 - tcdf(abs(slm.t),slm.df))

        tmp = zeros(20484,1);
        for i = 1:200
            tmp(find(parcels400==i+1)) = slm.t(i);
        end
        for i = 1:200
            tmp(find(parcels400==i+1001)) =slm.t(i+200);
        end
        
        f=figure,
        BoSurfStatViewData(tmp,SN,'')
        BoSurfStatColLim([-4 4])
        colormap(flipud(cbrewer('div','RdBu',11)))
        exportfigbo(f,[RPATH 'F1.Perspective.t.png'],'png', 10)
        
        slm      = SurfStatLinMod(ggg_400_last(keep,:),M);
        slm      = SurfStatT(slm,(GN.Affect))
        p        = 2*(1 - tcdf(abs(slm.t),slm.df))

        tmp = zeros(20484,1);
        for i = 1:200
            tmp(find(parcels400==i+1)) = slm.t(i);
        end
        for i = 1:200
            tmp(find(parcels400==i+1001)) =slm.t(i+200);
        end
        
        f=figure,
        BoSurfStatViewData(tmp,SN,'')
        BoSurfStatColLim([-4 4])
        colormap(flipud(cbrewer('div','RdBu',11)))
        exportfigbo(f,[RPATH 'F1.Affect.t.png'],'png', 10)
        
        slm      = SurfStatLinMod(ggg_400_last(keep,:),M);
        slm      = SurfStatT(slm,(GN.Control1))
        p        = 2*(1 - tcdf(abs(slm.t),slm.df))

        tmp = zeros(20484,1);
        for i = 1:200
            tmp(find(parcels400==i+1)) = slm.t(i);
        end
        for i = 1:200
            tmp(find(parcels400==i+1001)) =slm.t(i+200);
        end
        
        f=figure,
        BoSurfStatViewData(tmp,SN,'')
        BoSurfStatColLim([-4 4])
        colormap(flipud(cbrewer('div','RdBu',11)))
        exportfigbo(f,[RPATH 'F1.Cnotrol.t.png'],'png', 10)
 
        slm      = SurfStatLinMod(ggg_400_last(keep,:),M);
        slm      = SurfStatT(slm,(GN.Presence-GN.Perspective))
        p        = 2*(1 - tcdf(abs(slm.t),slm.df))
           
        hp = fdr_bh(p(1:400),0.05)
        h  = hp;
        
        f1_h  = h;

        tmp = zeros(20484,1);
        for i = 1:200
            tmp(find(parcels400==i+1)) = slm.t(i).*h(i);
        end
        for i = 1:200
            tmp(find(parcels400==i+1001)) =slm.t(i+200).*h(i+200);
        end
        
        f=figure,
        BoSurfStatViewData(tmp,SN,'')
        BoSurfStatColLim([-4 4])
        colormap(flipud(cbrewer('div','RdBu',11)))
        exportfigbo(f,[RPATH 'F1.Perspective-Presence.fdr.png'],'png', 10)
    
        tren = p<0.01;
            
        tmp = zeros(20484,1);
        for i = 1:200
            tmp(find(parcels400==i+1)) = slm.t(i).*tren(i);
        end
        for i = 1:200
            tmp(find(parcels400==i+1001)) =slm.t(i+200).*tren(i+200);
        end
        
        f=figure,
        BoSurfStatViewData(tmp,SN,'')
        BoSurfStatColLim([-4 4])
        colormap(flipud(cbrewer('div','RdBu',11)))
        exportfigbo(f,[RPATH 'F1.Perspective-Presence.trend.png'],'png', 10)
  
        Xn = ggg_400_last;
        neuros = zeros(5,3);
        for i = 1:5
            slm      = SurfStatLinMod(mean(Xn(keep,maps_ns(:,i)>0),2),M);
            slm      = SurfStatT(slm,(GN.Presence-GN.Perspective))
            neuros(i,1) = slm.t;
            neurop(i,1) = 2*(1 - tcdf(abs(slm.t),slm.df));
        end
        for i = 1:5
            slm      = SurfStatLinMod(mean(Xn(keep,maps_ns(:,i)>0),2),M);
            slm      = SurfStatT(slm,(GN.Presence-GN.Affect))
            neuros(i,2) = slm.t;
            neurop(i,2) = 2*(1 - tcdf(abs(slm.t),slm.df));
        end
        for i = 1:5
            slm      = SurfStatLinMod(mean(Xn(keep,maps_ns(:,i)>0),2),M);
            slm      = SurfStatT(slm,(GN.Perspective-GN.Affect))
            neuros(i,3) = slm.t;
            neurop(i,3) = 2*(1 - tcdf(abs(slm.t),slm.df));
        end
        
       csvwrite([RPATH 'F2.GGG.csv'],[neuros(:,1),neurop(:,1),neuros(:,2),neurop(:,2),neuros(:,3),neurop(:,3)])

       for key_descriptives = 1
           Xn = ggg_400_last;
           meansdata = zeros(5,4);
           keep_presence = (find(strcmp(groupN(keep),'Presence')))
           keep_affect = (find(strcmp(groupN(keep),'Affect')))
           keep_perspective = (find(strcmp(groupN(keep),'Perspective')))
           keep_control   = (find(strcmp(groupN(keep),'Control1')))
           
           for i = 1:5
               x = mean(Xn(keep(keep_control),maps_ns(:,i)>0),2);
               meansdata(i,4) = mean(x);
               stddata(i,4) = std(x);
               SEM = std(x)/sqrt(length(x)); % Standard Error
               ts = tinv([0.025  0.975],length(x)-1);      % T-Score
               cidata(i,4,:) = mean(x) + ts*SEM;
           end
           for i = 1:5
               x = mean(Xn(keep(keep_affect),maps_ns(:,i)>0),2);
               meansdata(i,2) = mean(x);
               stddata(i,2) = std(x);
               SEM = std(x)/sqrt(length(x)); % Standard Error
               ts = tinv([0.025  0.975],length(x)-1);      % T-Score
               cidata(i,2,:) = mean(x) + ts*SEM;
           end
           for i = 1:5
               x = mean(Xn(keep(keep_perspective),maps_ns(:,i)>0),2);
               meansdata(i,3) = mean(x);
               stddata(i,3) = std(x);
               SEM = std(x)/sqrt(length(x)); % Standard Error
               ts = tinv([0.025  0.975],length(x)-1);      % T-Score
               cidata(i,3,:) = mean(x) + ts*SEM;
           end
           for i = 1:5
               x = mean(Xn(keep(keep_presence),maps_ns(:,i)>0),2);
               meansdata(i,1) = mean(x);
               stddata(i,1) = std(x);
               SEM = std(x)/sqrt(length(x)); % Standard Error
               ts = tinv([0.025  0.975],length(x)-1);      % T-Score
               cidata(i,1,:) = mean(x) + ts*SEM;
           end
           
           csvwrite([RPATH 'sF2.GGGraws_network_presence.csv'],[meansdata(:,1),stddata(:,1),cidata(:,1,1),cidata(:,1,2)])
           csvwrite([RPATH 'sF2.GGGraws_network_affect.csv'],[meansdata(:,2),stddata(:,2),cidata(:,2,1),cidata(:,2,2)])
           csvwrite([RPATH 'sF2.GGGraws_network_perspective.csv'],[meansdata(:,3),stddata(:,3),cidata(:,3,1),cidata(:,3,2)])
           csvwrite([RPATH 'sF2.GGGraws_network_control.csv'],[meansdata(:,4),stddata(:,4),cidata(:,4,1),cidata(:,4,2)])

   
       end
        
        %% post-hoc scatter change 
        for run_plot =1
            bh = mean(GGG400(:,f1_h>0),2);
            keep = (find(keeproi2.*keeproi3))
            Ck = bh;
            gk = group4(keep);
            ik = id(keep);
            dk = daysfrombl(keep);
            tk = tpnum(keep);
            myaxis = [0 400 1300 1700];
            %% mean thickness in ROI
            mY  = Ck(keep);
            
            %% define range
            
            myaxis = [min(dk)-2 max(dk)+2 min(mY)-0.01 max(mY)+0.01];
       
            
            %% define colors for groups
                        
            c_1  = [0.4196    0.6824    0.8392]; % controls
            c_1b = [0.2 0.2 0.2]; % retest
            
            c_2 = [0.9412,    0.2314,    0.1255];
            c_3 = [0.1922,    0.6392,    0.3294];
            
            c_4 = [0.4 0.2 0];
            c_4b = [ 0.9844,    0.8828,    0.0195];
            
            
            %% Group 1
            f=figure,
            subplot(2,3,1)
            hold on
            % 1.  change in affective
            keep = find(strcmp(gk,'Group_1'));
            ikk  = ik(keep);
            dkk  = dk(keep);
            tkk   = tk(keep);
            X    = mY(keep);
            
            hold on
            scatter(dkk,X,20,'k','fill')
            axis(myaxis);
            Subjects = unique(ikk);
            for j = 1:length(Subjects)
                su = find(strcmp(ikk,Subjects(j)));
                for s = 1:(length(su)-1)
                    plot([dkk(su(s)),dkk(su(s+1))],[X(su(s)),X(su(s+1))],...
                        'Color',[0.5 0.5 0.5])
                end
                outliers.G1.subj_code{j} = Subjects{j};
                outliers.G1.subj_range(j,1:2) = [min(X(su)) max(X(su))];
            end
            
            keeps = intersect(find(strcmp(gk,'Group_1')),find(tk<2));
            ik2 = ik(keeps);
            dk2 = tk(keeps);
            model1  = 1 +  term(dk2) + random(term(ik2)) + I;
            slm     = SurfStatLinMod( mY(keeps), model1);
            slm1     = SurfStatT(slm, dk2)
            
            keeps = mintersect(find(strcmp(gk,'Group_1')),find(tk>0),find(tk<3));
            ik2 = ik(keeps);
            dk2 = tk(keeps);
            model1  = 1 +  term(dk2) + random(term(ik2)) + I;
            slm     = SurfStatLinMod( mY(keeps), model1);
            slm2     = SurfStatT(slm, dk2)
            
            keeps = mintersect(find(strcmp(gk,'Group_1')),find(tk>1),find(tk<4));
            ik2 = ik(keeps);
            dk2 = tk(keeps);
            model1  = 1 +  term(dk2) + random(term(ik2)) + I;
            slm     = SurfStatLinMod( mY(keeps), model1);
            slm3     = SurfStatT(slm, dk2)
            
            xrange  = linspace(0,1,100);
            xrange1  = linspace(0,100,100);
            ymod    = slm1.coef(1) + (xrange)*slm1.coef(2);
            plot(xrange1,ymod,'Color',c_4b,'LineWidth',4)
            
            xrange  = linspace(100,200,100);
            xrange1 = linspace(1,2,100);
            ymod    = slm2.coef(1) + (xrange1)*slm2.coef(2);
            plot(xrange,ymod,'Color',c_2,'LineWidth',4)
            
            xrange  = linspace(2,3,100);
            xrange1 = linspace(200,300,100);
            ymod    = slm3.coef(1) + (xrange)*slm3.coef(2);
            plot(xrange1,ymod,'Color',c_3,'LineWidth',4)
            title('Cohort 1');
            
            
            %% Group 2
            subplot(2,3,2)
            keep = find(strcmp(gk,'Group_2'));
            ikk  = ik(keep);
            dkk  = dk(keep);
            X    = mY(keep);
            tkk  = tk(keep);
            %xlabel('days from T0')
            
            scatter(dkk,X,20,'k','fill')
            axis(myaxis);
            hold on

            Subjects = unique(ikk);
            for j = 1:length(Subjects)
                su = find(strcmp(ikk,Subjects(j)));
                for s = 1:(length(su)-1)
                    plot([dkk(su(s)),dkk(su(s+1))],[X(su(s)),X(su(s+1))],...
                        'Color',[0.5 0.5 0.5])
                end
                outliers.G2.subj_code{j} = Subjects{j};
                outliers.G2.subj_range(j,1:2) = [min(X(su)) max(X(su))];
            end
            keeps = intersect(find(strcmp(gk,'Group_2')),find(tk<2));
            ik2 = ik(keeps);
            dk2 = tk(keeps);
            model1  = 1 +  term(dk2) + random(term(ik2)) + I;
            slm     = SurfStatLinMod( mY(keeps), model1);
            slm1     = SurfStatT(slm, dk2)
            keeps = mintersect(find(strcmp(gk,'Group_2')),find(tk>0),find(tk<3));
            ik2 = ik(keeps);
            dk2 = tk(keeps);
            model1  = 1 +  term(dk2) + random(term(ik2)) + I;
            slm     = SurfStatLinMod( mY(keeps), model1);
            slm2     = SurfStatT(slm, dk2)
            keeps = mintersect(find(strcmp(gk,'Group_2')),find(tk>1),find(tk<4));
            ik2 = ik(keeps);
            dk2 = tk(keeps);
            model1  = 1 +  term(dk2) + random(term(ik2)) + I;
            slm     = SurfStatLinMod( mY(keeps), model1);
            slm3     = SurfStatT(slm, dk2)
            
            xrange  = linspace(0,1,100);
            xrange1  = linspace(0,100,100);
            yrange  = linspace( slm1.coef(1), slm2.coef(1),100)
            ymod    = slm1.coef(1) + (xrange)*slm1.coef(2);
            plot(xrange1,ymod,'Color',c_4b,'LineWidth',4)
            
            xrange  = linspace(100,200,100);
            xrange1 = linspace(1,2,100);
            yrange  = linspace( slm2.coef(1), slm3.coef(1),100)
            ymod    = slm2.coef(1) + (xrange1)*slm2.coef(2);
            plot(xrange,ymod,'Color',c_3,'LineWidth',4)
            xrange  = linspace(2,3,100);
            xrange1 = linspace(200,300,100);
            ymod    = slm3.coef(1) + (xrange)*slm3.coef(2);
            plot(xrange1,ymod,'Color',c_2,'LineWidth',4)

            title('Cohort 2');
            
            %% Controls
            subplot(2,3,3)
            keep = find(strcmp(gk,'Control'));
            ikk  = ik(keep);
            dkk  = dk(keep);
            X    = mY(keep);
            tkk  = tk(keep);
                        
            scatter(dkk,X,20,'k','fill')
            axis(myaxis);
            hold on

            Subjects = unique(ikk);
            for j = 1:length(Subjects)
                su = find(strcmp(ikk,Subjects(j)));
                for s = 1:(length(su)-1)
                    plot([dkk(su(s)),dkk(su(s+1))],[X(su(s)),X(su(s+1))],...
                        'Color',[0.5 0.5 0.5])
                end
                outliers.C.subj_code{j} = Subjects{j};
                outliers.C.subj_range(j,1:2) = [min(X(su)) max(X(su))];
            end
            keeps = intersect(find(strcmp(gk,'Control')),find(tk<2));
            ik2 = ik(keeps);
            dk2 = tk(keeps);
            model1  = 1 +  term(dk2) + random(term(ik2)) +I;
            slm     = SurfStatLinMod( mY(keeps), model1);
            slm1     = SurfStatT(slm, dk2)
            keeps = mintersect(find(strcmp(gk,'Control')),find(tk>0),find(tk<3));
            ik2 = ik(keeps);
            dk2 = tk(keeps);
            model1  = 1 +  term(dk2) + random(term(ik2)) +I;
            slm     = SurfStatLinMod( mY(keeps), model1);
            slm2     = SurfStatT(slm, dk2)
            keeps = mintersect(find(strcmp(gk,'Control')),find(tk>1),find(tk<4));
            ik2 = ik(keeps);
            dk2 = tk(keeps);
            model1  = 1 +  term(dk2) + random(term(ik2)) +I;
            slm     = SurfStatLinMod( mY(keeps), model1);
            slm3     = SurfStatT(slm, dk2)
            
            xrange  = linspace(0,1,100);
            xrange1  = linspace(0,100,100);
            yrange  = linspace( slm1.coef(1), slm2.coef(1),100)
            ymod    = slm1.coef(1) + (xrange)*slm1.coef(2);
            plot(xrange1,ymod,'Color',c_1,'LineWidth',4)
            
            xrange  = linspace(100,200,100);
            xrange1 = linspace(1,2,100);
            yrange  = linspace( slm2.coef(1), slm3.coef(1),100)
            ymod    = slm2.coef(1) + (xrange1)*slm2.coef(2);
            plot(xrange,ymod,'Color',c_1,'LineWidth',4)
            xrange  = linspace(2,3,100);
            xrange1 = linspace(200,300,100);
            ymod    = slm3.coef(1) + (xrange)*slm3.coef(2);
            plot(xrange1,ymod,'Color',c_1,'LineWidth',4)
            title('Controls');
            
            keep = find(tp_t0_there.*tp_t1_there.*(bh>0));
            Ck = bh;
            gk = group4(keep);
            ik = id(keep);
            dk = daysfrombl(keep);
            tk = tpnum(keep);
            myaxes = [0 400 2.1 2.8];
            %% mean thickness in ROI
            mY  = Ck(keep);
            subplot(2,3,4)
            keep = find(strcmp(gk,'Group3'));
            ikk  = ik(keep);
            dkk  = dk(keep);
            X    = mY(keep);
            tkk  = tk(keep);
            
            scatter(dkk,X,20,'k','fill')
            xlim([0 400])
            hold on
            Subjects = unique(ikk);
            for j = 1:length(Subjects)
                su = find(strcmp(ikk,Subjects(j)));
                for s = 1:(length(su)-1)
                    plot([dkk(su(s)),dkk(su(s+1))],[X(su(s)),X(su(s+1))],...
                        'Color',[0.5 0.5 0.5])
                end
                outliers.G3.subj_code{j} = Subjects{j};
                outliers.G3.subj_range(j,1:2) = [min(X(su)) max(X(su))];
            end
            keeps = intersect(find(strcmp(gk,'Group3')),find(tk<2));
            ik2 = ik(keeps);
            dk2 = tk(keeps);
            model1  = 1 +  term(dk2) + term(ik2);
            slm     = SurfStatLinMod( mY(keeps), model1);
            slm1     = SurfStatT(slm, dk2)
            xrange  = linspace(0,1,100);
            xrange1  = linspace(0,100,100);
            yrange  = linspace( slm1.coef(1), slm2.coef(1),100)
            ymod    = slm1.coef(1) + (xrange)*slm1.coef(2);
            plot(xrange1,ymod,'Color',c_2,'LineWidth',4)            
            title('Group 3');
        end
        exportfigbo(f,[RPATH 'F2.meanGGG.pers-pres.png'],'png', 10)
     
        %% post-hoc which G1-G2-G3
                
        for tmaps = 1
            keep_presence = (find(strcmp(groupN(keep),'Presence')))
            keep_affect = (find(strcmp(groupN(keep),'Affect')))
            keep_perspective = (find(strcmp(groupN(keep),'Perspective')))
            keep_control   = (find(strcmp(groupN(keep),'Control1')))
         
            neuros = zeros(1,3);
            for i = 1:5
                slm      = SurfStatLinMod(mean(g1_400_last(keep,find(f1_h>0)),2),M);
                slm      = SurfStatT(slm,(GN.Presence-GN.Perspective))
                neuros(i,1) = slm.t;
                neurop(i,1) = 2*(1 - tcdf(abs(slm.t),slm.df));
            end
            for i = 1
                slm      = SurfStatLinMod(mean(g1_400_last(keep,find(f1_h>0)),2),M);
                slm      = SurfStatT(slm,(GN.Presence-GN.Affect))
                neuros(i,2) = slm.t;
                neurop(i,2) = 2*(1 - tcdf(abs(slm.t),slm.df));
            end
            for i = 1
                slm      = SurfStatLinMod(mean(g1_400_last(keep,find(f1_h>0)),2),M);
                slm      = SurfStatT(slm,(GN.Perspective-GN.Affect))
                neuros(i,3) = slm.t;
                neurop(i,3) = 2*(1 - tcdf(abs(slm.t),slm.df));
            end
            
  
        for human_scatter = 1
            to_brain_m = g1_400_last;
            change_p = [];
            change_p.c = mean(to_brain_m(keep(keep_control),find(hp>0)),2);
            change_p.v = mean(to_brain_m(keep(keep_presence),find(hp>0)),2);
            change_p.sn = mean(to_brain_m(keep(keep_affect),find(hp>0)),2);
            change_p.da = mean(to_brain_m(keep(keep_perspective),find(hp>0)),2);
            cl = [0.6196    0.7922    0.8824; 0.9844,    0.8828,    0.0195;0.9412,    0.2314,    0.1255;0.1922,    0.6392,    0.3294];
            fig_position = [200 200 600 400]; % coordinates for figures
            d{1} = change_p.c';
            d{2} = change_p.v';
            d{3} = change_p.sn';
            d{4} = change_p.da';
            means = cellfun(@mean, d);
            variances = cellfun(@std, d);
            f = figure('Position', fig_position);
            h1 = raincloud_plot(d{1}, 'box_on', 1, 'color', cl(1,:), 'alpha', 0.5,...
                'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15,...
                'box_col_match', 0);
            h2 = raincloud_plot(d{2}, 'box_on', 1, 'color', cl(2,:), 'alpha', 0.5,...
                'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0);
            h3 = raincloud_plot(d{3}, 'box_on', 1, 'color', cl(3,:), 'alpha', 0.5,...
                'box_dodge', 1, 'box_dodge_amount', .55, 'dot_dodge_amount', .55, 'box_col_match', 0);
            h4 = raincloud_plot(d{4}, 'box_on', 1, 'color', cl(4,:), 'alpha', 0.5,...
                'box_dodge', 1, 'box_dodge_amount', .75, 'dot_dodge_amount', .75,...
                'box_col_match', 0);
            box off
            exportfigbo(f,[RPATH 'F2.func.fdr.g1.png'],'png', 10)
        end
           
     
        end
        
        for tmaps = 1
            
            % mean in yeo networks for each..
            keep_presence = (find(strcmp(groupN(keep),'Presence')))
            keep_affect = (find(strcmp(groupN(keep),'Affect')))
            keep_perspective = (find(strcmp(groupN(keep),'Perspective')))
            keep_control   = (find(strcmp(groupN(keep),'Control1')))
         
            neuros = zeros(1,3);
            for i = 1
                slm      = SurfStatLinMod(mean(g2_400_last(keep,find(f1_h>0)),2),M);
                slm      = SurfStatT(slm,(GN.Presence-GN.Perspective))
                neuros(i,1) = slm.t;
                neurop(i,1) = 2*(1 - tcdf(abs(slm.t),slm.df));
            end
            for i = 1
                slm      = SurfStatLinMod(mean(g2_400_last(keep,find(f1_h>0)),2),M);
                slm      = SurfStatT(slm,(GN.Presence-GN.Affect))
                neuros(i,2) = slm.t;
                neurop(i,2) = 2*(1 - tcdf(abs(slm.t),slm.df));
            end
            for i = 1
                slm      = SurfStatLinMod(mean(g2_400_last(keep,find(f1_h>0)),2),M);
                slm      = SurfStatT(slm,(GN.Perspective-GN.Affect))
                neuros(i,3) = slm.t;
                neurop(i,3) = 2*(1 - tcdf(abs(slm.t),slm.df));
            end
            
  
        for human_scatter = 1
            to_brain_m = g2_400_last;
            change_p = [];
            change_p.c = mean(to_brain_m(keep(keep_control),find(hp>0)),2);
            change_p.v = mean(to_brain_m(keep(keep_presence),find(hp>0)),2);
            change_p.sn = mean(to_brain_m(keep(keep_affect),find(hp>0)),2);
            change_p.da = mean(to_brain_m(keep(keep_perspective),find(hp>0)),2);
            cl = [0.6196    0.7922    0.8824; 0.9844,    0.8828,    0.0195;0.9412,    0.2314,    0.1255;0.1922,    0.6392,    0.3294];
            fig_position = [200 200 600 400]; % coordinates for figures
            d{1} = change_p.c';
            d{2} = change_p.v';
            d{3} = change_p.sn';
            d{4} = change_p.da';
            means = cellfun(@mean, d);
            variances = cellfun(@std, d);
            f = figure('Position', fig_position);
            h1 = raincloud_plot(d{1}, 'box_on', 1, 'color', cl(1,:), 'alpha', 0.5,...
                'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15,...
                'box_col_match', 0);
            h2 = raincloud_plot(d{2}, 'box_on', 1, 'color', cl(2,:), 'alpha', 0.5,...
                'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0);
            h3 = raincloud_plot(d{3}, 'box_on', 1, 'color', cl(3,:), 'alpha', 0.5,...
                'box_dodge', 1, 'box_dodge_amount', .55, 'dot_dodge_amount', .55, 'box_col_match', 0);
            h4 = raincloud_plot(d{4}, 'box_on', 1, 'color', cl(4,:), 'alpha', 0.5,...
                'box_dodge', 1, 'box_dodge_amount', .75, 'dot_dodge_amount', .75,...
                'box_col_match', 0);
            %set(gca,'XLim', [-1 1], 'YLim', [-5 5]);
            box off
            exportfigbo(f,[RPATH 'F2.func.fdr.g2.png'],'png', 10)
        end
           
     
        end

        for tmaps = 1
            
               % mean in yeo networks for each..
            keep_presence = (find(strcmp(groupN(keep),'Presence')))
            keep_affect = (find(strcmp(groupN(keep),'Affect')))
            keep_perspective = (find(strcmp(groupN(keep),'Perspective')))
            keep_control   = (find(strcmp(groupN(keep),'Control1')))
         
            neuros = zeros(1,3);
            for i = 1
                slm      = SurfStatLinMod(mean(g3_400_last(keep,find(f1_h>0)),2),M);
                slm      = SurfStatT(slm,(GN.Presence-GN.Perspective))
                neuros(i,1) = slm.t;
                neurop(i,1) = 1 - tcdf(slm.t,slm.df);
            end
            for i = 1
                slm      = SurfStatLinMod(mean(g3_400_last(keep,find(f1_h>0)),2),M);
                slm      = SurfStatT(slm,(GN.Presence-GN.Affect))
                neuros(i,2) = slm.t;
                neurop(i,2) = 1 - tcdf(slm.t,slm.df);
            end
            for i = 1
                slm      = SurfStatLinMod(mean(g3_400_last(keep,find(f1_h>0)),2),M);
                slm      = SurfStatT(slm,(GN.Perspective-GN.Affect))
                neuros(i,3) = slm.t;
                neurop(i,3) = 1 - tcdf(slm.t,slm.df);
            end
            
  
        for human_scatter = 1
            to_brain_m = g3_400_last;
            change_p = [];
            change_p.c = mean(to_brain_m(keep(keep_control),find(hp>0)),2);
            change_p.v = mean(to_brain_m(keep(keep_presence),find(hp>0)),2);
            change_p.sn = mean(to_brain_m(keep(keep_affect),find(hp>0)),2);
            change_p.da = mean(to_brain_m(keep(keep_perspective),find(hp>0)),2);
            cl = [0.6196    0.7922    0.8824; 0.9844,    0.8828,    0.0195;0.9412,    0.2314,    0.1255;0.1922,    0.6392,    0.3294];
            fig_position = [200 200 600 400]; % coordinates for figures
            d{1} = change_p.c';
            d{2} = change_p.v';
            d{3} = change_p.sn';
            d{4} = change_p.da';
            means = cellfun(@mean, d);
            variances = cellfun(@std, d);
            f = figure('Position', fig_position);
            h1 = raincloud_plot(d{1}, 'box_on', 1, 'color', cl(1,:), 'alpha', 0.5,...
                'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15,...
                'box_col_match', 0);
            h2 = raincloud_plot(d{2}, 'box_on', 1, 'color', cl(2,:), 'alpha', 0.5,...
                'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0);
            h3 = raincloud_plot(d{3}, 'box_on', 1, 'color', cl(3,:), 'alpha', 0.5,...
                'box_dodge', 1, 'box_dodge_amount', .55, 'dot_dodge_amount', .55, 'box_col_match', 0);
            h4 = raincloud_plot(d{4}, 'box_on', 1, 'color', cl(4,:), 'alpha', 0.5,...
                'box_dodge', 1, 'box_dodge_amount', .75, 'dot_dodge_amount', .75,...
                'box_col_match', 0);
            %set(gca,'XLim', [-1 1], 'YLim', [-5 5]);
            box off
            exportfigbo(f,[RPATH 'F2.func.fdr.g3.png'],'png', 10)
        end
           
     
        end
       
        %% G1-G2-G3 - tabel supplements
        
        for PeA =1
            Xn = g1_400_last;
            neuros = zeros(5,3);
            for i = 1:5
                slm      = SurfStatLinMod(mean(Xn(keep,maps_ns(:,i)>0),2),M);
                slm      = SurfStatT(slm,(GN.Perspective-GN.Affect))
                neuros(i,1) = slm.t;
                neurop(i,1) = 2*(1 - tcdf(abs(slm.t),slm.df));
            end
            Xn = g2_400_last;
            for i = 1:5
                slm      = SurfStatLinMod(mean(Xn(keep,maps_ns(:,i)>0),2),M);
                slm      = SurfStatT(slm,(GN.Perspective-GN.Affect))
                neuros(i,2) = slm.t;
                neurop(i,2) = 2*(1 - tcdf(abs(slm.t),slm.df));
            end
            Xn = g3_400_last;
            for i = 1:5
                slm      = SurfStatLinMod(mean(Xn(keep,maps_ns(:,i)>0),2),M);
                slm      = SurfStatT(slm,(GN.Perspective-GN.Affect))
                neuros(i,3) = slm.t;
                neurop(i,3) = 2*(1 - tcdf(abs(slm.t),slm.df));
            end
   
            csvwrite([RPATH 'F2.Perspective-Affect_G1-G3.csv'],[neuros(:,1),neurop(:,1),neuros(:,2),neurop(:,2),neuros(:,3),neurop(:,3)])
            
        end
       
        for PrA =1
            Xn = g1_400_last;
            neuros = zeros(5,3);
            for i = 1:5
                slm      = SurfStatLinMod(mean(Xn(keep,maps_ns(:,i)>0),2),M);
                slm      = SurfStatT(slm,(GN.Presence-GN.Affect))
                neuros(i,1) = slm.t;
                neurop(i,1) = 2*(1 - tcdf(abs(slm.t),slm.df));
            end
            Xn = g2_400_last;
            for i = 1:5
                slm      = SurfStatLinMod(mean(Xn(keep,maps_ns(:,i)>0),2),M);
                slm      = SurfStatT(slm,(GN.Presence-GN.Affect))
                neuros(i,2) = slm.t;
                neurop(i,2) = 2*(1 - tcdf(abs(slm.t),slm.df));
            end
            Xn = g3_400_last;
            for i = 1:5
                slm      = SurfStatLinMod(mean(Xn(keep,maps_ns(:,i)>0),2),M);
                slm      = SurfStatT(slm,(GN.Presence-GN.Affect))
                neuros(i,3) = slm.t;
                neurop(i,3) = 2*(1 - tcdf(abs(slm.t),slm.df));
            end
            
        csvwrite([RPATH 'F2.Presence-Affect_G1-G3.csv'],[neuros(:,1),neurop(:,1),neuros(:,2),neurop(:,2),neuros(:,3),neurop(:,3)])
        
        end
      
        for PrP =1 
        Xn = g1_400_last;
        neuros = zeros(5,3);
        for i = 1:5
            slm      = SurfStatLinMod(mean(Xn(keep,maps_ns(:,i)>0),2),M);
            slm      = SurfStatT(slm,(GN.Presence-GN.Perspective))
            neuros(i,1) = slm.t;
            neurop(i,1) = 2*(1 - tcdf(abs(slm.t),slm.df));
        end
        Xn = g2_400_last;
        for i = 1:5
            slm      = SurfStatLinMod(mean(Xn(keep,maps_ns(:,i)>0),2),M);
            slm      = SurfStatT(slm,(GN.Presence-GN.Perspective))
            neuros(i,2) = slm.t;
            neurop(i,2) = 2*(1 - tcdf(abs(slm.t),slm.df));
        end
        Xn = g3_400_last;
        for i = 1:5
            slm      = SurfStatLinMod(mean(Xn(keep,maps_ns(:,i)>0),2),M);
            slm      = SurfStatT(slm,(GN.Presence-GN.Perspective))
            neuros(i,3) = slm.t;
            neurop(i,3) = 2*(1 - tcdf(abs(slm.t),slm.df));
        end
        
       csvwrite([RPATH 'F2.Presence-Perspective_G1-G3.csv'],[neuros(:,1),neurop(:,1),neuros(:,2),neurop(:,2),neuros(:,3),neurop(:,3)])

       end

    
end

% eccentriciy - tc1 tc2  - pres vs tc3 - all vs control
% cortical thicknness corrected Main findings
for f2_subs = 1

    for tc1_ggg = 1
        Xn            = squeeze(layer_400_last(:,1,:));
        Xn(isnan(Xn)) = 0;
        Xn(isinf(Xn)) = 0;
        Xm            = ggg_400_last;
        Xm(isnan(Xm)) = 0;
        Xm(isinf(Xm)) = 0;
        keep1   = union(find(strcmp(groupN,'Affect')),    find(strcmp(groupN,'Perspective')));
        keep1   = union(keep1,find(strcmp(groupN,'Presence')))
        keep1   = union(keep1,find(strcmp(groupN,'Control1')))
        keep2   = intersect(find(abs(mean(Xn,2)) <666), find(sum(Xn,2)~=0));
        keep3   = intersect(find(strcmp(group4,'Group_1')),find(tpnum>0));  
        keep4   = intersect(find(abs(mean(Xm,2)) <666), find(sum(Xm,2)~=0));
        keep    = mintersect(keep1, keep2,keep3,keep4);
        
        ik      = id(keep,:);
        ink     = idnum(keep);
        tk      = tpnum(keep,:);
        tnk     = tp(keep,:);
        gk      = group(keep);
        g4k     = group4(keep);
        gNk     = groupN(keep);
        ak      = age(keep);
        sk      = cellstr(sex(keep));   
        
        A       = term(ak);
        S       = term(sk);
        G       = term(gk);
        G4      = term(g4k);
        GN      = term(gNk);
        Sub     = term(ink);
        Tn      = term(tnk);
        
        M        = 1 + A + S + GN + random(Sub) + I;
        
        slm      = SurfStatLinMod(ggg_400_last(keep,:),M);
        slm      = SurfStatT(slm,(GN.Presence-GN.Perspective))
        p        = 2*(1 - tcdf(abs(slm.t),slm.df))
           
        hp = fdr_bh(p(1:400),0.025)
        h  = hp;
        
        tmp = zeros(20484,1);
        for i = 1:200
            tmp(find(parcels400==i+1)) = slm.t(i).*h(i);
        end
        for i = 1:200
            tmp(find(parcels400==i+1001)) =slm.t(i+200).*h(i+200);
        end
        
        f=figure,
        BoSurfStatViewData(tmp,SN,'')
        BoSurfStatColLim([-4 4])
        colormap(flipud(cbrewer('div','RdBu',11)))
        exportfigbo(f,[RPATH 'F1.Perspective-Presence.fdr.TC1.png'],'png', 10)
    
 
        Xn = ggg_400_last;
        neuros = zeros(5,3);
        neurop = zeros(5,3);
        for i = 1:5
            slm      = SurfStatLinMod(mean(Xn(keep,maps_ns(:,i)>0),2),M);
            slm      = SurfStatT(slm,(GN.Presence-GN.Perspective))
            neuros(i,1) = slm.t;
            neurop(i,1) = 2*(1 - tcdf(abs(slm.t),slm.df));
        end
        for i = 1:5
            slm      = SurfStatLinMod(mean(Xn(keep,maps_ns(:,i)>0),2),M);
            slm      = SurfStatT(slm,(GN.Presence-GN.Affect))
            neuros(i,2) = slm.t;
            neurop(i,2) = 2*(1 - tcdf(abs(slm.t),slm.df));
        end
        for i = 1:5
            slm      = SurfStatLinMod(mean(Xn(keep,maps_ns(:,i)>0),2),M);
            slm      = SurfStatT(slm,(GN.Perspective-GN.Affect))
            neuros(i,3) = slm.t;
            neurop(i,3) = 2*(1 - tcdf(abs(slm.t),slm.df));
        end
        
       csvwrite([RPATH 'F2s.GGG_TC1.csv'],[neuros(:,1),neurop(:,1),neuros(:,2),neurop(:,2),neuros(:,3),neurop(:,3)])
    end
    
    for tc2_ggg = 1
        Xn            = squeeze(layer_400_last(:,1,:));
        Xn(isnan(Xn)) = 0;
        Xn(isinf(Xn)) = 0;
        Xm            = ggg_400_last;
        Xm(isnan(Xm)) = 0;
        Xm(isinf(Xm)) = 0;
        keep1   = union(find(strcmp(groupN,'Affect')),    find(strcmp(groupN,'Perspective')));
        keep1   = union(keep1,find(strcmp(groupN,'Presence')))
        keep1   = union(keep1,find(strcmp(groupN,'Control1')))
        keep2   = intersect(find(abs(mean(Xn,2)) <666), find(sum(Xn,2)~=0));
        keep3   = intersect(find(strcmp(group4,'Group_2')),find(tpnum>0));  
        keep4   = intersect(find(abs(mean(Xm,2)) <666), find(sum(Xm,2)~=0));
        keep    = mintersect(keep1, keep2,keep3,keep4);
        
        ik      = id(keep,:);
        ink     = idnum(keep);
        tk      = tpnum(keep,:);
        tnk     = tp(keep,:);
        gk      = group(keep);
        g4k     = group4(keep);
        gNk     = groupN(keep);
        ak      = age(keep);
        sk      = cellstr(sex(keep));   
        
        A       = term(ak);
        S       = term(sk);
        G       = term(gk);
        G4      = term(g4k);
        GN      = term(gNk);
        Sub     = term(ink);
        Tn      = term(tnk);
        
        M        = 1 + A + S + GN + random(Sub) + I;
 
        slm      = SurfStatLinMod(ggg_400_last(keep,:),M);
        slm      = SurfStatT(slm,(GN.Presence-GN.Perspective))
        p        = 2*(1 - tcdf(abs(slm.t),slm.df))
           
        hp = fdr_bh(p(1:400),0.025)
        h  = hp;
 
        tmp = zeros(20484,1);
        for i = 1:200
            tmp(find(parcels400==i+1)) = slm.t(i).*h(i);
        end
        for i = 1:200
            tmp(find(parcels400==i+1001)) =slm.t(i+200).*h(i+200);
        end
        
        f=figure,
        BoSurfStatViewData(tmp,SN,'')
        BoSurfStatColLim([-4 4])
        colormap(flipud(cbrewer('div','RdBu',11)))
        exportfigbo(f,[RPATH 'F1.Perspective-Presence.fdr.TC2.png'],'png', 10)
    
        Xn = ggg_400_last;
        neuros = zeros(5,3);
        for i = 1:5
            slm      = SurfStatLinMod(mean(Xn(keep,maps_ns(:,i)>0),2),M);
            slm      = SurfStatT(slm,(GN.Presence-GN.Perspective))
            neuros(i,1) = slm.t;
            neurop(i,1) = 2*(1 - tcdf(abs(slm.t),slm.df));
        end
        for i = 1:5
            slm      = SurfStatLinMod(mean(Xn(keep,maps_ns(:,i)>0),2),M);
            slm      = SurfStatT(slm,(GN.Presence-GN.Affect))
            neuros(i,2) = slm.t;
            neurop(i,2) = 2*(1 - tcdf(abs(slm.t),slm.df));
        end
        for i = 1:5
            slm      = SurfStatLinMod(mean(Xn(keep,maps_ns(:,i)>0),2),M);
            slm      = SurfStatT(slm,(GN.Perspective-GN.Affect))
            neuros(i,3) = slm.t;
            neurop(i,3) = 2*(1 - tcdf(abs(slm.t),slm.df));
        end
        
       csvwrite([RPATH 'F2s.GGG_TC2.csv'],[neuros(:,1),neurop(:,1),neuros(:,2),neurop(:,2),neuros(:,3),neurop(:,3)])
    end

    for T0_T1_ggg = 1
        Xn            = squeeze(layer_400_last(:,1,:));
        Xn(isnan(Xn)) = 0;
        Xn(isinf(Xn)) = 0;
        Xm            = ggg_400_last;
        Xm(isnan(Xm)) = 0;
        Xm(isinf(Xm)) = 0;
        keep1   = union(find(strcmp(groupN,'Affect')),    find(strcmp(groupN,'Perspective')));
        keep1   = union(keep1,find(strcmp(groupN,'Presence')))
        keep1   = union(keep1,find(strcmp(groupN,'Control1')))
        keep2   = intersect(find(abs(mean(Xn,2)) <666), find(sum(Xn,2)~=0));
        keep3   = find(tpnum==1);  
        keep4   = intersect(find(abs(mean(Xm,2)) <666), find(sum(Xm,2)~=0));
        keep    = mintersect(keep1, keep2,keep3,keep4);
        
        ik      = id(keep,:);
        ink     = idnum(keep);
        tk      = tpnum(keep,:);
        tnk     = tp(keep,:);
        gk      = group(keep);
        g4k     = group4(keep);
        gNk     = groupN(keep);
        ak      = age(keep);
        sk      = cellstr(sex(keep));   
        
        A       = term(ak);
        S       = term(sk);
        G       = term(gk);
        G4      = term(g4k);
        GN      = term(gNk);
        Sub     = term(ink);
        Tn      = term(tnk);
        
        M        = 1 + A + S + GN + random(Sub) + I;
 
        Xn = ggg_400_last;
        neuros = zeros(5,3);
        for i = 1:5
            slm      = SurfStatLinMod(mean(Xn(keep,maps_ns(:,i)>0),2),M);
            slm      = SurfStatT(slm,(GN.Presence-GN.Affect))
            neuros(i,1) = slm.t;
            neurop(i,1) = 2*(1 - tcdf(abs(slm.t),slm.df));
        end
        for i = 1:5
            slm      = SurfStatLinMod(mean(Xn(keep,maps_ns(:,i)>0),2),M);
            slm      = SurfStatT(slm,(GN.Presence-GN.Control1))
            neuros(i,2) = slm.t;
            neurop(i,2) = 2*(1 - tcdf(abs(slm.t),slm.df));
        end
        for i = 1:5
            slm      = SurfStatLinMod(mean(Xn(keep,maps_ns(:,i)>0),2),M);
            slm      = SurfStatT(slm,(GN.Affect-GN.Control1))
            neuros(i,3) = slm.t;
            neurop(i,3) = 2*(1 - tcdf(abs(slm.t),slm.df));
        end
        
       csvwrite([RPATH 'F2s.GGG_T0-T1.csv'],[neuros(:,1),neurop(:,1),neuros(:,2),neurop(:,2),neuros(:,3),neurop(:,3)])
    end
    
    for T1_T3_ggg = 1
        Xn            = squeeze(layer_400_last(:,1,:));
        Xn(isnan(Xn)) = 0;
        Xn(isinf(Xn)) = 0;
        Xm            = ggg_400_last;
        Xm(isnan(Xm)) = 0;
        Xm(isinf(Xm)) = 0;
        keep1   = union(find(strcmp(groupN,'Affect')),    find(strcmp(groupN,'Perspective')));
        keep1   = union(keep1,find(strcmp(groupN,'Presence')))
        keep1   = union(keep1,find(strcmp(groupN,'Control1')))
        keep2   = intersect(find(abs(mean(Xn,2)) <666), find(sum(Xn,2)~=0));
        keep3   = find(tpnum>1);  
        keep4   = intersect(find(abs(mean(Xm,2)) <666), find(sum(Xm,2)~=0));
        keep    = mintersect(keep1, keep2,keep3,keep4);
        
        ik      = id(keep,:);
        ink     = idnum(keep);
        tk      = tpnum(keep,:);
        tnk     = tp(keep,:);
        gk      = group(keep);
        g4k     = group4(keep);
        gNk     = groupN(keep);
        ak      = age(keep);
        sk      = cellstr(sex(keep));   
        
        A       = term(ak);
        S       = term(sk);
        G       = term(gk);
        G4      = term(g4k);
        GN      = term(gNk);
        Sub     = term(ink);
        Tn      = term(tnk);
        
        M        = 1 + A + S + GN + random(Sub) + I;
 
        Xn = ggg_400_last;
        neuros = zeros(5,3);
        for i = 1:5
            slm      = SurfStatLinMod(mean(Xn(keep,maps_ns(:,i)>0),2),M);
            slm      = SurfStatT(slm,(GN.Perspective-GN.Control1))
            neuros(i,2) = slm.t;
            neurop(i,2) = 2*(1 - tcdf(abs(slm.t),slm.df));
        end
        for i = 1:5
            slm      = SurfStatLinMod(mean(Xn(keep,maps_ns(:,i)>0),2),M);
            slm      = SurfStatT(slm,(GN.Affect-GN.Control1))
            neuros(i,3) = slm.t;
            neurop(i,3) = 2*(1 - tcdf(abs(slm.t),slm.df));
        end
        
       csvwrite([RPATH 'F2s.GGG_T1-T3.csv'],[neuros(:,2),neurop(:,2),neuros(:,3),neurop(:,3)])
    end

    for cortical_thickness_control = 1
        Xn            = squeeze(layer_400_last(:,1,:));
        Xn(isnan(Xn)) = 0;
        Xn(isinf(Xn)) = 0;
        Xm            = ggg_400_last;
        Xm(isnan(Xm)) = 0;
        Xm(isinf(Xm)) = 0;
        keep1   = union(find(strcmp(groupN,'Affect')),    find(strcmp(groupN,'Perspective')));
        keep1   = union(keep1,find(strcmp(groupN,'Presence')))
        keep1   = union(keep1,find(strcmp(groupN,'Control1')))
        keep2   = intersect(find(abs(mean(Xn,2)) <666), find(sum(Xn,2)~=0));
        keep3   = intersect(find(~strcmp(groupN,'Group3')),find(tpnum>0));  
        keep4   = intersect(find(abs(mean(Xm,2)) <666), find(sum(Xm,2)~=0));
        keep    = mintersect(keep1, keep2,keep3,keep4);
        
        ik      = id(keep,:);
        ink     = idnum(keep);
        tk      = tpnum(keep,:);
        tnk     = tp(keep,:);
        gk      = group(keep);
        g4k     = group4(keep);
        gNk     = groupN(keep);
        ak      = age(keep);
        sk      = cellstr(sex(keep));   
        
        A       = term(ak);
        S       = term(sk);
        G       = term(gk);
        G4      = term(g4k);
        GN      = term(gNk);
        Sub     = term(ink);
        Tn      = term(tnk);
        
        M        = 1 + A + S + GN + random(Sub) + I;
 
        Xn = ggg_400_last_ctx;
        neuros = zeros(5,3);
        for i = 1:5
            slm      = SurfStatLinMod(mean(Xn(keep,maps_ns(:,i)>0),2),M);
            slm      = SurfStatT(slm,(GN.Presence-GN.Perspective))
            neuros(i,1) = slm.t;
            neurop(i,1) = 2*(1 - tcdf(abs(slm.t),slm.df));
        end
        for i = 1:5
            slm      = SurfStatLinMod(mean(Xn(keep,maps_ns(:,i)>0),2),M);
            slm      = SurfStatT(slm,(GN.Presence-GN.Affect))
            neuros(i,2) = slm.t;
            neurop(i,2) = 2*(1 - tcdf(abs(slm.t),slm.df));
        end
        for i = 1:5
            slm      = SurfStatLinMod(mean(Xn(keep,maps_ns(:,i)>0),2),M);
            slm      = SurfStatT(slm,(GN.Perspective-GN.Affect))
            neuros(i,3) = slm.t;
            neurop(i,3) = 2*(1 - tcdf(abs(slm.t),slm.df));
        end
        
       csvwrite([RPATH 'F2s.GGGctx.csv'],[neuros(:,1),neurop(:,1),neuros(:,2),neurop(:,2),neuros(:,3),neurop(:,3)])

    end
    
    for all_training_vs_control = 1
         Xn            = squeeze(layer_400_bl(:,1,:));
        Xn(isnan(Xn)) = 0;
        Xn(isinf(Xn)) = 0;
        Xm            = ggg_400_bl;
        Xm(isnan(Xm)) = 0;
        Xm(isinf(Xm)) = 0;
        keep1   = union(find(strcmp(groupN,'Affect')),    find(strcmp(groupN,'Perspective')));
        keep1   = union(keep1,find(strcmp(groupN,'Presence')))
        keep1   = union(keep1,find(strcmp(groupN,'Control1')))
        keep2   = intersect(find(abs(mean(Xn,2)) <666), find(sum(Xn,2)~=0));
        keep3   = intersect(find(~strcmp(groupN,'Group3')),find(tpnum>2));
        keep4   = intersect(find(abs(mean(Xm,2)) <666), find(sum(Xm,2)~=0));
        keep    = mintersect(keep1, keep2,keep3,keep4);
        
        ik      = id(keep,:);
        ink     = idnum(keep);
        tk      = tpnum(keep,:);
        tnk     = tp(keep,:);
        gk      = group(keep);
        g4k     = group4(keep);
        gNk     = groupN(keep);
        ak      = age(keep);
        sk      = cellstr(sex(keep));
        
        A       = term(ak);
        S       = term(sk);
        G       = term(gk);
        G4      = term(g4k);
        GN      = term(gNk);
        Sub     = term(ink);
        Tn      = term(tnk);
        g2k         = gNk;
        g2k(strcmp(gNk,'Control1')) = {'Control'};
        g2k(~strcmp(gNk,'Control1')) = {'Training'};
      
        G2  = term(g2k);   
        
        M        = 1 + A + S + G2;
   
        for i  = 1
            xxx = ggg_400_bl;
            slm      = SurfStatLinMod(xxx(keep,:),M);
            slm      = SurfStatT(slm,(G2.Training-G2.Control));
            p        = 2*(1 - tcdf(abs(slm.t),slm.df))
            
            hp = fdr_bh(p(1:400),0.05)
            h  = hp;
                        
            tmp = zeros(20484,1);
            for i = 1:200
                tmp(find(parcels400==i+1)) = slm.t(i).*h(i);
            end
            for i = 1:200
                tmp(find(parcels400==i+1001)) =slm.t(i+200).*h(i+200);
            end
            
            f=figure,
            BoSurfStatViewData(tmp,SN,'')
            BoSurfStatColLim([-4 4])
            colormap(flipud(cbrewer('div','RdBu',11)))
            exportfigbo(f,[RPATH 'F3s.training_control.ggg.fdr.png'],'png', 10)
            close(f)
        end
 
        neuros = zeros(5,3);
        for i  = 1
            xxx = ggg_400_bl;
            for j = 1:5
                slm      = SurfStatLinMod(mean(xxx(keep,maps_ns(:,j)>0),2),M);
                slm      = SurfStatT(slm,G2.Training-G2.Control)
                neuros(j,i) = slm.t;
                neurop(j,i) = 2*(1 - tcdf(abs(slm.t),slm.df));
            end
            
            csvwrite([RPATH 'F2s.GGG_all_control.csv'],[neuros(:,1),neurop(:,1)])

        end
         

    end

    for g1_g2_g3 = 1
        Xn            = squeeze(layer_400_last(:,1,:));
        Xn(isnan(Xn)) = 0;
        Xn(isinf(Xn)) = 0;
        Xm            = ggg_400_last;
        Xm(isnan(Xm)) = 0;
        Xm(isinf(Xm)) = 0;
        keep1   = union(find(strcmp(groupN,'Affect')),    find(strcmp(groupN,'Perspective')));
        keep1   = union(keep1,find(strcmp(groupN,'Presence')))
        keep1   = union(keep1,find(strcmp(groupN,'Control1')))
        keep2   = intersect(find(abs(mean(Xn,2)) <666), find(sum(Xn,2)~=0));
        keep3   = intersect(find(~strcmp(groupN,'Group3')),find(tpnum>0));
        keep4   = intersect(find(abs(mean(Xm,2)) <666), find(sum(Xm,2)~=0));
        keep    = mintersect(keep1, keep2,keep3,keep4);
        
        ik      = id(keep,:);
        ink     = idnum(keep);
        tk      = tpnum(keep,:);
        tnk     = tp(keep,:);
        gk      = group(keep);
        g4k     = group4(keep);
        gNk     = groupN(keep);
        ak      = age(keep);
        sk      = cellstr(sex(keep));
        
        A       = term(ak);
        S       = term(sk);
        G       = term(gk);
        G4      = term(g4k);
        GN      = term(gNk);
        Sub     = term(ink);
        Tn      = term(tnk);
        
        M        = 1 + A + S + GN + random(Sub) + I;
        
        for gradient1  =1
            slm      = SurfStatLinMod(g1_400_last(keep,:),M);
            slm      = SurfStatT(slm,(GN.Presence-GN.Perspective))
            p        = 2*(1 - tcdf(abs(slm.t),slm.df));
            
            hp = fdr_bh(p(1:400),0.05)
            h  = hp;
            
            tmp = zeros(20484,1);
            for i = 1:200
                tmp(find(parcels400==i+1)) = slm.t(i).*h(i);
            end
            for i = 1:200
                tmp(find(parcels400==i+1001)) =slm.t(i+200).*h(i+200);
            end
            
            f=figure,
            BoSurfStatViewData(tmp,SN,'')
            BoSurfStatColLim([-4 4])
            colormap(flipud(cbrewer('div','RdBu',11)))
            exportfigbo(f,[RPATH 'F2s.Presence-Perspective.G1.fdr.png'],'png', 10)
            
            tren = ((p<0.01))
            
            tmp = zeros(20484,1);
            for i = 1:200
                tmp(find(parcels400==i+1)) = slm.t(i).*tren(i);
            end
            for i = 1:200
                tmp(find(parcels400==i+1001)) =slm.t(i+200).*tren(i+200);
            end
            
            f=figure,
            BoSurfStatViewData(tmp,SN,'')
            BoSurfStatColLim([-4 4])
            colormap(flipud(cbrewer('div','RdBu',11)))
            exportfigbo(f,[RPATH 'F2s.Presence-Perspective.G1.p.png'],'png', 10)
            
        end
        for gradient2  =1
            slm      = SurfStatLinMod(g2_400_last(keep,:),M);
            slm      = SurfStatT(slm,(GN.Presence-GN.Perspective))
            p        = 2*(1 - tcdf(abs(slm.t),slm.df));
            
            hp = fdr_bh(p(1:400),0.05)
            h  = hp;
            
            tmp = zeros(20484,1);
            for i = 1:200
                tmp(find(parcels400==i+1)) = slm.t(i).*h(i);
            end
            for i = 1:200
                tmp(find(parcels400==i+1001)) =slm.t(i+200).*h(i+200);
            end
            
            f=figure,
            BoSurfStatViewData(tmp,SN,'')
            BoSurfStatColLim([-4 4])
            colormap(flipud(cbrewer('div','RdBu',11)))
            exportfigbo(f,[RPATH 'F2s.Presence-Perspective.G2.fdr.png'],'png', 10)
            
            
            tren = ((p<0.01) )
            
            tmp = zeros(20484,1);
            for i = 1:200
                tmp(find(parcels400==i+1)) = slm.t(i).*tren(i);
            end
            for i = 1:200
                tmp(find(parcels400==i+1001)) =slm.t(i+200).*tren(i+200);
            end
            
            f=figure,
            BoSurfStatViewData(tmp,SN,'')
            BoSurfStatColLim([-4 4])
            colormap(flipud(cbrewer('div','RdBu',11)))
            exportfigbo(f,[RPATH 'F2s.Presence-Perspective.G2.p.png'],'png', 10)
            
        end
        for gradient3  =1
            slm      = SurfStatLinMod(g3_400_last(keep,:),M);
            slm      = SurfStatT(slm,(GN.Presence-GN.Perspective))
            p        = 2*(1 - tcdf(abs(slm.t),slm.df));
            
            hp = fdr_bh(p(1:400),0.025)
            h  = hp;
            
            tmp = zeros(20484,1);
            for i = 1:200
                tmp(find(parcels400==i+1)) = slm.t(i).*h(i);
            end
            for i = 1:200
                tmp(find(parcels400==i+1001)) =slm.t(i+200).*h(i+200);
            end
            
            f=figure,
            BoSurfStatViewData(tmp,SN,'')
            BoSurfStatColLim([-4 4])
            colormap(flipud(cbrewer('div','RdBu',11)))
            exportfigbo(f,[RPATH 'F2s.Presence-Perspective.G3.fdr.png'],'png', 10)
            
            
            tren = ((p<0.01))
            
            tmp = zeros(20484,1);
            for i = 1:200
                tmp(find(parcels400==i+1)) = slm.t(i).*tren(i);
            end
            for i = 1:200
                tmp(find(parcels400==i+1001)) =slm.t(i+200).*tren(i+200);
            end
            
            f=figure,
            BoSurfStatViewData(tmp,SN,'')
            BoSurfStatColLim([-4 4])
            colormap(flipud(cbrewer('div','RdBu',11)))
            exportfigbo(f,[RPATH 'F2s.Presence-Perspective.G3.p.png'],'png', 10)
            
        end
        
        for gradient1  =1
            slm      = SurfStatLinMod(g1_400_last(keep,:),M);
            slm      = SurfStatT(slm,(GN.Presence-GN.Affect))
            p        = 2*(1 - tcdf(abs(slm.t),slm.df));
            
            hp = fdr_bh(p(1:400),0.05)
            h  = hp;
            
            tmp = zeros(20484,1);
            for i = 1:200
                tmp(find(parcels400==i+1)) = slm.t(i).*h(i);
            end
            for i = 1:200
                tmp(find(parcels400==i+1001)) =slm.t(i+200).*h(i+200);
            end
            
            f=figure,
            BoSurfStatViewData(tmp,SN,'')
            BoSurfStatColLim([-4 4])
            colormap(flipud(cbrewer('div','RdBu',11)))
            exportfigbo(f,[RPATH 'F2s.Presence-Affect.G1.fdr.png'],'png', 10)
            
            
            tren = ((p<0.01))
            
            tmp = zeros(20484,1);
            for i = 1:200
                tmp(find(parcels400==i+1)) = slm.t(i).*tren(i);
            end
            for i = 1:200
                tmp(find(parcels400==i+1001)) =slm.t(i+200).*tren(i+200);
            end
            
            f=figure,
            BoSurfStatViewData(tmp,SN,'')
            BoSurfStatColLim([-4 4])
            colormap(flipud(cbrewer('div','RdBu',11)))
            exportfigbo(f,[RPATH 'F2s.Presence-Affect.G1.p.png'],'png', 10)
            
        end
        for gradient2  =1
            slm      = SurfStatLinMod(g2_400_last(keep,:),M);
            slm      = SurfStatT(slm,(GN.Presence-GN.Affect))
            p        = 2*(1 - tcdf(abs(slm.t),slm.df))
            
            hp = fdr_bh(p(1:400),0.05)
            h  = hp;
            
            tmp = zeros(20484,1);
            for i = 1:200
                tmp(find(parcels400==i+1)) = slm.t(i).*h(i);
            end
            for i = 1:200
                tmp(find(parcels400==i+1001)) =slm.t(i+200).*h(i+200);
            end
            
            f=figure,
            BoSurfStatViewData(tmp,SN,'')
            BoSurfStatColLim([-4 4])
            colormap(flipud(cbrewer('div','RdBu',11)))
            exportfigbo(f,[RPATH 'F2s.Presence-Affect.G2.fdr.png'],'png', 10)
            
            tren = ((p<0.01))
            
            tmp = zeros(20484,1);
            for i = 1:200
                tmp(find(parcels400==i+1)) = slm.t(i).*tren(i);
            end
            for i = 1:200
                tmp(find(parcels400==i+1001)) =slm.t(i+200).*tren(i+200);
            end
            
            f=figure,
            BoSurfStatViewData(tmp,SN,'')
            BoSurfStatColLim([-4 4])
            colormap(flipud(cbrewer('div','RdBu',11)))
            exportfigbo(f,[RPATH 'F2s.Presence-Affect.G2.p.png'],'png', 10)
            
        end
        for gradient3  =1
            slm      = SurfStatLinMod(g3_400_last(keep,:),M);
            slm      = SurfStatT(slm,(GN.Presence-GN.Affect))
            p        = 2*(1 - tcdf(abs(slm.t),slm.df));
            
            hp = fdr_bh(p(1:400),0.025)
            h  = hp;
            
            tmp = zeros(20484,1);
            for i = 1:200
                tmp(find(parcels400==i+1)) = slm.t(i).*h(i);
            end
            for i = 1:200
                tmp(find(parcels400==i+1001)) =slm.t(i+200).*h(i+200);
            end
            
            f=figure,
            BoSurfStatViewData(tmp,SN,'')
            BoSurfStatColLim([-4 4])
            colormap(flipud(cbrewer('div','RdBu',11)))
            exportfigbo(f,[RPATH 'F2s.Presence-Affect.G3.fdr.png'],'png', 10)
            
            tren = ((p<0.01))
            
            tmp = zeros(20484,1);
            for i = 1:200
                tmp(find(parcels400==i+1)) = slm.t(i).*tren(i);
            end
            for i = 1:200
                tmp(find(parcels400==i+1001)) =slm.t(i+200).*tren(i+200);
            end
            
            f=figure,
            BoSurfStatViewData(tmp,SN,'')
            BoSurfStatColLim([-4 4])
            colormap(flipud(cbrewer('div','RdBu',11)))
            exportfigbo(f,[RPATH 'F2s.Presence-Affect.G3.p.png'],'png', 10)
            
        end
        
        for gradient1  =1
            slm      = SurfStatLinMod(g1_400_last(keep,:),M);
            slm      = SurfStatT(slm,(GN.Perspective-GN.Affect))
            p        = 2*(1 - tcdf(abs(slm.t),slm.df));
            
            hp = fdr_bh(p(1:400),0.05)
            h  = hp;
            
            tmp = zeros(20484,1);
            for i = 1:200
                tmp(find(parcels400==i+1)) = slm.t(i).*h(i);
            end
            for i = 1:200
                tmp(find(parcels400==i+1001)) =slm.t(i+200).*h(i+200);
            end
            
            f=figure,
            BoSurfStatViewData(tmp,SN,'')
            BoSurfStatColLim([-4 4])
            colormap(flipud(cbrewer('div','RdBu',11)))
            exportfigbo(f,[RPATH 'F2s.Perspective-Affect.G1.fdr.png'],'png', 10)
            
            tren = ((p<0.01))
            
            tmp = zeros(20484,1);
            for i = 1:200
                tmp(find(parcels400==i+1)) = slm.t(i).*tren(i);
            end
            for i = 1:200
                tmp(find(parcels400==i+1001)) =slm.t(i+200).*tren(i+200);
            end
            
            f=figure,
            BoSurfStatViewData(tmp,SN,'')
            BoSurfStatColLim([-4 4])
            colormap(flipud(cbrewer('div','RdBu',11)))
            exportfigbo(f,[RPATH 'F2s.Perspective-Affect.G1.p.png'],'png', 10)
            
        end
        for gradient2  =1
            slm      = SurfStatLinMod(g2_400_last(keep,:),M);
            slm      = SurfStatT(slm,(GN.Perspective-GN.Affect))
            p        = 2*(1 - tcdf(abs(slm.t),slm.df));
            
            hp = fdr_bh(p(1:400),0.05)
            h  = hp;
            
            tmp = zeros(20484,1);
            for i = 1:200
                tmp(find(parcels400==i+1)) = slm.t(i).*h(i);
            end
            for i = 1:200
                tmp(find(parcels400==i+1001)) =slm.t(i+200).*h(i+200);
            end
            
            f=figure,
            BoSurfStatViewData(tmp,SN,'')
            BoSurfStatColLim([-4 4])
            colormap(flipud(cbrewer('div','RdBu',11)))
            exportfigbo(f,[RPATH 'F2s.Perspective-Affect.G2.fdr.png'],'png', 10)
            
            tren = ((p<0.01) )
            
            tmp = zeros(20484,1);
            for i = 1:200
                tmp(find(parcels400==i+1)) = slm.t(i).*tren(i);
            end
            for i = 1:200
                tmp(find(parcels400==i+1001)) =slm.t(i+200).*tren(i+200);
            end
            
            f=figure,
            BoSurfStatViewData(tmp,SN,'')
            BoSurfStatColLim([-4 4])
            colormap(flipud(cbrewer('div','RdBu',11)))
            exportfigbo(f,[RPATH 'F2s.Perspective-Affect.G2.p.png'],'png', 10)
        end
        for gradient3  =1
            slm      = SurfStatLinMod(g3_400_last(keep,:),M);
            slm      = SurfStatT(slm,(GN.Perspective-GN.Affect))
            p        = 2*(1 - tcdf(abs(slm.t),slm.df))
            
            hp = fdr_bh(p(1:400),0.05)
            h  = hp;
            
            tmp = zeros(20484,1);
            for i = 1:200
                tmp(find(parcels400==i+1)) = slm.t(i).*h(i);
            end
            for i = 1:200
                tmp(find(parcels400==i+1001)) =slm.t(i+200).*h(i+200);
            end
            
            f=figure,
            BoSurfStatViewData(tmp,SN,'')
            BoSurfStatColLim([-4 4])
            colormap(flipud(cbrewer('div','RdBu',11)))
            exportfigbo(f,[RPATH 'F2s.Perspective-Affect.G3.fdr.png'],'png', 10)
            
            tren = ((p<0.01) )
            
            tmp = zeros(20484,1);
            for i = 1:200
                tmp(find(parcels400==i+1)) = slm.t(i).*tren(i);
            end
            for i = 1:200
                tmp(find(parcels400==i+1001)) =slm.t(i+200).*tren(i+200);
            end
            
            f=figure,
            BoSurfStatViewData(tmp,SN,'')
            BoSurfStatColLim([-4 4])
            colormap(flipud(cbrewer('div','RdBu',11)))
            exportfigbo(f,[RPATH 'F2s.Perspective-Affect.G3.p.png'],'png', 10)
        end
        
    end
    
    for gsr_change = 1
          Xn            = squeeze(layer_400_last(:,1,:));
        Xn(isnan(Xn)) = 0;
        Xn(isinf(Xn)) = 0;
        Xm            = ggg_400_gsr;
        Xm(isnan(Xm)) = 0;
        Xm(isinf(Xm)) = 0;
        keep1   = union(find(strcmp(groupN,'Affect')),    find(strcmp(groupN,'Perspective')));
        keep1   = union(keep1,find(strcmp(groupN,'Presence')))
        keep1   = union(keep1,find(strcmp(groupN,'Control1')))
        keep2   = intersect(find(abs(mean(Xn,2)) <666), find(sum(Xn,2)~=0));
        keep3   = intersect(find(~strcmp(groupN,'Group3')),find(tpnum>0));  
        keep4   = intersect(find(abs(mean(Xm,2)) <666), find(sum(Xm,2)~=0));
        keep    = mintersect(keep1, keep2,keep3,keep4);
        
        ik      = id(keep,:);
        ink     = idnum(keep);
        tk      = tpnum(keep,:);
        tnk     = tp(keep,:);
        gk      = group(keep);
        g4k     = group4(keep);
        gNk     = groupN(keep);
        ak      = age(keep);
        sk      = cellstr(sex(keep));   
        
        A       = term(ak);
        S       = term(sk);
        G       = term(gk);
        G4      = term(g4k);
        GN      = term(gNk);
        Sub     = term(ink);
        Tn      = term(tnk);
        
        M        = 1 + A + S + GN + random(Sub) + I;
        
        Xn = ggg_400_gsr;
        neuros = zeros(5,3);
        for i = 1:5
            slm      = SurfStatLinMod(mean(Xn(keep,maps_ns(:,i)>0),2),M);
            slm      = SurfStatT(slm,(GN.Presence-GN.Perspective))
            neuros(i,1) = slm.t;
            neurop(i,1) = 2*(1 - tcdf(abs(slm.t),slm.df));
        end
        for i = 1:5
            slm      = SurfStatLinMod(mean(Xn(keep,maps_ns(:,i)>0),2),M);
            slm      = SurfStatT(slm,(GN.Presence-GN.Affect))
            neuros(i,2) = slm.t;
            neurop(i,2) = 2*(1 - tcdf(abs(slm.t),slm.df));
        end
        for i = 1:5
            slm      = SurfStatLinMod(mean(Xn(keep,maps_ns(:,i)>0),2),M);
            slm      = SurfStatT(slm,(GN.Perspective-GN.Affect))
            neuros(i,3) = slm.t;
            neurop(i,3) = 2*(1 - tcdf(abs(slm.t),slm.df));
        end
    end
end

for f3 = 1
    % scatter
    keepbl = mintersect(isthere_fc,find(~isnan(mean(GGG400,2))),find(fmri_qc==1),...
        find(tpnum==0), find(~isnan(mean(measure_used,2)).* (nanmean(measure_used,2)~=0)));
    
    
    load('/Users/sofievalk/Documents/Tools/ScientificColourMaps6/bilbao/bilbao.mat')
    
    meanrs = squeeze(mean(MPC_layers(keepbl,1:12,:),2));
    tmp = zeros(20484,1);
    for i = 1:200
        tmp(find(parcels400==i+1)) = mean(meanrs(:,i));
    end
    for i = 1:200
        tmp(find(parcels400==i+1001)) =mean(meanrs(:,i+200));
    end
    
    f=figure,
    BoSurfStatViewData(tmp,SN,'')
    BoSurfStatColLim([1350 1600])
    colormap(flipud(bilbao))
    exportfigbo(f,[RPATH 'F1.t1q.png'],'png', 10)
   
    Xn            = squeeze(layer_400_last(:,1,:));
    Xn(isnan(Xn)) = 0;
    Xn(isinf(Xn)) = 0;
    Xm            = ggg_400_last;
    Xm(isnan(Xm)) = 0;
    Xm(isinf(Xm)) = 0;
    keep1   = union(find(strcmp(groupN,'Affect')),    find(strcmp(groupN,'Perspective')));
    keep1   = union(keep1,find(strcmp(groupN,'Presence')))
    keep1   = union(keep1,find(strcmp(groupN,'Control1')))
    keep2   = intersect(find(abs(mean(Xn,2)) <666), find(sum(Xn,2)~=0));
    keep3   = intersect(find(~strcmp(groupN,'Group3')),find(tpnum>0));
    keep4   = intersect(find(abs(mean(Xm,2)) <666), find(sum(Xm,2)~=0));
    keep    = mintersect(keep1, keep2,keep3,keep4);
    
    ik      = id(keep,:);
    ink     = idnum(keep);
    tk      = tpnum(keep,:);
    tnk     = tp(keep,:);
    gk      = group(keep);
    g4k     = group4(keep);
    gNk     = groupN(keep);
    ak      = age(keep);
    sk      = cellstr(sex(keep));
    
    A       = term(ak);
    S       = term(sk);
    G       = term(gk);
    G4      = term(g4k);
    GN      = term(gNk);
    Sub     = term(ink);
    Tn      = term(tnk);
    
    M        = 1 + A + S + GN + random(Sub) + I;
    
    for func_in_t1q = 1
        t_dept = zeros(12,3);
        p_dept = zeros(12,3);
        for k=1:12
            Xx =  squeeze(mean(layer_400_last(:,k,:),2));
            slm      = SurfStatLinMod(mean(Xx(keep,f1_h>0),2),M);
            slm      = SurfStatT(slm,(GN.Presence-GN.Perspective))
            slm2 = slm.t;
            p_dept(k,1)  = 2*(1 - tcdf(abs(slm.t),slm.df))
            t_dept(k,1) = slm.t;
            slm      = SurfStatLinMod(median(Xx(keep,f1_h>0),2),M);
            slm      = SurfStatT(slm,(GN.Presence-GN.Affect))
            slm2 = slm.t;
            p_dept(k,2)  = 2*(1 - tcdf(abs(slm.t),slm.df))
            t_dept(k,2) = slm.t;
            slm      = SurfStatLinMod(median(Xx(keep,f1_h>0),2),M);
            slm      = SurfStatT(slm,(GN.Perspective-GN.Affect))
            slm2 = slm.t;
            p_dept(k,3)  = 2*(1 - tcdf(abs(slm.t),slm.df))
            t_dept(k,3) = slm.t;
        end
        
        f=figure,
        imagesc(t_dept.*(fdr_bh(p_dept(:,1))),[-4 4])
        colormap(flipud(cbrewer('div','RdBu',99)))
        colorbar
        exportfigbo(f,[RPATH 'F3.layer2func.super.stat.png'],'png', 10)
        
        f=figure,
        imagesc(t_dept,[-4 4])
        colormap(flipud(cbrewer('div','RdBu',99)))
        exportfigbo(f,[RPATH 'F3.layer2func.super.raw.png'],'png', 10)
        
        
        for human_profile = 1
            f=figure,
            ax = axes;
            ax.ColorOrder = [0.9844,    0.8828,    0.0195; 0.1922,    0.6392,    0.3294; 0.9412,    0.2314,    0.1255; ...
                0.8275    0.8902    0.9529; 0 0 0];
            hold on
            keep_presence = (find(strcmp(groupN(keep),'Presence')))
            keep_perspective = (find(strcmp(groupN(keep),'Perspective')))
            keep_affect = (find(strcmp(groupN(keep),'Affect')))
            keep_control1 = (find(strcmp(groupN(keep),'Control1')))
            
            plot(mean(mean(squeeze(layer_400_last(keep(keep_presence),:,hp>0))),3),1:12,'LineWidth',10)
            hold on
            plot(mean(mean(squeeze(layer_400_last(keep(keep_perspective),:,hp>0))),3),1:12,'LineWidth',10)
            hold on
            plot(mean(mean(squeeze(layer_400_last(keep(keep_affect),:,hp>0))),3),1:12,'LineWidth',10)
            hold on
            plot(mean(mean(squeeze(layer_400_last(keep(keep_control1),:,hp>0))),3),1:12,'LineWidth',10)
            hold on
            exportfigbo(f,[RPATH 'F3.Presence-Perspective.func.layerT1q.fdr.png'],'png', 10)
        end
    end
    
    for profile_network = 1
        for i  =1:5
            for human_profile = 1
                f=figure,
                ax = axes;
                ax.ColorOrder = [0.9844,    0.8828,    0.0195; 0.1922,    0.6392,    0.3294; 0.9412,    0.2314,    0.1255; 0.6196    0.7922    0.8824];
                hold on
                keep_presence = (find(strcmp(groupN(keep),'Presence')))
                keep_perspective = (find(strcmp(groupN(keep),'Perspective')))
                keep_affect = (find(strcmp(groupN(keep),'Affect')))
                keep_control = (find(strcmp(groupN(keep),'Control1')))
                plot(mean(mean(squeeze(layer_400_last(keep(keep_presence),:,find(maps_ns(:,i)>0)))),3),1:12,'LineWidth',10)
                hold on
                plot(mean(mean(squeeze(layer_400_last(keep(keep_perspective),:,find(maps_ns(:,i)>0)))),3),1:12,'LineWidth',10)
                hold on
                plot(mean(mean(squeeze(layer_400_last(keep(keep_affect),:,find(maps_ns(:,i)>0)))),3),1:12,'LineWidth',10)
                hold on
                plot(mean(mean(squeeze(layer_400_last(keep(keep_control),:,find(maps_ns(:,i)>0)))),3),1:12,'LineWidth',10)
                exportfigbo(f,[RPATH 'F2.networks.meanT1q.',num2str(i) ,'.png'],'png', 10)
            end
        end
        for i = 1:5
            for k=1:12
                Xx =  squeeze(mean(layer_400_last(:,k,:),2));
                slm      = SurfStatLinMod(mean(Xx(keep,maps_ns(:,i)>0),2),M);
                slm      = SurfStatT(slm,(GN.Presence-GN.Perspective));
                slm2 = slm.t;
                p_dept(k,1,i)  = 2*(1 - tcdf(abs(slm.t),slm.df));
                t_dept(k,1,i) = slm.t;
                Xx =  squeeze(mean(layer_400_last(:,k,:),2));
                slm      = SurfStatLinMod(mean(Xx(keep,maps_ns(:,i)>0),2),M);
                slm      = SurfStatT(slm,(GN.Presence-GN.Affect));
                slm2 = slm.t;
                p_dept(k,2,i)  = 2*(1 - tcdf(abs(slm.t),slm.df));
                t_dept(k,2,i) = slm.t;
                Xx =  squeeze(mean(layer_400_last(:,k,:),2));
                slm      = SurfStatLinMod(mean(Xx(keep,maps_ns(:,i)>0),2),M);
                slm      = SurfStatT(slm,(GN.Perspective-GN.Affect));
                slm2 = slm.t;
                p_dept(k,3,i)  = 2*(1 - tcdf(abs(slm.t),slm.df));
                t_dept(k,3,i) = slm.t;
            end
        end
        
        fdrcor = (fdr_bh(p_dept,0.05));
        
        for i = 1:5
            f=figure,
            imagesc(t_dept(:,:,i).*fdrcor(:,:,i),[-5 5])
            colormap(flipud(cbrewer('div','RdBu',99)))
            exportfigbo(f,[RPATH 'F3.layer-network', num2str(i) ,'.png'],'png', 10)
        end
            %f=figure,
            %imagesc(t_dept(:,:,i),[-5 5])
            %colormap(flipud(cbrewer('div','RdBu',99)))
            %exportfigbo(f,[RPATH 'F3.layer-network', num2str(i) ,'.raw.png'],'png', 10)
      
        
        
        
        csvwrite([RPATH 'network', num2str(i),'layer_myelin_presence-perspective.csv'], [...
            t_dept(:,1,1), p_dept(:,1,1),...
            t_dept(:,1,2), p_dept(:,1,2),...
            t_dept(:,1,3), p_dept(:,1,3),...
            t_dept(:,1,4), p_dept(:,1,4),...
            t_dept(:,1,5), p_dept(:,1,5)])
        
        csvwrite([RPATH 'network', num2str(i),'layer_myelin_presence-affect.csv'], [...
            t_dept(:,2,1), p_dept(:,2,1),...
            t_dept(:,2,2), p_dept(:,2,2),...
            t_dept(:,2,3), p_dept(:,2,3),...
            t_dept(:,2,4), p_dept(:,2,4),...
            t_dept(:,2,5), p_dept(:,2,5)])
        
        csvwrite([RPATH 'network', num2str(i),'layer_myelin_perspective-affect.csv'], [...
            t_dept(:,3,1), p_dept(:,3,1),...
            t_dept(:,3,2), p_dept(:,3,2),...
            t_dept(:,3,3), p_dept(:,3,3),...
            t_dept(:,3,4), p_dept(:,3,4),...
            t_dept(:,3,5), p_dept(:,3,5)])
        
    end
  
    for change_attention = 1
        val             = Att;
        Xn              = ggg_400_last;
        Xn(isnan(Xn))   = 0;
        Xm              = squeeze(layer_400_last(:,1,:));
        Xm(isnan(Xm))   = 0;
        keep1   = (find(strcmp(groupN,'Presence')));
        keep2   = intersect(find(Xn(:,1) >-666), find(sum(Xn,2)~=0));
        keep2b   = intersect(find(Xm(:,1) >-666), find(sum(Xm,2)~=0));
        keep4   = find(~(strcmp(group4,'Group3')));
        keep5   = find(tpnum==1);
        keep    = mintersect(keep1,keep2,keep2b,keep4,keep5);
        
        %eccentricity spin
        for gradient_diff = 1
            cd('/Users/sofievalk/brainstat_data/surface_data/tpl-fsaverage/fsaverage5/surf/')
            SPHERE_L      = SurfStatReadSurf({'lh.sphere'});
            SPHERE_R      = SurfStatReadSurf({'lh.sphere'});
            
            gradient_mask = mean(ggg_400_last(keep,:))';
            
            gradient_mask_fsa5 = nan(20484,1);
            for i = 1:200
                gradient_mask_fsa5(find(parcels400==i+1)) = gradient_mask(i);
            end
            for i = 1:200
                gradient_mask_fsa5(find(parcels400==i+1001)) = gradient_mask(i+200);
            end
   
            n_permutations = 1000;
            y_rand = spin_permutations({[gradient_mask_fsa5(1:10242)],[gradient_mask_fsa5(10243:end)]}, ...
                {SPHERE_L SPHERE_R}, ...
                n_permutations);
            
            % Merge the rotated data into single vectors
            g_rotated = squeeze([y_rand{1}(:,1,:);y_rand{2}(:,1,:)]);
            
            for i = 1:200
                gradient_mask_400(i,:) = median(g_rotated(find(parcels400==1+i),:));
            end
            for i = 1:200
                gradient_mask_400(i+200,:) = median(g_rotated(find(parcels400==1001+i),:));
            end
            
            r_original_thick_1 = zeros(12,1);
            pval_thick_spin_1 = zeros(12,1);
            prctile_rank_thick_1 = zeros(12,1);
            for i =  1:12
                [r_original_thick_1(i), pval_thick_spin_1(i)] = corr(mean(squeeze(layer_400_last(keep,i,:)))',mean(ggg_400_last(keep,:))', ...
                    'rows','pairwise')
                
                r_rand_thick = corr(mean(squeeze(layer_400_last(keep,i,:)))',gradient_mask_400, ...
                    'rows','pairwise');
                
                nperm = n_permutations;
                % Compute percentile rank.
                if r_original_thick_1(i) > 0
                    prctile_rank_thick_1(i) = sum(r_rand_thick > r_original_thick_1(i))/nperm;
                else
                    prctile_rank_thick_1(i) = sum(r_rand_thick < r_original_thick_1(i))/nperm;
                end
            end
            
        end
        %G1 spin
        for gradient_diff = 1
            cd('/Users/sofievalk/brainstat_data/surface_data/tpl-fsaverage/fsaverage5/surf/')
            SPHERE_L      = SurfStatReadSurf({'lh.sphere'});
            SPHERE_R      = SurfStatReadSurf({'lh.sphere'});
            
            gradient_mask = mean(g1_400_last(keep,:))';
            
            gradient_mask_fsa5 = nan(20484,1);
            for i = 1:200
                gradient_mask_fsa5(find(parcels400==i+1)) = gradient_mask(i);
            end
            for i = 1:200
                gradient_mask_fsa5(find(parcels400==i+1001)) = gradient_mask(i+200);
            end
            
            n_permutations = 1000;
            y_rand = spin_permutations({[gradient_mask_fsa5(1:10242)],[gradient_mask_fsa5(10243:end)]}, ...
                {SPHERE_L SPHERE_R}, ...
                n_permutations);
            
            % Merge the rotated data into single vectors
            g_rotated = squeeze([y_rand{1}(:,1,:);y_rand{2}(:,1,:)]);
            
            for i = 1:200
                gradient_mask_400(i,:) = median(g_rotated(find(parcels400==1+i),:));
            end
            for i = 1:200
                gradient_mask_400(i+200,:) = median(g_rotated(find(parcels400==1001+i),:));
            end
            
            r_original_thick_2 = zeros(12,1);
            pval_thick_spin_2 = zeros(12,1);
            prctile_rank_thick_2 = zeros(12,1);
            for i =  1:12
                [r_original_thick_2(i), pval_thick_spin_2(i)] = corr(mean(squeeze(layer_400_last(keep,i,:)))',mean(g1_400_last(keep,:))', ...
                    'rows','pairwise')
                
                r_rand_thick = corr(mean(squeeze(layer_400_last(keep,i,:)))',gradient_mask_400, ...
                    'rows','pairwise');
                
                % Compute percentile rank.
                if r_original_thick_2(i) > 0
                    prctile_rank_thick_2(i) = sum(r_rand_thick > r_original_thick_2(i))/nperm;
                else
                    prctile_rank_thick_2(i) = sum(r_rand_thick < r_original_thick_2(i))/nperm;
                end
            end
            
        end
        %G2 spin
        for gradient_diff = 1
            cd('/Users/sofievalk/brainstat_data/surface_data/tpl-fsaverage/fsaverage5/surf/')
            SPHERE_L      = SurfStatReadSurf({'lh.sphere'});
            SPHERE_R      = SurfStatReadSurf({'lh.sphere'});
            
            gradient_mask = mean(g2_400_last(keep,:))';
            
            gradient_mask_fsa5 = nan(20484,1);
            for i = 1:200
                gradient_mask_fsa5(find(parcels400==i+1)) = gradient_mask(i);
            end
            for i = 1:200
                gradient_mask_fsa5(find(parcels400==i+1001)) = gradient_mask(i+200);
            end
            
            n_permutations = 1000;
            y_rand = spin_permutations({[gradient_mask_fsa5(1:10242)],[gradient_mask_fsa5(10243:end)]}, ...
                {SPHERE_L SPHERE_R}, ...
                n_permutations);
            
            % Merge the rotated data into single vectors
            g_rotated = squeeze([y_rand{1}(:,1,:);y_rand{2}(:,1,:)]);
            
            for i = 1:200
                gradient_mask_400(i,:) = median(g_rotated(find(parcels400==1+i),:));
            end
            for i = 1:200
                gradient_mask_400(i+200,:) = median(g_rotated(find(parcels400==1001+i),:));
            end
            
            r_original_thick_3 = zeros(12,1);
            pval_thick_spin_3 = zeros(12,1);
            prctile_rank_thick_3 = zeros(12,1);
            for i =  1:12
                [r_original_thick_3(i), pval_thick_spin_3(i)] = corr(mean(squeeze(layer_400_last(keep,i,:)))',mean(g2_400_last(keep,:))', ...
                    'rows','pairwise')
                
                r_rand_thick = corr(mean(squeeze(layer_400_last(keep,i,:)))',gradient_mask_400, ...
                    'rows','pairwise');
                
                % Compute percentile rank.
                if r_original_thick_3(i) > 0
                    prctile_rank_thick_3(i) = sum(r_rand_thick > r_original_thick_3(i))/nperm;
                else
                    prctile_rank_thick_3(i) = sum(r_rand_thick < r_original_thick_3(i))/nperm;
                end
            end
            
        end
        %G3 spin
        for gradient_diff = 1
            cd('/Users/sofievalk/brainstat_data/surface_data/tpl-fsaverage/fsaverage5/surf/')
            SPHERE_L      = SurfStatReadSurf({'lh.sphere'});
            SPHERE_R      = SurfStatReadSurf({'lh.sphere'});
            
            gradient_mask = mean(g3_400_last(keep,:))';
            
            gradient_mask_fsa5 = nan(20484,1);
            for i = 1:200
                gradient_mask_fsa5(find(parcels400==i+1)) = gradient_mask(i);
            end
            for i = 1:200
                gradient_mask_fsa5(find(parcels400==i+1001)) = gradient_mask(i+200);
            end
            
            n_permutations = 1000;
            y_rand = spin_permutations({[gradient_mask_fsa5(1:10242)],[gradient_mask_fsa5(10243:end)]}, ...
                {SPHERE_L SPHERE_R}, ...
                n_permutations);
            
            % Merge the rotated data into single vectors
            g_rotated = squeeze([y_rand{1}(:,1,:);y_rand{2}(:,1,:)]);
            
            for i = 1:200
                gradient_mask_400(i,:) = median(g_rotated(find(parcels400==1+i),:));
            end
            for i = 1:200
                gradient_mask_400(i+200,:) = median(g_rotated(find(parcels400==1001+i),:));
            end
            
            r_original_thick_4 = zeros(12,1);
            pval_thick_spin_4 = zeros(12,1);
            prctile_rank_thick_4 = zeros(12,1);
            for i =  1:12
                [r_original_thick_4(i), pval_thick_spin_4(i)] = corr(mean(squeeze(layer_400_last(keep,i,:)))',mean(g3_400_last(keep,:))', ...
                    'rows','pairwise')
                
                r_rand_thick = corr(mean(squeeze(layer_400_last(keep,i,:)))',gradient_mask_400, ...
                    'rows','pairwise');
                
                % Compute percentile rank.
                if r_original_thick_4(i) > 0
                    prctile_rank_thick_4(i) = sum(r_rand_thick > r_original_thick_4(i))/nperm;
                else
                    prctile_rank_thick_4(i) = sum(r_rand_thick < r_original_thick_4(i))/nperm;
                end
            end
            
        end
        
        f=figure,
        ax = axes;
        ax.ColorOrder = cbrewer('qual','Set2',4);
        hold on
        plot(r_original_thick_1,'LineWidth',10),
        hold on
        plot(r_original_thick_2,'LineWidth',10),
        hold on
        plot(r_original_thick_3,'LineWidth',10),
        hold on
        plot(r_original_thick_4,'LineWidth',10),
        ylim([-0.5 0.5])
        exportfigbo(f,[RPATH 'F3.layervsfunc.presence.raw.png'],'png', 10)
        
        r_original_thick_1((prctile_rank_thick_1>0.05) )  = nan
        r_original_thick_2((prctile_rank_thick_2>0.05) )  = nan
        r_original_thick_3((prctile_rank_thick_3>0.05) )  = nan
        r_original_thick_4((prctile_rank_thick_4>0.05) )  = nan
        f=figure,
        ax = axes;
        ax.ColorOrder = cbrewer('qual','Set2',4);
        hold on
        plot(r_original_thick_1,'LineWidth',10),
        hold on
        plot(r_original_thick_2,'LineWidth',10),
        hold on
        plot(r_original_thick_3,'LineWidth',10),
        hold on
        plot(r_original_thick_4,'LineWidth',10),
        ylim([-0.5 0.5])
        exportfigbo(f,[RPATH 'F3.layervsfunc.presence.p.png'],'png', 10)

        
    end
    
    for change_compassion = 1
        val             = Comp;
        Xn              = ggg_400_last;
        Xn(isnan(Xn))   = 0;
        Xm              = squeeze(layer_400_last(:,1,:));
        Xm(isnan(Xm))   = 0;
        keep1   = (find(strcmp(groupN,'Affect')));
        keep2   = intersect(find(Xn(:,1) >-666), find(sum(Xn,2)~=0));
        keep2b   = intersect(find(Xm(:,1) >-666), find(sum(Xm,2)~=0));
        keep4   = find(~(strcmp(group4,'Group3')));
        keep5   = find(tpnum>0);
        keep6   = find(abs(val)<666);
        keep    = mintersect(keep1,keep2,keep2b,keep4,keep5);
        
        %eccentricity spin
               %eccentricity spin
        for gradient_diff = 1
            cd('/Users/sofievalk/brainstat_data/surface_data/tpl-fsaverage/fsaverage5/surf/')
            SPHERE_L      = SurfStatReadSurf({'lh.sphere'});
            SPHERE_R      = SurfStatReadSurf({'lh.sphere'});
            
            gradient_mask = mean(ggg_400_last(keep,:))';
            
            gradient_mask_fsa5 = nan(20484,1);
            for i = 1:200
                gradient_mask_fsa5(find(parcels400==i+1)) = gradient_mask(i);
            end
            for i = 1:200
                gradient_mask_fsa5(find(parcels400==i+1001)) = gradient_mask(i+200);
            end
            
            n_permutations = 1000;
            y_rand = spin_permutations({[gradient_mask_fsa5(1:10242)],[gradient_mask_fsa5(10243:end)]}, ...
                {SPHERE_L SPHERE_R}, ...
                n_permutations);
            
            % Merge the rotated data into single vectors
            g_rotated = squeeze([y_rand{1}(:,1,:);y_rand{2}(:,1,:)]);
            
            for i = 1:200
                gradient_mask_400(i,:) = median(g_rotated(find(parcels400==1+i),:));
            end
            for i = 1:200
                gradient_mask_400(i+200,:) = median(g_rotated(find(parcels400==1001+i),:));
            end
            
            r_original_thick_1 = zeros(12,1);
            pval_thick_spin_1 = zeros(12,1);
            prctile_rank_thick_1 = zeros(12,1);
            for i =  1:12
                [r_original_thick_1(i), pval_thick_spin_1(i)] = corr(mean(squeeze(layer_400_last(keep,i,:)))',mean(ggg_400_last(keep,:))', ...
                    'rows','pairwise')
                
                r_rand_thick = corr(mean(squeeze(layer_400_last(keep,i,:)))',gradient_mask_400, ...
                    'rows','pairwise');
                
                % Compute percentile rank.
                if r_original_thick_1(i) > 0
                    prctile_rank_thick_1(i) = sum(r_rand_thick > r_original_thick_1(i))/nperm;
                else
                    prctile_rank_thick_1(i) = sum(r_rand_thick < r_original_thick_1(i))/nperm;
                end
            end
            
        end
        %G1 spin
        for gradient_diff = 1
            cd('/Users/sofievalk/brainstat_data/surface_data/tpl-fsaverage/fsaverage5/surf/')
            SPHERE_L      = SurfStatReadSurf({'lh.sphere'});
            SPHERE_R      = SurfStatReadSurf({'lh.sphere'});
            
            gradient_mask = mean(g1_400_last(keep,:))';
            
            gradient_mask_fsa5 = nan(20484,1);
            for i = 1:200
                gradient_mask_fsa5(find(parcels400==i+1)) = gradient_mask(i);
            end
            for i = 1:200
                gradient_mask_fsa5(find(parcels400==i+1001)) = gradient_mask(i+200);
            end
            
            n_permutations = 1000;
            y_rand = spin_permutations({[gradient_mask_fsa5(1:10242)],[gradient_mask_fsa5(10243:end)]}, ...
                {SPHERE_L SPHERE_R}, ...
                n_permutations);
            
            % Merge the rotated data into single vectors
            g_rotated = squeeze([y_rand{1}(:,1,:);y_rand{2}(:,1,:)]);
            
            for i = 1:200
                gradient_mask_400(i,:) = median(g_rotated(find(parcels400==1+i),:));
            end
            for i = 1:200
                gradient_mask_400(i+200,:) = median(g_rotated(find(parcels400==1001+i),:));
            end
            
            r_original_thick_2 = zeros(12,1);
            pval_thick_spin_2 = zeros(12,1);
            prctile_rank_thick_2 = zeros(12,1);
            for i =  1:12
                [r_original_thick_2(i), pval_thick_spin_2(i)] = corr(mean(squeeze(layer_400_last(keep,i,:)))',mean(g1_400_last(keep,:))', ...
                    'rows','pairwise')
                
                r_rand_thick = corr(mean(squeeze(layer_400_last(keep,i,:)))',gradient_mask_400, ...
                    'rows','pairwise');
                
                % Compute percentile rank.
                if r_original_thick_2(i) > 0
                    prctile_rank_thick_2(i) = sum(r_rand_thick > r_original_thick_2(i))/nperm;
                else
                    prctile_rank_thick_2(i) = sum(r_rand_thick < r_original_thick_2(i))/nperm;
                end
            end
            
        end
        %G2 spin
        for gradient_diff = 1
            cd('/Users/sofievalk/brainstat_data/surface_data/tpl-fsaverage/fsaverage5/surf/')
            SPHERE_L      = SurfStatReadSurf({'lh.sphere'});
            SPHERE_R      = SurfStatReadSurf({'lh.sphere'});
            
            gradient_mask = mean(g2_400_last(keep,:))';
            
            gradient_mask_fsa5 = nan(20484,1);
            for i = 1:200
                gradient_mask_fsa5(find(parcels400==i+1)) = gradient_mask(i);
            end
            for i = 1:200
                gradient_mask_fsa5(find(parcels400==i+1001)) = gradient_mask(i+200);
            end
            
            n_permutations = 1000;
            y_rand = spin_permutations({[gradient_mask_fsa5(1:10242)],[gradient_mask_fsa5(10243:end)]}, ...
                {SPHERE_L SPHERE_R}, ...
                n_permutations);
            
            % Merge the rotated data into single vectors
            g_rotated = squeeze([y_rand{1}(:,1,:);y_rand{2}(:,1,:)]);
            
            for i = 1:200
                gradient_mask_400(i,:) = median(g_rotated(find(parcels400==1+i),:));
            end
            for i = 1:200
                gradient_mask_400(i+200,:) = median(g_rotated(find(parcels400==1001+i),:));
            end
            
            r_original_thick_3 = zeros(12,1);
            pval_thick_spin_3 = zeros(12,1);
            prctile_rank_thick_3 = zeros(12,1);
            for i =  1:12
                [r_original_thick_3(i), pval_thick_spin_3(i)] = corr(mean(squeeze(layer_400_last(keep,i,:)))',mean(g2_400_last(keep,:))', ...
                    'rows','pairwise')
                
                r_rand_thick = corr(mean(squeeze(layer_400_last(keep,i,:)))',gradient_mask_400, ...
                    'rows','pairwise');
                
                % Compute percentile rank.
                if r_original_thick_3(i) > 0
                    prctile_rank_thick_3(i) = sum(r_rand_thick > r_original_thick_3(i))/nperm;
                else
                    prctile_rank_thick_3(i) = sum(r_rand_thick < r_original_thick_3(i))/nperm;
                end
            end
            
        end
        %G3 spin
        for gradient_diff = 1
            cd('/Users/sofievalk/brainstat_data/surface_data/tpl-fsaverage/fsaverage5/surf/')
            SPHERE_L      = SurfStatReadSurf({'lh.sphere'});
            SPHERE_R      = SurfStatReadSurf({'lh.sphere'});
            
            gradient_mask = mean(g3_400_last(keep,:))';
            
            gradient_mask_fsa5 = nan(20484,1);
            for i = 1:200
                gradient_mask_fsa5(find(parcels400==i+1)) = gradient_mask(i);
            end
            for i = 1:200
                gradient_mask_fsa5(find(parcels400==i+1001)) = gradient_mask(i+200);
            end
            
            n_permutations = 1000;
            y_rand = spin_permutations({[gradient_mask_fsa5(1:10242)],[gradient_mask_fsa5(10243:end)]}, ...
                {SPHERE_L SPHERE_R}, ...
                n_permutations);
            
            % Merge the rotated data into single vectors
            g_rotated = squeeze([y_rand{1}(:,1,:);y_rand{2}(:,1,:)]);
            
            for i = 1:200
                gradient_mask_400(i,:) = median(g_rotated(find(parcels400==1+i),:));
            end
            for i = 1:200
                gradient_mask_400(i+200,:) = median(g_rotated(find(parcels400==1001+i),:));
            end
            
            r_original_thick_4 = zeros(12,1);
            pval_thick_spin_4 = zeros(12,1);
            prctile_rank_thick_4 = zeros(12,1);
            for i =  1:12
                [r_original_thick_4(i), pval_thick_spin_4(i)] = corr(mean(squeeze(layer_400_last(keep,i,:)))',mean(g3_400_last(keep,:))', ...
                    'rows','pairwise')
                
                r_rand_thick = corr(mean(squeeze(layer_400_last(keep,i,:)))',gradient_mask_400, ...
                    'rows','pairwise');
                
                % Compute percentile rank.
                if r_original_thick_4(i) > 0
                    prctile_rank_thick_4(i) = sum(r_rand_thick > r_original_thick_4(i))/nperm;
                else
                    prctile_rank_thick_4(i) = sum(r_rand_thick < r_original_thick_4(i))/nperm;
                end
            end
            
        end
        
        f=figure,
        ax = axes;
        ax.ColorOrder = cbrewer('qual','Set2',4);
        hold on
        plot(r_original_thick_1,'LineWidth',10),
        hold on
        plot(r_original_thick_2,'LineWidth',10),
        hold on
        plot(r_original_thick_3,'LineWidth',10),
        hold on
        plot(r_original_thick_4,'LineWidth',10),
        ylim([-0.5 0.5])
        exportfigbo(f,[RPATH 'F3.layervsfunc.affect.raw.png'],'png', 10)
        
        r_original_thick_1((prctile_rank_thick_1>0.05) )  = nan
        r_original_thick_2((prctile_rank_thick_2>0.05) )  = nan
        r_original_thick_3((prctile_rank_thick_3>0.05) )  = nan
        r_original_thick_4((prctile_rank_thick_4>0.05) )  = nan
        f=figure,
        ax = axes;
        ax.ColorOrder = cbrewer('qual','Set2',4);
        hold on
        plot(r_original_thick_1,'LineWidth',10),
        hold on
        plot(r_original_thick_2,'LineWidth',10),
        hold on
        plot(r_original_thick_3,'LineWidth',10),
        hold on
        plot(r_original_thick_4,'LineWidth',10),
        ylim([-0.5 0.5])
        exportfigbo(f,[RPATH 'F3.layervsfunc.affect.p.png'],'png', 10)
 
    end
    
    for change_tom = 1
        val             = Tom;
        Xn              = ggg_400_last;
        Xn(isnan(Xn))   = 0;
        Xm              = squeeze(layer_400_last(:,1,:));
        Xm(isnan(Xm))   = 0;
        keep1   = (find(strcmp(groupN,'Perspective')));
        keep2   = intersect(find(Xn(:,1) >-666), find(sum(Xn,2)~=0));
        keep2b   = intersect(find(Xm(:,1) >-666), find(sum(Xm,2)~=0));
        keep4   = find(~(strcmp(group4,'Group3')));
        keep5   = find(tpnum>0);
        keep6   = find(abs(val)<666);
        keep    = mintersect(keep1,keep2,keep2b,keep4,keep5);
        
        %eccentricity spin
        for gradient_diff = 1
            cd('/Users/sofievalk/brainstat_data/surface_data/tpl-fsaverage/fsaverage5/surf/')
            SPHERE_L      = SurfStatReadSurf({'lh.sphere'});
            SPHERE_R      = SurfStatReadSurf({'lh.sphere'});
            
            gradient_mask = mean(ggg_400_last(keep,:))';
            
            gradient_mask_fsa5 = nan(20484,1);
            for i = 1:200
                gradient_mask_fsa5(find(parcels400==i+1)) = gradient_mask(i);
            end
            for i = 1:200
                gradient_mask_fsa5(find(parcels400==i+1001)) = gradient_mask(i+200);
            end
            
            n_permutations = 1000;
            nperm = n_permutations;
            y_rand = spin_permutations({[gradient_mask_fsa5(1:10242)],[gradient_mask_fsa5(10243:end)]}, ...
                {SPHERE_L SPHERE_R}, ...
                n_permutations);
            
            % Merge the rotated data into single vectors
            g_rotated = squeeze([y_rand{1}(:,1,:);y_rand{2}(:,1,:)]);
            
            for i = 1:200
                gradient_mask_400(i,:) = median(g_rotated(find(parcels400==1+i),:));
            end
            for i = 1:200
                gradient_mask_400(i+200,:) = median(g_rotated(find(parcels400==1001+i),:));
            end
            
            r_original_thick_1 = zeros(12,1);
            pval_thick_spin_1 = zeros(12,1);
            prctile_rank_thick_1 = zeros(12,1);
            for i =  1:12
                [r_original_thick_1(i), pval_thick_spin_1(i)] = corr(mean(squeeze(layer_400_last(keep,i,:)))',mean(ggg_400_last(keep,:))', ...
                    'rows','pairwise')
                
                r_rand_thick = corr(mean(squeeze(layer_400_last(keep,i,:)))',gradient_mask_400, ...
                    'rows','pairwise');
                
                if r_original_thick_1(i) > 0
                    prctile_rank_thick_1(i) = sum(r_rand_thick > r_original_thick_1(i))/nperm;
                else
                    prctile_rank_thick_1(i) = sum(r_rand_thick < r_original_thick_1(i))/nperm;
                end
            end
            
        end
        %G1 spin
        for gradient_diff = 1
            cd('/Users/sofievalk/brainstat_data/surface_data/tpl-fsaverage/fsaverage5/surf/')
            SPHERE_L      = SurfStatReadSurf({'lh.sphere'});
            SPHERE_R      = SurfStatReadSurf({'lh.sphere'});
            
            gradient_mask = mean(g1_400_last(keep,:))';
            
            gradient_mask_fsa5 = nan(20484,1);
            for i = 1:200
                gradient_mask_fsa5(find(parcels400==i+1)) = gradient_mask(i);
            end
            for i = 1:200
                gradient_mask_fsa5(find(parcels400==i+1001)) = gradient_mask(i+200);
            end
            
            n_permutations = 1000;
            y_rand = spin_permutations({[gradient_mask_fsa5(1:10242)],[gradient_mask_fsa5(10243:end)]}, ...
                {SPHERE_L SPHERE_R}, ...
                n_permutations);
            
            % Merge the rotated data into single vectors
            g_rotated = squeeze([y_rand{1}(:,1,:);y_rand{2}(:,1,:)]);
            
            for i = 1:200
                gradient_mask_400(i,:) = median(g_rotated(find(parcels400==1+i),:));
            end
            for i = 1:200
                gradient_mask_400(i+200,:) = median(g_rotated(find(parcels400==1001+i),:));
            end
            
            r_original_thick_2 = zeros(12,1);
            pval_thick_spin_2 = zeros(12,1);
            prctile_rank_thick_2 = zeros(12,1);
            for i =  1:12
                [r_original_thick_2(i), pval_thick_spin_2(i)] = corr(mean(squeeze(layer_400_last(keep,i,:)))',mean(g1_400_last(keep,:))', ...
                    'rows','pairwise')
                
                r_rand_thick = corr(mean(squeeze(layer_400_last(keep,i,:)))',gradient_mask_400, ...
                    'rows','pairwise');
                
                % Compute percentile rank.
                if r_original_thick_2(i) > 0
                    prctile_rank_thick_2(i) = sum(r_rand_thick > r_original_thick_2(i))/nperm;
                else
                    prctile_rank_thick_2(i) = sum(r_rand_thick < r_original_thick_2(i))/nperm;
                end
            end
            
        end
        %G2 spin
        for gradient_diff = 1
            cd('/Users/sofievalk/brainstat_data/surface_data/tpl-fsaverage/fsaverage5/surf/')
            SPHERE_L      = SurfStatReadSurf({'lh.sphere'});
            SPHERE_R      = SurfStatReadSurf({'lh.sphere'});
            
            gradient_mask = mean(g2_400_last(keep,:))';
            
            gradient_mask_fsa5 = nan(20484,1);
            for i = 1:200
                gradient_mask_fsa5(find(parcels400==i+1)) = gradient_mask(i);
            end
            for i = 1:200
                gradient_mask_fsa5(find(parcels400==i+1001)) = gradient_mask(i+200);
            end
            
            n_permutations = 1000;
            y_rand = spin_permutations({[gradient_mask_fsa5(1:10242)],[gradient_mask_fsa5(10243:end)]}, ...
                {SPHERE_L SPHERE_R}, ...
                n_permutations);
            
            % Merge the rotated data into single vectors
            g_rotated = squeeze([y_rand{1}(:,1,:);y_rand{2}(:,1,:)]);
            
            for i = 1:200
                gradient_mask_400(i,:) = median(g_rotated(find(parcels400==1+i),:));
            end
            for i = 1:200
                gradient_mask_400(i+200,:) = median(g_rotated(find(parcels400==1001+i),:));
            end
            
            r_original_thick_3 = zeros(12,1);
            pval_thick_spin_3 = zeros(12,1);
            prctile_rank_thick_3 = zeros(12,1);
            for i =  1:12
                [r_original_thick_3(i), pval_thick_spin_3(i)] = corr(mean(squeeze(layer_400_last(keep,i,:)))',mean(g2_400_last(keep,:))', ...
                    'rows','pairwise')
                
                r_rand_thick = corr(mean(squeeze(layer_400_last(keep,i,:)))',gradient_mask_400, ...
                    'rows','pairwise');
                
                % Compute percentile rank.
                if r_original_thick_3(i) > 0
                    prctile_rank_thick_3(i) = sum(r_rand_thick > r_original_thick_3(i))/nperm;
                else
                    prctile_rank_thick_3(i) = sum(r_rand_thick < r_original_thick_3(i))/nperm;
                end
            end
            
        end
        %G3 spin
        for gradient_diff = 1
            cd('/Users/sofievalk/brainstat_data/surface_data/tpl-fsaverage/fsaverage5/surf/')
            SPHERE_L      = SurfStatReadSurf({'lh.sphere'});
            SPHERE_R      = SurfStatReadSurf({'lh.sphere'});
            
            gradient_mask = mean(g3_400_last(keep,:))';
            
            gradient_mask_fsa5 = nan(20484,1);
            for i = 1:200
                gradient_mask_fsa5(find(parcels400==i+1)) = gradient_mask(i);
            end
            for i = 1:200
                gradient_mask_fsa5(find(parcels400==i+1001)) = gradient_mask(i+200);
            end
            
            n_permutations = 1000;
            y_rand = spin_permutations({[gradient_mask_fsa5(1:10242)],[gradient_mask_fsa5(10243:end)]}, ...
                {SPHERE_L SPHERE_R}, ...
                n_permutations);
            
            % Merge the rotated data into single vectors
            g_rotated = squeeze([y_rand{1}(:,1,:);y_rand{2}(:,1,:)]);
            
            for i = 1:200
                gradient_mask_400(i,:) = median(g_rotated(find(parcels400==1+i),:));
            end
            for i = 1:200
                gradient_mask_400(i+200,:) = median(g_rotated(find(parcels400==1001+i),:));
            end
            
            r_original_thick_4 = zeros(12,1);
            pval_thick_spin_4 = zeros(12,1);
            prctile_rank_thick_4 = zeros(12,1);
            for i =  1:12
                [r_original_thick_4(i), pval_thick_spin_4(i)] = corr(mean(squeeze(layer_400_last(keep,i,:)))',mean(g3_400_last(keep,:))', ...
                    'rows','pairwise')
                
                r_rand_thick = corr(mean(squeeze(layer_400_last(keep,i,:)))',gradient_mask_400, ...
                    'rows','pairwise');
                
                % Compute percentile rank.
                if r_original_thick_4(i) > 0
                    prctile_rank_thick_4(i) = sum(r_rand_thick > r_original_thick_4(i))/nperm;
                else
                    prctile_rank_thick_4(i) = sum(r_rand_thick < r_original_thick_4(i))/nperm;
                end
            end
            
        end
        
        f=figure,
        ax = axes;
        ax.ColorOrder = cbrewer('qual','Set2',4);
        hold on
        plot(r_original_thick_1,'LineWidth',10),
        hold on
        plot(r_original_thick_2,'LineWidth',10),
        hold on
        plot(r_original_thick_3,'LineWidth',10),
        hold on
        plot(r_original_thick_4,'LineWidth',10),
        ylim([-0.5 0.5])
        exportfigbo(f,[RPATH 'F3.layervsfunc.perspective.raw.png'],'png', 10)
        
        r_original_thick_1((prctile_rank_thick_1>0.05) )  = nan
        r_original_thick_2((prctile_rank_thick_2>0.05) )  = nan
        r_original_thick_3((prctile_rank_thick_3>0.05) )  = nan
        r_original_thick_4((prctile_rank_thick_4>0.05) )  = nan
        f=figure,
        ax = axes;
        ax.ColorOrder = cbrewer('qual','Set2',4);
        hold on
        plot(r_original_thick_1,'LineWidth',10),
        hold on
        plot(r_original_thick_2,'LineWidth',10),
        hold on
        plot(r_original_thick_3,'LineWidth',10),
        hold on
        plot(r_original_thick_4,'LineWidth',10),
        ylim([-0.5 0.5])
        exportfigbo(f,[RPATH 'F3.layervsfunc.perspective.p.png'],'png', 10)
   
    end
        
end

for f3_subs =  1
    
    for all_layers = 1
        Xn            = squeeze(layer_400_last(:,1,:));
        Xn(isnan(Xn)) = 0;
        Xn(isinf(Xn)) = 0;
        Xm            = ggg_400_last;
        Xm(isnan(Xm)) = 0;
        Xm(isinf(Xm)) = 0;
        keep1   = union(find(strcmp(groupN,'Affect')),    find(strcmp(groupN,'Perspective')));
        keep1   = union(keep1,find(strcmp(groupN,'Presence')))
        keep1   = union(keep1,find(strcmp(groupN,'Control1')))
        keep2   = intersect(find(abs(mean(Xn,2)) <666), find(sum(Xn,2)~=0));
        keep3   = intersect(find(~strcmp(group4,'Group3')),find(tpnum>0));
        keep4   = intersect(find(abs(mean(Xm,2)) <666), find(sum(Xm,2)~=0));
        keep    = mintersect(keep1, keep2,keep3,keep4);
        
        ik      = id(keep,:);
        ink     = idnum(keep);
        tk      = tpnum(keep,:);
        tnk     = tp(keep,:);
        gk      = group(keep);
        g4k     = group4(keep);
        gNk     = groupN(keep);
        ak      = age(keep);
        sk      = cellstr(sex(keep));
        
        A       = term(ak);
        S       = term(sk);
        G       = term(gk);
        G4      = term(g4k);
        GN      = term(gNk);
        Sub     = term(ink);
        Tn      = term(tnk);
        
        M        = 1 + A + S + GN + random(Sub) + I;
        
        for i = 1:5
            for k=1:12
                Xx =  squeeze(mean(layer_400_last(:,k,:),2));
                slm      = SurfStatLinMod(mean(Xx(keep,maps_ns(:,i)>0),2),M);
                slm      = SurfStatT(slm,(GN.Presence-GN.Perspective));
                slm2 = slm.t;
                p_dept(k,1,i)  = 1 - tcdf(slm.t,slm.df);
                t_dept(k,1,i) = slm.t;
                Xx =  squeeze(mean(layer_400_last(:,k,:),2));
                slm      = SurfStatLinMod(mean(Xx(keep,maps_ns(:,i)>0),2),M);
                slm      = SurfStatT(slm,(GN.Presence-GN.Affect));
                slm2 = slm.t;
                p_dept(k,2,i)  = 1 - tcdf(slm.t,slm.df);
                t_dept(k,2,i) = slm.t;
                Xx =  squeeze(mean(layer_400_last(:,k,:),2));
                slm      = SurfStatLinMod(mean(Xx(keep,maps_ns(:,i)>0),2),M);
                slm      = SurfStatT(slm,(GN.Perspective-GN.Affect));
                slm2 = slm.t;
                p_dept(k,3,i)  = 1 - tcdf(slm.t,slm.df);
                t_dept(k,3,i) = slm.t;
            end
        end
        
        csvwrite([RPATH 'network', num2str(i),'layer_myelin_presence-perspective_all.csv'], [...
            t_dept(:,1,1), p_dept(:,1,1),...
            t_dept(:,1,2), p_dept(:,1,2),...
            t_dept(:,1,3), p_dept(:,1,3),...
            t_dept(:,1,4), p_dept(:,1,4),...
            t_dept(:,1,5), p_dept(:,1,5)])
        
        csvwrite([RPATH 'network', num2str(i),'layer_myelin_presence-affect_all.csv'], [...
            t_dept(:,2,1), p_dept(:,2,1),...
            t_dept(:,2,2), p_dept(:,2,2),...
            t_dept(:,2,3), p_dept(:,2,3),...
            t_dept(:,2,4), p_dept(:,2,4),...
            t_dept(:,2,5), p_dept(:,2,5)])
        
        csvwrite([RPATH 'network', num2str(i),'layer_myelin_perspective-affect_all.csv'], [...
            t_dept(:,3,1), p_dept(:,3,1),...
            t_dept(:,3,2), p_dept(:,3,2),...
            t_dept(:,3,3), p_dept(:,3,3),...
            t_dept(:,3,4), p_dept(:,3,4),...
            t_dept(:,3,5), p_dept(:,3,5)])
        
        for key_descriptives = 1
            for k = 1:12
                Xn = squeeze(mean(layer_400_last(:,k,:),2));
                keep_presence = (find(strcmp(groupN(keep),'Presence')))
                keep_affect = (find(strcmp(groupN(keep),'Affect')))
                keep_perspective = (find(strcmp(groupN(keep),'Perspective')))
                keep_control   = (find(strcmp(groupN(keep),'Control1')))
                
                for i = 1:5
                    x = mean(Xn(keep(keep_control),maps_ns(:,i)>0),2);
                    meansdata(k,4,i) = mean(x);
                    stddata(k,4,i) = std(x);
                    SEM = std(x)/sqrt(length(x)); % Standard Error
                    ts = tinv([0.025  0.975],length(x)-1);      % T-Score
                    cidata = mean(x) + ts*SEM;
                    cidatamin(k,4,i) = cidata(1);
                    cidatamax(k,4,i) = cidata(2);
                end
                for i = 1:5
                    x = mean(Xn(keep(keep_affect),maps_ns(:,i)>0),2);
                    meansdata(k,2,i) = mean(x);
                    stddata(k,2,i) = std(x);
                    SEM = std(x)/sqrt(length(x)); % Standard Error
                    ts = tinv([0.025  0.975],length(x)-1);      % T-Score
                    cidata = mean(x) + ts*SEM;
                    cidatamin(k,2,i) = cidata(1);
                    cidatamax(k,2,i) = cidata(2);
                end
                for i = 1:5
                    x = mean(Xn(keep(keep_perspective),maps_ns(:,i)>0),2);
                    meansdata(k,3,i) = mean(x);
                    stddata(k,3,i) = std(x);
                    SEM = std(x)/sqrt(length(x)); % Standard Error
                    ts = tinv([0.025  0.975],length(x)-1);      % T-Score
                    cidata = mean(x) + ts*SEM;
                    cidatamin(k,3,i) = cidata(1);
                    cidatamax(k,3,i) = cidata(2);
                end
                for i = 1:5
                    x = mean(Xn(keep(keep_presence),maps_ns(:,i)>0),2);
                    meansdata(k,1,i) = mean(x);
                    stddata(k,1,i) = std(x);
                    SEM = std(x)/sqrt(length(x)); % Standard Error
                    ts = tinv([0.025  0.975],length(x)-1);      % T-Score
                    cidata = mean(x) + ts*SEM;
                    cidatamin(k,1,i) = cidata(1);
                    cidatamax(k,1,i) = cidata(2);
                end
            end
        end
        
        csvwrite([RPATH 'stats_layer_myelin_controlall.csv'], [...
            meansdata(:,4,1), stddata(:,4,1),cidatamin(:,4,1),cidatamax(:,4,1),...
            meansdata(:,4,2), stddata(:,4,2),cidatamin(:,4,2),cidatamax(:,4,2),...
            meansdata(:,4,3), stddata(:,4,3),cidatamin(:,4,3),cidatamax(:,4,3),...
            meansdata(:,4,4), stddata(:,4,4),cidatamin(:,4,4),cidatamax(:,4,4),...
            meansdata(:,4,5), stddata(:,4,5),cidatamin(:,4,5),cidatamax(:,4,5)])
        
        for j  = 1:4:12
            xxx      = squeeze(mean(layer_400_last(:,j:(j+3),:),2));
            slm      = SurfStatLinMod(xxx(keep,:),M);
            slm      = SurfStatT(slm,(GN.Presence));
            p        = 2*(1 - tcdf(abs(slm.t),slm.df));
            
            hp = fdr_bh(p(1:400),0.05)
            h  = hp;
            
            tmp = zeros(20484,1);
            for i = 1:200
                tmp(find(parcels400==i+1)) = slm.t(i).*h(i);
            end
            for i = 1:200
                tmp(find(parcels400==i+1001)) =slm.t(i+200).*h(i+200);
            end
            
            f=figure,
            BoSurfStatViewData(tmp,SN,'')
            BoSurfStatColLim([-4 4])
            colormap(flipud(cbrewer('div','RdBu',11)))
            exportfigbo(f,[RPATH 'F3s.presence.layer', num2str(j) ,'.fdr.png'],'png', 10)
            close(f)
            
            tmp = zeros(20484,1);
            for i = 1:200
                tmp(find(parcels400==i+1)) = slm.t(i);
            end
            for i = 1:200
                tmp(find(parcels400==i+1001)) =slm.t(i+200);
            end
            
            f=figure,
            BoSurfStatViewData(tmp,SN,'')
            BoSurfStatColLim([-4 4])
            colormap(flipud(cbrewer('div','RdBu',11)))
            exportfigbo(f,[RPATH 'F3s.presence.layer', num2str(j) ,'.t.png'],'png', 10)
            close(f)
        end
        
        for j  = 1:4:12
            xxx      = squeeze(mean(layer_400_last(:,j:(j+3),:),2));
            slm      = SurfStatLinMod(xxx(keep,:),M);
            slm      = SurfStatT(slm,(GN.Perspective));
            p        = 2*(1 - tcdf(abs(slm.t),slm.df));
            
            hp = fdr_bh(p(1:400),0.05)
            h  = hp;
            
            tmp = zeros(20484,1);
            for i = 1:200
                tmp(find(parcels400==i+1)) = slm.t(i).*h(i);
            end
            for i = 1:200
                tmp(find(parcels400==i+1001)) =slm.t(i+200).*h(i+200);
            end
            
            f=figure,
            BoSurfStatViewData(tmp,SN,'')
            BoSurfStatColLim([-4 4])
            colormap(flipud(cbrewer('div','RdBu',11)))
            exportfigbo(f,[RPATH 'F3s.perspective.layer', num2str(j) ,'.fdr.png'],'png', 10)
            close(f)
            
            tmp = zeros(20484,1);
            for i = 1:200
                tmp(find(parcels400==i+1)) = slm.t(i);
            end
            for i = 1:200
                tmp(find(parcels400==i+1001)) =slm.t(i+200);
            end
            
            f=figure,
            BoSurfStatViewData(tmp,SN,'')
            BoSurfStatColLim([-4 4])
            colormap(flipud(cbrewer('div','RdBu',11)))
            exportfigbo(f,[RPATH 'F3s.perspective.layer', num2str(j) ,'.t.png'],'png', 10)
            close(f)
            
        end
        
        for j  = 1:4:12
            xxx      = squeeze(mean(layer_400_last(:,j:(j+3),:),2));
            slm      = SurfStatLinMod(xxx(keep,:),M);
            slm      = SurfStatT(slm,(GN.Affect));
            p        = 2*(1 - tcdf(abs(slm.t),slm.df));
            
            hp = fdr_bh(p(1:400),0.05)
            h  = hp;
            
            tmp = zeros(20484,1);
            for i = 1:200
                tmp(find(parcels400==i+1)) = slm.t(i).*h(i);
            end
            for i = 1:200
                tmp(find(parcels400==i+1001)) =slm.t(i+200).*h(i+200);
            end
            
            f=figure,
            BoSurfStatViewData(tmp,SN,'')
            BoSurfStatColLim([-4 4])
            colormap(flipud(cbrewer('div','RdBu',11)))
            exportfigbo(f,[RPATH 'F3s.affect.layer', num2str(j) ,'.fdr.png'],'png', 10)
            close(f)
            
            tmp = zeros(20484,1);
            for i = 1:200
                tmp(find(parcels400==i+1)) = slm.t(i);
            end
            for i = 1:200
                tmp(find(parcels400==i+1001)) =slm.t(i+200);
            end
            
            f=figure,
            BoSurfStatViewData(tmp,SN,'')
            BoSurfStatColLim([-4 4])
            colormap(flipud(cbrewer('div','RdBu',11)))
            exportfigbo(f,[RPATH 'F3s.affect.layer', num2str(j) ,'.t.png'],'png', 10)
            close(f)
            
        end
        
        for j  = 1:4:12
            xxx      = squeeze(mean(layer_400_last(:,j:(j+3),:),2));
            slm      = SurfStatLinMod(xxx(keep,:),M);
            slm      = SurfStatT(slm,(GN.Control1));
            p        = 2*(1 - tcdf(abs(slm.t),slm.df));
            
            hp = fdr_bh(p(1:400),0.05)
            h  = hp;
            
            tmp = zeros(20484,1);
            for i = 1:200
                tmp(find(parcels400==i+1)) = slm.t(i).*h(i);
            end
            for i = 1:200
                tmp(find(parcels400==i+1001)) =slm.t(i+200).*h(i+200);
            end
            
            f=figure,
            BoSurfStatViewData(tmp,SN,'')
            BoSurfStatColLim([-4 4])
            colormap(flipud(cbrewer('div','RdBu',11)))
            exportfigbo(f,[RPATH 'F3s.control.layer', num2str(j) ,'.fdr.png'],'png', 10)
            close(f)
            
            tmp = zeros(20484,1);
            for i = 1:200
                tmp(find(parcels400==i+1)) = slm.t(i);
            end
            for i = 1:200
                tmp(find(parcels400==i+1001)) =slm.t(i+200);
            end
            
            f=figure,
            BoSurfStatViewData(tmp,SN,'')
            BoSurfStatColLim([-4 4])
            colormap(flipud(cbrewer('div','RdBu',11)))
            exportfigbo(f,[RPATH 'F3s.control.layer', num2str(j) ,'.t.png'],'png', 10)
            close(f)
            
        end
        
        for j  = 1:4:12
            xxx      = squeeze(mean(layer_400_last(:,j:(j+3),:),2));
            slm      = SurfStatLinMod(xxx(keep,:),M);
            slm      = SurfStatT(slm,(GN.Perspective)-(GN.Affect));
            p        = 2*(1 - tcdf(abs(slm.t),slm.df));
            
            hp = fdr_bh(p(1:400),0.05)
            h  = hp;
            
            tmp = zeros(20484,1);
            for i = 1:200
                tmp(find(parcels400==i+1)) = slm.t(i).*h(i);
            end
            for i = 1:200
                tmp(find(parcels400==i+1001)) =slm.t(i+200).*h(i+200);
            end
            
            f=figure,
            BoSurfStatViewData(tmp,SN,'')
            BoSurfStatColLim([-4 4])
            colormap(flipud(cbrewer('div','RdBu',11)))
            exportfigbo(f,[RPATH 'F3s.pa.layer', num2str(j) ,'.fdr.png'],'png', 10)
            close(f)
            
            h = p<0.01
            tmp = zeros(20484,1);
            for i = 1:200
                tmp(find(parcels400==i+1)) = slm.t(i).*h(i);
            end
            for i = 1:200
                tmp(find(parcels400==i+1001)) =slm.t(i+200).*h(i+200);
            end
            
            f=figure,
            BoSurfStatViewData(tmp,SN,'')
            BoSurfStatColLim([-4 4])
            colormap(flipud(cbrewer('div','RdBu',11)))
            exportfigbo(f,[RPATH 'F3s.pa.layer', num2str(j) ,'.t.png'],'png', 10)
            close(f)
            
        end
        
        for j  = 1:4:12
            xxx      = squeeze(mean(layer_400_last(:,j:(j+3),:),2));
            slm      = SurfStatLinMod(xxx(keep,:),M);
            slm      = SurfStatT(slm,(GN.Presence)-(GN.Affect));
            p        = 2*(1 - tcdf(abs(slm.t),slm.df));
            
            hp = fdr_bh(p(1:400),0.05)
            h  = hp;
            
            tmp = zeros(20484,1);
            for i = 1:200
                tmp(find(parcels400==i+1)) = slm.t(i).*h(i);
            end
            for i = 1:200
                tmp(find(parcels400==i+1001)) =slm.t(i+200).*h(i+200);
            end
            
            f=figure,
            BoSurfStatViewData(tmp,SN,'')
            BoSurfStatColLim([-4 4])
            colormap(flipud(cbrewer('div','RdBu',11)))
            exportfigbo(f,[RPATH 'F3s.pra.layer', num2str(j) ,'.fdr.png'],'png', 10)
            close(f)
            
            h = p<0.01
            tmp = zeros(20484,1);
            for i = 1:200
                tmp(find(parcels400==i+1)) = slm.t(i).*h(i);
            end
            for i = 1:200
                tmp(find(parcels400==i+1001)) =slm.t(i+200).*h(i+200);
            end
            
            f=figure,
            BoSurfStatViewData(tmp,SN,'')
            BoSurfStatColLim([-4 4])
            colormap(flipud(cbrewer('div','RdBu',11)))
            exportfigbo(f,[RPATH 'F3s.pra.layer', num2str(j) ,'.t.png'],'png', 10)
            close(f)
            
        end
        
        for j  = 1:4:12
            xxx      = squeeze(mean(layer_400_last(:,j:(j+3),:),2));
            slm      = SurfStatLinMod(xxx(keep,:),M);
            slm      = SurfStatT(slm,(GN.Presence)-(GN.Perspective));
            p        = 2*(1 - tcdf(abs(slm.t),slm.df));
            
            hp = fdr_bh(p(1:400),0.05)
            h  = hp;
            
            tmp = zeros(20484,1);
            for i = 1:200
                tmp(find(parcels400==i+1)) = slm.t(i).*h(i);
            end
            for i = 1:200
                tmp(find(parcels400==i+1001)) =slm.t(i+200).*h(i+200);
            end
            
            f=figure,
            BoSurfStatViewData(tmp,SN,'')
            BoSurfStatColLim([-4 4])
            colormap(flipud(cbrewer('div','RdBu',11)))
            exportfigbo(f,[RPATH 'F3s.prp.layer', num2str(j) ,'.fdr.png'],'png', 10)
            close(f)
            
            h = p<0.01
            tmp = zeros(20484,1);
            for i = 1:200
                tmp(find(parcels400==i+1)) = slm.t(i).*h(i);
            end
            for i = 1:200
                tmp(find(parcels400==i+1001)) =slm.t(i+200).*h(i+200);
            end
            
            f=figure,
            BoSurfStatViewData(tmp,SN,'')
            BoSurfStatColLim([-4 4])
            colormap(flipud(cbrewer('div','RdBu',11)))
            exportfigbo(f,[RPATH 'F3s.prp.layer', num2str(j) ,'.t.png'],'png', 10)
            close(f)
            
        end
        
    end
     
    for tc1_layers = 1
        Xn            = squeeze(layer_400_last(:,1,:));
        Xn(isnan(Xn)) = 0;
        Xn(isinf(Xn)) = 0;
        Xm            = ggg_400_last;
        Xm(isnan(Xm)) = 0;
        Xm(isinf(Xm)) = 0;
        keep1   = union(find(strcmp(groupN,'Affect')),    find(strcmp(groupN,'Perspective')));
        keep1   = union(keep1,find(strcmp(groupN,'Presence')))
        keep1   = union(keep1,find(strcmp(groupN,'Control1')))
        keep2   = intersect(find(abs(mean(Xn,2)) <666), find(sum(Xn,2)~=0));
        keep3   = intersect(find(strcmp(group4,'Group_1')),find(tpnum>0));  
        keep4   = intersect(find(abs(mean(Xm,2)) <666), find(sum(Xm,2)~=0));
        keep    = mintersect(keep1, keep2,keep3,keep4);
        
        ik      = id(keep,:);
        ink     = idnum(keep);
        tk      = tpnum(keep,:);
        tnk     = tp(keep,:);
        gk      = group(keep);
        g4k     = group4(keep);
        gNk     = groupN(keep);
        ak      = age(keep);
        sk      = cellstr(sex(keep));   
        
        A       = term(ak);
        S       = term(sk);
        G       = term(gk);
        G4      = term(g4k);
        GN      = term(gNk);
        Sub     = term(ink);
        Tn      = term(tnk);
        
        M        = 1 + A + S + GN + random(Sub) + I;
 
        for i = 1:5
                for k=1:12
                    Xx =  squeeze(mean(layer_400_last(:,k,:),2));
                    slm      = SurfStatLinMod(mean(Xx(keep,maps_ns(:,i)>0),2),M);
                    slm      = SurfStatT(slm,(GN.Presence-GN.Perspective));
                    slm2 = slm.t;
                    p_dept(k,1,i)  = 2*(1 - tcdf(abs(slm.t),slm.df));
                    t_dept(k,1,i) = slm.t;
                    Xx =  squeeze(mean(layer_400_last(:,k,:),2));
                    slm      = SurfStatLinMod(mean(Xx(keep,maps_ns(:,i)>0),2),M);
                    slm      = SurfStatT(slm,(GN.Presence-GN.Affect));
                    slm2 = slm.t;
                    p_dept(k,2,i)  = 2*(1 - tcdf(abs(slm.t),slm.df));
                    t_dept(k,2,i) = slm.t;
                    Xx =  squeeze(mean(layer_400_last(:,k,:),2));
                    slm      = SurfStatLinMod(mean(Xx(keep,maps_ns(:,i)>0),2),M);
                    slm      = SurfStatT(slm,(GN.Perspective-GN.Affect));
                    slm2 = slm.t;
                    p_dept(k,3,i)  = 2*(1 - tcdf(abs(slm.t),slm.df));
                    t_dept(k,3,i) = slm.t;
                end
        end
        
          csvwrite([RPATH 'network', num2str(i),'layer_myelin_presence-perspective_TC1.csv'], [...
                    t_dept(:,1,1), p_dept(:,1,1),...
                    t_dept(:,1,2), p_dept(:,1,2),...
                    t_dept(:,1,3), p_dept(:,1,3),...
                    t_dept(:,1,4), p_dept(:,1,4),...
                    t_dept(:,1,5), p_dept(:,1,5)]) 
                
           csvwrite([RPATH 'network', num2str(i),'layer_myelin_presence-affect_TC1.csv'], [...
                    t_dept(:,2,1), p_dept(:,2,1),...
                    t_dept(:,2,2), p_dept(:,2,2),...
                    t_dept(:,2,3), p_dept(:,2,3),...
                    t_dept(:,2,4), p_dept(:,2,4),...
                    t_dept(:,2,5), p_dept(:,2,5)])     
                
                csvwrite([RPATH 'network', num2str(i),'layer_myelin_perspective-affect_TC1.csv'], [...
                    t_dept(:,3,1), p_dept(:,3,1),...
                    t_dept(:,3,2), p_dept(:,3,2),...
                    t_dept(:,3,3), p_dept(:,3,3),...
                    t_dept(:,3,4), p_dept(:,3,4),...
                    t_dept(:,3,5), p_dept(:,3,5)])  
            
    end
    
    for tc2_layers = 1
        Xn            = squeeze(layer_400_last(:,1,:));
        Xn(isnan(Xn)) = 0;
        Xn(isinf(Xn)) = 0;
        Xm            = ggg_400_last;
        Xm(isnan(Xm)) = 0;
        Xm(isinf(Xm)) = 0;
        keep1   = union(find(strcmp(groupN,'Affect')),    find(strcmp(groupN,'Perspective')));
        keep1   = union(keep1,find(strcmp(groupN,'Presence')))
        keep1   = union(keep1,find(strcmp(groupN,'Control1')))
        keep2   = intersect(find(abs(mean(Xn,2)) <666), find(sum(Xn,2)~=0));
        keep3   = intersect(find(strcmp(group4,'Group_2')),find(tpnum>0));  
        keep4   = intersect(find(abs(mean(Xm,2)) <666), find(sum(Xm,2)~=0));
        keep    = mintersect(keep1, keep2,keep3,keep4);
        
        ik      = id(keep,:);
        ink     = idnum(keep);
        tk      = tpnum(keep,:);
        tnk     = tp(keep,:);
        gk      = group(keep);
        g4k     = group4(keep);
        gNk     = groupN(keep);
        ak      = age(keep);
        sk      = cellstr(sex(keep));   
        
        A       = term(ak);
        S       = term(sk);
        G       = term(gk);
        G4      = term(g4k);
        GN      = term(gNk);
        Sub     = term(ink);
        Tn      = term(tnk);
        
        M        = 1 + A + S + GN + random(Sub) + I;
 
        for i = 1:5
                for k=1:12
                    Xx =  squeeze(mean(layer_400_last(:,k,:),2));
                    slm      = SurfStatLinMod(mean(Xx(keep,maps_ns(:,i)>0),2),M);
                    slm      = SurfStatT(slm,(GN.Presence-GN.Perspective));
                    slm2 = slm.t;
                    p_dept(k,1,i)  = 2*(1 - tcdf(abs(slm.t),slm.df));
                    t_dept(k,1,i) = slm.t;
                    Xx =  squeeze(mean(layer_400_last(:,k,:),2));
                    slm      = SurfStatLinMod(mean(Xx(keep,maps_ns(:,i)>0),2),M);
                    slm      = SurfStatT(slm,(GN.Presence-GN.Affect));
                    slm2 = slm.t;
                    p_dept(k,2,i)  = 2*(1 - tcdf(abs(slm.t),slm.df));
                    t_dept(k,2,i) = slm.t;
                    Xx =  squeeze(mean(layer_400_last(:,k,:),2));
                    slm      = SurfStatLinMod(mean(Xx(keep,maps_ns(:,i)>0),2),M);
                    slm      = SurfStatT(slm,(GN.Perspective-GN.Affect));
                    slm2 = slm.t;
                    p_dept(k,3,i)  = 2*(1 - tcdf(abs(slm.t),slm.df));
                    t_dept(k,3,i) = slm.t;
                end
        end
        
          csvwrite([RPATH 'network', num2str(i),'layer_myelin_presence-perspective_TC2.csv'], [...
                    t_dept(:,1,1), p_dept(:,1,1),...
                    t_dept(:,1,2), p_dept(:,1,2),...
                    t_dept(:,1,3), p_dept(:,1,3),...
                    t_dept(:,1,4), p_dept(:,1,4),...
                    t_dept(:,1,5), p_dept(:,1,5)]) 
                
           csvwrite([RPATH 'network', num2str(i),'layer_myelin_presence-affect_TC2.csv'], [...
                    t_dept(:,2,1), p_dept(:,2,1),...
                    t_dept(:,2,2), p_dept(:,2,2),...
                    t_dept(:,2,3), p_dept(:,2,3),...
                    t_dept(:,2,4), p_dept(:,2,4),...
                    t_dept(:,2,5), p_dept(:,2,5)])     
                
                csvwrite([RPATH 'network', num2str(i),'layer_myelin_perspective-affect_TC2.csv'], [...
                    t_dept(:,3,1), p_dept(:,3,1),...
                    t_dept(:,3,2), p_dept(:,3,2),...
                    t_dept(:,3,3), p_dept(:,3,3),...
                    t_dept(:,3,4), p_dept(:,3,4),...
                    t_dept(:,3,5), p_dept(:,3,5)])  
            
    end
    
    for t0_t1_layers = 1
        Xn            = squeeze(layer_400_last(:,1,:));
        Xn(isnan(Xn)) = 0;
        Xn(isinf(Xn)) = 0;
        Xm            = ggg_400_last;
        Xm(isnan(Xm)) = 0;
        Xm(isinf(Xm)) = 0;
        keep1   = union(find(strcmp(group4,'Group3')),    find(strcmp(groupN,'Perspective')));
        keep1   = union(keep1,find(strcmp(groupN,'Presence')))
        keep1   = union(keep1,find(strcmp(groupN,'Control1')))
        keep2   = intersect(find(abs(mean(Xn,2)) <666), find(sum(Xn,2)~=0));
        keep3   = find(tpnum==1);  
        keep4   = intersect(find(abs(mean(Xm,2)) <666), find(sum(Xm,2)~=0));
        keep    = mintersect(keep1, keep2,keep3,keep4);
        
        ik      = id(keep,:);
        ink     = idnum(keep);
        tk      = tpnum(keep,:);
        tnk     = tp(keep,:);
        gk      = group(keep);
        g4k     = group4(keep);
        gNk     = groupN(keep);
        ak      = age(keep);
        sk      = cellstr(sex(keep));   
        
        A       = term(ak);
        S       = term(sk);
        G       = term(gk);
        G4      = term(g4k);
        GN      = term(gNk);
        Sub     = term(ink);
        Tn      = term(tnk);
        
        M        = 1 + A + S + GN + random(Sub) + I;
 
        for i = 1:5
                for k=1:12
                    Xx =  squeeze(mean(layer_400_last(:,k,:),2));
                    slm      = SurfStatLinMod(mean(Xx(keep,maps_ns(:,i)>0),2),M);
                    slm      = SurfStatT(slm,(GN.Presence-GN.Affect));
                    slm2 = slm.t;
                    p_dept(k,1,i)  = 2*(1 - tcdf(abs(slm.t),slm.df));
                    t_dept(k,1,i) = slm.t;
                    Xx =  squeeze(mean(layer_400_last(:,k,:),2));
                    slm      = SurfStatLinMod(mean(Xx(keep,maps_ns(:,i)>0),2),M);
                    slm      = SurfStatT(slm,(GN.Presence-GN.Control));
                    slm2 = slm.t;
                    p_dept(k,2,i)  = 2*(1 - tcdf(abs(slm.t),slm.df));
                    t_dept(k,2,i) = slm.t;
                    Xx =  squeeze(mean(layer_400_last(:,k,:),2));
                    slm      = SurfStatLinMod(mean(Xx(keep,maps_ns(:,i)>0),2),M);
                    slm      = SurfStatT(slm,(GN.Affect-GN.Control1));
                    slm2 = slm.t;
                    p_dept(k,3,i)  = 2*(1 - tcdf(abs(slm.t),slm.df));
                    t_dept(k,3,i) = slm.t;
                end
        end
        
         csvwrite([RPATH 'network', num2str(i),'layer_myelin_presence-affectTC3.csv'], [...
                    t_dept(:,1,1), p_dept(:,1,1),...
                    t_dept(:,1,2), p_dept(:,1,2),...
                    t_dept(:,1,3), p_dept(:,1,3),...
                    t_dept(:,1,4), p_dept(:,1,4),...
                    t_dept(:,1,5), p_dept(:,1,5)])  
        
                csvwrite([RPATH 'network', num2str(i),'layer_myelin_presence-control.csv'], [...
                    t_dept(:,2,1), p_dept(:,2,1),...
                    t_dept(:,2,2), p_dept(:,2,2),...
                    t_dept(:,2,3), p_dept(:,2,3),...
                    t_dept(:,2,4), p_dept(:,2,4),...
                    t_dept(:,2,5), p_dept(:,2,5)])     
                
                csvwrite([RPATH 'network', num2str(i),'layer_myelin_affectTC3-control.csv'], [...
                    t_dept(:,3,1), p_dept(:,3,1),...
                    t_dept(:,3,2), p_dept(:,3,2),...
                    t_dept(:,3,3), p_dept(:,3,3),...
                    t_dept(:,3,4), p_dept(:,3,4),...
                    t_dept(:,3,5), p_dept(:,3,5)])  
                
                for key_descriptives = 1
                    for k = 1:12
                        Xn = squeeze(mean(layer_400_last(:,k,:),2));
                        keep_presence = (find(strcmp(groupN(keep),'Presence')))
                        keep_affect = (find(strcmp(groupN(keep),'Affect')))
                        keep_perspective = (find(strcmp(groupN(keep),'Perspective')))
                        keep_control   = (find(strcmp(groupN(keep),'Control1')))
                        
                        for i = 1:5
                            x = mean(Xn(keep(keep_control),maps_ns(:,i)>0),2);
                            meansdata(k,4,i) = mean(x);
                            stddata(k,4,i) = std(x);
                            SEM = std(x)/sqrt(length(x)); % Standard Error
                            ts = tinv([0.025  0.975],length(x)-1);      % T-Score
                            cidata = mean(x) + ts*SEM;
                            cidatamin(k,4,i) = cidata(1);
                            cidatamax(k,4,i) = cidata(2);
                        end
                        for i = 1:5
                            x = mean(Xn(keep(keep_affect),maps_ns(:,i)>0),2);
                            meansdata(k,2,i) = mean(x);
                            stddata(k,2,i) = std(x);
                            SEM = std(x)/sqrt(length(x)); % Standard Error
                            ts = tinv([0.025  0.975],length(x)-1);      % T-Score
                            cidata = mean(x) + ts*SEM;
                            cidatamin(k,2,i) = cidata(1);
                            cidatamax(k,2,i) = cidata(2);
                        end
%                         for i = 1:5
%                             x = mean(Xn(keep(keep_perspective),maps_ns(:,i)>0),2);
%                             meansdata(k,3,i) = mean(x);
%                             stddata(k,3,i) = std(x);
%                             SEM = std(x)/sqrt(length(x)); % Standard Error
%                             ts = tinv([0.025  0.975],length(x)-1);      % T-Score
%                             cidata = mean(x) + ts*SEM;
%                             cidatamin(k,3,i) = cidata(1);
%                             cidatamax(k,3,i) = cidata(2);
%                         end
                        for i = 1:5
                            x = mean(Xn(keep(keep_presence),maps_ns(:,i)>0),2);
                            meansdata(k,1,i) = mean(x);
                            stddata(k,1,i) = std(x);
                            SEM = std(x)/sqrt(length(x)); % Standard Error
                            ts = tinv([0.025  0.975],length(x)-1);      % T-Score
                            cidata = mean(x) + ts*SEM;
                            cidatamin(k,1,i) = cidata(1);
                            cidatamax(k,1,i) = cidata(2);
                        end
                    end
                end
                
                   csvwrite([RPATH 'stats_layer_myelin_control.csv'], [...
                    meansdata(:,4,1), stddata(:,4,1),cidatamin(:,4,1),cidatamax(:,4,1),...
                    meansdata(:,4,2), stddata(:,4,2),cidatamin(:,4,2),cidatamax(:,4,2),...
                    meansdata(:,4,3), stddata(:,4,3),cidatamin(:,4,3),cidatamax(:,4,3),...
                    meansdata(:,4,4), stddata(:,4,4),cidatamin(:,4,4),cidatamax(:,4,4),...
                    meansdata(:,4,5), stddata(:,4,5),cidatamin(:,4,5),cidatamax(:,4,5)])  
                
                csvwrite([RPATH 'stats_layer_myelin_presence.csv'], [...
                    meansdata(:,1,1), stddata(:,1,1),cidatamin(:,1,1),cidatamax(:,1,1),...
                    meansdata(:,1,2), stddata(:,1,2),cidatamin(:,1,2),cidatamax(:,1,2),...
                    meansdata(:,1,3), stddata(:,1,3),cidatamin(:,1,3),cidatamax(:,1,3),...
                    meansdata(:,1,4), stddata(:,1,4),cidatamin(:,1,4),cidatamax(:,1,4),...
                    meansdata(:,1,5), stddata(:,1,5),cidatamin(:,1,5),cidatamax(:,1,5)])  
                
                csvwrite([RPATH 'stats_layer_myelin_affect3.csv'], [...
                    meansdata(:,2,1), stddata(:,2,1),cidatamin(:,2,1),cidatamax(:,2,1),...
                    meansdata(:,2,2), stddata(:,2,2),cidatamin(:,2,2),cidatamax(:,2,2),...
                    meansdata(:,2,3), stddata(:,2,3),cidatamin(:,2,3),cidatamax(:,2,3),...
                    meansdata(:,2,4), stddata(:,2,4),cidatamin(:,2,4),cidatamax(:,2,4),...
                    meansdata(:,2,5), stddata(:,2,5),cidatamin(:,2,5),cidatamax(:,2,5)])  
                
                
         
   
    end 

    for t1_t3_layers = 1
        Xn            = squeeze(layer_400_last(:,1,:));
        Xn(isnan(Xn)) = 0;
        Xn(isinf(Xn)) = 0;
        Xm            = ggg_400_last;
        Xm(isnan(Xm)) = 0;
        Xm(isinf(Xm)) = 0;
        keep1   = union(find(strcmp(groupN,'Affect')),    find(strcmp(groupN,'Perspective')));
        keep1   = union(keep1,find(strcmp(groupN,'Presence')))
        keep1   = union(keep1,find(strcmp(groupN,'Control1')))
        keep2   = intersect(find(abs(mean(Xn,2)) <666), find(sum(Xn,2)~=0));
        keep3   = find(tpnum>1);  
        keep4   = intersect(find(abs(mean(Xm,2)) <666), find(sum(Xm,2)~=0));
        keep    = mintersect(keep1, keep2,keep3,keep4);
        
        ik      = id(keep,:);
        ink     = idnum(keep);
        tk      = tpnum(keep,:);
        tnk     = tp(keep,:);
        gk      = group(keep);
        g4k     = group4(keep);
        gNk     = groupN(keep);
        ak      = age(keep);
        sk      = cellstr(sex(keep));   
        
        A       = term(ak);
        S       = term(sk);
        G       = term(gk);
        G4      = term(g4k);
        GN      = term(gNk);
        Sub     = term(ink);
        Tn      = term(tnk);
        
        M        = 1 + A + S + GN + random(Sub) + I;
 
        for i = 1:5
                for k=1:12
                    Xx =  squeeze(mean(layer_400_last(:,k,:),2));
                    slm      = SurfStatLinMod(mean(Xx(keep,maps_ns(:,i)>0),2),M);
                    slm      = SurfStatT(slm,(GN.Affect-GN.Control1));
                    slm2 = slm.t;
                    p_dept(k,2,i)  = 2*(1 - tcdf(abs(slm.t),slm.df));
                    t_dept(k,2,i) = slm.t;
                    Xx =  squeeze(mean(layer_400_last(:,k,:),2));
                    slm      = SurfStatLinMod(mean(Xx(keep,maps_ns(:,i)>0),2),M);
                    slm      = SurfStatT(slm,(GN.Perspective-GN.Control1));
                    slm2 = slm.t;
                    p_dept(k,3,i)  = 2*(1 - tcdf(abs(slm.t),slm.df));
                    t_dept(k,3,i) = slm.t;
                end
        end
        
                csvwrite([RPATH 'network', num2str(i),'layer_myelin_affect-control.csv'], [...
                    t_dept(:,2,1), p_dept(:,2,1),...
                    t_dept(:,2,2), p_dept(:,2,2),...
                    t_dept(:,2,3), p_dept(:,2,3),...
                    t_dept(:,2,4), p_dept(:,2,4),...
                    t_dept(:,2,5), p_dept(:,2,5)])     
                
                csvwrite([RPATH 'network', num2str(i),'layer_myelin_persepctive-control.csv'], [...
                    t_dept(:,3,1), p_dept(:,3,1),...
                    t_dept(:,3,2), p_dept(:,3,2),...
                    t_dept(:,3,3), p_dept(:,3,3),...
                    t_dept(:,3,4), p_dept(:,3,4),...
                    t_dept(:,3,5), p_dept(:,3,5)])  
                
           for key_descriptives = 1
                    for k = 1:12
                        Xn = squeeze(mean(layer_400_last(:,k,:),2));
                        keep_presence = (find(strcmp(groupN(keep),'Presence')))
                        keep_affect = (find(strcmp(groupN(keep),'Affect')))
                        keep_perspective = (find(strcmp(groupN(keep),'Perspective')))
                        keep_control   = (find(strcmp(groupN(keep),'Control1')))
                        
                        for i = 1:5
                            x = mean(Xn(keep(keep_control),maps_ns(:,i)>0),2);
                            meansdata(k,4,i) = mean(x);
                            stddata(k,4,i) = std(x);
                            SEM = std(x)/sqrt(length(x)); % Standard Error
                            ts = tinv([0.025  0.975],length(x)-1);      % T-Score
                            cidata = mean(x) + ts*SEM;
                            cidatamin(k,4,i) = cidata(1);
                            cidatamax(k,4,i) = cidata(2);
                        end
                        for i = 1:5
                            x = mean(Xn(keep(keep_affect),maps_ns(:,i)>0),2);
                            meansdata(k,2,i) = mean(x);
                            stddata(k,2,i) = std(x);
                            SEM = std(x)/sqrt(length(x)); % Standard Error
                            ts = tinv([0.025  0.975],length(x)-1);      % T-Score
                            cidata = mean(x) + ts*SEM;
                            cidatamin(k,2,i) = cidata(1);
                            cidatamax(k,2,i) = cidata(2);
                        end
                        for i = 1:5
                            x = mean(Xn(keep(keep_perspective),maps_ns(:,i)>0),2);
                            meansdata(k,3,i) = mean(x);
                            stddata(k,3,i) = std(x);
                            SEM = std(x)/sqrt(length(x)); % Standard Error
                            ts = tinv([0.025  0.975],length(x)-1);      % T-Score
                            cidata = mean(x) + ts*SEM;
                            cidatamin(k,3,i) = cidata(1);
                            cidatamax(k,3,i) = cidata(2);
                        end
                        for i = 1:5
                            x = mean(Xn(keep(keep_presence),maps_ns(:,i)>0),2);
                            meansdata(k,1,i) = mean(x);
                            stddata(k,1,i) = std(x);
                            SEM = std(x)/sqrt(length(x)); % Standard Error
                            ts = tinv([0.025  0.975],length(x)-1);      % T-Score
                            cidata = mean(x) + ts*SEM;
                            cidatamin(k,1,i) = cidata(1);
                            cidatamax(k,1,i) = cidata(2);
                        end
                    end
                end
                
                   csvwrite([RPATH 'stats_layer_myelin_controlT1-T3.csv'], [...
                    meansdata(:,4,1), stddata(:,4,1),cidatamin(:,4,1),cidatamax(:,4,1),...
                    meansdata(:,4,2), stddata(:,4,2),cidatamin(:,4,2),cidatamax(:,4,2),...
                    meansdata(:,4,3), stddata(:,4,3),cidatamin(:,4,3),cidatamax(:,4,3),...
                    meansdata(:,4,4), stddata(:,4,4),cidatamin(:,4,4),cidatamax(:,4,4),...
                    meansdata(:,4,5), stddata(:,4,5),cidatamin(:,4,5),cidatamax(:,4,5)])  
                
                csvwrite([RPATH 'stats_layer_myelin_perspective.csv'], [...
                    meansdata(:,3,1), stddata(:,3,1),cidatamin(:,3,1),cidatamax(:,3,1),...
                    meansdata(:,3,2), stddata(:,3,2),cidatamin(:,3,2),cidatamax(:,3,2),...
                    meansdata(:,3,3), stddata(:,3,3),cidatamin(:,3,3),cidatamax(:,3,3),...
                    meansdata(:,3,4), stddata(:,3,4),cidatamin(:,3,4),cidatamax(:,3,4),...
                    meansdata(:,3,5), stddata(:,3,5),cidatamin(:,3,5),cidatamax(:,3,5)])  
                
                csvwrite([RPATH 'stats_layer_myelin_affect.csv'], [...
                    meansdata(:,2,1), stddata(:,2,1),cidatamin(:,2,1),cidatamax(:,2,1),...
                    meansdata(:,2,2), stddata(:,2,2),cidatamin(:,2,2),cidatamax(:,2,2),...
                    meansdata(:,2,3), stddata(:,2,3),cidatamin(:,2,3),cidatamax(:,2,3),...
                    meansdata(:,2,4), stddata(:,2,4),cidatamin(:,2,4),cidatamax(:,2,4),...
                    meansdata(:,2,5), stddata(:,2,5),cidatamin(:,2,5),cidatamax(:,2,5)])  
                
      
    end

    for cortical_thickness_correct = 1
        Xn            = squeeze(layer_400_last(:,1,:));
        Xn(isnan(Xn)) = 0;
        Xn(isinf(Xn)) = 0;
        Xm            = ggg_400_last;
        Xm(isnan(Xm)) = 0;
        Xm(isinf(Xm)) = 0;
        keep1   = union(find(strcmp(groupN,'Affect')),    find(strcmp(groupN,'Perspective')));
        keep1   = union(keep1,find(strcmp(groupN,'Presence')))
        keep1   = union(keep1,find(strcmp(groupN,'Control1')))
        keep2   = intersect(find(abs(mean(Xn,2)) <666), find(sum(Xn,2)~=0));
        keep3   = intersect(find(~strcmp(groupN,'Group3')),find(tpnum>0));
        keep4   = intersect(find(abs(mean(Xm,2)) <666), find(sum(Xm,2)~=0));
        keep    = mintersect(keep1, keep2,keep3,keep4);
        
        ik      = id(keep,:);
        ink     = idnum(keep);
        tk      = tpnum(keep,:);
        tnk     = tp(keep,:);
        gk      = group(keep);
        g4k     = group4(keep);
        gNk     = groupN(keep);
        ak      = age(keep);
        sk      = cellstr(sex(keep));
        
        A       = term(ak);
        S       = term(sk);
        G       = term(gk);
        G4      = term(g4k);
        GN      = term(gNk);
        Sub     = term(ink);
        Tn      = term(tnk);
        
        M        = 1 + A + S + GN + random(Sub) + I;
        
        for i = 1:5
            for k=1:12
                Xx =  squeeze(mean(layer_400_last_ctx(:,k,:),2));
                slm      = SurfStatLinMod(mean(Xx(keep,maps_ns(:,i)>0),2),M);
                slm      = SurfStatT(slm,(GN.Presence-GN.Perspective));
                slm2 = slm.t;
                p_dept(k,1,i)  = 2*(1 - tcdf(abs(slm.t),slm.df));
                t_dept(k,1,i) = slm.t;
                Xx =  squeeze(mean(layer_400_last_ctx(:,k,:),2));
                slm      = SurfStatLinMod(mean(Xx(keep,maps_ns(:,i)>0),2),M);
                slm      = SurfStatT(slm,(GN.Presence-GN.Affect));
                slm2 = slm.t;
                p_dept(k,2,i)  = 2*(1 - tcdf(abs(slm.t),slm.df));
                t_dept(k,2,i) = slm.t;
                Xx =  squeeze(mean(layer_400_last_ctx(:,k,:),2));
                slm      = SurfStatLinMod(mean(Xx(keep,maps_ns(:,i)>0),2),M);
                slm      = SurfStatT(slm,(GN.Perspective-GN.Affect));
                slm2 = slm.t;
                p_dept(k,3,i)  = 2*(1 - tcdf(abs(slm.t),slm.df));
                t_dept(k,3,i) = slm.t;
            end
        end
        
        csvwrite([RPATH 'network', num2str(i),'layer_myelin_presence-perspective_ctx.csv'], [...
            t_dept(:,1,1), p_dept(:,1,1),...
            t_dept(:,1,2), p_dept(:,1,2),...
            t_dept(:,1,3), p_dept(:,1,3),...
            t_dept(:,1,4), p_dept(:,1,4),...
            t_dept(:,1,5), p_dept(:,1,5)])
        
        csvwrite([RPATH 'network', num2str(i),'layer_myelin_presence-affect_ctx.csv'], [...
            t_dept(:,2,1), p_dept(:,2,1),...
            t_dept(:,2,2), p_dept(:,2,2),...
            t_dept(:,2,3), p_dept(:,2,3),...
            t_dept(:,2,4), p_dept(:,2,4),...
            t_dept(:,2,5), p_dept(:,2,5)])
        
        csvwrite([RPATH 'network', num2str(i),'layer_myelin_perspective-affect_ctx.csv'], [...
            t_dept(:,3,1), p_dept(:,3,1),...
            t_dept(:,3,2), p_dept(:,3,2),...
            t_dept(:,3,3), p_dept(:,3,3),...
            t_dept(:,3,4), p_dept(:,3,4),...
            t_dept(:,3,5), p_dept(:,3,5)])
        
        fdr_bh(p_dept)
    end
 
    for all_layers_contrasts = 1   
        Xn            = squeeze(layer_400_last(:,1,:));
        Xn(isnan(Xn)) = 0;
        Xn(isinf(Xn)) = 0;
        Xm            = ggg_400_last;
        Xm(isnan(Xm)) = 0;
        Xm(isinf(Xm)) = 0;
        keep1   = union(find(strcmp(groupN,'Affect')),    find(strcmp(groupN,'Perspective')));
        keep1   = union(keep1,find(strcmp(groupN,'Presence')))
        keep1   = union(keep1,find(strcmp(groupN,'Control1')))
        keep2   = intersect(find(abs(mean(Xn,2)) <666), find(sum(Xn,2)~=0));
        keep3   = intersect(find(~strcmp(groupN,'Group3')),find(tpnum>0));
        keep4   = intersect(find(abs(mean(Xm,2)) <666), find(sum(Xm,2)~=0));
        keep    = mintersect(keep1, keep2,keep3,keep4);
        
        ik      = id(keep,:);
        ink     = idnum(keep);
        tk      = tpnum(keep,:);
        tnk     = tp(keep,:);
        gk      = group(keep);
        g4k     = group4(keep);
        gNk     = groupN(keep);
        ak      = age(keep);
        sk      = cellstr(sex(keep));
        
        A       = term(ak);
        S       = term(sk);
        G       = term(gk);
        G4      = term(g4k);
        GN      = term(gNk);
        Sub     = term(ink);
        Tn      = term(tnk);
        
        M        = 1 + A + S + GN + random(Sub) + I;
        
        for i  = 1:12
            xxx = squeeze(layer_400_last(:,i,:));
            slm      = SurfStatLinMod(xxx(keep,:),M);
            slm      = SurfStatT(slm,(GN.Presence-GN.Perspective));
            p        = 2*(1 - tcdf(abs(slm.t),slm.df))
            
            hp = fdr_bh(p(1:400),0.05)
            h  = hp;
                        
            tmp = zeros(20484,1);
            for i = 1:200
                tmp(find(parcels400==i+1)) = slm.t(i).*h(i);
            end
            for i = 1:200
                tmp(find(parcels400==i+1001)) =slm.t(i+200).*h(i+200);
            end
            
            f=figure,
            BoSurfStatViewData(tmp,SN,'')
            BoSurfStatColLim([-4 4])
            colormap(flipud(cbrewer('div','RdBu',11)))
            exportfigbo(f,[RPATH 'F3s.Perspective-Presence.layer', num2str(i) ,'.fdr.png'],'png', 10)
            close(f)
        end
        for i  = 1:12
            xxx = squeeze(layer_400_last(:,i,:));
            slm      = SurfStatLinMod(xxx(keep,:),M);
            slm      = SurfStatT(slm,(GN.Presence-GN.Affect));
            p        = 1 - tcdf(slm.t,slm.df)
            
            hp = fdr_bh(p(1:400),0.025)
            hn = fdr_bh(1-p(1:400),0.025)
            h  = hp+hn;
                        
            tmp = zeros(20484,1);
            for i = 1:200
                tmp(find(parcels400==i+1)) = slm.t(i).*h(i);
            end
            for i = 1:200
                tmp(find(parcels400==i+1001)) =slm.t(i+200).*h(i+200);
            end
            
            f=figure,
            BoSurfStatViewData(tmp,SN,'')
            BoSurfStatColLim([-4 4])
            colormap(flipud(cbrewer('div','RdBu',11)))
            exportfigbo(f,[RPATH 'F3s.Presence-Affect.layer', num2str(i) ,'.fdr.png'],'png', 10)
            close(f)
        end
        for i  = 1:12
            xxx = squeeze(layer_400_last(:,i,:));
            slm      = SurfStatLinMod(xxx(keep,:),M);
            slm      = SurfStatT(slm,(GN.Perspective-GN.Affect));
            p        = 1 - tcdf(slm.t,slm.df)
            
            hp = fdr_bh(p(1:400),0.025)
            hn = fdr_bh(1-p(1:400),0.025)
            h  = hp+hn;
                        
            tmp = zeros(20484,1);
            for i = 1:200
                tmp(find(parcels400==i+1)) = slm.t(i).*h(i);
            end
            for i = 1:200
                tmp(find(parcels400==i+1001)) =slm.t(i+200).*h(i+200);
            end
            
            f=figure,
            BoSurfStatViewData(tmp,SN,'')
            BoSurfStatColLim([-4 4])
            colormap(flipud(cbrewer('div','RdBu',11)))
            exportfigbo(f,[RPATH 'F3s.Perspective-Affect.layer', num2str(i) ,'.fdr.png'],'png', 10)
            close(f)
        end

    end
    
    for all_training_vs_control = 1
        Xn            = squeeze(layer_400_bl(:,1,:));
        Xn(isnan(Xn)) = 0;
        Xn(isinf(Xn)) = 0;
        Xm            = ggg_400_bl;
        Xm(isnan(Xm)) = 0;
        Xm(isinf(Xm)) = 0;
        keep1   = union(find(strcmp(groupN,'Affect')),    find(strcmp(groupN,'Perspective')));
        keep1   = union(keep1,find(strcmp(groupN,'Presence')))
        keep1   = union(keep1,find(strcmp(groupN,'Control1')))
        keep2   = intersect(find(abs(mean(Xn,2)) <666), find(sum(Xn,2)~=0));
        keep3   = intersect(find(~strcmp(groupN,'Group3')),find(tpnum>2));
        keep4   = intersect(find(abs(mean(Xm,2)) <666), find(sum(Xm,2)~=0));
        keep    = mintersect(keep1, keep2,keep3,keep4);
        
        ik      = id(keep,:);
        ink     = idnum(keep);
        tk      = tpnum(keep,:);
        tnk     = tp(keep,:);
        gk      = group(keep);
        g4k     = group4(keep);
        gNk     = groupN(keep);
        ak      = age(keep);
        sk      = cellstr(sex(keep));
        
        A       = term(ak);
        S       = term(sk);
        G       = term(gk);
        G4      = term(g4k);
        GN      = term(gNk);
        Sub     = term(ink);
        Tn      = term(tnk);
        g2k         = gNk;
        g2k(strcmp(gNk,'Control1')) = {'Control'};
        g2k(~strcmp(gNk,'Control1')) = {'Training'};
      
        G2  = term(g2k);   
        
        M        = 1 + A + S + G2;
   
        for j  = 1:4:12
            xxx      = squeeze(mean(layer_400_bl(:,j:(j+3),:),2));
            slm      = SurfStatLinMod(xxx(keep,:),M);
            slm      = SurfStatT(slm,(G2.Training-G2.Control));
            p        = 2*(1 - tcdf(abs(slm.t),slm.df));
            
            hp = fdr_bh(p(1:400),0.05)
            h  = hp;
                        
            tmp = zeros(20484,1);
            for i = 1:200
                tmp(find(parcels400==i+1)) = slm.t(i).*h(i);
            end
            for i = 1:200
                tmp(find(parcels400==i+1001)) =slm.t(i+200).*h(i+200);
            end
            
            f=figure,
            BoSurfStatViewData(tmp,SN,'')
            BoSurfStatColLim([-4 4])
            colormap(flipud(cbrewer('div','RdBu',11)))
            exportfigbo(f,[RPATH 'F3s.training_control.layer', num2str(j) ,'.fdr.png'],'png', 10)
            close(f)
        end
        
        for j  = 1:4:12
            xxx      = squeeze(mean(layer_400_bl(:,j:(j+3),:),2));
            slm      = SurfStatLinMod(xxx(keep,:),M);
            slm      = SurfStatT(slm,(G2.Control));
            p        = 2*(1 - tcdf(abs(slm.t),slm.df));
            
            hp = fdr_bh(p(1:400),0.05)
            h  = hp;
                        
            tmp = zeros(20484,1);
            for i = 1:200
                tmp(find(parcels400==i+1)) = slm.t(i).*h(i);
            end
            for i = 1:200
                tmp(find(parcels400==i+1001)) =slm.t(i+200).*h(i+200);
            end
            
            f=figure,
            BoSurfStatViewData(tmp,SN,'')
            BoSurfStatColLim([-4 4])
            colormap(flipud(cbrewer('div','RdBu',11)))
            exportfigbo(f,[RPATH 'F3s.control.layer', num2str(j) ,'.fdr.png'],'png', 10)
            close(f)
            
            h  = p(1:400)<0.01;
                        
            tmp = zeros(20484,1);
            for i = 1:200
                tmp(find(parcels400==i+1)) = slm.t(i).*h(i);
            end
            for i = 1:200
                tmp(find(parcels400==i+1001)) =slm.t(i+200).*h(i+200);
            end
            
            f=figure,
            BoSurfStatViewData(tmp,SN,'')
            BoSurfStatColLim([-4 4])
            colormap(flipud(cbrewer('div','RdBu',11)))
            exportfigbo(f,[RPATH 'F3s.control.layer', num2str(j) ,'.p.png'],'png', 10)
            close(f)
        end
 
        neuros = zeros(5,3);
        for i  = 1:12
            xxx = squeeze(layer_400_bl(:,i,:));
            for j = 1:5
                slm      = SurfStatLinMod(mean(xxx(keep,maps_ns(:,j)>0),2),M);
                slm      = SurfStatT(slm,G2.Training-G2.Control)
                t_dept(i,1,j) = slm.t;
                p_dept(i,1,j) = 2*(1 - tcdf(abs(slm.t),slm.df));
            end
        end
        
        fdr_bh(p_dept)
           csvwrite([RPATH 'network', num2str(i),'layer_myelin_all_training-control.csv'], [...
                    t_dept(:,1,1), p_dept(:,1,1),...
                    t_dept(:,1,2), p_dept(:,1,2),...
                    t_dept(:,1,3), p_dept(:,1,3),...
                    t_dept(:,1,4), p_dept(:,1,4),...
                    t_dept(:,1,5), p_dept(:,1,5)])  
    
        
        for  i = 1:5
        for human_profile = 1
           f=figure,
           ax = axes;
           ax.ColorOrder =  [0.8275    0.8902    0.9529;     0.5294    0.5294    0.5294;      0.8784    0.8784    0.8784];
           hold on
           keep_control1 = intersect(find(strcmp(groupN(keep),'Control1')),find(strcmp(tp(keep),'T3')))
           keep_aff = intersect(find(strcmp(groupN(keep),'Affect')),find(strcmp(tp(keep),'T3')))
           keep_pers = intersect(find(strcmp(groupN(keep),'Perspective')),find(strcmp(tp(keep),'T3')))

           plot(mean(mean(squeeze(layer_400_bl(keep(keep_control1),:,maps_ns(:,i)>0))),3),1:12,'LineWidth',10)
           hold on
           plot(mean(mean(squeeze(layer_400_bl(keep(keep_aff),:,maps_ns(:,i)>0))),3),1:12,'LineWidth',10)
           hold on
           plot(mean(mean(squeeze(layer_400_bl(keep(keep_pers),:,maps_ns(:,i)>0))),3),1:12,'LineWidth',10)
           exportfigbo(f,[RPATH 'S.ince_baseline.network',num2str(i),'.layerz.png'],'png', 10)
        end
        end
    end
end


for f4 = 1
    for make_attention = 1
        val             = Att;
        Xn            = squeeze(layer_400_last(:,1,:));
        Xn(isnan(Xn)) = 0;
        Xn(isinf(Xn)) = 0;
        Xm              = ggg_400_last;
        Xm(isnan(Xm))   = 0;
        keep1   = (find(strcmp(groupN,'Presence')));
        keep2   = intersect(find(Xn(:,1) >-666), find(sum(Xn,2)~=0));
        keep2b   = intersect(find(Xm(:,1) >-666), find(sum(Xm,2)~=0));
        keep4   = find(~(strcmp(group4,'Group3')));
        keep5   = find(tpnum>0);
        keep6   = find(abs(val)<666);
        keep    = mintersect(keep1,keep2,keep2b,keep4,keep5,keep6);
        
        ak      = age(keep);
        sk      = sex(keep);
        
        A       = term(ak);
        S       = term(sk);
        
        Ck = []
        for k=1:5
            Xx =  mean(squeeze(mean(layer_400_last(keep,1:4,maps_ns(:,k)>0),2)),2);
            Ck(:,k)  =  Xx;
            Xx =  mean(squeeze(mean(layer_400_last(keep,5:8,maps_ns(:,k)>0),2)),2);
            Ck(:,k+5)  =  Xx;
            Xx =  mean(squeeze(mean(layer_400_last(keep,9:12,maps_ns(:,k)>0),2)),2);
            Ck(:,k+10)  =  Xx;
            Xx =  mean(ggg_400_last(keep,maps_ns(:,k)>0),2);
            Ck(:,k+15)  =  Xx;
            Xx =  mean(g1_400_last(keep,maps_ns(:,k)>0),2);
            Ck(:,k+20)  =  Xx;
            Xx =  mean(g2_400_last(keep,maps_ns(:,k)>0),2);
            Ck(:,k+25)  =  Xx;
            Xx =  mean(g3_400_last(keep,maps_ns(:,k)>0),2);
            Ck(:,k+30)  =  Xx;
        end
        
        filename = ['attention_data.csv']
        age1     = age(keep);
        sex1     = sex(keep);
        node     = Ck;
        val1     = val(keep);
        T = table(val1, age1, sex1, node)
        writetable(T, [filename]);
        
    end
    
    for make_compassion = 1
        val             = Comp;
        Xn            = squeeze(layer_400_last(:,1,:));
        Xn(isnan(Xn)) = 0;
        Xn(isinf(Xn)) = 0;
        Xm              = ggg_400_last;
        Xm(isnan(Xm))   = 0;
        keep1   = (find(strcmp(groupN,'Affect')));
        keep2   = intersect(find(Xn(:,1) >-666), find(sum(Xn,2)~=0));
        keep2b   = intersect(find(Xm(:,1) >-666), find(sum(Xm,2)~=0));
        keep4   = find(~(strcmp(group4,'Group3')));
        keep5   = find(tpnum>0);
        keep6   = find(abs(val)<666);
        keep    = mintersect(keep1,keep2,keep2b,keep4,keep5,keep6);    
        
        Ck = []
        for k=1:5
            Xx =  mean(squeeze(mean(layer_400_last(keep,1:4,maps_ns(:,k)>0),2)),2);
            Ck(:,k)  =  Xx;
            Xx =  mean(squeeze(mean(layer_400_last(keep,5:8,maps_ns(:,k)>0),2)),2);
            Ck(:,k+5)  =  Xx;
            Xx =  mean(squeeze(mean(layer_400_last(keep,9:12,maps_ns(:,k)>0),2)),2);
            Ck(:,k+10)  =  Xx;
            Xx =  mean(ggg_400_last(keep,maps_ns(:,k)>0),2);
            Ck(:,k+15)  =  Xx;
            Xx =  mean(g1_400_last(keep,maps_ns(:,k)>0),2);
            Ck(:,k+20)  =  Xx;
            Xx =  mean(g2_400_last(keep,maps_ns(:,k)>0),2);
            Ck(:,k+25)  =  Xx;
            Xx =  mean(g3_400_last(keep,maps_ns(:,k)>0),2);
            Ck(:,k+30)  =  Xx;
        end
        
        filename = ['compassion_data.csv']
        age1     = age(keep);
        sex1     = sex(keep);
        node     = Ck;
        val1     = val(keep);
        T = table(val1, age1, sex1, node)
        writetable(T, [filename]);    
    end
    
    for make_tom = 1
        val             = Tom;
        Xn            = squeeze(layer_400_last(:,1,:));
        Xn(isnan(Xn)) = 0;
        Xn(isinf(Xn)) = 0;
        Xm              = ggg_400_last;
        Xm(isnan(Xm))   = 0;
        keep1   = (find(strcmp(groupN,'Perspective')));
        keep2   = intersect(find(Xn(:,1) >-666), find(sum(Xn,2)~=0));
        keep2b   = intersect(find(Xm(:,1) >-666), find(sum(Xm,2)~=0));
        keep4   = find(~(strcmp(group4,'Group3')));
        keep5   = find(tpnum>0);
        keep6   = find(abs(val)<666);
        keep    = mintersect(keep1,keep2,keep2b,keep4,keep5,keep6);
       
        Ck = []
        for k=1:5
            Xx =  mean(squeeze(mean(layer_400_last(keep,1:4,maps_ns(:,k)>0),2)),2);
            Ck(:,k)  =  Xx;
            Xx =  mean(squeeze(mean(layer_400_last(keep,5:8,maps_ns(:,k)>0),2)),2);
            Ck(:,k+5)  =  Xx;
            Xx =  mean(squeeze(mean(layer_400_last(keep,9:12,maps_ns(:,k)>0),2)),2);
            Ck(:,k+10)  =  Xx;
            Xx =  mean(ggg_400_last(keep,maps_ns(:,k)>0),2);
            Ck(:,k+15)  =  Xx;
            Xx =  mean(g1_400_last(keep,maps_ns(:,k)>0),2);
            Ck(:,k+20)  =  Xx;
            Xx =  mean(g2_400_last(keep,maps_ns(:,k)>0),2);
            Ck(:,k+25)  =  Xx;
             Xx =  mean(g3_400_last(keep,maps_ns(:,k)>0),2);
            Ck(:,k+30)  =  Xx;
        end
  
        filename = ['tom_data.csv']
        age1     = age(keep);
        sex1     = sex(keep);
        node     = Ck;
        val1     = val(keep);
        T = table(val1, age1, sex1, node)
        writetable(T, [filename]);
  
    end
end









     
  


 
        
        

        
        

 

 
 