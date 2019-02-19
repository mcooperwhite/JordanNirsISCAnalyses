%This script organizes Jordan data and computes inter-subject correlations
%for all channels

%MCW - first analysis

%This script will only run in 2016b onward because of function "contains"

raw_dir=pwd;
dataprefix='JOR';
cd('Life');
currdir=dir(strcat(dataprefix,'*'));
frames_life=2531;
frames_choice=2195;

lifePs_alldata_life=zeros(frames_life,20,length(currdir));
lifePs_alldata_choice=zeros(frames_choice,20,length(currdir));
for i=1:length(currdir)
    cd(currdir(i).name);
    scannames=dir(strcat(currdir(i).name,'*'));
    for j=1:length(scannames)
        if ~isempty(strfind(scannames(j).name,'life')) && isempty(strfind(scannames(j).name,'NO_EVENTS'))
            file=dir(strcat(scannames(j).name,'/','JOR_*','.mat'));
            filename=file.name;
            load(strcat(scannames(j).name,'/',filename));
            lifePs_alldata_life(:,:,i)=z_oxy1;
        end
        if ~isempty(strfind(scannames(j).name,'choice')) && isempty(strfind(scannames(j).name,'NO_EVENTS'))
            file=dir(strcat(scannames(j).name,'/','JOR_*','.mat'));
            filename=file.name;
            load(strcat(scannames(j).name,'/',filename));
            lifePs_alldata_choice(:,:,i)=z_oxy1;
        end
    end
    cd ..
    
  

%first, compute averages
%d=zeros(2531,20);


end
% choicesize=size(lifePs_alldata_choice);
% lifesize=size(lifePs_alldata_life);
% lifePs_avgtimecourse_life=zeros(lifesize(1),lifesize(2));
% lifePs_avgtimecourse_choice=zeros(choicesize(1),choicesize(2));
%  for y=1:choicesize(2)
%      for x=1:choicesize(1) 
%          lifePs_avgtimecourse_choice(x,y)=mean(lifePs_alldata_choice(x,y,~isnan(lifePs_alldata_choice(x,y,:)))); 
%      end
% end
% for y=1:lifesize(2)
%      for x=1:lifesize(1) 
%          lifePs_avgtimecourse_life(x,y)=mean(lifePs_alldata_life(x,y,~isnan(lifePs_alldata_life(x,y,:)))); 
%      end
%  end

cd ..
cd('Choice')
clear currdir
currdir=dir(strcat(dataprefix,'*'));


choicePs_alldata_life=zeros(frames_life,20,length(currdir));
choicePs_alldata_choice=zeros(frames_choice,20,length(currdir));
for i=1:length(currdir)
    msg = sprintf('\n\t P number %d/%d ...',i,length(currdir));
    fprintf(msg);
    cd(currdir(i).name);
    scannames=dir(strcat(currdir(i).name,'*'));
    for j=1:length(scannames)
        if ~isempty(strfind(scannames(j).name,'life')) && isempty(strfind(scannames(j).name,'NO_EVENTS'))
            file=dir(strcat(scannames(j).name,'/','JOR_*','.mat'));
            filename=file.name;
            load(strcat(scannames(j).name,'/',filename));
            choicePs_alldata_life(:,:,i)=z_oxy1;
        end
        if ~isempty(strfind(scannames(j).name,'choice')) && isempty(strfind(scannames(j).name,'NO_EVENTS'))
            file=dir(strcat(scannames(j).name,'/','JOR_*','.mat'));
            filename=file.name;
            load(strcat(scannames(j).name,'/',filename));
            choicePs_alldata_choice(:,:,i)=z_oxy1;
        end
    end
    cd ..
    
  

%first, compute averages
%d=zeros(2531,20);


end
% choicesize=size(choicePs_alldata_choice);
% lifesize=size(choicePs_alldata_life);
% choicePs_avgtimecourse_life=zeros(lifesize(1),lifesize(2));
% choicePs_avgtimecourse_choice=zeros(choicesize(1),choicesize(2));
% for y=1:choicesize(2)
%     for x=1:choicesize(1) 
%         choicePs_avgtimecourse_choice(x,y)=mean(choicePs_alldata_choice(x,y,~isnan(choicePs_alldata_choice(x,y,:)))); 
%     end
% end
% for y=1:lifesize(2)
%     for x=1:lifesize(1) 
%         choicePs_avgtimecourse_life(x,y)=mean(choicePs_alldata_life(x,y,~isnan(choicePs_alldata_life(x,y,:)))); 
%     end
% end    
clear z_deoxy1 z_oxy1 z_totaloxy1 totaloxy1 oxy1 deoxy1
cd ..
cd('Life')
currdir=dir(strcat(dataprefix,'*'));
choicesize=size(choicePs_alldata_choice);
lifesize=size(choicePs_alldata_life);
choicePs_avgtimecourse_life=zeros(lifesize(1),lifesize(2));
choicePs_avgtimecourse_choice=zeros(choicesize(1),choicesize(2));
for y=1:choicesize(2)
    for x=1:choicesize(1) 
        choicePs_avgtimecourse_choice(x,y)=mean(choicePs_alldata_choice(x,y,~isnan(choicePs_alldata_choice(x,y,:)))); 
    end
end
for y=1:lifesize(2)
    for x=1:lifesize(1) 
        choicePs_avgtimecourse_life(x,y)=mean(choicePs_alldata_life(x,y,~isnan(choicePs_alldata_life(x,y,:)))); 
    end
end 
for i=1:length(currdir)
    msg = sprintf('\n\t P number %d/%d ...',i,length(currdir));
    fprintf(msg);
    cd(currdir(i).name);
    scannames=dir(strcat(currdir(i).name,'*'));
    
    
        
        choicesize=size(lifePs_alldata_choice);
        lifesize=size(lifePs_alldata_life);
        lifePs_avgtimecourse_life=zeros(lifesize(1),lifesize(2));
        lifePs_avgtimecourse_choice=zeros(choicesize(1),choicesize(2));
        for y=1:choicesize(2)
            for x=1:choicesize(1) 
                mylifePs_alldata_choice=lifePs_alldata_choice;
                mylifePs_alldata_choice(:,:,j)=[]; %delete P's timecourse from whole dataset (temporarily)
                lifePs_avgtimecourse_choice(x,y)=mean(mylifePs_alldata_choice(x,y,~isnan(mylifePs_alldata_choice(x,y,:)))); 
            end
        end
        for y=1:lifesize(2)
            for x=1:lifesize(1) 
                mylifePs_alldata_life=lifePs_alldata_life;
                mylifePs_alldata_life(:,:,j)=[];
                lifePs_avgtimecourse_life(x,y)=mean(mylifePs_alldata_life(x,y,~isnan(mylifePs_alldata_life(x,y,:)))); 
            end
        end
    
    for j=1:length(scannames)
        
        

        if ~isempty(strfind(scannames(j).name,'life')) && isempty(strfind(scannames(j).name,'NO_EVENTS'))
            file=dir(strcat(scannames(j).name,'/','JOR_*','.mat'));
            filename=file.name;
            load(strcat(scannames(j).name,'/',filename));
            for z=1:lifesize(2)
                
                all_life_corrs_withlife_forlife(i,z)=corr(z_oxy1(:,z),lifePs_avgtimecourse_life(:,z));
                all_life_corrs_withchoice_forlife(i,z)=corr(z_oxy1(:,z),choicePs_avgtimecourse_life(:,z));
            end
        end
        if ~isempty(strfind(scannames(j).name,'choice')) && isempty(strfind(scannames(j).name,'NO_EVENTS'))
            file=dir(strcat(scannames(j).name,'/','JOR_*','.mat'));
            filename=file.name;
            load(strcat(scannames(j).name,'/',filename));
            for z=1:choicesize(2)
                all_life_corrs_withlife_forchoice(i,z)=corr(z_oxy1(:,z),lifePs_avgtimecourse_choice(:,z));
                all_life_corrs_withchoice_forchoice(i,z)=corr(z_oxy1(:,z),choicePs_avgtimecourse_choice(:,z));
            end
        end
        

    end  
    cd ..
end
for i=1:lifesize(2)
    avg_life_withlife_forlife(i)=mean(all_life_corrs_withlife_forlife(~isnan(all_life_corrs_withlife_forlife(:,i)),i));
    avg_life_withchoice_forlife(i)=mean(all_life_corrs_withchoice_forlife(~isnan(all_life_corrs_withchoice_forlife(:,i)),i));
    avg_life_withlife_forchoice(i)=mean(all_life_corrs_withlife_forchoice(~isnan(all_life_corrs_withlife_forchoice(:,i)),i));
    avg_life_withchoice_forchoice(i)=mean(all_life_corrs_withchoice_forchoice(~isnan(all_life_corrs_withchoice_forchoice(:,i)),i));
end

clear z_deoxy1 z_oxy1 z_totaloxy1 totaloxy1 oxy1 deoxy1 lifePs_avgtimecourse_life lifePs_avgtimecourse_choice choicePs_avgtimecourse_life choicePs_avgtimecourse_choice
cd ..
cd('Choice')
currdir=dir(strcat(dataprefix,'*'));
choicesize=size(lifePs_alldata_choice);
lifesize=size(lifePs_alldata_life);
clear 
lifePs_avgtimecourse_life=zeros(lifesize(1),lifesize(2));
lifePs_avgtimecourse_choice=zeros(choicesize(1),choicesize(2));

 for y=1:choicesize(2)
     for x=1:choicesize(1) 
         lifePs_avgtimecourse_choice(x,y)=mean(lifePs_alldata_choice(x,y,~isnan(lifePs_alldata_choice(x,y,:)))); 
     end
end
for y=1:lifesize(2)
     for x=1:lifesize(1) 
         lifePs_avgtimecourse_life(x,y)=mean(lifePs_alldata_life(x,y,~isnan(lifePs_alldata_life(x,y,:)))); 
     end
 end
        
for i=1:length(currdir)
    cd(currdir(i).name);
    scannames=dir(strcat(currdir(i).name,'*'));
    msg = sprintf('\n\t P Choice number %d/%d ...',i,length(currdir));
    fprintf(msg);
    
      
        choicesize=size(choicePs_alldata_choice);
        lifesize=size(choicePs_alldata_life);
        mychoicePs_alldata_life=choicePs_alldata_life;
        mychoicePs_alldata_life(:,:,j)=[];
        choicePs_avgtimecourse_life(x,y)=mean(mychoicePs_alldata_life(x,y,~isnan(mylifePs_alldata_life(x,y,:)))); 
        
        mychoicePs_alldata_choice=choicePs_alldata_choice;
        mychoicePs_alldata_choice(:,:,j)=[];
        choicePs_avgtimecourse_choice(x,y)=mean(mychoicePs_alldata_choice(x,y,~isnan(mychoicePs_alldata_choice(x,y,:)))); 
        
        
        
    
    for j=1:length(scannames)
        
        if ~isempty(strfind(scannames(j).name,'life')) && isempty(strfind(scannames(j).name,'NO_EVENTS'))
            file=dir(strcat(scannames(j).name,'/','JOR_*','.mat'));
            filename=file.name;
            load(strcat(scannames(j).name,'/',filename));
            for z=1:lifesize(2)
                all_choice_corrs_withlife_forlife(i,z)=corr(z_oxy1(:,z),lifePs_avgtimecourse_life(:,z));
                all_choice_corrs_withchoice_forlife(i,z)=corr(z_oxy1(:,z),choicePs_avgtimecourse_life(:,z));
            end
        end
        if ~isempty(strfind(scannames(j).name,'choice')) && isempty(strfind(scannames(j).name,'NO_EVENTS'))
            file=dir(strcat(scannames(j).name,'/','JOR_*','.mat'));
            filename=file.name;
            load(strcat(scannames(j).name,'/',filename));
            for z=1:choicesize(2)
                all_choice_corrs_withlife_forchoice(i,z)=corr(z_oxy1(:,z),lifePs_avgtimecourse_choice(:,z));
                all_choice_corrs_withchoice_forchoice(i,z)=corr(z_oxy1(:,z),choicePs_avgtimecourse_choice(:,z));
            end
        end
        

    end  
    cd ..
end

for i=1:lifesize(2)%note: lifesize(2) and choicesize(2) are equivalent, which is why these are all together
    avg_choice_withlife_forlife(i)=mean(all_choice_corrs_withlife_forlife(~isnan(all_choice_corrs_withlife_forlife(:,i)),i));
    avg_choice_withchoice_forlife(i)=mean(all_choice_corrs_withchoice_forlife(~isnan(all_choice_corrs_withchoice_forlife(:,i)),i));
    avg_choice_withlife_forchoice(i)=mean(all_choice_corrs_withlife_forchoice(~isnan(all_choice_corrs_withlife_forchoice(:,i)),i));
    avg_choice_withchoice_forchoice(i)=mean(all_choice_corrs_withchoice_forchoice(~isnan(all_choice_corrs_withchoice_forchoice(:,i)),i));
end

%fisherz transforms
allchoice_size=size(all_choice_corrs_withchoice_forchoice);
for x=1:allchoice_size(1)
    for y=1:allchoice_size(2)
        z_all_choice_withlife_forlife(x,y)=0.5*log((1+all_choice_corrs_withlife_forlife(x,y))/(1-all_choice_corrs_withlife_forlife(x,y)));
        z_all_choice_withchoice_forlife(x,y)=0.5*log((1+all_choice_corrs_withchoice_forlife(x,y))/(1-all_choice_corrs_withchoice_forlife(x,y)));
        z_all_choice_withlife_forchoice(x,y)=0.5*log((1+all_choice_corrs_withlife_forchoice(x,y))/(1-all_choice_corrs_withlife_forchoice(x,y)));
        z_all_choice_withchoice_forchoice(x,y)=0.5*log((1+all_choice_corrs_withchoice_forchoice(x,y))/(1-all_choice_corrs_withchoice_forchoice(x,y)));
    end
end
for i=1:allchoice_size(2)
    z_avg_choice_withlife_forlife(i)=mean(z_all_choice_withlife_forlife(~isnan(z_all_choice_withlife_forlife(:,i)),i));
    z_avg_choice_withchoice_forlife(i)=mean(z_all_choice_withchoice_forlife(~isnan(z_all_choice_withchoice_forlife(:,i)),i));
    z_avg_choice_withlife_forchoice(i)=mean(z_all_choice_withlife_forchoice(~isnan(z_all_choice_withlife_forchoice(:,i)),i));
    z_avg_choice_withchoice_forchoice(i)=mean(z_all_choice_withchoice_forchoice(~isnan(z_all_choice_withchoice_forchoice(:,i)),i));
end

alllife_size=size(all_life_corrs_withchoice_forchoice);
for x=1:alllife_size(1)
    for y=1:alllife_size(2)
        z_all_life_withlife_forlife(x,y)=0.5*log((1+all_life_corrs_withlife_forlife(x,y))/(1-all_life_corrs_withlife_forlife(x,y)));
        z_all_life_withchoice_forlife(x,y)=0.5*log((1+all_life_corrs_withchoice_forlife(x,y))/(1-all_life_corrs_withchoice_forlife(x,y)));
        z_all_life_withlife_forchoice(x,y)=0.5*log((1+all_life_corrs_withlife_forchoice(x,y))/(1-all_life_corrs_withlife_forchoice(x,y)));
        z_all_life_withchoice_forchoice(x,y)=0.5*log((1+all_life_corrs_withchoice_forchoice(x,y))/(1-all_life_corrs_withchoice_forchoice(x,y)));
    end
end
for i=1:alllife_size(2)
    z_avg_life_withlife_forlife(i)=mean(z_all_life_withlife_forlife(~isnan(z_all_life_withlife_forlife(:,i)),i));
    z_avg_life_withchoice_forlife(i)=mean(z_all_life_withchoice_forlife(~isnan(z_all_life_withchoice_forlife(:,i)),i));
    z_avg_life_withlife_forchoice(i)=mean(z_all_life_withlife_forchoice(~isnan(z_all_life_withlife_forchoice(:,i)),i));
    z_avg_life_withchoice_forchoice(i)=mean(z_all_life_withchoice_forchoice(~isnan(z_all_life_withchoice_forchoice(:,i)),i));
end

% for x=1:alllife_size(1)
%     for y=1:alllife_size(2)
%         if isnan(z_all_life_withlife_forlife(x,y))
%             z_all_life_withlife_forlife(x,y)=999;
%         end
%         if isnan(z_all_life_withchoice_forlife(x,y))
%             z_all_life_withchoice_forlife(x,y)=999;
%         end
%         if isnan(z_all_life_withlife_forchoice(x,y))
%             z_all_life_withlife_forchoice(x,y)=999;
%         end
%         if isnan(z_all_life_withchoice_forchoice(x,y))
%             z_all_life_withchoice_forchoice(x,y)=999;
%         end
%     end
%     
% end
% 
% for x=1:allchoice_size(1)
%     for y=1:allchoice_size(2)
%         if isnan(z_all_choice_withlife_forlife(x,y))
%             z_all_choice_withlife_forlife(x,y)=999;
%         end
%         if isnan(z_all_choice_withchoice_forlife(x,y))
%             z_all_choice_withchoice_forlife(x,y)=999;
%         end
%         if isnan(z_all_choice_withlife_forchoice(x,y))
%             z_all_choice_withlife_forchoice(x,y)=999;
%         end
%         if isnan(z_all_choice_withchoice_forchoice(x,y))
%             z_all_choice_withchoice_forchoice(x,y)=999;
%         end
%     end
%     
% end
% 
% 
