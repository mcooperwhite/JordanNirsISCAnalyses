#import the dataset from matlab file
library(R.matlab)
data<-readMat('/Users/macrinacooper-white/Dropbox/MINERVA3/Data/Jordan/PreProcessedFiles_cut/variables_for_analysis_updated.mat') #FYI I think these are z-scores of oxygenated Hb. They were Fisher-z-transformed in order to be compared

#set variables from matlab variables
life_withchoice_forchoice<-data[["z.all.life.withchoice.forchoice"]]
life_withchoice_forlife<-data[["z.all.life.withchoice.forlife"]]
life_withlife_forchoice<-data[["z.all.life.withlife.forchoice"]]
life_withlife_forlife<-data[["z.all.life.withlife.forlife"]]
choice_withchoice_forchoice<-data[["z.all.choice.withchoice.forchoice"]]
choice_withchoice_forlife<-data[["z.all.choice.withchoice.forlife"]]
choice_withlife_forchoice<-data[["z.all.choice.withlife.forchoice"]]
choice_withlife_forlife<-data[["z.all.choice.withlife.forlife"]]

#reorganize the data so easily comparable in paired t-tests
l_withc_forc_compare <- lapply (1 : 20, function (x)life_withchoice_forchoice[,x])
l_withc_forl_compare <- lapply (1 : 20, function (x)life_withchoice_forlife[,x])
l_withl_forc_compare <- lapply (1 : 20, function (x)life_withlife_forchoice[,x])
l_withl_forl_compare <- lapply (1 : 20, function (x)life_withlife_forlife[,x])
c_withc_forc_compare <- lapply (1 : 20, function (x)choice_withchoice_forchoice[,x])
c_withc_forl_compare <- lapply (1 : 20, function (x)choice_withchoice_forlife[,x])
c_withl_forc_compare <- lapply (1 : 20, function (x)choice_withlife_forchoice[,x])
c_withl_forl_compare <- lapply (1 : 20, function (x)choice_withlife_forlife[,x])

#conduct paired t-tests in all of the channels and store results
life_forlife_compare_ltoc <- lapply (1 : 20, function (x)t.test(l_withl_forl_compare[[x]],l_withc_forl_compare[[x]],paired=TRUE))
life_forchoice_compare_ltoc <- lapply (1 : 20, function (x)t.test(l_withl_forc_compare[[x]],l_withc_forc_compare[[x]],paired=TRUE))
choice_forlife_compare_ctol <- lapply (1 : 20, function (x)t.test(c_withc_forl_compare[[x]],c_withl_forl_compare[[x]],paired=TRUE))
choice_forchoice_compare_ctol <- lapply (1 : 20, function (x)t.test(c_withc_forc_compare[[x]],c_withl_forc_compare[[x]],paired=TRUE))

#Print out p-value results for each of the comparison (organized with one p-value corresponding to each channel for chans 1-20)
print("Current comparison: Life for choice - life to choice")
for (i in 1:20){print(life_forchoice_compare_ltoc[[i]]$p.value)}

print("Current comparison: Life for life - life to choice")
for (i in 1:20){print(life_forlife_compare_ltoc[[i]]$p.value)}

print("Current comparison: Choice for life - choice to life")
for (i in 1:20){print(choice_forlife_compare_ctol[[i]]$p.value)}

print("Current comparison: Choice for choice - choice to life")
for (i in 1:20){print(choice_forchoice_compare_ctol[[i]]$p.value)}

#store all in one dataset
allpvalues<-matrix(,nrow=20,ncol=4)
alltvalues<-matrix(,nrow=20,ncol=4)
alldiffs<-matrix(,nrow=20,ncol=4)
allpvalues_ostt<-matrix(,nrow=20,ncol=8)
alltvalues_ostt<-matrix(,nrow=20,ncol=8)
allzvalues_ostt<-matrix(,nrow=20,ncol=8)
#life_withlife_forchoice<-matrix(,nrow=20,ncol=4)
#life_withlife_forlife<-matrix(,nrow=20,ncol=4)
#life_withchoice_forchoice<-matrix(,nrow=20,ncol=4)
#life_withchoice_forlife<-matrix(,nrow=20,ncol=4)
#choice_withlife_forchoice<-matrix(,nrow=20,ncol=4)
#choice_withlife_forlife<-matrix(,nrow=20,ncol=4)
#choice_withchoice_forchoice<-matrix(,nrow=20,ncol=4)
#choice_withchoice_forlife<-matrix(,nrow=20,ncol=4)

for (i in 1:20){allpvalues[[i,1]]<-(life_forchoice_compare_ltoc[[i]]$p.value);allpvalues[[i,2]]<-(life_forlife_compare_ltoc[[i]]$p.value);allpvalues[[i,3]]<-(choice_forlife_compare_ctol[[i]]$p.value);allpvalues[[i,4]]<-(choice_forchoice_compare_ctol[[i]]$p.value)}
for (i in 1:20){alldiffs[[i,1]]<-(life_forchoice_compare_ltoc[[i]]$estimate[[1]]);alldiffs[[i,2]]<-(life_forlife_compare_ltoc[[i]]$estimate[[1]]);alldiffs[[i,3]]<-(choice_forlife_compare_ctol[[i]]$estimate[[1]]);alldiffs[[i,4]]<-(choice_forchoice_compare_ctol[[i]]$estimate[[1]])}
for (i in 1:20){alltvalues[[i,1]]<-(life_forchoice_compare_ltoc[[i]]$statistic[[1]]);alltvalues[[i,2]]<-(life_forlife_compare_ltoc[[i]]$statistic[[1]]);alltvalues[[i,3]]<-(choice_forlife_compare_ctol[[i]]$statistic[[1]]);alltvalues[[i,4]]<-(choice_forchoice_compare_ctol[[i]]$statistic[[1]])}

#for (i in 1:20){alldiffs[[i,1]]<-(life_forchoice_compare_ltoc[[i]]$estimate[[1]]-life_forchoice_compare_ltoc[[i]]$estimate[[2]]);alldiffs[[i,2]]<-(life_forlife_compare_ltoc[[i]]$estimate[[1]]-life_forlife_compare_ltoc[[i]]$estimate[[2]]);alldiffs[[i,3]]<-(choice_forlife_compare_ctol[[i]]$estimate[[1]]-choice_forlife_compare_ctol[[i]]$estimate[[2]]);alldiffs[[i,4]]<-(choice_forchoice_compare_ctol[[i]]$estimate[[1]]-choice_forchoice_compare_ctol[[i]]$estimate[[2]])}

life_withlife_life_ostt <- lapply (1 : 20, function (x)t.test(l_withl_forl_compare[[x]]))
life_withlife_choice_ostt <- lapply (1 : 20, function (x)t.test(l_withl_forc_compare[[x]]))
life_withchoice_life_ostt <- lapply (1 : 20, function (x)t.test(l_withc_forl_compare[[x]]))
life_withchoice_choice_ostt <- lapply (1 : 20, function (x)t.test(l_withc_forc_compare[[x]]))
choice_withlife_life_ostt <- lapply (1 : 20, function (x)t.test(c_withl_forl_compare[[x]]))
choice_withlife_choice_ostt <- lapply (1 : 20, function (x)t.test(c_withl_forc_compare[[x]]))
choice_withchoice_life_ostt <- lapply (1 : 20, function (x)t.test(c_withc_forl_compare[[x]]))
choice_withchoice_choice_ostt <- lapply (1 : 20, function (x)t.test(c_withc_forc_compare[[x]]))

for (i in 1:20){allpvalues_ostt[[i,1]]<-(life_withlife_life_ostt[[i]]$p.value);allpvalues_ostt[[i,2]]<-(life_withlife_choice_ostt[[i]]$p.value);allpvalues_ostt[[i,3]]<-(life_withchoice_life_ostt[[i]]$p.value);allpvalues_ostt[[i,4]]<-(life_withchoice_choice_ostt[[i]]$p.value);allpvalues_ostt[[i,5]]<-(choice_withlife_life_ostt[[i]]$p.value);allpvalues_ostt[[i,6]]<-(choice_withlife_choice_ostt[[i]]$p.value);allpvalues_ostt[[i,7]]<-(choice_withchoice_life_ostt[[i]]$p.value);allpvalues_ostt[[i,8]]<-(choice_withchoice_choice_ostt[[i]]$p.value)}
for (i in 1:20){allzvalues_ostt[[i,1]]<-(life_withlife_life_ostt[[i]]$estimate[[1]]);allzvalues_ostt[[i,2]]<-(life_withlife_choice_ostt[[i]]$estimate[[1]]);allzvalues_ostt[[i,3]]<-(life_withchoice_life_ostt[[i]]$estimate[[1]]);allzvalues_ostt[[i,4]]<-(life_withchoice_choice_ostt[[i]]$estimate[[1]]);allzvalues_ostt[[i,5]]<-(choice_withlife_life_ostt[[i]]$estimate[[1]]);allzvalues_ostt[[i,6]]<-(choice_withlife_choice_ostt[[i]]$estimate[[1]]);allzvalues_ostt[[i,7]]<-(choice_withchoice_life_ostt[[i]]$estimate[[1]]);allzvalues_ostt[[i,8]]<-(choice_withchoice_choice_ostt[[i]]$estimate[[1]])}
for (i in 1:20){alltvalues_ostt[[i,1]]<-(life_withlife_life_ostt[[i]]$statistic[[1]]);alltvalues_ostt[[i,2]]<-(life_withlife_choice_ostt[[i]]$statistic[[1]]);alltvalues_ostt[[i,3]]<-(life_withchoice_life_ostt[[i]]$statistic[[1]]);alltvalues_ostt[[i,4]]<-(life_withchoice_choice_ostt[[i]]$statistic[[1]]);alltvalues_ostt[[i,5]]<-(choice_withlife_life_ostt[[i]]$statistic[[1]]);alltvalues_ostt[[i,6]]<-(choice_withlife_choice_ostt[[i]]$statistic[[1]]);alltvalues_ostt[[i,7]]<-(choice_withchoice_life_ostt[[i]]$statistic[[1]]);alltvalues_ostt[[i,8]]<-(choice_withchoice_choice_ostt[[i]]$statistic[[1]])}



#Use Benjamini-Hochman procedure to correct for multiple comparisons (FDR correction)
allpvalues_corrected<-matrix(,nrow=20,ncol=4)
allpvalues_corrected[,1]<-p.adjust(allpvalues[,1],"BH")
allpvalues_corrected[,2]<-p.adjust(allpvalues[,2],"BH")
allpvalues_corrected[,3]<-p.adjust(allpvalues[,3],"BH")
allpvalues_corrected[,4]<-p.adjust(allpvalues[,4],"BH")

allpvalues_corrected_ostt<-matrix(,nrow=20,ncol=8)
allpvalues_corrected_ostt[,1]<-p.adjust(allpvalues_ostt[,1],"BH")
allpvalues_corrected_ostt[,2]<-p.adjust(allpvalues_ostt[,2],"BH")
allpvalues_corrected_ostt[,3]<-p.adjust(allpvalues_ostt[,3],"BH")
allpvalues_corrected_ostt[,4]<-p.adjust(allpvalues_ostt[,4],"BH")
allpvalues_corrected_ostt[,5]<-p.adjust(allpvalues_ostt[,5],"BH")
allpvalues_corrected_ostt[,6]<-p.adjust(allpvalues_ostt[,6],"BH")
allpvalues_corrected_ostt[,7]<-p.adjust(allpvalues_ostt[,7],"BH")
allpvalues_corrected_ostt[,8]<-p.adjust(allpvalues_ostt[,8],"BH")



#combine Ps (across groups) based on whether they were watching same or opposite video
combinedsame_watchingsame<-matrix(nrow=63,ncol=20)
for (i in 1:20){combinedsame_watchingsame[,i]<-c(choice_withchoice_forchoice[,i],life_withlife_forlife[,i])}

combinedother_watchingsame<-matrix(nrow=63,ncol=20)
for (i in 1:20){combinedother_watchingsame[,i]<-c(choice_withlife_forchoice[,i],life_withchoice_forlife[,i])}

combinedsame_watchingother<-matrix(nrow=63,ncol=20)
for (i in 1:20){combinedsame_watchingother[,i]<-c(choice_withchoice_forlife[,i],life_withlife_forchoice[,i])}

combinedother_watchingother<-matrix(nrow=63,ncol=20)
for (i in 1:20){combinedother_watchingother[,i]<-c(choice_withlife_forlife[,i],life_withchoice_forchoice[,i])}


watchingsame_same_v_other <- lapply (1 : 20, function (x)t.test(combinedsame_watchingsame[,x],combinedother_watchingsame[,x],paired=TRUE))
watchingother_same_v_other <- lapply (1 : 20, function (x)t.test(combinedsame_watchingother[,x],combinedother_watchingother[,x],paired=TRUE))
same_watchingsame_v_watchingother<-lapply (1 : 20, function (x)t.test(combinedsame_watchingsame[,x],combinedsame_watchingother[,x],paired=TRUE))
other_watchingsame_v_watchingother<-lapply (1 : 20, function (x)t.test(combinedother_watchingsame[,x],combinedother_watchingother[,x],paired=TRUE))

alltvalues_combinedcons<-matrix(,nrow=20,ncol=2)
alldiffs_combinedcons<-matrix(,nrow=20,ncol=2)
allpvalues_combinedcons<-matrix(,nrow=20,ncol=2)

for (i in 1:20){alltvalues_combinedcons[[i,1]]<-(watchingsame_same_v_other[[i]]$statistic[[1]]);alltvalues_combinedcons[[i,2]]<-(watchingother_same_v_other[[i]]$statistic[[1]]);alltvalues_combinedcons[[i,3]]<-(same_watchingsame_v_watchingother[[i]]$statistic[[1]]);alltvalues_combinedcons[[i,4]]<-(other_watchingsame_v_watchingother[[i]]$statistic[[1]])}
for (i in 1:20){alldiffs_combinedcons[[i,1]]<-(watchingsame_same_v_other[[i]]$estimate[[1]]);alldiffs_combinedcons[[i,2]]<-(watchingother_same_v_other[[i]]$estimate[[1]]);alltvalues_combinedcons[[i,3]]<-(same_watchingsame_v_watchingother[[i]]$statistic[[1]]);alltvalues_combinedcons[[i,4]]<-(other_watchingsame_v_watchingother[[i]]$statistic[[1]])}
for (i in 1:20){allpvalues_combinedcons[[i,1]]<-(watchingsame_same_v_other[[i]]$p.value[[1]]);allpvalues_combinedcons[[i,2]]<-(watchingother_same_v_other[[i]]$p.value[[1]]);alltvalues_combinedcons[[i,3]]<-(same_watchingsame_v_watchingother[[i]]$statistic[[1]]);alltvalues_combinedcons[[i,4]]<-(other_watchingsame_v_watchingother[[i]]$statistic[[1]])}

combinedsame_watchingsame_ostt <- lapply (1 : 20, function (x)t.test(combinedsame_watchingsame[[x]]))
combinedsame_watchingother_ostt <- lapply (1 : 20, function (x)t.test(combinedsame_watchingother[[x]]))
combinedother_watchingsame_ostt <- lapply (1 : 20, function (x)t.test(combinedother_watchingsame[[x]]))
combinedother_watchingother_ostt <- lapply (1 : 20, function (x)t.test(combinedother_watchingother[[x]]))

for (i in 1:20){allpvalues_combined_ostt[[i,1]]<-(combinedsame_watchingsame_ostt[[i]]$p.value);allpvalues_combined_ostt[[i,2]]<-(combinedsame_watchingother_ostt[[i]]$p.value);allpvalues_combined_ostt[[i,3]]<-(combinedother_watchingsame_ostt[[i]]$p.value);allpvalues_combined_ostt[[i,4]]<-(combinedother_watchingother_ostt[[i]]$p.value)}
for (i in 1:20){allzvalues_combined_ostt[[i,1]]<-(combinedsame_watchingsame_ostt[[i]]$estimate[[1]]);allzvalues_combined_ostt[[i,2]]<-(combinedsame_watchingother_ostt[[i]]$estimate[[1]]);allzvalues_combined_ostt[[i,3]]<-(combinedother_watchingsame_ostt[[i]]$estimate[[1]]);allzvalues_combined_ostt[[i,4]]<-(combinedother_watchingother_ostt[[i]]$estimate[[1]])}
for (i in 1:20){alltvalues_combined_ostt[[i,1]]<-(combinedsame_watchingsame_ostt[[i]]$statistic[[1]]);alltvalues_combined_ostt[[i,2]]<-(combinedsame_watchingother_ostt[[i]]$statistic[[1]]);alltvalues_combined_ostt[[i,3]]<-(combinedother_watchingsame_ostt[[i]]$statistic[[1]]);alltvalues_combined_ostt[[i,4]]<-(combinedother_watchingother_ostt[[i]]$statistic[[1]])}

allpvalues_corrected_comb_ostt<-matrix(,nrow=20,ncol=4)
allpvalues_corrected_comb_ostt[,1]<-p.adjust(allpvalues_ostt[,1],"BH")
allpvalues_corrected_comb_ostt[,2]<-p.adjust(allpvalues_ostt[,2],"BH")
allpvalues_corrected_comb_ostt[,3]<-p.adjust(allpvalues_ostt[,3],"BH")
allpvalues_corrected_comb_ostt[,4]<-p.adjust(allpvalues_ostt[,4],"BH")

allpvalues_corrected_c<-matrix(,nrow=20,ncol=4)
allpvalues_corrected_c[,1]<-p.adjust(allpvalues_combinedcons[,1],"BH")
allpvalues_corrected_c[,2]<-p.adjust(allpvalues_combinedcons[,2],"BH")
allpvalues_corrected_c[,3]<-p.adjust(allpvalues_combinedcons[,3],"BH")
allpvalues_corrected_c[,4]<-p.adjust(allpvalues_combinedcons[,4],"BH")

write.csv(alltvalues_combinedcons, "/Users/macrinacooper-white/Dropbox/MINERVA3/Data/Jordan/PreProcessedFiles_cut/ISC analyses/combined_tstt-test-t-values_020319.csv")
write.csv(alldiffs_combinedcons, "/Users/macrinacooper-white/Dropbox/MINERVA3/Data/Jordan/PreProcessedFiles_cut/ISC analyses/combined_tstt-test-all-difference-scores_020319.csv")
write.csv(allpvalues_corrected_c, "/Users/macrinacooper-white/Dropbox/MINERVA3/Data/Jordan/PreProcessedFiles_cut/ISC analyses/combined_tstt-test-p-values_020319.csv")
write.csv(allpvalues_corrected_comb_ostt, "/Users/macrinacooper-white/Dropbox/MINERVA3/Data/Jordan/PreProcessedFiles_cut/ISC analyses/combined_ostt-test-p-values_020319.csv")
write.csv(allzvalues_combined_ostt, "/Users/macrinacooper-white/Dropbox/MINERVA3/Data/Jordan/PreProcessedFiles_cut/ISC analyses/ostt-combined_test-z-values_020319.csv")
write.csv(alltvalues_combined_ostt, "/Users/macrinacooper-white/Dropbox/MINERVA3/Data/Jordan/PreProcessedFiles_cut/ISC analyses/ostt-combined_test-t-values_020319.csv")


#save the p-values to an excel file
write.csv(allpvalues_corrected, "/Users/macrinacooper-white/Dropbox/MINERVA3/Data/Jordan/PreProcessedFiles_cut/ISC analyses/tstt-test-p-values_010319.csv")
write.csv(alldiffs, "/Users/macrinacooper-white/Dropbox/MINERVA3/Data/Jordan/PreProcessedFiles_cut/ISC analyses/tstt-test-all-difference-scores_020319.csv")
write.csv(alltvalues, "/Users/macrinacooper-white/Dropbox/MINERVA3/Data/Jordan/PreProcessedFiles_cut/ISC analyses/tstt-test-t-values_020319.csv")

write.csv(alltvalues_ostt, "/Users/macrinacooper-white/Dropbox/MINERVA3/Data/Jordan/PreProcessedFiles_cut/ISC analyses/ostt-test-t-values_020319.csv")
write.csv(allzvalues_ostt, "/Users/macrinacooper-white/Dropbox/MINERVA3/Data/Jordan/PreProcessedFiles_cut/ISC analyses/ostt-test-z-values_020319.csv")
write.csv(allpvalues_corrected_ostt, "/Users/macrinacooper-white/Dropbox/MINERVA3/Data/Jordan/PreProcessedFiles_cut/ISC analyses/ostt-test-p-values_020319.csv")


#making graphs
for (i in 1:20){assign(paste0("life_same_choice",i),c(l_withl_forc_compare[[i]]))}
for (i in 1:20){assign(paste0("life_other_choice",i),c(l_withc_forc_compare[[i]]))}
for (i in 1:20){assign(paste0("choice_same_choice",i),c(c_withc_forc_compare[[i]]))}
for (i in 1:20){assign(paste0("choice_other_choice",i),c(c_withl_forc_compare[[i]]))}
for (i in 1:20){assign(paste0("life_same_life",i),c(l_withl_forl_compare[[i]]))}
for (i in 1:20){assign(paste0("life_other_life",i),c(l_withc_forl_compare[[i]]))}
for (i in 1:20){assign(paste0("choice_same_life",i),c(c_withc_forl_compare[[i]]))}
for (i in 1:20){assign(paste0("choice_other_life",i),c(c_withl_forl_compare[[i]]))}

for (i in 1:20){assign(paste0("life_forchoice",i),data.frame(group=rep(c("same group","other group"),each=31),ISC=c(get(paste0("life_same_choice",i)),get(paste0("life_other_choice",i)))))}
for (i in 1:20){assign(paste0("life_forlife",i),data.frame(group=rep(c("same group","other group"),each=30),ISC=c(get(paste0("life_same_life",i)),get(paste0("life_other_life",i)))))}
for (i in 1:20){assign(paste0("choice_forlife",i),data.frame(group=rep(c("same group","other group"),each=30),ISC=c(get(paste0("choice_same_choice",i)),get(paste0("choice_other_choice",i)))))}
for (i in 1:20){assign(paste0("choice_forchoice",i),data.frame(group=rep(c("same group","other group"),each=31),ISC=c(get(paste0("choice_same_life",i)),get(paste0("choice_other_life",i)))))}

#library("ggpubr")
#ggboxplot(life_forchoice1, x = "group", y = "ISC", 
#          color = "group", palette = c("#00AFBB", "#E7B800"),
#          order = c("before", "after"),
#          ylab = "ISC", xlab = "Group")

samegroup <- c(life_same_choice9)
othergroup <- c(life_other_choice9)
m        = c(mean(samegroup,na.rm = TRUE), mean(othergroup,na.rm = TRUE))
names(m) = c("same group", "other group")
CIlower       = c(life_withlife_choice_ostt[[9]]$conf.int[[1]], 
                  life_withchoice_choice_ostt[[9]]$conf.int[[1]])
CIupper       = c(life_withlife_choice_ostt[[9]]$conf.int[[2]], 
                  life_withchoice_choice_ostt[[9]]$conf.int[[2]])
windows()
bp = barplot(m, ylim=c(0, 0.5), xpd=FALSE, col=c("dark blue","red"))
box()
arrows(x0=bp, y0=CIlower, y1=CIupper, code=3, angle=90)