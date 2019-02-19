%This script makes nirs images for viewing in surfICE

%You need to download the full version of xjview2 and add it to your path

%set mymni to whatever your mni coordinates are

%set a variable that contains your values to be put into the coordinates
%(in order)

mymni=[   -46    42    21;-30    44    39;-46    48     0; -32    62    -8;-24    55    31
   -25    66     4
   -24    65    21
   -10    44    48
     2    54    38
    12    44    48
   -19    70    -5
     3    70    11
    14    70    -5
    16    64    20
    26    55    31
    28    66     4
    32    43    40
    48    42    22
    35    63    -8
    49    48     1]; % may need to nudged if ROIs are being clipped...

combinedwatchingsame_same_v_other=[4.56884902500000;4.29958553600000;3.57699099700000;3.95554116100000;4.56072042400000;6.55615833900000;5.52111535300000;6.35065886300000;5.09096148100000;5.56921270900000;4.12962625800000;5.72825455600000;3.10360305700000;5.25786024300000;5.66279644700000;4.91344741300000;7.12625211300000;6.24841246000000;4.48116257500000;4.82628499700000];

combinedwatchingother_same_v_other=[   5.0792
    6.1610
    4.8746
    5.2384
    3.9252
    6.7479
    4.9217
    6.6501
    5.3814
    6.8606
    5.6660
    6.1855
    4.8688
    5.0212
    4.0823
    6.2053
    5.1798
    4.3244
    4.8097
    5.5291];


nirs2img('combined_watchingsame_same_v_other_trial.img',mymni,combined_watchingsame_same_v_other,0,0,0)
nirs2img('combined_watchingother_same_v_other_trial.img',mymni,combinedwatchingother_same_v_other,0,0,0)