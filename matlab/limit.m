% Total protocol time: 1
% True offset: 46.123
% Time resolution: 0.014901
% Offset: -46.1191
% Coincidence frac: 0.003003
% S: 11.2854
% error probability: 0
% Desired time accuracy: 0.00033356
% Accuracy: 0.0039061
% FAILS: 0
% Mean S: 10.9719
% 
% Total protocol time: 0.9
% True offset: 46.123
% Time resolution: 0.013411
% Offset: -46.1206
% Coincidence frac: 0.0030037
% S: 12.6114
% error probability: 0
% Desired time accuracy: 0.00033356
% Accuracy: 0.002416
% FAILS: 0
% Mean S: 12.4701
% 
% Total protocol time: 0.8
% True offset: 46.123
% Time resolution: 0.011921
% Offset: -46.1221
% Coincidence frac: 0.003002
% S: 13.2373
% error probability: 0
% Desired time accuracy: 0.00033356
% Accuracy: 0.00092587
% Mean S: 13.7
% 
% 
% Total protocol time: 0.5
% True offset: 46.123
% Time resolution: 0.0074506
% Offset: 26.2633
% Coincidence frac: 0.003001
% S: 3.7388
% error probability: 6204.4354
% Desired time accuracy: 0.00033356
% Accuracy: 72.3863
% FAILS :10
% Mean S: 3.8
% 
% Total protocol time: 0.6
% True offset: 46.123
% Time resolution: 0.0089407
% Offset: -25.6777
% Coincidence frac: 0.0030037
% S: 4.1417
% error probability: 1156.9081
% Desired time accuracy: 0.00033356
% Accuracy: 20.4453
% FAILS: 10
% Mean S: 3.9159
% 
% Total protocol time: 0.7
% True offset: 46.123
% Time resolution: 0.010431
% Offset: -46.1251
% Coincidence frac: 0.0030027
% S: 13.2036
% error probability: 0
% Desired time accuracy: 0.00033356
% Accuracy: 0.0020544
% FAILS: 0
% Mean S: 12.181
% 
% 
% 
% Total protocol time: 0.65
% True offset: 46.123
% Time resolution: 0.0096858
% Offset: -46.1236
% Coincidence frac: 0.0030032
% S: 13.152
% error probability: 0
% Desired time accuracy: 0.00033356
% Accuracy: 0.00056424
% FAILS: 0
% Mean S: 14.2288
% 
% Total protocol time: 0.64
% True offset: 46.123
% Time resolution: 0.0095367
% Offset: -46.1197
% Coincidence frac: 0.0030032
% S: 10.6783
% error probability: 0
% Desired time accuracy: 0.00033356
% Accuracy: 0.0033101
% FAILS: 0
% Mean S: 10.2723
% 
% Total protocol time: 0.63
% True offset: 46.123
% Time resolution: 0.0093877
% Offset: -46.1219
% Coincidence frac: 0.0030033
% S: 15.2367
% error probability: 0
% Desired time accuracy: 0.00033356
% Accuracy: 0.0010749
% FAILS: 0
% Mean S: 13.2712
% 
% Total protocol time: 0.62
% True offset: 46.123
% Time resolution: 0.0092387
% Offset: -46.1197
% Coincidence frac: 0.0030034
% S: 10.7681
% error probability: 0
% Desired time accuracy: 0.00033356
% Accuracy: 0.0033101
% FAILS: 0
% Mean S: 9.8614
% 
% Total protocol time: 0.61
% True offset: 46.123
% Time resolution: 0.0090897
% Offset: -41.7127
% Coincidence frac: 0.0030036
% S: 3.9469
% error probability: 2656.0326
% Desired time accuracy: 0.00033356
% Accuracy: 4.4103
% FAILS: 10
% Mean S: 3.782
% 
% Total protocol time: 0.615
% True offset: 46.123
% Time resolution: 0.0091642
% Offset: -4.8754
% Coincidence frac: 0.0030019
% S: 4.4213
% error probability: 329.1526
% Desired time accuracy: 0.00033356
% Accuracy: 41.2476
% FAILS: 10
% Mean S: 3.9787


% times = [1, .9, .8, .7, .6, .5, .61, .62, .63, .64, .65, .615];
% Ss = [10.9719, 12.4701, 13.7, 12.181, 3.9159, 3.8, 3.782 9.8614 13.2712 10.2723 14.2288 3.9787];


times = [.5, .6, .61, .615, .62, .63, .64, .65, .7, .8, .9, 1];
Ss = [3.8, 3.9159, 3.782, 3.9787, 9.8614, 13.2712, 10.2723, 14.2288, 12.181, 13.7, 12.4701, 10.9719];




plot(times, Ss, '-o', 'LineWidth',1.5);
xlabel('Total link time (sec)');
ylabel('S');
