data_root = 'E:\LungIMCData';
% load communities
raw_community_path = fullfile(data_root, 'BatchCorrection', 'RData', 'RawCommunities.txt');
raw_list = importdata(raw_community_path);
correct_community_path = fullfile(data_root, 'BatchCorrection', 'RData', 'CorrectCommunities.txt');
correct_list = importdata(correct_community_path);
cc = spx.cluster.ClusterComparison(raw_list, correct_list);
fm = cc.fMeasure();
spx.cluster.ClusterComparison.printF1MeasureResult(fm);
