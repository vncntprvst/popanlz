function [normData,normFactor]=pop_a_normalization(data, epochs, sigma) %,stdata,conn)
%get normalization factor:
% For each neuron, find peak response for each designed epoch.
% e.g.,  visual epoch 0 to 250 ms from stimulus onset
%         motor epoch -300 to 300 ms around saccade onset
% The largest value is kept as the neuron's normalization factor.

% input data is a n by p cell array, with n neurons and p alignments
% each cell is structure with the following fields:
% rast: n by p array of n trials and p time points, with binary  spike train data (0=no spike, 1=spike)
% alignt: alignment time

% epochs is a n by 1 cell array with n epochs boundary times (default [-300 300])
% if n is inferior to data's p alignements, extra alignments will not be
% considered for calulating the normalization factor. 
% However, normalization is applied accross all provided data.

% global directory slash;
% if isempty(directory)
%     [directory,slash]=SetUserDir;
% end

switch nargin
    case 1
    epochs={[300 300]};
    sigma=15;
    case 2
    sigma=15;
end

% cellnum=sum(~cellfun('isempty',gsdata.allsacdelay));
% st.cellnum=sum(~cellfun('isempty',stdata.allsacdelay));
% 
% disp([num2str(cellnum) ' cells for gs, '  num2str(st.cellnum) ' cells for st'])
% 
% queries{1} = ['SELECT a_file FROM sorts s INNER JOIN recordings r on s.recording_id_fk = r.recording_id WHERE sort_id IN (' ...
%     sprintf('%.0f,' ,cellfun(@(x) x.sort_id,gsdata.alldb(1:end-1,1))) num2str(gsdata.alldb{end,1}.sort_id) ')'];
% % recnames=fetch(conn,query);
% 
% queries{2} = ['SELECT a_file FROM sorts s INNER JOIN recordings r on s.recording_id_fk = r.recording_id WHERE sort_id IN (' ...
%     sprintf('%.0f,' ,cellfun(@(x) x.sort_id,stdata.alldb(1:end-1,1))) num2str(stdata.alldb{end,1}.sort_id) ')'];
% % st.recnames=fetch(conn,query);
% 
% [~,fileidx]=compare_db_filelists(queries,conn);
% 
% commoncells=find(fileidx{1});

%% find recordings with one trial
properLengthRecs=logical(cellfun(@(x) size(x(1,1).rast,1)>1, data(:,1)));

%% Find normalization factor 
% Convolve traces around epochs and get. 
% sacActivity=cellfun(@(x) conv_raster(x(1,1).rast,sigma,x(1,1).alignt-(300+sigma*3),x(1,1).alignt+(299+sigma*3)), data.allndata(:,1), 'UniformOutput',false); %600ms epoch
% tgtActivity=cellfun(@(x) conv_raster(x(1,1).rast,sigma,x(1,1).alignt-sigma*3,x(1,1).alignt+(249+sigma*3)), data.allndata(:,2), 'UniformOutput',false); % 250ms epoch
normFactor=nan(size(properLengthRecs,1),size(epochs,1)+1);
normFactor(properLengthRecs,:)=FindNormFactor(data(properLengthRecs,1:size(epochs,1)), epochs);

%% Normalize data
normData=RespNormalization(data, normFactor);




