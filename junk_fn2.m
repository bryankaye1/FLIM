    




            pth_imap = 'C:\Users\Bryan\Documents\MATLAB\Bayes_2014\data\2015-2-14\';
           fn_imap = 'atto565-10um.sdt';

    %pth_imap = varargin{1};
   % fn_imap = varargin{2};
    data=dir([pth_imap fn_imap]); %does this line do anything?
    block=1;
    sdt = bh_readsetup([pth_imap data.name]);
    ch = bh_getdatablock(sdt,block); %raw lifetime data
    single_imap = ch; %raw lifetime data in single form (normally stored as 4 bit number)
    imaphis = double(single_imap);
    imap = squeeze(sum(imaphis,1));
    
    figure(3);clf; imshow(imap,[]);


figure(4); clf; imshow(smoothn(imap,150),[]);
