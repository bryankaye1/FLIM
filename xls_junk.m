h1 = 'header 1';
h2 = 'header 2';

d1 = 1;
d2 = 22;
data = {'Category 1', 'Category 2'; 14,34; 12,78; 456, 97;};

file_name = 'FLIM Analysis.xls';
xlswrite(file_name, data);





[~,~,alldata] = xlsread(file_name);
newdata = {'it worked', 'it really worked'};
concat = {alldata ; newdata};
file_name2 = 'testfile.xls';
xlswrite(file_name2, concat);

%fseek(fileID, offset, origin)

