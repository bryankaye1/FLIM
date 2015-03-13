fid = fopen('Z:\bkaye\cluster\junk.txt','r');     %# Open the file as a binary
lastLine = '';                   %# Initialize to empty
offset = 1;                      %# Offset from the end of file
fseek(fid,-offset,'eof');        %# Seek to the file end, minus the offset
newChar = fread(fid,1,'*char');  %# Read one character
while (~strcmp(newChar,char(10))) || (offset == 1)
  lastLine = [newChar lastLine];   %# Add the character to a string
  offset = offset+1;
  fseek(fid,-offset,'eof');        %# Seek to the file end, minus the offset
  newChar = fread(fid,1,'*char');  %# Read one character
end
fclose(fid);  %# Close the file
if lastLine(end)=='r'%checks if the last line has an "r" or nor
oldrun = str2num(lastLine(1:end-1)); % Reads the number, not including the r
else
oldrun = str2num(lastLine(1:end)); % Reads the number
end
newrun = oldrun + 1;
fileID = fopen('Z:\bkaye\cluster\junk.txt','a');
fprintf(fileID,'hello world');%,num2str(newrun));
fclose(fileID);  %# Close the file