function fullfilename = EELS_writeraw(III,varargin)
% function fullfilename = writeraw(II,filename)
% writeraw, writeraw(II) or writeraw(II,filename)
% Based on writeraw MATLAB code by Huolin Xin
% Modified by Barnaby Levin

if nargin==0
    varlist = evalin('base','who');
    varselnum = listdlg('ListString',varlist);
    varname = varlist{varselnum};
    II = evalin('base',varname); 
    filename = inputdlg('Input the name for this file (w/o suffix):');
    filename = filename{1};    
elseif nargin==1
    filename = inputdlg('Input the name for this file (w/o suffix):');
    filename = filename{1};
else
    filename = varargin{1};
end

II=permute(III,[2 1 3]);

[pathstr, name] = fileparts(filename);
if isempty(pathstr)
    dirname = pwd;
else
    dirname = '';
end
suffix = [];
matsize = size(II);
for i=1:length(matsize)
    if i==1
        suffix = num2str(matsize(i));
    else
        suffix = [suffix,'x',num2str(matsize(i))];
    end
end

fullfilename=fullfile(dirname,[filename,'_',suffix,'_32BitReal.dat']);

fp=fopen(fullfilename,'w');
fwrite(fp,II,'float'); % This will output 32 bit data.
fclose(fp);

