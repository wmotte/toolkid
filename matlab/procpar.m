% Name : procpar
%
% Purpose : To read parameters from the procpar file belonging to an fid
%           obtained using Varian spectrometer.
%
% Calling sequence : param=read_procpar(filename,parameter);
%
% Input : filename : Name of directory containing both the fid and procpar
%               data (without the .fid extension)
%         parameter: The name of the parameter of which you want to know
%         the value(s).
%
% Output : param : The parameter values.
%
% Essential subroutines : none
%
% Modification history : Author : Annette van der Toorn 26-10-2004

function param=procpar(filename,parameter);

procparname=strcat(filename,'.fid/procpar');

% Open procpar file

fileid=fopen(procparname);
procpar=fscanf(fileid,'%c',[inf]);
fclose(fileid);
addwhite=' ';
parameter=[parameter addwhite];
position=findstr(parameter,procpar);
if length(position) == 0;
    param=[];
else;
    szpos=size(position);
    if szpos(2)>1
        check=isspace(procpar(position-1));
        index=find(check ==1);
        position=position(index);
    end;
    sz=size(procpar);
    procpar=procpar(position:sz(2));
    letters=isletter(procpar);
    spaces=isspace(procpar);
    posspaces=find(spaces==1);
    basictype=str2num(procpar(posspaces(2)+1:posspaces(3)-1));
    
    if basictype==1;
        nvalues=str2num(procpar(posspaces(11)+1:posspaces(12)-1));
        param=[];
        for icount=1:nvalues;
            param1=str2num(procpar(posspaces(11+icount)+1:posspaces(12+icount)-1));
            param=[param param1];
        end;
    else;
        nvalues=str2num(procpar(posspaces(11)+1:posspaces(12)-1));
        param=[];
        for icount=1:nvalues;
            param1=procpar(posspaces(11+icount)+2:posspaces(12+icount)-2);
            param=[param param1];
        end;
    end;
end;   

