function data=slurp_dat(filepath,format)
%---------------------------------------!
% struct = slurp_dat(filename, format); :
%---------------------------------------!
% reads in a single data file from a    :
% colony sizer, in the format specified :
% by 'format'.                          :
%---------------------------------------!
% Currently supported format strings:   :
% collins_v1:  row|col|sze|crc          :
% collins_v2:  row|col|sze|crc|int      :
% iris_v0_ecogrowth:  row|col|sze       :
% iris_v0_ecoopacity:row|col|sze|crc|opc:
% iris_v0_ecobiofilm:row|col|sze|crc|int|bsze|bint|rsze|nint|opc:
% iris_kritikos: row|col|opc            :
% krit_dat: opc | crc                   :
%---------------------------------------!
% Anthony Shiver (2013)                 :
%---------------------------------------!
fid=fopen(filepath,'r');
switch format
    case 'collins_v1'
        cell=textscan(fid,'%*u%*u%f%f','HeaderLines',1);
        [data.sze,data.crc]=cell{:};
    case 'collins_v2'
        cell=textscan(fid,'%*u%*u%f%f%f','HeaderLines',1);
        [data.sze,data.crc,data.int]=cell{:};
    case 'iris_v0_ecogrowth'
        cell=textscan(fid,'%*u%*u%f','HeaderLines',7);
        [data.sze]=cell{:};
    case 'iris_v0_ecoopacity'
        cell=textscan(fid,'%*u%*u%f%f%f','HeaderLines',7);
        [data.sze,data.crc,data.opc]=cell{:};
    case 'iris_kritikos'
        cell=textscan(fid,'%*u%*u%f','HeaderLines',7,'TreatAsEmpty','NA');
        [data.opc]=cell{:};
    case 'krit_dat'
        cell=textscan(fid,'%*u%*u%f%f','HeaderLines',1);
        [data.opc,data.crc]=cell{:};
    case 'iris_v0_ecobiofilm'
        cell=textscan(fid,'%*u%*u%f%f%f%f%f%f%f%f','HeaderLines',7);
        [data.sze,data.crc,data.int,data.bsze,data.bint,data.rsze,data.nint,data.opc]=cell{:};
    otherwise
        disp([filepath ': Proper format not specified!']);
end
fclose(fid);
end
