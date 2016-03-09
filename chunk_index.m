function linind = chunk_index(plt_sze,format_string)
%--------------------------------------------------!
% linind=chunk_index(plt_sze,format_string)        :
%--------------------------------------------------!
if ~isempty(format_string)
    %---
    %generate key between plate position and linear index
    plt2lin=NaN*ones(plt_sze);
    for i = 1 : plt_sze(1)
        for j = 1 : plt_sze(2)
            plt2lin(i,j)=(i-1)*plt_sze(2)+j;
        end
    end

    %---
    %break up the string according to chunks (rectangles, split by ",").
    chunks=regexp(format_string,',','split');
    chunks=regexprep(chunks,{'[',']'},'');

    %---
    %interpret individual chunks, add to linind().
    linind=[];
    for i = 1:length(chunks)
        %break up according to row(Y),column(X)
        XY=regexp(chunks{i},':','split');
        Yrange=regexp(XY{1},'-','split');
        Xrange=regexp(XY{2},'-','split');
        for i = 1 : length(Yrange)
            Y(i)=str2num(Yrange{i});
        end
        for j = 1 : length(Xrange)
            X(j)=str2num(Xrange{j});
        end
        range_size=[length(Y),length(X)];
        %look up linear index
        new_ind=[];
        if range_size==[1,1]
                new_ind=plt2lin(Y,X);
        elseif range_size==[1,2]
                new_ind=plt2lin(Y,X(1):X(2));
        elseif range_size==[2,1]
                new_ind=plt2lin(Y(1):Y(2),X);
        elseif range_size==[2,2]
                new_ind=plt2lin(Y(1):Y(2),X(1):X(2));
        end
        linind=[linind;new_ind(:)];
    end
else
    lindind=[];
end
end
