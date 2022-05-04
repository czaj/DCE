function ResultsOut = genOutput_LCMXL(EstimOpt, Results, Head, Tail, Names, Template1, Template2, Heads, ST)

if ~isfield(EstimOpt,'xlsOverwrite')
    EstimOpt.xlsOverwrite = 1;
end

head1 = {'var.', 'dist.', 'coef.','sign.' ,'st.err.' , 'p-value'};
head2 = {'dist.','coef.','sign.' ,'st.err.' , 'p-value'};
head3 = {'var.', '', 'coef.','sign.' ,'st.err.' , 'p-value'};
head4 = {'','coef.','sign.' ,'st.err.' , 'p-value'};
Coords = struct;
Dim1 = size(Template1,1);
Dim2 = size(Template1,2);
DimA = size(Template2,1);

Block = Results.(Template1{1,1});
RowOut = num2cell(Block);
for l = 1:(size(Block,2)/4)
    RowOut(:,(l-1)*4+2) = star_sig_cell(Block(:,l*4));
end
if ismember(Template1{1,1},ST)
    fixed = 1;
else
    fixed = 0;
end
RowOut = [Names.(Template1{1,1}), distType(Results.Dist(1,:), fixed, size(Block,1)), RowOut];


if fixed == 0
    RowOut = [[head1(1:2) ,repmat(head2(2:end),1,size(Block,2)/4)];RowOut];
else
    RowOut = [[head3(1:2) ,repmat(head2(2:end),1,size(Block,2)/4)];RowOut];
end

headssize = size(Heads.(Template1{1,1}),2);
HeadsTmp = cell(headssize,4*size(Heads.(Template1{1,1}),1)-2);
for s=1:headssize
    indxh = find(~cellfun(@isempty,Heads.(Template1{1,1})(:,s)));
    indxh = indxh(end);
    for m=1:indxh-1
        HeadsTmp(s,4*m-1) = Heads.(Template1{1,1})(m,s);
    end
end
%HeadsTmp(1,3) = Heads.(Template1{1,1});
RowOut = [HeadsTmp; RowOut];


ResultsOut = [];               
for i = 1:Dim1
    for j = 2:Dim2
        if strcmp(Template1{i,j},'NULL')
           Block = Results.(Template1{i,j-1});
           headssize = size(Heads.(Template1{i,j-1}),2);
           ResultsTmp = cell(size(Block,1)+1+headssize,4);
           RowOut = [RowOut, ResultsTmp];
        elseif ~isempty(Template1{i,j})
           Block = Results.(Template1{i,j});
           HeadsTmp = [];
           if ismember(Template1{i,j},ST)
                fixed = 1;
           else
                fixed = 0;
           end
           if size(Block,2) == 4
               ResultsTmp = num2cell(Block);
               ResultsTmp(:,2) = star_sig_cell(Block(:,4));
               if strcmp(Template1{i,j-1},'NULL')
                ResultsTmp = [distType(Results.Dist((j+1)/2,:),fixed,size(Block,1)),ResultsTmp];
               else
                   ResultsTmp = [distType(Results.Dist(j,:),fixed,size(Block,1)),ResultsTmp];
               end
               if fixed == 0
                   ResultsTmp = [head2;ResultsTmp];
               else
                   ResultsTmp = [head4;ResultsTmp];
               end
               headssize = size(Heads.(Template1{i,j}),2);%zakaz dawania roznej wielkosci headsow do rzedu blokow
               HeadsTmp = cell(headssize,5);
               for s=1:headssize
                   HeadsTmp(s,2) = Heads.(Template1{i,j})(1,s);
               end
               ResultsTmp = [HeadsTmp; ResultsTmp];
               %HeadsTmp = [HeadsTmp, Heads.(Template1{i,j})];
           else % This is for Xm, and other similar
               ResultsTmp = num2cell(Block);
               for l = 1:(size(Block,2)/4)
                   ResultsTmp(:,(l-1)*4+2) = star_sig_cell(Block(:,l*4));
               end
               ResultsTmp = [distType(Results.Dist(j,:), fixed, size(Block,1)),ResultsTmp];
               headssize = size(Heads.(Template1{i,j}),2);
               HeadsTmp = cell(headssize,size(ResultsTmp,2));
               for s=1:headssize
                   indxh = find(~cellfun(@isempty,Heads.(Template1{i,j})(:,s)));
                   indxh = indxh(end);
                   for l = 1:indxh-1
                       if l == 1
                           HeadsTmp(s,2) = Heads.(Template1{i,j})(l,s);
                       else
                           HeadsTmp(s,4*l-2) = Heads.(Template1{i,j})(l,s);
                       end
                   end
               end
               if fixed == 0
                   ResultsTmp = [['dist.',repmat(head2(2:end),1,size(Block,2)/4)];ResultsTmp];
               else
                   ResultsTmp = [['',repmat(head2(2:end),1,size(Block,2)/4)];ResultsTmp];
               end
               ResultsTmp = [HeadsTmp; ResultsTmp];
               %HeadsTmp = [HeadsTmp, Heads.(Template1{i,j})'];
           end
           if size(RowOut,1) < size(ResultsTmp,1)
               RowOuttmp = cell(size(ResultsTmp,1),size(RowOut,2));
               RowOuttmp(size(ResultsTmp,1) - size(RowOut,1)+1:end,:) = RowOut;
               RowOut = [RowOuttmp, ResultsTmp];
               %Heads_tmp = cell(size(Heads.(Template1{i,1}),1),max(size(Heads.(Template1{i,1}),2),size(Heads.(Template1{i,j}),2)));
               %Heads_tmp(:,1:size(Heads.(Template1{i,1}),2)) = Heads.(Template1{i,1});
               %Heads.(Template1{i,1}) = Heads_tmp;
               Coords.(Template1{i,1})(1) = 1; Coords.(Template1{i,1})(2) = 0;
           elseif size(RowOut,1) > size(ResultsTmp,1)
               Results_tmp = cell(size(RowOut,1),size(ResultsTmp,2));
               Results_tmp(size(RowOut,1)-size(ResultsTmp,1)+1:end,:) = ResultsTmp;
               RowOut = [RowOut, Results_tmp];
           else
               RowOut = [RowOut, ResultsTmp];
           end   
           %HeadsOut = [HeadsOut, HeadsTmp];          

        end
    end
    
    MinVal = min(size(ResultsOut,2), size(RowOut,2));
    MaxVal = max(size(ResultsOut,2), size(RowOut,2));
    if MinVal == size(ResultsOut,2)
        ResTemp = cell(size(ResultsOut,1), MaxVal);
        ResTemp(:,1:MinVal) = ResultsOut;
        ResultsOut = [ResTemp; RowOut];
    else
        ResTemp = cell(size(RowOut,1), MaxVal);
        ResTemp(:,1:MinVal) = RowOut;
        ResultsOut = [ResultsOut; ResTemp];
    end
    if i ~= Dim1
        Block = Results.(Template1{i+1,1});
        if ismember(Template1{i+1,1},ST)
            fixed = 1;
        else
            fixed = 0;
        end
        RowOut = num2cell(Block);
        for s=1:size(Block,2)/4
            RowOut(:,4*s-2) = star_sig_cell(Block(:,s*4));
        end
        RowOut = [Names.(Template1{i+1,1}), distType(Results.Dist, fixed, size(Block,1)), RowOut]; %it will crash if size of the block and number of variables will differ
        if fixed == 0
            headn1 = [head1(1:2) ,repmat(head2(2:end),1,size(Block,2)/4)];
        else
            headn1 = [head3(1:2) ,repmat(head2(2:end),1,size(Block,2)/4)];
        end
        RowOut = [headn1;RowOut];

        headssize = size(Heads.(Template1{i+1,1}),2);
        HeadsTmp = cell(headssize,2+size(Block,2));
        for s=1:headssize
            indxh = find(~cellfun(@isempty,Heads.(Template1{i+1,1})(:,s)));
            indxh = indxh(end);
            for m=1:indxh-1
                HeadsTmp(s,4*m-1) = Heads.(Template1{i+1,1})(m,s);
            end
        end
        RowOut = [HeadsTmp; RowOut];
%         if length(Head{1,1}) > 5 && strcmp(Head{1,1}(1:5),'LCMXL')
%             len = (size(RowOut,2)-2)/8 - 1;
%             for c=1:len
%                 Cdist = cellstr([repmat(" ",[headssize, 1]); 'dist.'; distType(Results.Dist(c+1,:), fixed,size(Block,1))]);
%                 RowOut = [RowOut(:,1:10+(c-1)*9), Cdist, RowOut(:,11+(c-1)*9:end)];
%             end
%         end
    end
end


Changed = struct;

for i = 1:Dim1
    for j = 1:Dim2
        if ~isempty(Template1{i,j})
            if j>=2
                if ~strcmp(Template1{i,j-1},'NULL')
                    Blockh = Results.(Template1{i,j-1});
                    Coords.(Template1{i,j}) = [Coords.(Template1{i,j-1})(1),Coords.(Template1{i,j-1})(2) + size(Blockh,2)+1];
                else
                    Blockh = Results.(Template1{i,j-2});
                    Coords.(Template1{i,j}) = [Coords.(Template1{i,j-2})(1),Coords.(Template1{i,j-2})(2) + size(Blockh,2)+1+4];
                end
            else
                if i>=2
                    Blockv = Results.(Template1{i-1,j});
                    Coords.(Template1{i,j}) = [Coords.(Template1{i-1,j})(1)+size(Blockv,1)+2,4*j-1];
                    if size(Heads.(Template1{i,j}),2) > 1
                        Coords.(Template1{i,j})(1) = Coords.(Template1{i,j})(1) + size(Heads.(Template1{i,j}),2) - 1;
                    end
                else
                    if isfield(Coords, Template1{i,j}) && ~isempty(Coords.(Template1{i,j}))
                        Coords.(Template1{i,j}) = [Coords.(Template1{i,j})(1)+2+i,Coords.(Template1{i,j})(2)+4*j-1];
                    else
                        Coords.(Template1{i,j}) = [2+i,4*j-1];
                    end
                    
                    if size(Heads.(Template1{i,j}),2) > 1
                        Coords.(Template1{i,j})(1) = Coords.(Template1{i,j})(1) + size(Heads.(Template1{i,j}),2) - 1;
                    end
                end 
            end
        end
    end
end

if EstimOpt.Display~=0
    spacing = 2;
    precision = 4;
    
    
fprintf('\n')
fprintf('__________________________________________________________________________________________________________________')
fprintf('\n')
fprintf('\n')
cprintf('*Black',[Head{1,1},' ']);
cprintf('*Black',strcat(Head{1,2}, '\n'));

for i=1:DimA
    FirstBlock = Coords.(Template2{i,1})(1);
    [~,CWt] = CellColumnWidth(ResultsOut(FirstBlock:FirstBlock+size(Results.(Template2{i,1}),1)-1,:));
    if i>1
        CW(CWt >= CW) = CWt(CWt >= CW);
    else
        CW = CWt;
    end
    indx = find(~cellfun(@isempty,Template2(i,:)));
    indx = indx(end);
    for c =1:indx
        Y = Coords.(Template2{i,c})(2);
        if ~strcmp(Template2{i,c}, 'NULL')
            indxh = find(~cellfun(@isempty,Heads.(Template2{i,c})(:,1)));
            if ~isempty(indxh)
            	indxh = indxh(end);
            else
                    indxh = 1;
            end
        end
        if ~strcmp(Template2{i,c}, 'NULL')
            method = Heads.(Template2{i,c}){indxh,1};
        else
            method = 'lc';
        end
        if strcmp(method,'lc') | strcmp(method,'lb')
            if strcmp(Template2{i,c},'NULL')
                name = ' '; 
            else
                name = Heads.(Template2{i,c}){1,1};
            end
        else
            name = ' '; 
        end
        if length(name) > (CW(1)+spacing+CW(2)+4+CW(Y))
            CW(1) = length(name)-(CW(1)+spacing+CW(2)+4+CW(Y))+spacing;
        end
    end
end

for i=1:DimA
    indx = find(~cellfun(@isempty,Template2(i,:)));
    indx = indx(end);
    %UPPERHEADER
    headssize = size(Heads.(Template2{i,1}),2);
    for s = 1:headssize 
        for c =1:indx
            if strcmp(Template2{i,c}, 'NULL')
                Y = Coords.(Template2{i,c-1})(2)+4;
            else            
                Y = Coords.(Template2{i,c})(2);
            end
            if ~strcmp(Template2{i,c}, 'NULL')
                indxh = find(~cellfun(@isempty,Heads.(Template2{i,c})(:,s)));
                if ~isempty(indxh)
                    indxh = indxh(end);
                else
                    indxh = 1;
                end
            end
            if ~strcmp(Template2{i,c}, 'NULL')
            	method = Heads.(Template2{i,c}){indxh,s};
            else
                method = 'lc';
            end
          
            if strcmp(method,'lc')
                if strcmp(Template2{i,c},'NULL')
                   name = ' '; 
                else
                    name = Heads.(Template2{i,c}){1,s};
                end
                if c == 1                
                    fprintf('%-*s',CW(1)+spacing*2+CW(2)+4+CW(Y),name)
                else
                    if ~(strcmp(Template2{i,c},'NULL'))
                        fprintf('%*s',5,' ')
                    end
                     fprintf('%-*s',sum(CW(Y+2:Y+3)) + precision*3 + 16 + CW(Y), name)
                end
                for m = 2:indxh-1
                    name = Heads.(Template2{i,c}){m,s};
                    fprintf('%-*s',sum(CW(Y+(m-1)*4+2:Y+(m-1)*4+3)) + precision*3 + 16 + CW(Y+m*4), name)
               end
            elseif strcmp(method,'tb') || strcmp(method,'lb') || strcmp(method,'tc')
                if method(1) == 't'
                    fprintf('%*s',CW(1)+spacing*2+CW(2)+4+CW(Y),' ')
                end
                if method(1) == 'l' & ~(strcmp(Template2{i,c},'NULL')) & c ~= 1
                    fprintf('%*s',5,' ')
                end
                
                for m = 1:indxh-1
                    if strcmp(Template2{i,c},'NULL')
                       name = ' '; 
                    else
                        name = Heads.(Template2{i,c}){m,s};
                    end
                    if m~=(indxh-1)
                        fprintf('%-*s',sum(CW(Y+(m-1)*4+2:Y+(m-1)*4+3)) + precision*3 + 16 + CW(Y+m*4), name)
                    else
                        fprintf('%-*s',sum(CW(Y+(m-1)*4+2:Y+(m-1)*4+3)) + precision*3 + 16+1, name)
                    end
                end
            end
            if method(2) == 'b'
                fprintf('\n')
            end
            if strcmp(Template2{i,c}, 'NULL') & c==indx
                fprintf('\n')
            end
        end
    end
    %\UPPERHEADER
    %HEADER
        fprintf('%-*s%-*s',CW(1)+spacing,ResultsOut{Coords.(Template2{i,1})(1) - 1,1},CW(2)+ 4, ResultsOut{Coords.(Template2{i,1})(1) - 1,2})
        for c =1:indx
            if strcmp(Template2{i,c}, 'NULL')
                X = Coords.(Template2{i,c-1})(1);
                Y = Coords.(Template2{i,c-1})(2)+4;
            else            
                X = Coords.(Template2{i,c})(1);
                Y = Coords.(Template2{i,c})(2);
            end
            if strcmp(Template2{i,c}, 'NULL')
                fprintf('%1s%*s%*s%*s%s',' ',CW(Y)+spacing+precision,"     ", CW(Y+2)+spacing+precision+4,"       ", CW(Y+3)+spacing+precision+2,"       ",'   ')           
            else
                    if c>1                
                        fprintf('%-*s', 5, ResultsOut{Coords.(Template2{i,c})(1) - 1,Coords.(Template2{i,c})(2)-1})
                    end
                    for m=1:size(Results.(Template2{i,c}),2)/4
                        fprintf('%1s%*s%*s%*s%s',' ',CW(Y+(m-1)*4)+spacing+precision,ResultsOut{X-1,Y+(m-1)*4}, CW(Y+(m-1)*4+2)+spacing+precision+4,ResultsOut{X-1,Y+(m-1)*4+2}, CW(Y+(m-1)*4+3)+spacing+precision+2,ResultsOut{X-1,Y+(m-1)*4+3},'   ')
                    end
           end
        end
        fprintf('\n')
    %\HEADER
    %VALUES
    
            for d=1:size(Results.(Template2{i,1}),1)
                fprintf('%-*s%-*s', CW(1)+spacing+1,ResultsOut{Coords.(Template2{i,1})(1) + d - 1,1}, CW(2)+3, ResultsOut{Coords.(Template2{i,1})(1) + d - 1,2})
                for c =1:indx
                    if strcmp(Template2{i,c}, 'NULL')
                        X = Coords.(Template2{i,c-1})(1);
                        Y = Coords.(Template2{i,c-1})(2)+4;
                    else            
                        X = Coords.(Template2{i,c})(1);
                        Y = Coords.(Template2{i,c})(2);
                    end
                    if c~=1 & ~strcmp(Template2{i,c}, 'NULL')
                        fprintf('%-*s', 5, ResultsOut{X,Y-1})
                    end
                    if strcmp(Template2{i,c}, 'NULL')
                        fprintf('%1s%*s%*s%*s%s',' ',CW(Y+(m-1)*4)+spacing+precision,"     ", CW(Y+(m-1)*4+2)+spacing+precision+4,"       ", CW(Y+(m-1)*4+3)+spacing+precision+2,"       ",'   ')           
                    else
                        for m=1:size(Results.(Template2{i,c}),2)/4
                            fprintf('% *.*f%-3s% *.*f% *.*f%s', CW(Y+(m-1)*4)+spacing+precision+1,precision,ResultsOut{X+d-1,Y+(m-1)*4}, ResultsOut{X+d-1,Y+1+(m-1)*4}, CW(Y+(m-1)*4+2)+spacing+precision+1,precision,ResultsOut{X+d-1,Y+2+(m-1)*4}, CW(Y+3+(m-1)*4)+spacing+precision+2,precision,ResultsOut{X+d-1,Y+3+(m-1)*4},'   ')
                        end
                    end
                end
                fprintf('\n')
            end
        disp(' ');        
end
    %\VALUES
    [~,CWm] = CellColumnWidth(num2cell(Results.stats));
    cprintf('*Black', 'Model diagnostics: \n')
    fprintf('%-29s%*.*f\n', 'LL at convergence:',CWm(1)+spacing+precision+1,precision, Results.stats(1))
    fprintf('%-29s%*.*f\n', 'LL at constant(s) only:',CWm(1)+spacing+precision+1,precision, Results.stats(2))
    fprintf('%-29s%*.*f\n', strcat('McFadden''s pseudo-R',char(178),':'),CWm(1)+spacing+precision+1,precision, Results.stats(3))
    fprintf('%-29s%*.*f\n', strcat('Ben-Akiva-Lerman''s pseudo-R',char(178),':'),CWm(1)+spacing+precision+1,precision, Results.stats(4))
    fprintf('%-29s%*.*f\n','AIC/n:',CWm(1)+spacing+precision+1,precision, Results.stats(5))
    fprintf('%-29s%*.*f\n','BIC/n:',CWm(1)+spacing+precision+1,precision, Results.stats(6))
    fprintf('%-29s%*.*f\n','n (observations):',CWm(1)+spacing,0, Results.stats(7))
    fprintf('%-29s%*.*f\n','r (respondents):',CWm(1)+spacing,0, Results.stats(8))
    fprintf('%-29s%*.*f\n','k (parameters):',CWm(1)+spacing,0, Results.stats(9))
    disp(' ')
    for i = 13:size(Tail,1)
        fprintf('%-23s%-s\n', Tail{i,1}, Tail{i,2})
    end
    
    disp(' ')
        clocknote = clock;
        tocnote = toc;
        [~,DayName] = weekday(now,'long');
        disp(['Estimation completed on ' DayName ', ' num2str(clocknote(1)) '-' sprintf('%02.0f',clocknote(2)) '-' sprintf('%02.0f',clocknote(3)) ' at ' sprintf('%02.0f',clocknote(4)) ':' sprintf('%02.0f',clocknote(5)) ':' sprintf('%02.0f',clocknote(6))])
        disp(['Estimation took ' num2str(tocnote) ' seconds ('  num2str(floor(tocnote/(60*60))) ' hours ' num2str(floor(rem(tocnote,60*60)/60)) ' minutes ' num2str(rem(tocnote,60)) ' seconds).']);
    
end


% Adding head and tail

Indx = size(ResultsOut,2);
HeadOut = cell(size(Head,1), Indx);
HeadOut(:,1:2) = Head;
ResultsOut = [HeadOut; ResultsOut];
TailOut = cell(size(Tail,1), Indx);
TailOut(:,1:2) = Tail;
ResultsOut = [ResultsOut; TailOut];

% excel

fullOrgTemplate = which('template.xls');
currFld = pwd;

if isfield(EstimOpt,'ProjectName')
    fullSaveName = strcat(currFld,'\', Head(1,1), '_results_',EstimOpt.ProjectName,'.xls');
else
    fullSaveName = strcat(currFld,'\', Head(1,1), '_results.xls');
end

copyfile(fullOrgTemplate,'templateTMP.xls')
fullTMPTemplate = which('templateTMP.xls');
try
    ex = actxGetRunningServer('Excel.Application');
catch
end

excel = actxserver('Excel.Application');
excelWorkbook = excel.Workbooks.Open(fullTMPTemplate);
excel.Visible = 1;
excel.DisplayAlerts = 0;
excelSheets = excel.ActiveWorkbook.Sheets;
excelSheet1 = excelSheets.get('Item',1);
excelSheet1.Activate;
column = size(ResultsOut,2);
columnName = [];
while column > 0
    modulo = mod(column - 1,26);
    columnName = [char(65 + modulo) , columnName]; %#ok<AGROW>
    column = floor(((column - modulo) / 26));
end
rangeE = strcat('A1:',columnName,num2str(size(ResultsOut,1)));
excelActivesheetRange = get(excel.Activesheet,'Range',rangeE);
excelActivesheetRange.Value = ResultsOut;
fullSaveName = strjoin(fullSaveName);
if isfield(EstimOpt,'xlsOverwrite') && EstimOpt.xlsOverwrite == 0
    i = 1;
    while exist(fullSaveName,'file') == 2
        if ~contains(fullSaveName, '(')
            pos = strfind(fullSaveName, '.xls');
            fullSaveName = strcat(fullSaveName(1:pos-1),'(',num2str(i),').xls');
        else
            pos = strfind(fullSaveName, '(');
            fullSaveName = strcat(fullSaveName(1:pos),num2str(i),').xls');
        end
        i = i+1;
    end
elseif isfield(EstimOpt,'xlsOverwrite') && EstimOpt.xlsOverwrite == 1
    if exist('ex','var')
        wbs = ex.Workbooks;
        for j = 1:wbs.Count
            if strcmp(wbs.Item(j).FullName,fullSaveName)
                if ~contains(fullSaveName, '(')
                    pos = strfind(fullSaveName, '.xls');
                    fullSaveName = strcat(fullSaveName(1:pos-1),'(',num2str(1),').xls');
                else
                    pos = strfind(fullSaveName, '(');
                    pos2 = strfind(fullSaveName, ')');
                    num = str2double(fullSaveName(pos+1:pos2-1)) + 1;
                    fullSaveName = strcat(fullSaveName(1:pos),num2str(num),').xls');
                end
            end
        end
    end
end

excelWorkbook.ConflictResolution = 2;
SaveAs(excelWorkbook,fullSaveName);
excel.DisplayAlerts = 0;
excelWorkbook.Saved = 1;
Close(excelWorkbook)
Quit(excel)
delete(excel)
delete(fullTMPTemplate)