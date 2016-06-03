function [xPositions, yPositions, Label, RangeCut] = UnivarScatter(data, varargin)
%--------------------------------------------------------------------------
% UnivariateScatter - Draw an Univariate Scatter plot out of a nx2 table with a
%                     categorical/string and a numerical variable, or out of a
%                     numerical array, with groups/categories corresponding
%                     to the columns of the numerical array.
%                     Many custom options are available as Name,Value
%                     pairs. For optimal visualization of your data I
%                     recommend to play with 'RangeCut', and also with pbaspect of the plot, it
%                     really changes the appearance.
%                     You need the following functions for this function to
%                     work:
%                                              
%                       +CatTable2Array: included in the file
%
%                     Also, a very simple function to assign colors is
%                     provided, ColorCoder, you can see an example in the
%                     script attached
%
%  Input: data       - Can take two types of input:
%
%                       + a nx2 table with a categorical/string and a numerical
%                       variable, the groups will be made according to the
%                       categories
%
%                       + a numerical array in which the different columns
%                       correspond to the groups/categories of the Univariate 
%                       Scatter plot in x axis. When plotting data with
%                       different number of points for the different
%                       groups, the empty elements of the array should be
%                       filled with nan. In order to do this much faster I
%                       recommend the function padcat(available in
%                       mathworks, Copyright (c) 2009, Jos van der Geest), which already
%                       concatenates row/column vectors of different
%                       lengths and fills with nans the missing values.
%
%         Name,Value  - Name-Value pairs arguments to customize the plot:
%
%                       +'Label' - Only when data is a numerical array
%                       A cell array of strings containing the
%                       labels that you want to give to categories/groups
%                       in the plot, by default they are 1,2,3...n The
%                       number of columns of this array must be the same as
%                       the number of columns in data, example:
%                       UniVarScatter(data,'Label',{'a','b','c'})
%
%                       +'Whiskers' - A string which can take 3 values
%                           *'lines': std, 95% SEM, mean are represented as
%                           lines
%                           *'boxes':std, 95% SEM, mean are represented as
%                           boxes. Default value
%                           *'none':none of these are shown
%
%                       +'PointStyle' - A String containing the different
%                       values that can be given as markertype for scatter,
%                       such as 'o','*''.', etc. Default is 'o'
%
%                       +'Width' - A number that determines the spread of
%                       each subset of points in x axis. Default is 0.4
%
%                       +'Compression' - A number that determines how
%                       separated from each other the points are in x
%                       dimension. Default is 5
%
%                       +'RangeCut' - A number or column/row vector with
%                       length= size(data,2) that determines how broad
%                       groups of points are made in y dimension, groups
%                       width in y dimension corresponds to 1/RangeCut of the
%                       95% confidence interval of a normal distribution
%                       with std and mean of each subset of points.
%                       Default is 2+2*sqrt(number of points) for each subset
%                       of points, this has shown to be a good value, but
%                       it does not perform optimally for every dataset, so
%                       I strongly recommend not to discard the graph and
%                       try different combinations of values if you don't
%                       obtained a nice result.
%
%                       +'WhiskersWidthRatio' - Lines/boxes that represent
%                       the mean,SEM,std will be Width/WhiskersWidthRatio
%                       long in x dimension. Default is 4
%
%                       +'PointSize' - A number indicating the point size.
%                       Default is 36
%
%                       +'LineWidth' - A number indicating the width of the
%                       edge of the points. Default is 0.5
%
%                       +'MarkerEdgeColor' color(s) of the edges of the points; def. 'k'
%                       +'MarkerFaceColor' color(s) of the points; def. 'flat'
%                       +'MeanColor' color(s) of the mean line; def. black
%                       +'SEMColor'  color(s) of the SEM box; def. dark gray
%                       +'StdColor'  color(s) of the Std box; def. light gray
%
%                       All this last Color Properties behave similarly,
%                       they can take several types of values
%
%                           *RGB vector 1x3: all the representations
%                           are made in the same colour represented by that vector
%                     
%                           *A string: it may contain the name of the color
%                           for instance 'k' or 'black', also it can take
%                           the values 'none', or 'flat'
%
%                           *RGB array nx3 where n>=size(data,2): each row
%                           of this array corresponds to the color of the n
%                           group/category in x. You can use the function
%                           ColorCoder included in the zip to Easily create
%                           a custom RGB array by chosing the colors from
%                           the matlab palette
%
%  Output: the figure
%
%          xPositions - An Array with the x  positions of the points, with
%                       the same size as yPositions
%                           
%          
%          yPositions - a numerical array in which the different columns
%                       correspond to the groups/categories of the Univariate 
%                       Scatter plot in x axis, its the same as data if data was
%                       already a numerical array, and if it was a table,
%                       it generates a numerical array, and the columns
%                       correspond to the Label groups.
%       
%          Label      - A string cell array with the names of the groups
%                       specially useful when using a table as input, so
%                       you know the order in which the scatter plots are
%                       showed.
%
%          RangeCut   - The value of RangeCut(explained at the Name,Values
%                       section), specially useful to optimize the
%                       representation of your data if the RangeCut value
%                       generated by the function does not perform well for
%                       your data.
%
% version: <1.00> from 30/11/2015
% 
%       Manuel Lera Ramírez: manulera14@gmail.com
%
%       Developed during Master Thesis Internship in Marcos González-Gaitán's Lab
%       Departments of Biochemistry and Molecular Biology, Sciences II, 
%       30 Quai Ernest-Ansermet, CH-1211 Geneva 4, Switzerland
%
% Modified to eliminate figure output: GAR - University of Glasgow - 25 May 2016
%--------------------------------------------------------------------------


%% Variable handling and default values generation

if nargin==0
    error(sprintf('There\''s no input'))
end
if mod(length(varargin),2)~=0
    error(sprintf('The arguments should be in pairs,in which the first one is a string, as in (data, \''Label\'',{\''A\'',\''B\''})'))
end

stringVars=varargin(1:2:end);
valueVars=varargin(2:2:end);
if any(~cellfun(@isstr,stringVars))
    error(sprintf('The arguments should be in pairs,in which the first one is a string, as in (data, \''Label\'',{\''A\'',\''B\''})'))
end

%data_istable is a logical value used to suppress the default Label value
%assignation

if istable(data)
    data_istable=true;
    [data,Label]=CatTable2Array(data);
    warning('The label will be taken from the table, any label input will be ignored');
else
    data_istable=false;
end

%% String Variables
if ~data_istable
LabelInd=strmatch('Label',stringVars);
if ~isempty(LabelInd) && data_istable==false
    %Check that the Label is a string cell array and has a propper length
    if iscellstr(valueVars{LabelInd}) && length(valueVars{LabelInd})==size(data,2) 
        Label=valueVars{LabelInd};
    else
        error('wrong Label, Label should be a linear cell array of strings, its length should be the same as the number of columns of data')
    end
else
    Label=[];
end
end

Ind=strmatch('Whiskers',stringVars);
if ~isempty(Ind)
    %Check that Whiskers is a string with a valid value
    if ischar(valueVars{Ind}) && any([strcmp('box',valueVars{Ind}),strcmp('none',valueVars{Ind}),strcmp('lines',valueVars{Ind})]) 
        Whiskers=valueVars{Ind};
    else
        error(sprintf('wrong Whiskers, it should be a string saying \''none\'', \''box\'' or \''lines\'''))
    end
else
    Whiskers='box';
end


Ind=strmatch('PointStyle',stringVars);
if ~isempty(Ind)
    %Check that PointStyle is a string
    if ischar(valueVars{Ind}) 
        PointStyle=valueVars{Ind};
    else
        error(sprintf('wrong PointStyle, it should be a string defining the shape of the point \''o\'', \''.\'', \''*\'', etc.'))
    end
else
    PointStyle='o';
end

%% Numeric Variables
Ind=strmatch('Width',stringVars);
if ~isempty(Ind)
    %Check that Width is a number
    if isnumeric(valueVars{Ind}) && all(size(valueVars{Ind})== [1,1])
        Width=valueVars{Ind};
    else
        error(sprintf('wrong Width, it should be a number'))
    end
else
    Width=0.4;
end

Ind=strmatch('Compression',stringVars);
if ~isempty(Ind)
    %Check that Compression is a number
    if isnumeric(valueVars{Ind}) && all(size(valueVars{Ind})== [1,1])
        Compression=valueVars{Ind};
    else
        error(sprintf('wrong Compression, it should be a number'))
    end
else
    Compression=5;
end

Ind=strmatch('RangeCut',stringVars);
if ~isempty(Ind)
    %Check that RangeCut is a number
    if isnumeric(valueVars{Ind})
        if size(valueVars{Ind})==1
            RangeCut=valueVars{Ind}*ones(size(data,2),1);
        elseif all(size(valueVars{Ind})== [size(data,2),1]) || all(size(valueVars{Ind}) == [1,size(data,2)])
            RangeCut=valueVars{Ind};
        else
            error(sprintf('wrong RangeCut, it should be a number, or a column/row vector with length=size(data,2)'))
        end
    else
        error(sprintf('wrong RangeCut, it should be a number, or a column/row vector with length=size(data,2)'))
    end
else
    numPoints=sum(~isnan(data));
    RangeCut=(log(numPoints)/log(150)+3*(numPoints).^(1/3))';
end

Ind=strmatch('WhiskersWidthRatio',stringVars);
if ~isempty(Ind)
    %Check that WhiskersWidthRatio is a number
    if isnumeric(valueVars{Ind}) && all(size(valueVars{Ind})== [1,1])
        WhiskersWidthRatio=valueVars{Ind};
    else
        error(sprintf('wrong WhiskersWidthRatio, it should be a number'))
    end
else
   WhiskersWidthRatio=4;
end

Ind=strmatch('PointSize',stringVars);
if ~isempty(Ind)
    %Check that PointSize is a number
    if isnumeric(valueVars{Ind}) && all(size(valueVars{Ind})== [1,1])
        PointSize=valueVars{Ind};
    else
        error(sprintf('wrong PointSize, it should be a positive number'))
    end
else
   PointSize=36;
end

Ind=strmatch('LineWidth',stringVars);
if ~isempty(Ind)
    %Check that LineWidth is a number
    if isnumeric(valueVars{Ind}) && all(size(valueVars{Ind})== [1,1])
        LineWidth=valueVars{Ind};
    else
        error(sprintf('wrong LineWidth, it should be a positive number'))
    end
else
   LineWidth=0.5;
end


%% Mixed Variables
%Some of the following variables are the classic scatter personalizing tools,
%please revisit the scatter documentation for deeper understanding.

Ind=strmatch('MarkerEdgeColor',stringVars);
EdgeIsArray=false;
if ~isempty(Ind)
    %Check that MarkerEdgeColor takes valid arguments
    if isnumeric(valueVars{Ind}) && all(size(valueVars{Ind}) >= [size(data,2),3])
        MarkerEdgeColor=valueVars{Ind};
        EdgeIsArray=true;
    elseif (isnumeric(valueVars{Ind}) && all(size(valueVars{Ind}) == [1,3])) || ischar(valueVars{Ind})
        MarkerEdgeColor=valueVars{Ind};
    else
        error(sprintf('wrong MarkerEdgeColor input'))
    end
else
    MarkerEdgeColor='k';
end

Ind=strmatch('MarkerFaceColor',stringVars);
%This will be used later when plotting
FaceIsArray=false;

if ~isempty(Ind)
    %Check that MarkerFaceColor takes valid arguments
    if isnumeric(valueVars{Ind}) && all(size(valueVars{Ind}) >= [size(data,2),3])
        MarkerFaceColor=valueVars{Ind};
        FaceIsArray=true;
    elseif (isnumeric(valueVars{Ind}) && all(size(valueVars{Ind}) == [1,3])) || ischar(valueVars{Ind})
        MarkerFaceColor=valueVars{Ind};
    else
        error(sprintf('wrong MarkerFaceColor input')) 
    end
else
    MarkerFaceColor='flat';
end

Ind=strmatch('MeanColor',stringVars);
MeanColorIsArray=false;
if ~isempty(Ind)
    %Check that MeanColor takes valid arguments
    if isnumeric(valueVars{Ind}) && all(size(valueVars{Ind}) >= [size(data,2),3])
        MeanColor=valueVars{Ind};
        MeanColorIsArray=true;
    elseif (isnumeric(valueVars{Ind}) && all(size(valueVars{Ind}) == [1,3])) || ischar(valueVars{Ind})
        MeanColor=valueVars{Ind};
    else
        error(sprintf('wrong MeanColor input'))
    end
else
    MeanColor='k';
end

Ind=strmatch('SEMColor',stringVars);
SEMColorIsArray=false;
if ~isempty(Ind)
    %Check that SEMColor takes valid arguments
    if isnumeric(valueVars{Ind}) && all(size(valueVars{Ind}) >= [size(data,2),3])
        SEMColor=valueVars{Ind};
        SEMColorIsArray=true;
    elseif (isnumeric(valueVars{Ind}) && all(size(valueVars{Ind}) == [1,3])) || ischar(valueVars{Ind})
        SEMColor=valueVars{Ind};
    else
        error(sprintf('wrong SEMColor input'))
    end
elseif strcmp(Whiskers,'box')
    SEMColor=[1 1 1]*0.95;
else
    SEMColor='k';
end

Ind=strmatch('StdColor',stringVars);
StdColorIsArray=false;
if ~isempty(Ind)
    %Check that StdColor takes valid arguments
    if isnumeric(valueVars{Ind}) && all(size(valueVars{Ind}) >= [size(data,2),3])
        StdColor=valueVars{Ind};
        StdColorIsArray=true;
    elseif (isnumeric(valueVars{Ind}) && all(size(valueVars{Ind}) == [1,3])) || ischar(valueVars{Ind})
        StdColor=valueVars{Ind};
    else
        error(sprintf('wrong StdColor input'))
    end
elseif strcmp(Whiskers,'box')
    StdColor=[1 1 1]*0.85;
else
    StdColor='k';
end



%% Plotting data
HoldWasOn=ishold;
%This is just a trick to not change anything if there was a plot with hold on already,
%to create a new figure if there was none, and to clear the figure if it
%had hold off
plot(nan,nan)
hold on

xPositions=nan(size(data));
yPositions=data;
for i=1:size(data,2)

%% Sort the data and determine the limits of the different groups in which we divide the points 

% yValues is a column vector containing one of the columns of the matrix
% data, it contains the y values of one of the subset of points we want to
% represent.
% We sort it in order to make groups according to its value, however, we
% want to keep the information of the order of the yvalues, just in case
% its needed for something, we used SortingIndex for that

% eliminate the nans and keep the information to later put it in
% xPositions, this info is kept in yValues_nanIndex
yValues_nanIndex=~isnan(data(:,i));
yValues=data(yValues_nanIndex,i);

[yValues,y_SortingIndex]=sort(yValues);
[~,x_SortingIndex]=sort(y_SortingIndex);

%If one of the columns is empty, we don't plot anything
if isempty(yValues)
    continue
end
% xValues is a column vector as big as yValues, but with all values equal to i. If we
% represented plot(xValues, yValues), we would get all the points with the
% same x,a dot blot, and therefore there would be a massive overlapping
% depending on how many points we had. Thats why the while loops that we
% find afterwards change the x of each point, so they spread and do not
% overlap.
xValues=ones(size(yValues))*i;

%% Plotting the error bars and the standard deviation bars

%Color type managing

if MeanColorIsArray
    MeanIndex=i;
else
    MeanIndex=1;
end

if SEMColorIsArray
    SEMIndex=i;
else
    SEMIndex=1;
end

if StdColorIsArray
    StdIndex=i;
else
    StdIndex=1;
end

%If we want boxes for the mean and std
if strcmp(Whiskers,'box')
    
    yMean=mean(yValues);
    yStd=std(yValues);
    ySem=yStd/sqrt(size(yValues,1));
    yCI=ySem*1.96;
    %plot the standard deviation box
    rectangle('Position',[i-Width/WhiskersWidthRatio,yMean-yStd,2*Width/WhiskersWidthRatio,2*yStd ],'FaceColor',StdColor(StdIndex,:),'EdgeColor', StdColor(StdIndex,:),'LineWidth',0.1);
    %plot the mean+-SEM box
    rectangle('Position',[i-Width/WhiskersWidthRatio,yMean-yCI,2*Width/WhiskersWidthRatio,2*yCI ],'FaceColor',SEMColor(SEMIndex,:),'EdgeColor', SEMColor(SEMIndex,:),'LineWidth',0.1);
    %plot the mean line 
    plot([i-Width/WhiskersWidthRatio i+Width/WhiskersWidthRatio],[yMean yMean],'Color',MeanColor(MeanIndex,:),'LineWidth',2)

end
%% Changing the xValue of each point
% This is done the following way: The points are sorted in equidistant
% groups according to yValues, this distance is stablished as 1/RangeCut of the
% width of the 95% interval of the data. Of course this makes sense for
% more or less normally distributed data, but it can be changed easily, and also it can
% be changed easily changing RangeCut value

range_val=norminv([0.05 0.95],mean(yValues),std(yValues));

%range_val=quantile(yValues,[0.05 0.95]);
cuts=abs(range_val(1)-range_val(2))/RangeCut(i);

% In case one of the variables has equal values for all its points, the
% norminv will return a nan, since std(yValues) will be 0, in that case we
% arbitrarily give cuts the value of cuts=mean(yValues)/2, this value does not
% really matter because there will only be one group anyway
if isnan(cuts)
    cuts=mean(yValues)/2;
end
% The "seed" group contains the values which are in the range of the
% mean+-1/2cuts. cutUp is the higher border of the interval. So we
% take this value, and from it we move up and down in this two loops(steps is a variable 
% that allows us to use the same while loop to move in both directions,
% selecting the subset of yValues that are in each group range, and modifying
% their xValues, so that they are spread and do not overlap.

for steps=[-1,1]
    cutUp=mean(yValues)+cuts/2;
    subsetind= yValues<=cutUp & yValues>cutUp-cuts; 
    keep_going=true;
while keep_going
    % subsetind is a logical variable that represents how many points are
    % in that group
    if all(sum(subsetind)~=[1 0]) %If there's only one point, we represent it in the middle, and xValues remains equal to i
        
        % distmaker is a variable that is equal to Width when sum(subsetind)=1 , and gets smaller when the number of points
        % gets bigger in a group, exponentially, and it tends to zero in
        % the infinite. I don't have strong arguments for the choice of
        % this particular expression, but for me it works very well, and
        % you can tune with Width and Compression very well the appereance of
        % the plot
        distmaker=Width./(exp((sum(subsetind).^2-1)./(sum(subsetind)*Compression)));

        % xSubset is a column vector with all values equal to i, and as
        % long as the subset of the number of values in yValues that belong
        % to this particular group range
        xSubset=ones(sum(subsetind),1)*i;
        
        %oddie is a variable that indicates whether the number of points in
        %the group is odd or even, it's very important for the indexing
        %afterwards
        oddie=mod(size(xSubset,1),2);
        
        % xb is a symmetrical vector with mean=0, and the values in it have
        % the displacements in x that we want to apply to xValues, but if
        % we just applied it increasingly, the smaller yValue would get the most negative 
        % xValues displacement,etc. and what we would have is a series of lines of points, and
        % in publications, what you usually get is that the points with
        % either the highest or
        % lowest values get the positions that are more outside, and then the lowest or 
        % the highest get the central position. In this function, the lowest values of yValues
        % get the external positions, and the highest values get the central positions, making
        % some sort of eyebrow or sad mouth shape)
        % The range of xb gets bigger as the number of points inside the
        % group increases due to the action of distmaker, which gets
        % smaller as number of points increases.
        xb=linspace(1-Width+distmaker,1+Width-distmaker,size(yValues(subsetind),1))-1;
        
        % Since xb is symmetrical it's easy to add the more extreme values
        % of xb to the xValues with the lower yValues. oddie is needed to
        % index properly depending on wether the number of elements is even
        % or odd. If we sorted the values of xb by their absolute value,
        % their index would be (1,end,2,end-1,3,end-2,etc.) If you think of
        % an example or you debug the function and observe the values of xb
        % it's very easy to understand the indexing done here. Also, it's
        % here where oddie is useful
        
        xSubset(1:2:end-oddie)=xSubset(1:2:end-oddie)+xb(1:round(end/2)-oddie)';
        xSubset(2:2:end)=xSubset(2:2:end)-xb(1:round(end/2)-oddie)';
        xValues(subsetind)= xSubset;    
        
        %This following code line is very useful to see how the groups are made,
        %for watching it, comment the scatter lines in the end, and
        %uncomment this line, however groups with one point won't be shown
        %scatter(xValues(subsetind),yValues(subsetind),PointStyle)
    end
    % advance in the while loop as long as the cutUp value is out of the
    % range of the points
    keep_going=~(cutUp>max(yValues) || cutUp<min(yValues));
    cutUp=cutUp+steps*cuts;
    subsetind= yValues<cutUp & yValues>cutUp-cuts;
end
end

%% Drawing each subset of points
% This is to overcome the fact that MarkerEdgeColor and MarkerFaceColor can
% be either strings,1x3 vectors or nx3 arrays

if EdgeIsArray
    EdgeIndex=i;
else
    EdgeIndex=1;
end

if FaceIsArray
    FaceIndex=i;
else
    FaceIndex=1;
end

scatter(xValues,yValues,PointSize,PointStyle,'MarkerEdgeColor', MarkerEdgeColor(EdgeIndex,:),'MarkerFaceColor',MarkerFaceColor(FaceIndex,:),'LineWidth',LineWidth)

%If we want lines to represent the SEM and the std
if strcmp(Whiskers,'lines')
    %plot the mean line
    yMean=mean(yValues);
    plot([i-Width/WhiskersWidthRatio i+Width/WhiskersWidthRatio],[yMean yMean],'Color',MeanColor(MeanIndex),'LineWidth',1.5)
    
    %plot the sd 
    yStd=std(yValues);
    plot([i-Width/WhiskersWidthRatio*0.5 i+Width/WhiskersWidthRatio*0.5],[yMean+yStd yMean+yStd],'Color',StdColor(StdIndex,:),'LineWidth',1.5)
    plot([i-Width/WhiskersWidthRatio*0.5 i+Width/WhiskersWidthRatio*0.5],[yMean-yStd yMean-yStd],'Color',StdColor(StdIndex,:),'LineWidth',1.5)
    plot([i i],[yMean-yStd yMean+yStd],'Color',StdColor(StdIndex,:))
    %plot the conf. interval of the mean
    ySem=yStd/sqrt(size(yValues,1));
    yCI=ySem*1.96;
    plot([i-Width/WhiskersWidthRatio*0.7 i+Width/WhiskersWidthRatio*0.7],[yMean+yCI yMean+yCI],'Color',SEMColor(SEMIndex,:),'LineWidth',1.5)
    plot([i-Width/WhiskersWidthRatio*0.7 i+Width/WhiskersWidthRatio*0.7],[yMean-yCI yMean-yCI],'Color',SEMColor(SEMIndex,:),'LineWidth',1.5)
    plot([i i],[yMean-yCI yMean+yCI],'Color',SEMColor(SEMIndex,:))

end

xPositions(yValues_nanIndex,i)=xValues(x_SortingIndex);
end
set(gca,'xtick',1:size(data,2))

if ~isempty(Label)
set(gca,'xticklabel',Label)
pbaspect([size(data,2)+1,4,1])
end
if ~HoldWasOn
hold off
end
end









