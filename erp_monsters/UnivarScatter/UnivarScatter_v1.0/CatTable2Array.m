function [ output_array, Label] = CatTable2Array( input_table )
%--------------------------------------------------------------------------
%   This function converts a nx2 table (input_table) with one
%   string/categorical variable and a numerical variable into a numerical
%   array in which the columns correspond to the categories described by
%   the categorical/string variable of the table, and the values are the
%   ones of the numerical variable. If the different columns have different
%   number of observations, empty positions of the array are filled with
%   nans.
%   version: <1.00> from 30/11/2015
% 
%       Manuel Lera Ramírez: manulera14@gmail.com
%
%       Developed during Master Thesis Internship in Marcos González-Gaitán's Lab
%       Departments of Biochemistry and Molecular Biology, Sciences II, 
%       30 Quai Ernest-Ansermet, CH-1211 Geneva 4, Switzerland
%--------------------------------------------------------------------------
if istable(input_table) && size(input_table,2)==2
    Str_var_ind= [iscellstr(input_table{:,1}) iscellstr(input_table{:,2})];
    if sum([iscellstr(input_table{:,1}) iscellstr(input_table{:,2})])~=0;
        Cat_var=categorical(input_table{:,Str_var_ind});
        Cat_var_ind=Str_var_ind;
    else
        Cat_var_ind= [iscategorical(input_table{:,1}) iscategorical(input_table{:,2})];
        Cat_var=input_table{:,Cat_var_ind};
    end
    
    Num_var_ind= [isnumeric(input_table{:,1}) isnumeric(input_table{:,2})];
    if Cat_var_ind+Num_var_ind ~=[1 1]
        error('Wrong table format, it should be a nx2 table with one numerical variable and one string/categorical variable')
    end
    Num_var=input_table{:,Num_var_ind};
    output_array=nan(max(countcats(Cat_var)),length(unique(Cat_var)));
    j=1;
    for cat_i=unique(Cat_var)'
        output_array(1:sum(Cat_var==cat_i),j)=Num_var(Cat_var==cat_i);
        j=j+1;
    end
    Label=cellstr(char(unique(Cat_var)))';
else
    error('Wrong input format, it should be a nx2 table with one numerical variable and one string/categorical variable') 
end

end

