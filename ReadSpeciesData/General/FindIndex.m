function [Index] = FindIndex(A,B)

% Routine finds index of components of A in vector B

if (length(A)>length(B))
   fprintf('\nERROR: vector A larger than vector B!\n'); 
   Index = [ ];
   return
else
    Index = zeros(1,length(A));
    for i = 1:length(A)
        for j = 1:length(B)
            if strcmpi(A(i).Name,B(j).Name) 
                Index(i) = j; 
            end
        end
        if Index(i) == 0
            fprintf('\nERROR: component %s not found in reaction mechanism!\n',A(i).Name);
            return
        end
    end
end
