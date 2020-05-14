function [T,Q] = Data_interp(run,Tmin,Tmax,pts)
% Takes a run and interpolates data to make new set at set temperatures.
% Tmin, Tmax, pts - the inputs to the linspace function

T = transpose(linspace(Tmin,Tmax,pts));
Q = zeros(length(T),1);

for i = 1:length(T)
    check = find(T(i) == run(:,1),1);
    if isempty(check)
        j = find(run(:,1)<T(i),1,'last');
        Q(i) = run(j,2)+(T(i)-run(j,1))*...
            (run(j+1,2)-run(j,2))/(run(j+1,1)-run(j,1));    
    else
        Q(i) = run(check,2);
    end
end

end