%OMPL algorithm
function [x_sharp,k] = OMPL(A,y,s,M,m,N)
residual=y;S=[];
nb_iterations=sum(s,'all');
for k=1:nb_iterations+nb_iterations/2
    vals=A'*(residual);
    [s,j]=argmax(vals,s,M,S);
    S=union(S,j);
    S=unique(S,'rows',"sorted");
    A_s=restrict_matrix(A,S',m);
    z=A_s\y;
    x_sharp=unrestrict_vector(z,S',N);
    %remove smallest entry of z here when sum(s,'all')=0, to decide which
    %is the smaller, do A*resid, remoce that one
    residual=y-A*x_sharp;
end
end