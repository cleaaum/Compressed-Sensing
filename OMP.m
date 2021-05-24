%OMP algorithm
function [x_sharp, residual_plot,k] = OMP(A,y,N,m,nb_iterations)
residual_plot=zeros(nb_iterations,1);
residual=y;S=[];
for k=1:nb_iterations
    vals=A'*(residual);
    [M,j] = max(abs(vals));
    S(length(S)+1,1)=j;
    S=unique(S,'rows',"sorted");
    A_s=restrict_matrix(A,S,m);
    z=A_s\y;
    residual=y-A_s*z;
    residual_plot(k,1)=norm(residual);
end

x_sharp=unrestrict_vector(z,S,N);

end