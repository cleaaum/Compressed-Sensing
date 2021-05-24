%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%  CoSaMP ALGORITHM  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INPUTS: Sensing Matrix A, linear measurments y, dimensions of A (mxN),
%sparsity s, solution x (to calculate rel err), nb_terations we want to
%perform
%OUTPUTS: X^# (approxiamtion vector), and current_iteration
function [x_sharp, residual_plot, current_iteration]=COSAMP(A,y,m,N,s,nb_iterations,x)
x_sharp=zeros(N,1);
residual_plot=zeros(nb_iterations,1);
% for k=1:nb_iterations
current_iteration = 0;
error=10;
while (error > 10e-6) && (current_iteration<120)
    
    %CoSaMP 1: U^(n+1) = supp(x^n) Union L_2s(A*(y-Ax^n))
    residual=y-A*x_sharp;
    U = union(find(x_sharp),L_s(A'*residual,2*s));
    U = unique(U,'rows',"sorted");
    
    %CoSaMP 2: u^(n+1) in argmin{|y-Az|_2 : supp(z) subset U^(n+1)}
    A_U = restrict_matrix(A,U,m);
    
    z = A_U\y;
    %     z = lsqr(A_U,y);
    u = unrestrict_vector(z,U,N);
    
    %CoSaMP 3: x^(n+1) = H_s(u^(n+1))
    H = L_s(u,s);
    x_sharp = H_s(H,u,N);
    %     residual_plot(k,1) = norm(A*x_sharp-y);
    error=norm(x-x_sharp);
    current_iteration=current_iteration+1;
    %     residual_plot(k,1) = norm(x-x_sharp)/norm(x);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%  CoSaMPL ALGORITHM  %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INPUTS: Sensing Matrix A, linear measurments y, dimensions of A (mxN),
%sparsity vector s, levels M, and solution x (to calculate rel err)
%OUTPUTS: X^# (approxiamtion vector), and current_iteration
function [x_sharp,current_iteration]= COSAMPL(A,y,m,N,s,M,x)
x_sharp=zeros(N,1);
current_iteration = 0;
error=10;
while (error > 10e-6) && (current_iteration<120)
    
    %CoSaMP 1: U^(n+1) = supp(x^n) Union L_2s(A*(y-Ax^n))
    residual=y-A*x_sharp;
    U = union(find(x_sharp),L_s_M(A'*residual,2*s,M));
    U = unique(U,'rows',"sorted");
    
    %CoSaMP 2: u^(n+1) in argmin{|y-Az|_2 : supp(z) subset U^(n+1)}
    A_U = restrict_matrix(A,U,m);
    
    z = A_U\y;
    %     z = lsqr(A_U,y);
    u = unrestrict_vector(z,U,N);
    
    %CoSaMP 3: x^(n+1) = H_s(u^(n+1))
    H = L_s_M(u,s,M);
    x_sharp = H_s_M(H,u,N);
    %     residual_plot(k,1) = norm(A*x_sharp-y);
    error=norm(x-x_sharp);
    current_iteration=current_iteration+1;
    %     residual_plot(k,1) = norm(x-x_sharp)/norm(x);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%  OMP ALGORITHM  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INPUTS: Sensing Matrix A, linear measurments y, dimensions of A (mxN),nb_terations we want to
%perform
%OUTPUTS: X^# (approxiamtion vector), and current iteration k, and a plot
%of the residual

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%  OMPL ALGORITHM  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INPUTS: Sensing Matrix A, linear measurments y, dimensions of A (mxN),
%sparsity vector s, and levels M.
%OUTPUTS: X^# (approxiamtion vector), and current iteration k
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%  AUX FUNCTIONS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_residuals_vs_iterations(residual,type)
semilogy(residual)
ylim([10e-16 10e-6])
title('Residual as a function of k using')
legend('2-norm of residual')
xlabel('k iterations')
ylabel('$\|y-Ax^{(k)}\|_2$','Interpreter','latex')
end

%L_{s,M} operator, selects s largest entries, of each level
function max_index_set = L_s_M(z,s,M)
M=union(0,M);
max_index_set =zeros(sum(s,'all'),1);
counter = 1;
for i=2:length(M)
    lower_bound = (M(i-1)+1);
    upper_bound= M(i);
    current_level=abs(z(lower_bound:upper_bound));
    current_sparsity_level=s(i-1);
    for k=1:current_sparsity_level
        [Max,j]=max(current_level);
        max_index_set(counter,1)=j+M(i-1);
        current_level(j)=0;
        counter = counter+1;
    end
    
end
max_index_set=sort(max_index_set, 'ascend');
end

%argmax function for in levels, used in OMP only
function [s,j] = argmax(z,s,M,S)
M=union(0,M);
for i=2:length(M)
    lower_bound = (M(i-1)+1);
    upper_bound = M(i);
    if(sum(s,'all')==0)
        z_restricted=zeros(length(S),1);
        for a=1:length(S)
            z_restricted(a)=z(S(a));
        end
        %         z_removed_index=H_s(S,z,256);
        [m,index]=min(abs(z_restricted));
        %         z_restricted(index)=0;
        index_to_remove=S(index);
        for n=1:length(M)
            if(index_to_remove>M(n))
            else
                s(n-1)=s(n-1)+1;
                break
            end
        end
    end
    if(s(i-1) == 0)
        z(lower_bound:upper_bound)=0;
    end
end
[maximum,j]=max(abs(z));
for i=1:length(M)
    if(j>M(i))
    else
        s(i-1)=s(i-1)-1;
        break
    end
end
end

%checks coherance of the matrix A
function coherance = check_coherance(A)
dimension_A=size(A);
j=1;
total_inner_products=nchoosek(dimension_A(2),2);
inner_products=zeros(total_inner_products,1);
for i=1:dimension_A(2)-1
    for k=i+1:dimension_A(2)
        inner_products(j,1)=dot(A(:,i),A(:,k));
        j=j+1;
    end
end
[M,index]=max(abs(inner_products));
coherance=M;
end

%check if matrix A staisfies the RIP of order s with a given RIC
function RIC=check_RIP(RIC,A,x)
if (1-RIC)*norm(x)^2 <=norm(A*x)^2 && norm(A*x)^2 <=  (1+RIC)*norm(x)^2
    RIC=true;   %(1)
else
    RIC=false;  %(0)
end
end
%L_{s} operator, selects s largest entries
function max_index_set = L_s(z,s)
max_index_set =[];
for i=1:s
    z=abs(z);
    [M,j]=max(z);
    max_index_set(length(max_index_set)+1,1)=j;
    z(j)=0;
end
% max_index_set=unique(max_index_set,'rows',"sorted");
max_index_set=sort(max_index_set, 'ascend');
end

%Hard thresholding operator
function x = H_s(H,u,N)
x=zeros(N,1);
for i=1:length(H)
    x(H(i))=u(H(i));
end
end

%fn to restrict a matrix
function matrix_to_restrict = restrict_matrix(matrix, index_set,m)
% index_set=index_set';
matrix_to_restrict = zeros(m,length(index_set));
for i=1:length(index_set)
    matrix_to_restrict(:,i)=matrix(:,index_set(i,1));
end
end

%function to unrestrict a vector
function vector_to_unrestrict =unrestrict_vector(vector,index_set,N)
vector_to_unrestrict=zeros(N,1);
for i=1:length(index_set)
    vector_to_unrestrict(index_set(i,1))=vector(i,1);
end
end

%fn to generate a random signal
function y=generate_signal(m,N)
y=zeros(m,1);
for i=1:m
    y(i,1)=sin((12*pi)*((i-1)/N- 1.1)) + 2*max(1-abs(i/6 -7),0);
end

end

%Hard thresholding operator for in level model
function x = H_s_M(H,u,N)
x=zeros(N,1);
for i=1:length(H)
    x(H(i))=u(H(i));
end
end


