%Phase Transition Plot - CoSaMPL 4 levels
clear all;
close all;
N=256;
probability=zeros(15,23); %matrix with number of iterations until signal recovered
sparsities=4:4:60;
m_dim=20:10:240;
for z=1:50 
    z
    for k=1:23 %m
        m=m_dim(k);
        for j=1:15%s
            sparsity=sparsities(j);
            A = randn(m,N)/sqrt(m); %normally distributed random (Gaussian) matrix
            x = zeros(N,1);
            while (nnz(x(1:N/4))<(sparsity/2))
                x(randi([1,N/4],1)) = randi([1,10])/sqrt(m);
            end
            while (nnz(x(N/2+1:3*N/4))<(sparsity/2))
                x(randi([N/2+1,3*N/4],1)) = randi([1,10])/sqrt(m);
            end
            
            y=A*x;
            [x_sharp_cosampl,iteration]=COSAMPL(A,y,m,N,[sparsity/2,0,sparsity/2,0],[N/4,N/2,3*N/4,N],x);
            rel_error=norm(x-x_sharp_cosampl)/norm(x);
            if rel_error < 10e-4
                probability(16-j,k)=probability(16-j,k)+1;
            end
        end
    end
end
%phase transition plot
probability=probability/50;
imagesc(1-probability)
title('CoSaMPL-4 levels')
xlabel('size of m');
ylabel('total sparsity');
set(gca,'Xtick',1:12,'XTickLabel',{'20', '40', '60', '80','100', '120', '140', '160','180','200','220','240'})
set(gca,'Ytick',1:15,'YTickLabel',{'60', '56', '52', '48','44', '40', '36', '32','28','24','20','16','12','8','4'})
colorbar
colormap(gray)

%Phase Transition Plot - CoSaMPL 2 levels
clear all;
close all;

N=256;
probability=zeros(15,23); %matrix with number of iterations until signal recovered
sparsities=4:4:60;
m_dim=20:10:240;
for z=1:50
    z
    for k=1:23%m
        m=m_dim(k);
        for j=1:15%s
            sparsity=sparsities(j);
            A = randn(m,N)/sqrt(m); %normally distributed random (Gaussian) matrix
            x = zeros(N,1);
            while (nnz(x(1:N/2))<(3*sparsity/4))
                x(randi([1,N/2],1)) = randi([1,10])/sqrt(m);
            end
            while (nnz(x(N/2+1:N))<(sparsity/4))
                x(randi([N/2+1,N],1)) = randi([1,10])/sqrt(m);
            end
            
            y=A*x;
            [x_sharp_cosampl,iteration]=COSAMPL(A,y,m,N,[3*sparsity/4,sparsity/4],[N/2,N],x);
            rel_error=norm(x-x_sharp_cosampl)/norm(x);
            if rel_error < 10e-4
                probability(16-j,k)=probability(16-j,k)+1;
            end
        end
    end
end
%phase transition plot
probability=probability/50;
imagesc(1-probability)
title('CoSaMPL-2 levels')
xlabel('size of m');
ylabel('total sparsity');
set(gca,'Xtick',1:12,'XTickLabel',{'20', '40', '60', '80','100', '120', '140', '160','180','200','220','240'})
set(gca,'Ytick',1:15,'YTickLabel',{'60', '56', '52', '48','44', '40', '36', '32','28','24','20','16','12','8','4'})
colorbar
colormap(gray)