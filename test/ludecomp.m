% Decomposicao LU sem pivoteamento
% Autor: Afonso Paiva
% Data: 07/02/2010
% Verao 2010 -- ICMC-USP
% Input: Matriz quadrada A(nxn).
% Output: Matrizes triangulares inferior L e superior U de A = LU.

function [fillin]=ludecomp(A)
n=size(A,1);
L=sparse(eye(n)); U=sparse(zeros(n));

for i=1:n
	if (A(i,i)!=0)
		A(i,i) = n+1;
	end
end

for k=1:n
    for j=k:n
        U(k,j)= A(k,j)-L(k,1:k-1)*U(1:k-1,j);
    end
    for i=k+1:n
        L(i,k)=(A(i,k)-L(i,1:k-1)*U(1:k-1,k))/U(k,k);
    end
end

L(L!=0) = 1;
U(U!=0) = 1;

L = L+U;
L(L>1) = 1;

A(A>1)=1;
A = abs(A-L);

fillin = nnz(A);
