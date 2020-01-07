#### ______GROUP-MEMBERS_______ ####
####_______   Name________________________id
####_______1 yonas yifter ________________209/09
####_______2 tensay g/tsadik______________111/09
####_______3 Kibrom Takelle______________ 043/09
####_______4 weldegebreal bilad___________208/09
####_______5 Zeberhe Hagos ________________210/09
#Computational Methods Group Assignment one
#diverge for diagonally dominant
#strict diagonal dominance is a necessary condition for convergence of the
#Jacobi method is get updated after the entire calculation is done
# Xi=(bi-aijXj)/aii;
#finally update the new Xi
#in this function either you can apply all the parameters to the function or you can feed 
# as soon as it asks you to enter the desired matrix and other parameters

function [x,error] = gause_seidel_func(A,B,X,iter,err)
  # if you want parametrs to be feed from the keyboard
  [row,col]=size(A);
  if row ~= col
     disp('no unique solution for this system of linearEquation');
     disp('the matrice is not square matrice');
     return;
  endif
  [rowB,colB]=size(B);
  if ((rowB ~=row) || (colB ~= 1))
    disp('there is problem with B matrices');
    disp('either it havenot the same row as the matrices before or not column vector.');
    return;      
  endif
  [rowX,colX]=size(X);
  if ((rowX ~=row) || (colX ~= 1))
    disp('the initial guess does not fit with the matrices');
    return;
  endif
  if isnumeric(iter) == 0
    disp('the number of iteration is not number');
  endif
  if isnumeric(err) == 0
    disp('the error tolerance provided is not valid number');
  endif  
  % check for the diagonal dominance
 n=length(A);
 dom=0;
 diago=0;
 for i=1:n
   for j=1:n
     if i~=j
       dom=dom+abs(A(i,j));
     else
       diago=diago+abs(A(i,j));
     endif
   endfor
    if diago<dom
   diago=0;dom=0;
   disp('adust your matrices to fit convergence interchange the the rows to have diagonally dominance..');
   return;
 endif
endfor

 
 # lets iterate and find the solution
 #lets use initial guess for X
new=zeros(row,1);
#now it is time to compute the linear equation
#using the jacobi-algorithm
forError=zeros(row,1);
for k = 1:iter
    for i = 1:n
        sigma = 0;
        for j = 1:n
            if i~=j
                sigma = sigma + A(i,j)*X(j);
            endif  
        endfor
        X(i) = (B(i)-sigma)/A(i,i);  
    endfor
    for i=1:n
      for j=1:n
        summation(j)=(A(i,j)*X(j));
      endfor
        forError(i)=sum(summation)-B(i);  
    endfor
    err1=norm(forError);
    #compute the error in every iteration and stop if we reach 
    # desired error
  if err1<err
     disp('+================+=============+');
    fprintf('the iteration level is in %d\n',k);
    disp('');
    disp('the solution for this Equation is:');
    x=X
    disp("error is: ");
    error=err1
    disp("========================");
    return;
  endif
endfor
if err1>err
 fprintf('the iteration level is in %d\n',k);
 fprintf("we need more iteration to carry out in order to reach the optimum error.");
endif
 x=X
 disp("error is: ");
 error=err1
 disp("========================");
endfunction





