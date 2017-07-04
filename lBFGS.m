function dir= lBFGS(met, MEM, grad, mod)
% Limited memory quasi-newton method for the approximation of
% the Hessian based on the hystory of the gradient...
% Taken from D. Dagnino (CSI), originally from Nocedal's pseudo-code.
% it_num = ITERATION NUMBER
% grad = MATRIX CONTAINING THE MEM LAST GRADIENTS, WHERE MEM IS THE MEMORY OF THE SYSTEM
% mod = MATRIX CONTAINING THE MEM LAST MODELS.
% dir = SEARCH DIRECTION OUTPUT

% In the first loop, Q, initially equal to the gradient at the current iteration,
% is modified so that the gradient variations are subtracted with a weigth
% inversely proportional to the scalar product deltaGrad*deltaModel. It
% is basically a damping which acts in a sharper way if the model updates are
% less 'parallel' to the gradient changes.   

% Secondly, an initial guess of the Hessian is computed as the ratio bewtween
% the scalar product deltaModel*deltaGrad and the norm of deltaGrad. The size
% of the Hessian is basically proportional to the scalar product of the model
% variation and the gradient variation normalised by the gradient variation lenght.

% The second loop updates q adding a value proportional to deltaModel, scaled
% by the difference between (deltaModel*Q/deltaGrad*deltaModel)-(deltaGrad*Q/deltaGrad*deltaM)
% Basically, q is corrected when the variation of the model is non-parallel to the variation
% of the gradient, privilenging coherency in the model variation throughout the iterations...

% 	MEM=size(grad,2); % How many previous iterations we use?
	MEM_TOT=size(grad,2); % How many previous iterations we use?
	
	if (MEM~=1) 
    alpha=zeros(MEM-1,1);
  end
	
	if (MEM==1 || met==0) % That's equivalent to a steepest descent...

		dir=-grad(:,MEM_TOT);

	else

	%..First Loop

	q=grad(:,MEM_TOT); % The current gradient is stored at the end of the matrix

% 		for K=1:MEM-1 % We compute the summation and inner products to adjust the descent direction
		for K=(MEM-1):-1:1 % We compute the summation and inner products to adjust the descent direction
			% according to the gradient and model hystory...		

% 			dg=-grad(:,MEM-K)+grad(:,MEM-K+1); % Variation of the gradient

% 			dm=-mod(:,MEM-K)+mod(:,MEM-K+1); % Variation of the model.

			dg=-grad(:,MEM_TOT-(MEM-K))+grad(:,MEM_TOT-(MEM-K)+1); % Variation of the gradient

			dm=-mod(:,MEM_TOT-(MEM-K))+mod(:,MEM_TOT-(MEM-K)+1); % Variation of the model.

			num=dm'*q;
			
			den=dg'*dm;

% 			if den==0
% 
% 			alpha(K)=0;	
% 
% 			else

			alpha(K)=num/den;

% 			end

			q = q - alpha(K)*dg;
		
		end

		dg=-grad(:,MEM_TOT-1)+grad(:,MEM_TOT); % Local gradient variation rate estimate
		
		dm=-mod(:,MEM_TOT-1)+mod(:,MEM_TOT);

		H_0=(dm'*dg)/(dg'*dg); % First estimate of the Hessian
		
		q=H_0*q;

	%..Second loop


		for K=1:(MEM-1)

% 			dg=-grad(:,MEM-K)+grad(:,MEM-K+1);
% 			
% 			dm=-mod(:,MEM-K)+mod(:,MEM-K+1);

			dg=-grad(:,MEM_TOT-(MEM-K))+grad(:,MEM_TOT-(MEM-K)+1);
			
			dm=-mod(:,MEM_TOT-(MEM-K))+mod(:,MEM_TOT-(MEM-K)+1);
      
			num=dg'*q;
		
			den=dg'*dm;

% 			if den==0
% 			
% 			beta=0;
% 			
% 			else
% 	
% 			beta=num/(den+eps);
% 		
% 			end

			beta = num/den;

			q = q + (alpha(K)-beta)*dm;

		end
	
		dir=-q;

% 		dir=dir-mean(dir);% Check!

	end

end

