function [FeatsMatrix, scaling] = scale_data(FeatsMatrix, type, scaling_in, DIM) 
  DefaultVal('*DIM', '2');
  FeatsMatrix = full(FeatsMatrix);
  scaling = []; 
  
  if(any(any(any(isinf(FeatsMatrix)))))
      display('Found infinity values in Features!! Setting them to 0' )
      FeatsMatrix(isinf(FeatsMatrix)) = 0;
  end
  
  if(any(any(any(isnan(FeatsMatrix)))))
      display('Found NaN values in Features!! Setting them to 0' )
      FeatsMatrix(isnan(FeatsMatrix)) = 0;
  end
  
  if(isempty(type) || strcmp(type, 'none'))
      scaling = [];
      return;
  end 

  if(strcmp(type, 'zero_one'))      
      if(nargin == 2)
          % from libsvm's faq
          scaling = zero_one_scaling(FeatsMatrix);
      else
          scaling.to_subtract = scaling_in.to_subtract;
          scaling.to_divide = scaling_in.to_divide;
      end
      [FeatsMatrix] = normalize(FeatsMatrix, scaling);
  elseif(strcmp(type,'norm_1'))
      %%%% buggy previous version
      %FeatsMatrix = FeatsMatrix./(repmat(sum(FeatsMatrix), size(FeatsMatrix,1), 1)+eps);
      
      %%% correct version %%%
      FeatsMatrix = FeatsMatrix./(repmat(sum(abs(FeatsMatrix)), size(FeatsMatrix,1), 1)+eps);
      scaling = [];
  elseif(strcmp(type,'norm_2'))
      FeatsMatrix = FeatsMatrix./(repmat(vnorm(FeatsMatrix,1), size(FeatsMatrix,1), 1)+eps);
      scaling = [];
  elseif(strcmp(type, 'zscore'))
      if(nargin == 3)
          [FeatsMatrix] = normalize(FeatsMatrix, scaling_in);
          scaling = scaling_in;
      else
          [FeatsMatrix, m, sigma] = zscore(FeatsMatrix,0, DIM);
          scaling.to_subtract = m;
          scaling.to_divide = 1./(sigma');
      end
  else
      error('no such scaling type');
  end
end

