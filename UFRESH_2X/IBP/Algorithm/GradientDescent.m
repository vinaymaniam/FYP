function x = GradientDescent(lhs, rhs, initialGuess)
   maxIter = 100;
   iter = 0;
   eps = 0.01;
   
   x = initialGuess;
   res = lhs' * (rhs - lhs * x);
   mse = res' * res;
   mse0 = mse;
   while (iter < maxIter && mse > eps^2 * mse0)
       res = lhs' * (rhs - lhs * x);
       x = x + res;
       mse = res' * res;
       iter = iter + 1;
   end

end