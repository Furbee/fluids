function [ p ] = project2( rhs, w, h, dt, rho, dxy, limit)
%PROJECT2 Summary of this function goes here
%   Detailed explanation goes here

p = zeros(w*h,1);

        scale = dt/(rho*dxy*dxy);
        
        maxDelta = 0.0;
        for iter = 1:limit
            maxDelta = 0.0;
            for y = 1:h
                for x = 1:w
                    idx = getIdx(x,y,w);
                    
                    diag = 0.0;
                    offDiag = 0.0;
                    
                    if x > 1
                        diag    = diag + scale;
                        offDiag = offDiag - scale*p(idx - 1);
                    end
                    if y > 1
                        diag    = diag + scale;
                        offDiag = offDiag - scale*p(idx - w);
                    end
                    if x < w
                        diag    = diag + scale;
                        offDiag = offDiag - scale*p(idx + 1);
                    end
                    if y < h
                        diag    = diag + scale;
                        offDiag = offDiag - scale*p(idx + w);
                    end

                    newP = (rhs(idx) - offDiag)/diag;
                    
                    maxDelta = max(maxDelta, abs(p(idx) - newP));
                    
                    p(idx) = newP;
                end
            end

            if maxDelta < 1e-5
                fprintf('Exiting solver after %d iterations, maximum change is %f\n', iter, maxDelta);
                return;
            end
        end
        
        fprintf('Exceeded budget of %d iterations, maximum change was %f\n', limit, maxDelta);


end

