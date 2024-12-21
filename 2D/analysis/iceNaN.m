%%% set ice region to NaNs

function B = iceNaN(A,Nx,Ny,Nr,icetopo,zz)
    if ndims(A) ==4
        for i=1:Nx
            for j =1:Ny
                for k = 1:Nr
                    if zz(k)>icetopo(j)
                        if k>1
                            A(i,j,k-1,:) = NaN;
                        else
                            A(i,j,k,:) = NaN;
                        end
                    end
                end
            end
        end
        B = A;
    end
    if ndims(A) == 3
        for j =1:Ny
            for k = 1:Nr
                if zz(k)>icetopo(j)
                    if k>1
                        A(j,k-1,:) = NaN;
                    else
                        A(j,k,:) = NaN;
                    end
                end
            end
        end
        B = A;
    end

end % function
