# A recursive implementation of Householder transformations using Julia

using LinearAlgebra  # using norm function


function Householder(A, trunc=true, i=1, HH=[])
    """recursively forms Householder matrices H_i using A and i and 
    adds to array of Householder matrices HH returns R and HH
    trunc means to truncate floating points near zero for clean matrix
        """
    if i < length(A[1,:])
        M = A[i:end,i:end]  # take the elements we are working on in a new matrix
        m_ = M[:,1]  # the first column of M
        for j in 1:i  # add padding where required
            if j < i
                pushfirst!(m_,0.)
            end
        end
        θ = norm(m_)
        size = length(A[:,1])
        H = Matrix(I, size, size)  # initialize H as identity matrix
        e_ = H[:,i]  # creating base vector e_ by taking first column of H
        v_ = m_ - θ*e_  # define v_
        β = 2 / norm(v_)^2  # find β
        H = H - β*v_*transpose(v_)  # create Householder reflector
        A = H * A  # apply reflector
        pushfirst!(HH, H)  # add matrix to our array
        if trunc
            A[abs.(A) .≤ 1e-15] .= 0 end  # tidy up zero floats
        Householder(A, trunc, i+1, HH)  # rerusion baby!
    else
        return A, HH
    end
end


A = [0.8 4 3 4; 0 -3 0 3; -0.6 2 2 0; 0 0 0 4]  # A test matrix

Householder(A)[1]  # outputs R